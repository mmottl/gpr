open Lacaml.Impl.D

open Interfaces
open Utils

(* TODO: dimension sanity checks *)

module type Spec = sig
  module Kernel : Kernel

  val get_sigma2 : Kernel.t -> float
  val jitter : float
end

module type Loc = sig val loc : string end

module type Sig = functor (Spec : Spec) ->
  Inducing_input_gpr with module Kernel = Spec.Kernel

module Make_common (Spec : Spec) = struct
  include Spec

  open Kernel

  module Inducing = struct
    type t =
      {
        kernel : Kernel.t;
        inputs : inputs;
        km : mat;
        km_chol : mat;
        log_det_km : float;
      }

    let calc kernel inputs =
      let km = Kernel.upper kernel inputs in
      (* TODO: copy upper triangle only *)
      let km_chol = Mat.copy km in
      potrf ~jitter km_chol;
      let log_det_km = log_det km_chol in
      {
        kernel = kernel;
        inputs = inputs;
        km = km;
        km_chol = km_chol;
        log_det_km = log_det_km;
      }
  end

  module Induced = struct
    type t =
      {
        inducing : Inducing.t;
        input : input;
        k_m : vec;
      }

    let calc inducing input =
      let kernel = inducing.Inducing.kernel in
      {
        inducing = inducing;
        input = input;
        k_m = Kernel.evals kernel input inducing.Inducing.inputs;
      }

    let get_kernel t = t.inducing.Inducing.kernel
  end

  module Induceds = struct
    type t =
      {
        inducing : Inducing.t;
        inputs : inputs;
        kmn : mat;
      }

    let calc inducing inputs =
      let kernel = inducing.Inducing.kernel in
      let kmn = Kernel.cross kernel inducing.Inducing.inputs inputs in
      {
        inducing = inducing;
        inputs = inputs;
        kmn = kmn;
      }

    let get_kernel t = t.inducing.Inducing.kernel
  end

  let calc_basis_2_k basis_chol ~k =
    (* TODO: consider symmetric matrices *)
    let b_2_k = Mat.copy k in
    trtrs ~trans:`T basis_chol b_2_k;
    b_2_k

  (* Compute square root of Nystrom approximation, and diagonal of
     marginal variances *)
  let nystrom_2_marginals induceds =
    let inducing = induceds.Induceds.inducing in
    let kernel = Induceds.get_kernel induceds in
    let kn_diag = Kernel.diag_vec kernel induceds.Induceds.inputs in
    let kmn = induceds.Induceds.kmn in
    let km_2_kmn = calc_basis_2_k inducing.Inducing.km_chol ~k:kmn in
    kn_diag, km_2_kmn, kmn

  module Common_model = struct
    type t =
      {
        induceds : Induceds.t;
        kn_diag : vec;
        km_2_kmn : mat;
        lam_diag : vec;
        inv_lam_sigma2_diag : vec;
        b_chol : mat;
        neg_log_likelihood : float;
      }

    let calc induceds =
      let inducing = induceds.Induceds.inducing in
      let kernel = Induceds.get_kernel induceds in
      let sigma2 = get_sigma2 kernel in
      let kn_diag, km_2_kmn, kmn = nystrom_2_marginals induceds in
      let kmn_= Mat.copy kmn in
      let n_inputs = Vec.dim kn_diag in
      let lam_diag = Vec.create n_inputs in
      let inv_lam_sigma2_diag = Vec.create n_inputs in
      let rec loop log_det_lam_sigma2 i =
        if i = 0 then log_det_lam_sigma2
        else
          let kn_diag_i = kn_diag.{i} in
          (* TODO: optimize ssqr and col *)
          let qn_diag_i = Vec.ssqr (Mat.col km_2_kmn i) in
          let lam_i = kn_diag_i -. qn_diag_i in
          lam_diag.{i} <- lam_i;
          let lam_sigma2_i = lam_i +. sigma2 in
          let inv_lam_sigma2_i = 1. /. lam_sigma2_i in
          inv_lam_sigma2_diag.{i} <- inv_lam_sigma2_i;
          (* TODO: optimize scal col *)
          scal (sqrt inv_lam_sigma2_i) (Mat.col kmn_ i);
          loop (log_det_lam_sigma2 +. log lam_sigma2_i) (i - 1)
      in
      let log_det_lam_sigma2 = loop 0. n_inputs in
      (* TODO: copy upper triangle only *)
      let b_chol = syrk kmn_ ~beta:1. ~c:(Mat.copy inducing.Inducing.km) in
      potrf ~jitter b_chol;
      let log_det_km = inducing.Inducing.log_det_km in
      let log_det_b = log_det b_chol in
      let l1_2 = log_det_b -. log_det_km +. log_det_lam_sigma2 in
      let neg_log_likelihood = (l1_2 +. float n_inputs *. log_2pi) /. 2. in
      {
        induceds = induceds;
        kn_diag = kn_diag;
        km_2_kmn = km_2_kmn;
        lam_diag = lam_diag;
        inv_lam_sigma2_diag = inv_lam_sigma2_diag;
        b_chol = b_chol;
        neg_log_likelihood = neg_log_likelihood;
      }

    let neg_log_likelihood model = model.neg_log_likelihood
    let get_inducing model = model.induceds.Induceds.inducing
    let get_inducing_inputs model = (get_inducing model).Inducing.inputs
    let get_inputs model = model.induceds.Induceds.inputs
    let get_kernel model = (get_inducing model).Inducing.kernel
    let get_sigma2 model = get_sigma2 (get_kernel model)
    let get_kmn model = model.induceds.Induceds.kmn
    let get_lam_diag model = model.lam_diag
    let get_km model = (get_inducing model).Inducing.km
  end

  module Variational_model = struct
    include Common_model

    let calc induceds =
      let
        {
          lam_diag = lam_diag;
          inv_lam_sigma2_diag = inv_lam_sigma2_diag;
          neg_log_likelihood = neg_log_likelihood;
        } as model = calc induceds
      in
      let n = Vec.dim lam_diag in
      let trace_term = ref 0. in
      for i = 1 to n do
        trace_term := !trace_term +. lam_diag.{i} *. inv_lam_sigma2_diag.{i}
      done;
      {
        model with
        neg_log_likelihood = neg_log_likelihood +. 0.5 *. !trace_term;
      }
  end

  module Trained = struct
    type t =
      {
        kernel : Kernel.t;
        inducing_inputs : inputs;
        b_chol : mat;
        inv_b_chol_kmn_y__ : vec;
        neg_log_likelihood : float;
      }

    let calc model ~targets =
      let n_targets = Vec.dim targets in
      let y__ = Vec.create n_targets in
      let inv_lam_sigma2_diag = model.Common_model.inv_lam_sigma2_diag in
      for i = 1 to n_targets do
        y__.{i} <- targets.{i} *. inv_lam_sigma2_diag.{i}
      done;
      let inv_b_chol_kmn_y__ = gemv (Common_model.get_kmn model) y__ in
      let ssqr_y__ = dot ~x:targets y__ in
      let b_chol = model.Common_model.b_chol in
      trsv ~trans:`T b_chol inv_b_chol_kmn_y__;
      let fit_neg_log_likelihood =
        (ssqr_y__ -. Vec.ssqr inv_b_chol_kmn_y__) /. 2.
      in
      {
        kernel = Common_model.get_kernel model;
        inducing_inputs = (Common_model.get_inducing model).Inducing.inputs;
        b_chol = b_chol;
        inv_b_chol_kmn_y__ = inv_b_chol_kmn_y__;
        neg_log_likelihood =
          model.Common_model.neg_log_likelihood +. fit_neg_log_likelihood;
      }

    let neg_log_likelihood trained = trained.neg_log_likelihood
  end

  module Weights = struct
    type t =
      {
        kernel : Kernel.t;
        inducing_inputs : inputs;
        weights : vec;
      }

    let calc trained =
      let weights = copy trained.Trained.inv_b_chol_kmn_y__ in
      trsv trained.Trained.b_chol weights;
      {
        kernel = trained.Trained.kernel;
        inducing_inputs = trained.Trained.inducing_inputs;
        weights = weights;
      }

    let copy t = copy t.weights
  end

  module Mean = struct
    type t = { input : input; mean : float }

    let make ~input ~mean = { input = input; mean = mean }


    let calc w input =
      let module W = Weights in
      let mean =
        Kernel.weighted_eval
          w.W.kernel ~weights:w.W.weights input w.W.inducing_inputs
      in
      make ~input ~mean

    let calc_induced w { Induced.input = input; k_m = k_m } =
      make ~input ~mean:(dot ~x:k_m w.Weights.weights)

    let copy m = m.mean
  end

  module Means = struct
    type t = { inputs : inputs; means : vec }

    let make ~inputs ~means = { inputs = inputs; means = means }

    let calc_inducing { Weights.weights = weights } model =
      make
        ~inputs:(Common_model.get_inputs model)
        ~means:(symv (Common_model.get_km model) weights)

    let calc_inputs { Weights.weights = weights } model =
      make
        ~inputs:(Common_model.get_inputs model)
        ~means:(gemv ~trans:`T (Common_model.get_kmn model) weights)

    let calc w inputs =
      let module W = Weights in
      let means =
        Kernel.weighted_evals
          w.W.kernel ~weights:w.W.weights inputs w.W.inducing_inputs
      in
      make ~inputs ~means

    let calc_induceds w { Induceds.inputs = inputs; kmn = kmn } =
      make ~inputs ~means:(gemv ~trans:`T kmn w.Weights.weights)

    let copy m = copy m.means
  end

  module Variance = struct
    type t = { input : input; variance : float; sigma2 : float }

    let calc_induced model induced =
      let { Induced.input = input; k_m = k_m } = induced in
      let kernel = Common_model.get_kernel model in
      let prior_variance = Kernel.eval_one kernel input in
      let inv_km_chol_k_m = copy k_m in
      let inv_km_chol_k_m_mat = Mat.from_col_vec inv_km_chol_k_m in
      let inducing = induced.Induced.inducing in
      potrs ~factorize:false inducing.Inducing.km_chol inv_km_chol_k_m_mat;
      let km_arg = dot ~x:k_m inv_km_chol_k_m in
      let inv_b_chol_k_m = copy k_m ~y:inv_km_chol_k_m in
      let inv_b_chol_k_m_mat = inv_km_chol_k_m_mat in
      potrs ~factorize:false model.Common_model.b_chol inv_b_chol_k_m_mat;
      let b_arg = dot ~x:k_m inv_b_chol_k_m in
      let explained_variance = km_arg -. b_arg in
      let variance = prior_variance -. explained_variance in
      {
        input = input;
        variance = variance;
        sigma2 = Common_model.get_sigma2 model;
      }

    let copy ?predictive t =
      match predictive with
      | None | Some true -> t.variance +. t.sigma2
      | Some false -> t.variance
  end

  let calc_b_2_k model ~k = calc_basis_2_k model.Common_model.b_chol ~k

  module Variances = struct
    type t = { inputs : inputs; variances : vec; sigma2 : float }

    let make ~inputs ~variances ~model =
      let sigma2 = Common_model.get_sigma2 model in
      { inputs = inputs; variances = variances; sigma2 = sigma2 }

    let calc_inducing model =
      let b_2_km = calc_b_2_k model ~k:(Common_model.get_km model) in
      let m = Mat.dim2 b_2_km in
      let variances = Vec.create m in
      for i = 1 to m do
        (* TODO: optimize ssqr and col *)
        variances.{i} <- Vec.ssqr (Mat.col b_2_km i)
      done;
      make ~inputs:(Common_model.get_inducing_inputs model) ~variances ~model

    let calc_inputs model =
      let variances = copy (Common_model.get_lam_diag model) in
      let b_2_kmn = calc_b_2_k model ~k:(Common_model.get_kmn model) in
      let n = Mat.dim2 b_2_kmn in
      for i = 1 to n do
        (* TODO: optimize ssqr and col *)
        variances.{i} <- variances.{i} +. Vec.ssqr (Mat.col b_2_kmn i)
      done;
      make ~inputs:(Common_model.get_inputs model) ~variances ~model

    let calc_induceds model induceds =
      let kt_diag, km_2_kmt, kmt = nystrom_2_marginals induceds in
      let variances = copy kt_diag in
      let b_2_kmt = calc_b_2_k model ~k:kmt in
      let n = Mat.dim2 b_2_kmt in
      for i = 1 to n do
        let explained_variance =
          (* TODO: optimize ssqr and col *)
          Vec.ssqr (Mat.col km_2_kmt i) -. Vec.ssqr (Mat.col b_2_kmt i)
        in
        variances.{i} <- variances.{i} -. explained_variance
      done;
      make ~inputs:induceds.Induceds.inputs ~variances ~model

    let copy ?predictive { variances = variances; sigma2 = sigma2 } =
      match predictive with
      | None | Some true ->
          let predictive_variances = Vec.make (Vec.dim variances) sigma2 in
          axpy ~x:variances predictive_variances;
          predictive_variances
      | Some false -> copy variances
  end

  module Common_covariances = struct
    type t = { inputs : inputs; covariances : mat; sigma2 : float }

    let make ~inputs ~covariances ~model =
      let sigma2 = Common_model.get_sigma2 model in
      { inputs = inputs; covariances = covariances; sigma2 = sigma2 }

    let make_b_only ~inputs ~b_2_k ~model =
      make ~inputs ~covariances:(syrk ~trans:`T b_2_k) ~model

    let calc_inducing model =
      let inputs = Common_model.get_inducing_inputs model in
      let b_2_k = calc_b_2_k model ~k:(Common_model.get_km model) in
      make_b_only ~inputs ~b_2_k ~model

    let variances { inputs = inputs; covariances = covs; sigma2 = sigma2 } =
      { Variances.inputs = inputs; variances = Mat.diag covs; sigma2 = sigma2 }

    let copy ?predictive { covariances = covariances; sigma2 = sigma2 } =
      (* TODO: copy upper triangle only *)
      let res = Mat.copy covariances in
      match predictive with
      | None | Some true ->
          for i = 1 to Mat.dim1 res do res.{i, i} <- res.{i, i} +. sigma2 done;
          res
      | Some false -> res
  end

  module FITC_covariances = struct
    include Common_covariances

    let calc_common ~kn_diag ~kmn ~km_2_kmn ~inputs ~model =
      let kernel = Common_model.get_kernel model in
      let covariances = Kernel.upper_no_diag kernel inputs in
      for i = 1 to Vec.dim kn_diag do
        covariances.{i, i} <- kn_diag.{i}
      done;
      ignore (syrk ~trans:`T ~alpha:(-1.) km_2_kmn ~beta:1. ~c:covariances);
      let b_2_kmn = calc_b_2_k model ~k:kmn in
      ignore (syrk ~trans:`T ~alpha:1. b_2_kmn ~beta:1. ~c:covariances);
      make ~inputs ~covariances ~model

    let calc_inputs model =
      let kn_diag = model.Common_model.kn_diag in
      let kmn = model.Common_model.induceds.Induceds.kmn in
      let km_2_kmn = model.Common_model.km_2_kmn in
      let inputs = Common_model.get_inputs model in
      calc_common ~kn_diag ~kmn ~km_2_kmn ~inputs ~model

    let calc_induceds model induceds =
      let inputs = induceds.Induceds.inputs in
      let kn_diag, km_2_kmn, kmn = nystrom_2_marginals induceds in
      calc_common ~kn_diag ~kmn ~km_2_kmn ~inputs ~model
  end

  module FIC_covariances = struct
    include Common_covariances

    let calc_inputs model =
      let inputs = Common_model.get_inputs model in
      let lam_diag = model.Common_model.lam_diag in
      let kmn = model.Common_model.induceds.Induceds.kmn in
      let b_2_kmn = calc_b_2_k model ~k:kmn in
      let covariances = syrk ~trans:`T ~alpha:1. b_2_kmn in
      for i = 1 to Vec.dim lam_diag do
        covariances.{i, i} <- lam_diag.{i} +. covariances.{i, i}
      done;
      make ~inputs ~covariances ~model

    let calc_induceds model { Induceds.inputs = inputs; kmn = kmt } =
      make_b_only ~inputs ~b_2_k:(calc_b_2_k model ~k:kmt) ~model
  end

  module Common_sampler = struct
    type t = { mean : float; stddev : float }

    let calc ~loc ?predictive mean variance =
      if mean.Mean.input != variance.Variance.input then
        failwith (loc ^ ".Sampler: mean and variance disagree about input");
      let used_variance =
        match predictive with
        | None | Some true ->
            variance.Variance.variance +. variance.Variance.sigma2
        | Some false -> variance.Variance.variance
      in
      { mean = mean.Mean.mean; stddev = sqrt used_variance }

    let sample ?(rng = default_rng) sampler =
      let noise = Gsl_randist.gaussian rng ~sigma:sampler.stddev in
      sampler.mean +. noise
  end

  module Common_samplers = struct
    type t = { means : vec; cov_chol : mat }

    let calc ~loc ?predictive means covariances =
      let module Covariances = Common_covariances in
      if means.Means.inputs != covariances.Covariances.inputs then
        failwith
          (loc ^ ".Samplers: means and covariances disagree about inputs");
      (* TODO: copy upper triangle only *)
      let cov_chol = Mat.copy covariances.Covariances.covariances in
      begin
        match predictive with
        | None | Some true ->
            let sigma2 = covariances.Covariances.sigma2 in
            for i = 1 to Mat.dim1 cov_chol do
              cov_chol.{i, i} <- cov_chol.{i, i} +. sigma2
            done
        | Some false -> ()
      end;
      potrf ~jitter cov_chol;
      { means = means.Means.means; cov_chol = cov_chol }

    let sample ?(rng = default_rng) samplers =
      let n = Vec.dim samplers.means in
      let noises = Vec.init n (fun _ -> Gsl_randist.gaussian rng ~sigma:1.) in
      trmv ~trans:`T samplers.cov_chol noises;
      axpy ~x:samplers.means noises;
      noises
  end
end

module Make_traditional (Spec : Spec) = struct
  include Make_common (Spec)
  module Model = Common_model
end

module Make_variational (Spec : Spec) = struct
  include Make_common (Spec)
  module Model = Variational_model
end

module FIC_Loc = struct let loc = "FIC" end
module Variational_FITC_Loc = struct let loc = "Variational_FITC" end
module Variational_FIC_Loc = struct let loc = "Variational_FIC" end

let fitc_loc = "FITC"
let fic_loc = "FIC"
let variational_fitc_loc = "Variational_FITC"
let variational_fic_loc = "Variational_FIC"

module Make_FITC (Spec : Spec) = struct
  include Make_traditional (Spec)
  module Covariances = FITC_covariances

  module Sampler = struct
    include Common_sampler
    let calc = calc ~loc:fitc_loc
  end

  module Samplers = struct
    include Common_samplers
    let calc = calc ~loc:fitc_loc
  end
end

module Make_FIC (Spec : Spec) = struct
  include Make_traditional (Spec)
  module Covariances = FIC_covariances

  module Sampler = struct
    include Common_sampler
    let calc = calc ~loc:fic_loc
  end

  module Samplers = struct
    include Common_samplers
    let calc = calc ~loc:fic_loc
  end
end

module Make_variational_FITC (Spec : Spec) = struct
  include Make_variational (Spec)
  module Covariances = FITC_covariances

  module Sampler = struct
    include Common_sampler
    let calc = calc ~loc:variational_fitc_loc
  end

  module Samplers = struct
    include Common_samplers
    let calc = calc ~loc:variational_fitc_loc
  end
end

module Make_variational_FIC (Spec : Spec) = struct
  include Make_variational (Spec)
  module Covariances = FIC_covariances

  module Sampler = struct
    include Common_sampler
    let calc = calc ~loc:variational_fic_loc
  end

  module Samplers = struct
    include Common_samplers
    let calc = calc ~loc:variational_fic_loc
  end
end

module Make (Spec : Spec) = struct
  module type Sig = Inducing_input_gpr with module Kernel = Spec.Kernel

  module Common = Make_common (Spec)

  module FITC = struct
    include Common
    module Model = Common_model
    module Covariances = FITC_covariances

    module Sampler = struct
      include Common_sampler
      let calc = calc ~loc:fitc_loc
    end

    module Samplers = struct
      include Common_samplers
      let calc = calc ~loc:fitc_loc
    end
  end

  module FIC = struct
    include Common
    module Model = Common_model
    module Covariances = FIC_covariances

    module Sampler = struct
      include Common_sampler
      let calc = calc ~loc:fic_loc
    end

    module Samplers = struct
      include Common_samplers
      let calc = calc ~loc:fic_loc
    end
  end

  module Variational_FITC = struct
    include Common
    module Model = Variational_model
    module Covariances = FITC_covariances

    module Sampler = struct
      include Common_sampler
      let calc = calc ~loc:variational_fitc_loc
    end

    module Samplers = struct
      include Common_samplers
      let calc = calc ~loc:variational_fitc_loc
    end
  end

  module Variational_FIC = struct
    include Common
    module Model = Variational_model
    module Covariances = FIC_covariances

    module Sampler = struct
      include Common_sampler
      let calc = calc ~loc:variational_fic_loc
    end

    module Samplers = struct
      include Common_samplers
      let calc = calc ~loc:variational_fic_loc
    end
  end
end
