open Lacaml.Impl.D

open Interfaces
open Utils

(* TODO: dimension sanity checks *)
(* TODO: consistency checks; also finite differences *)

module type Spec = sig
  module Eval_spec : Inducing_input_gpr.Specs.Eval

  val get_sigma2 : Eval_spec.kernel -> float
  val jitter : float
end

module type Sig = functor (FITC_spec : Spec) ->
  Inducing_input_gpr.Sigs.Eval with module Spec = FITC_spec.Eval_spec

module Make_common (FITC_spec : Spec) = struct
  open FITC_spec

  module Spec = FITC_spec.Eval_spec

  open Spec

  module Inducing = struct
    type t = {
      kernel : kernel;
      points : Eval_spec.Inducing.t;
      km : mat;
      km_chol : mat;
      log_det_km : float;
    }

    let calc kernel points =
      let km = Eval_spec.Inducing.upper kernel points in
      (* TODO: copy upper triangle only *)
      let km_chol = Mat.copy km in
      potrf ~jitter km_chol;
      let log_det_km = log_det km_chol in
      {
        kernel = kernel;
        points = points;
        km = km;
        km_chol = km_chol;
        log_det_km = log_det_km;
      }
  end

  module Input = struct
    type t = {
      inducing : Inducing.t;
      point : Eval_spec.Input.t;
      k_m : vec;
    }

    let calc inducing point =
      let kernel = inducing.Inducing.kernel in
      {
        inducing = inducing;
        point = point;
        k_m =
          Eval_spec.Input.eval
            kernel ~inducing:inducing.Inducing.points ~input:point;
      }

    let get_kernel t = t.inducing.Inducing.kernel
  end

  let calc_basis_2_k basis_chol ~k =
    (* TODO: consider symmetric matrices *)
    let b_2_k = Mat.copy k in
    trtrs ~trans:`T basis_chol b_2_k;
    b_2_k

  module Inputs = struct
    type t = {
      inducing : Inducing.t;
      points : Inputs.t;
      kmn : mat;
    }

    let calc inducing points =
      let kernel = inducing.Inducing.kernel in
      let kmn =
        Inputs.cross kernel ~inducing:inducing.Inducing.points ~inputs:points
      in
      {
        inducing = inducing;
        points = points;
        kmn = kmn;
      }

    let get_kernel t = t.inducing.Inducing.kernel
    let get_inducing_points t = t.inducing.Inducing.points

    (* Compute square root of Nystrom approximation, and diagonal of
       marginal variances *)
    let nystrom_2_marginals inputs =
      let inducing = inputs.inducing in
      let kernel = get_kernel inputs in
      let kn_diag = Inputs.diag kernel inputs.points in
      let kmn = inputs.kmn in
      let km_2_kmn = calc_basis_2_k inducing.Inducing.km_chol ~k:kmn in
      km_2_kmn, kmn, kn_diag
  end

  module Common_model = struct
    type t = {
      inputs : Inputs.t;
      kn_diag : vec;
      km_2_kmn : mat;
      lam_diag : vec;
      inv_lam_sigma2_diag : vec;
      b_chol : mat;
      neg_log_marginal_likelihood : float;
    }

    let calc inputs =
      let inducing = inputs.Inputs.inducing in
      let kernel = Inputs.get_kernel inputs in
      let sigma2 = get_sigma2 kernel in
      let km_2_kmn, kmn, kn_diag = Inputs.nystrom_2_marginals inputs in
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
      {
        inputs = inputs;
        kn_diag = kn_diag;
        km_2_kmn = km_2_kmn;
        lam_diag = lam_diag;
        inv_lam_sigma2_diag = inv_lam_sigma2_diag;
        b_chol = b_chol;
        neg_log_marginal_likelihood = (l1_2 +. float n_inputs *. log_2pi) /. 2.;
      }

    let neg_log_marginal_likelihood model = model.neg_log_marginal_likelihood
    let get_inducing model = model.inputs.Inputs.inducing
    let get_inducing_points model = (get_inducing model).Inducing.points
    let get_input_points model = model.inputs.Inputs.points
    let get_kernel model = (get_inducing model).Inducing.kernel
    let get_sigma2 model = get_sigma2 (get_kernel model)
    let get_kmn model = model.inputs.Inputs.kmn
    let get_lam_diag model = model.lam_diag
    let get_km model = (get_inducing model).Inducing.km
  end

  module Variational_model = struct
    include Common_model

    let calc inputs =
      let
        {
          lam_diag = lam_diag;
          inv_lam_sigma2_diag = inv_lam_sigma2_diag;
          neg_log_marginal_likelihood = neg_log_marginal_likelihood;
        } as model = calc inputs
      in
      let n = Vec.dim lam_diag in
      let trace_term = ref 0. in
      for i = 1 to n do
        trace_term := !trace_term +. lam_diag.{i} *. inv_lam_sigma2_diag.{i}
      done;
      {
        model with
        neg_log_marginal_likelihood =
          neg_log_marginal_likelihood +. 0.5 *. !trace_term;
      }
  end

  module Trained = struct
    type t = {
      model : Common_model.t;
      inv_b_chol_kmn_y__ : vec;
      neg_log_marginal_likelihood : float;
    }

    let calc model ~targets =
      let n_targets = Vec.dim targets in
      let y__ = Vec.create n_targets in
      let inv_lam_sigma2_diag = model.Common_model.inv_lam_sigma2_diag in
      for i = 1 to n_targets do
        y__.{i} <- targets.{i} *. inv_lam_sigma2_diag.{i}
      done;
      let kmn = Common_model.get_kmn model in
      let inv_b_chol_kmn_y__ = gemv kmn y__ in
      let ssqr_y__ = dot ~x:targets y__ in
      let b_chol = model.Common_model.b_chol in
      trsv ~trans:`T b_chol inv_b_chol_kmn_y__;
      let fit_neg_log_marginal_likelihood =
        (ssqr_y__ -. Vec.ssqr inv_b_chol_kmn_y__) /. 2.
      in
      {
        model = model;
        inv_b_chol_kmn_y__ = inv_b_chol_kmn_y__;
        neg_log_marginal_likelihood =
          model.Common_model.neg_log_marginal_likelihood +.
          fit_neg_log_marginal_likelihood;
      }

    let neg_log_marginal_likelihood trained =
      trained.neg_log_marginal_likelihood
  end

  module Weights = struct
    type t = { model : Common_model.t; coeffs : vec }

    let get_kernel weights = Common_model.get_kernel weights.model
    let get_inducing weights = Common_model.get_inducing weights.model

    let get_inducing_points weights =
      Common_model.get_inducing_points weights.model

    let get_coeffs weights = weights.coeffs

    let calc trained =
      let coeffs = copy trained.Trained.inv_b_chol_kmn_y__ in
      trsv trained.Trained.model.Common_model.b_chol coeffs;
      {
        model = trained.Trained.model;
        coeffs = coeffs;
      }
  end

  module Mean = struct
    type t = { point : Eval_spec.Input.t; value : float }

    let make ~point ~value = { point = point; value = value }

    let calc_input weights point =
      let inducing_points = Weights.get_inducing_points weights in
      let kernel = Weights.get_kernel weights in
      let coeffs = Weights.get_coeffs weights in
      let value =
        Eval_spec.Input.weighted_eval
          kernel ~coeffs ~inducing:inducing_points ~input:point
      in
      make ~point ~value

    let calc_induced weights input =
      if Weights.get_inducing weights <> input.Input.inducing then
        failwith
          "Fitc.Make_common.Mean.calc_induced: \
          weights and input disagree about inducing points";
      let value = dot ~x:input.Input.k_m weights.Weights.coeffs in
      make ~point:input.Input.point ~value

    let get m = m.value
  end

  module Means = struct
    type t = { points : Eval_spec.Inputs.t; values : vec }

    let make ~points ~values = { points = points; values = values }

    let calc_model_inputs { Weights.coeffs = coeffs; model = model } =
      make
        ~points:(Common_model.get_input_points model)
        ~values:(gemv ~trans:`T (Common_model.get_kmn model) coeffs)

    let calc_inputs weights points =
      let kernel = Weights.get_kernel weights in
      let coeffs = Weights.get_coeffs weights in
      let inducing_points = Weights.get_inducing_points weights in
      let values =
        Eval_spec.Inputs.weighted_eval
          kernel ~coeffs ~inducing:inducing_points ~inputs:points
      in
      make ~points ~values

    let calc_induced weights inputs =
      let { Inputs.points = points; kmn = kmn } = inputs in
      if Weights.get_inducing weights <> inputs.Inputs.inducing then
        failwith
          "Fitc.Make_common.Means.calc_induced: \
          weights and inputs disagree about inducing points";
      make ~points ~values:(gemv ~trans:`T kmn weights.Weights.coeffs)

    let get means = means.values

    module Inducing = struct
      type t = { points : Eval_spec.Inducing.t; values : vec }

      let make ~points ~values = { points = points; values = values }

      let calc { Weights.coeffs = coeffs; model = model } =
        make
          ~points:(Common_model.get_inducing_points model)
          ~values:(symv (Common_model.get_km model) coeffs)

      let get means = means.values
    end
  end

  module Variance = struct
    type t = { point : Eval_spec.Input.t; variance : float; sigma2 : float }

    let calc_induced model induced =
      let { Input.point = point; k_m = k_m } = induced in
      let kernel = Common_model.get_kernel model in
      let prior_variance = Eval_spec.Input.eval_one kernel point in
      let inv_km_chol_k_m = copy k_m in
      let inv_km_chol_k_m_mat = Mat.from_col_vec inv_km_chol_k_m in
      let inducing = induced.Input.inducing in
      potrs ~factorize:false inducing.Inducing.km_chol inv_km_chol_k_m_mat;
      let km_arg = dot ~x:k_m inv_km_chol_k_m in
      let inv_b_chol_k_m = copy k_m ~y:inv_km_chol_k_m in
      let inv_b_chol_k_m_mat = inv_km_chol_k_m_mat in
      potrs ~factorize:false model.Common_model.b_chol inv_b_chol_k_m_mat;
      let b_arg = dot ~x:k_m inv_b_chol_k_m in
      let explained_variance = km_arg -. b_arg in
      let variance = prior_variance -. explained_variance in
      {
        point = point;
        variance = variance;
        sigma2 = Common_model.get_sigma2 model;
      }

    let get ?predictive t =
      match predictive with
      | None | Some true -> t.variance +. t.sigma2
      | Some false -> t.variance
  end

  let calc_b_2_k model ~k = calc_basis_2_k model.Common_model.b_chol ~k

  module Variances = struct
    type t = { points : Eval_spec.Inputs.t; variances : vec; sigma2 : float }

    let make ~points ~variances ~model =
      let sigma2 = Common_model.get_sigma2 model in
      { points = points; variances = variances; sigma2 = sigma2 }

    let calc_model_inputs model =
      let variances = copy (Common_model.get_lam_diag model) in
      let b_2_kmn = calc_b_2_k model ~k:(Common_model.get_kmn model) in
      let n = Mat.dim2 b_2_kmn in
      for i = 1 to n do
        (* TODO: optimize ssqr and col *)
        variances.{i} <- variances.{i} +. Vec.ssqr (Mat.col b_2_kmn i)
      done;
      make ~points:(Common_model.get_input_points model) ~variances ~model

    let calc_induced model inputs =
      if Common_model.get_inducing model <> inputs.Inputs.inducing then
        failwith
          "Fitc.Make_common.Variances.calc_induced: \
          model and inputs disagree about inducing points";
      let km_2_kmt, kmt, kt_diag = Inputs.nystrom_2_marginals inputs in
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
      make ~points:inputs.Inputs.points ~variances ~model

    let get_common ?predictive ~variances ~sigma2 =
      match predictive with
      | None | Some true ->
          let predictive_variances = Vec.make (Vec.dim variances) sigma2 in
          axpy ~x:variances predictive_variances;
          predictive_variances
      | Some false -> variances

    let get ?predictive { variances = variances; sigma2 = sigma2 } =
      get_common ?predictive ~variances ~sigma2

    module Inducing = struct
      type t = {
        points : Eval_spec.Inducing.t;
        variances : vec;
        sigma2 : float;
      }

      let make ~points ~variances ~model =
        let sigma2 = Common_model.get_sigma2 model in
        { points = points; variances = variances; sigma2 = sigma2 }

      let calc model =
        let b_2_km = calc_b_2_k model ~k:(Common_model.get_km model) in
        let m = Mat.dim2 b_2_km in
        let variances = Vec.create m in
        for i = 1 to m do
          (* TODO: optimize ssqr and col *)
          variances.{i} <- Vec.ssqr (Mat.col b_2_km i)
        done;
        make ~points:(Common_model.get_inducing_points model) ~variances ~model

      let get ?predictive { variances = variances; sigma2 = sigma2 } =
        get_common ?predictive ~variances ~sigma2
    end
  end

  module Common_covariances = struct
    type t = { points : Eval_spec.Inputs.t; covariances : mat; sigma2 : float }

    let make ~points ~covariances ~model =
      let sigma2 = Common_model.get_sigma2 model in
      { points = points; covariances = covariances; sigma2 = sigma2 }

    let make_b_only ~points ~b_2_k ~model =
      make ~points ~covariances:(syrk ~trans:`T b_2_k) ~model

    let get_common ?predictive ~covariances ~sigma2 =
      match predictive with
      | None | Some true ->
          (* TODO: copy upper triangle only *)
          let res = Mat.copy covariances in
          for i = 1 to Mat.dim1 res do res.{i, i} <- res.{i, i} +. sigma2 done;
          res
      | Some false -> covariances

    let get ?predictive { covariances = covariances; sigma2 = sigma2 } =
      get_common ?predictive ~covariances ~sigma2

    let variances { points = points; covariances = covs; sigma2 = sigma2 } =
      { Variances.points = points; variances = Mat.diag covs; sigma2 = sigma2 }

    module Inducing = struct
      type t = {
        points : Eval_spec.Inducing.t;
        covariances : mat;
        sigma2 : float;
      }

      let calc model =
        let points = Common_model.get_inducing_points model in
        let b_2_k = calc_b_2_k model ~k:(Common_model.get_km model) in
        let covariances = syrk ~trans:`T b_2_k in
        let sigma2 = Common_model.get_sigma2 model in
        { points = points; covariances = covariances; sigma2 = sigma2 }

      let get ?predictive { covariances = covariances; sigma2 = sigma2 } =
        get_common ?predictive ~covariances ~sigma2

      let variances { points = points; covariances = covs; sigma2 = sigma2 } =
        {
          Variances.Inducing.
          points = points;
          variances = Mat.diag covs;
          sigma2 = sigma2
        }
    end
  end

  module FITC_covariances = struct
    include Common_covariances

    let calc_common ~kn_diag ~kmn ~km_2_kmn ~points ~model =
      let kernel = Common_model.get_kernel model in
      let covariances = Eval_spec.Inputs.upper_no_diag kernel points in
      for i = 1 to Vec.dim kn_diag do
        covariances.{i, i} <- kn_diag.{i}
      done;
      ignore (syrk ~trans:`T ~alpha:(-1.) km_2_kmn ~beta:1. ~c:covariances);
      let b_2_kmn = calc_b_2_k model ~k:kmn in
      ignore (syrk ~trans:`T ~alpha:1. b_2_kmn ~beta:1. ~c:covariances);
      make ~points ~covariances ~model

    let calc_model_inputs model =
      let kn_diag = model.Common_model.kn_diag in
      let kmn = model.Common_model.inputs.Inputs.kmn in
      let km_2_kmn = model.Common_model.km_2_kmn in
      let points = Common_model.get_input_points model in
      calc_common ~kn_diag ~kmn ~km_2_kmn ~points ~model

    let calc_induced model inputs =
      if Common_model.get_inducing model <> inputs.Inputs.inducing then
        failwith (
          "Make_common.FITC_covariances.calc_induced: \
          model and inputs disagree about inducing points");
      let km_2_kmn, kmn, kn_diag = Inputs.nystrom_2_marginals inputs in
      let points = inputs.Inputs.points in
      calc_common ~kn_diag ~kmn ~km_2_kmn ~points ~model
  end

  module FIC_covariances = struct
    include Common_covariances

    let calc_model_inputs model =
      let points = Common_model.get_input_points model in
      let lam_diag = model.Common_model.lam_diag in
      let kmn = model.Common_model.inputs.Inputs.kmn in
      let b_2_kmn = calc_b_2_k model ~k:kmn in
      let covariances = syrk ~trans:`T ~alpha:1. b_2_kmn in
      for i = 1 to Vec.dim lam_diag do
        covariances.{i, i} <- lam_diag.{i} +. covariances.{i, i}
      done;
      make ~points ~covariances ~model

    let calc_induced model inputs =
      if Common_model.get_inducing model <> inputs.Inputs.inducing then
        failwith (
          "Make_common.FIC_covariances.calc_induced: \
          model and inputs disagree about inducing points");
      let kmt = inputs.Inputs.kmn in
      let points = inputs.Inputs.points in
      make_b_only ~points ~b_2_k:(calc_b_2_k model ~k:kmt) ~model
  end

  module Common_sampler = struct
    type t = { mean : float; stddev : float }

    let calc ~loc ?predictive mean variance =
      if mean.Mean.point <> variance.Variance.point then
        failwith (
          loc ^ ".Sampler: mean and variance disagree about input point");
      let used_variance =
        match predictive with
        | None | Some true ->
            variance.Variance.variance +. variance.Variance.sigma2
        | Some false -> variance.Variance.variance
      in
      { mean = mean.Mean.value; stddev = sqrt used_variance }

    let sample ?(rng = default_rng) sampler =
      let noise = Gsl_randist.gaussian rng ~sigma:sampler.stddev in
      sampler.mean +. noise

    let samples ?(rng = default_rng) sampler ~n =
      Vec.init n (fun _ -> sample ~rng sampler)
  end

  module Common_cov_sampler = struct
    type t = { means : vec; cov_chol : mat }

    let calc ~loc ?predictive means covariances =
      let module Covariances = Common_covariances in
      if means.Means.points <> covariances.Covariances.points then
        failwith (
          loc ^
          ".Cov_sampler: means and covariances disagree about input points");
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
      { means = means.Means.values; cov_chol = cov_chol }

    let sample ?(rng = default_rng) samplers =
      let n = Vec.dim samplers.means in
      let sample = Vec.init n (fun _ -> Gsl_randist.gaussian rng ~sigma:1.) in
      trmv ~trans:`T samplers.cov_chol sample;
      axpy ~x:samplers.means sample;
      sample

    let samples ?(rng = default_rng) { means = means; cov_chol = cov_chol } ~n =
      let n_means = Vec.dim means in
      let samples =
        Mat.init_cols n_means n (fun _ _ -> Gsl_randist.gaussian rng ~sigma:1.)
      in
      trmm ~trans:`T cov_chol ~b:samples;
      for col = 1 to n do
        for row = 1 to n_means do
          let mean = means.{row} in
          samples.{row, col} <- samples.{row, col} +. mean
        done
      done;
      samples
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

  module Cov_sampler = struct
    include Common_cov_sampler
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

  module Cov_sampler = struct
    include Common_cov_sampler
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

  module Cov_sampler = struct
    include Common_cov_sampler
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

  module Cov_sampler = struct
    include Common_cov_sampler
    let calc = calc ~loc:variational_fic_loc
  end
end

module Make (Spec : Spec) = struct
  module type Sig =
    Inducing_input_gpr.Sigs.Eval with module Spec = Spec.Eval_spec

  module Common = Make_common (Spec)

  module FITC = struct
    include Common
    module Model = Common_model
    module Covariances = FITC_covariances

    module Sampler = struct
      include Common_sampler
      let calc = calc ~loc:fitc_loc
    end

    module Cov_sampler = struct
      include Common_cov_sampler
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

    module Cov_sampler = struct
      include Common_cov_sampler
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

    module Cov_sampler = struct
      include Common_cov_sampler
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

    module Cov_sampler = struct
      include Common_cov_sampler
      let calc = calc ~loc:variational_fic_loc
    end
  end
end


(* Derivable *)

(*
module type Deriv_spec = sig
  module Eval : Spec
  module Deriv_kernel : Deriv_kernel
  module Inducing_deriv :
end

module Make_common_deriv (Spec : Deriv_spec) = struct
end
*)
