open Printf

open Lacaml.Impl.D

open Utils
open Interfaces
open Inducing_input_gpr

(* TODO: dimension sanity checks *)
(* TODO: consistency checks; also finite differences *)

module type Sig = functor (Spec : Specs.Eval) ->
  Sigs.Eval with module Spec = Spec

module Make_common (Spec : Specs.Eval) = struct
  module Spec = Spec

  open Spec

  let jitter = !cholesky_jitter

  module Inducing = struct
    type t = {
      kernel : Kernel.t;
      points : Spec.Inducing.t;
      km : mat;
      km_chol : mat;
      log_det_km : float;
    }

    let calc_internal kernel points km =
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

    let calc kernel points =
      let km = Spec.Inducing.upper kernel points in
      calc_internal kernel points km

    let get_kernel inducing = inducing.kernel
  end

  module Input = struct
    type t = {
      inducing : Inducing.t;
      point : Spec.Input.t;
      k_m : vec;
    }

    let calc inducing point =
      let kernel = inducing.Inducing.kernel in
      {
        inducing = inducing;
        point = point;
        k_m =
          Spec.Input.eval
            kernel ~inducing:inducing.Inducing.points ~input:point;
      }

    let get_kernel t = t.inducing.Inducing.kernel
  end

  module Inputs = struct
    type t = {
      inducing : Inducing.t;
      points : Inputs.t;
      kmn : mat;
    }

    let calc_internal inducing points kmn =
      {
        inducing = inducing;
        points = points;
        kmn = kmn;
      }

    let calc inducing points =
      let kernel = inducing.Inducing.kernel in
      let kmn =
        Inputs.cross kernel ~inducing:inducing.Inducing.points ~inputs:points
      in
      calc_internal inducing points kmn

    let get_kernel t = t.inducing.Inducing.kernel
    let get_sigma2 t = Kernel.get_sigma2 t.inducing.Inducing.kernel
    let get_inducing_points t = t.inducing.Inducing.points

    (* TODO: move *)
    let nystrom_chol_marginals inputs =
      let inducing = inputs.inducing in
      let kernel = get_kernel inputs in
      let kn_diag = Inputs.diag kernel inputs.points in
      let inv_km_chol_kmn =
        solve_triangular ~trans:`T inducing.Inducing.km_chol ~k:inputs.kmn
      in
      inv_km_chol_kmn, kn_diag
  end

  module Common_model = struct
    type t = {
      inputs : Inputs.t;
      kn_diag : vec;
      inv_km_chol_kmn : mat;
      lam_diag : vec;
      inv_lam_sigma2_diag : vec;
      b_chol : mat;
      evidence : float;
    }

    let calc_internal inputs kn_diag =
      let inducing = inputs.Inputs.inducing in
      let sigma2 = Inputs.get_sigma2 inputs in
      let kmn = inputs.Inputs.kmn in
      let inv_km_chol_kmn =
        solve_triangular ~trans:`T inducing.Inducing.km_chol ~k:kmn
      in
      let kmn_= Mat.copy kmn in
      let n_inputs = Vec.dim kn_diag in
      let lam_diag = Vec.create n_inputs in
      let inv_lam_sigma2_diag = Vec.create n_inputs in
      let rec loop log_det_lam_sigma2 i =
        if i = 0 then log_det_lam_sigma2
        else
          let kn_diag_i = kn_diag.{i} in
          (* TODO: optimize ssqr and col *)
          let qn_diag_i = Vec.ssqr (Mat.col inv_km_chol_kmn i) in
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
      let l1_2 = log_det_km -. log_det_b -. log_det_lam_sigma2 in
      {
        inputs = inputs;
        kn_diag = kn_diag;
        inv_km_chol_kmn = inv_km_chol_kmn;
        lam_diag = lam_diag;
        inv_lam_sigma2_diag = inv_lam_sigma2_diag;
        b_chol = b_chol;
        evidence = 0.5 *. (l1_2 -. float n_inputs *. log_2pi);
      }

    let calc inputs =
      let kernel = Inputs.get_kernel inputs in
      let kn_diag = Spec.Inputs.diag kernel inputs.Inputs.points in
      calc_internal inputs kn_diag

    let calc_evidence model = model.evidence
    let get_inducing model = model.inputs.Inputs.inducing
    let get_inducing_points model = (get_inducing model).Inducing.points
    let get_input_points model = model.inputs.Inputs.points
    let get_kernel model = (get_inducing model).Inducing.kernel
    let get_sigma2 model = Kernel.get_sigma2 (get_kernel model)
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
          inv_lam_sigma2_diag = x;
          evidence = evidence;
        } as model = calc inputs
      in
      { model with evidence = evidence +. 0.5 *. dot ~x lam_diag }
  end

  module Trained = struct
    type t = {
      model : Common_model.t;
      inv_b_chol_kmn_y__ : vec;
      evidence : float;
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
      let fit_evidence = 0.5 *. (Vec.ssqr inv_b_chol_kmn_y__ -. ssqr_y__) in
      {
        model = model;
        inv_b_chol_kmn_y__ = inv_b_chol_kmn_y__;
        evidence = model.Common_model.evidence +. fit_evidence;
      }

    let calc_evidence trained = trained.evidence
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
    type t = { point : Spec.Input.t; value : float }

    let make ~point ~value = { point = point; value = value }

    let calc_input weights point =
      let inducing_points = Weights.get_inducing_points weights in
      let kernel = Weights.get_kernel weights in
      let coeffs = Weights.get_coeffs weights in
      let value =
        Spec.Input.weighted_eval
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
    type t = { points : Spec.Inputs.t; values : vec }

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
        Spec.Inputs.weighted_eval
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
      type t = { points : Spec.Inducing.t; values : vec }

      let make ~points ~values = { points = points; values = values }

      let calc { Weights.coeffs = coeffs; model = model } =
        make
          ~points:(Common_model.get_inducing_points model)
          ~values:(symv (Common_model.get_km model) coeffs)

      let get means = means.values
    end
  end

  module Variance = struct
    type t = { point : Spec.Input.t; variance : float; sigma2 : float }

    let calc_induced model induced =
      let { Input.point = point; k_m = k_m } = induced in
      let kernel = Common_model.get_kernel model in
      let prior_variance = Spec.Input.eval_one kernel point in
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

  let solve_b_chol model ~k =
    solve_triangular ~trans:`T model.Common_model.b_chol ~k

  module Variances = struct
    type t = { points : Spec.Inputs.t; variances : vec; sigma2 : float }

    let make ~points ~variances ~model =
      let sigma2 = Common_model.get_sigma2 model in
      { points = points; variances = variances; sigma2 = sigma2 }

    let calc_model_inputs model =
      let variances = copy (Common_model.get_lam_diag model) in
      let inv_b_chol_kmn = solve_b_chol model ~k:(Common_model.get_kmn model) in
      let n = Mat.dim2 inv_b_chol_kmn in
      for i = 1 to n do
        (* TODO: optimize ssqr and col *)
        variances.{i} <- variances.{i} +. Vec.ssqr (Mat.col inv_b_chol_kmn i)
      done;
      make ~points:(Common_model.get_input_points model) ~variances ~model

    let calc_induced model inputs =
      if Common_model.get_inducing model <> inputs.Inputs.inducing then
        failwith
          "Fitc.Make_common.Variances.calc_induced: \
          model and inputs disagree about inducing points";
      let kmt = inputs.Inputs.kmn in
      let inv_km_chol_kmt, kt_diag = Inputs.nystrom_chol_marginals inputs in
      let variances = copy kt_diag in
      let inv_b_chol_kmt = solve_b_chol model ~k:kmt in
      let n = Mat.dim2 inv_b_chol_kmt in
      for i = 1 to n do
        let explained_variance =
          (* TODO: optimize ssqr and col *)
          Vec.ssqr (Mat.col inv_km_chol_kmt i) -.
            Vec.ssqr (Mat.col inv_b_chol_kmt i)
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
        points : Spec.Inducing.t;
        variances : vec;
        sigma2 : float;
      }

      let make ~points ~variances ~model =
        let sigma2 = Common_model.get_sigma2 model in
        { points = points; variances = variances; sigma2 = sigma2 }

      let calc model =
        let inv_b_chol_km = solve_b_chol model ~k:(Common_model.get_km model) in
        let m = Mat.dim2 inv_b_chol_km in
        let variances = Vec.create m in
        for i = 1 to m do
          (* TODO: optimize ssqr and col *)
          variances.{i} <- Vec.ssqr (Mat.col inv_b_chol_km i)
        done;
        make ~points:(Common_model.get_inducing_points model) ~variances ~model

      let get ?predictive { variances = variances; sigma2 = sigma2 } =
        get_common ?predictive ~variances ~sigma2
    end
  end

  module Common_covariances = struct
    type t = { points : Spec.Inputs.t; covariances : mat; sigma2 : float }

    let make ~points ~covariances ~model =
      let sigma2 = Common_model.get_sigma2 model in
      { points = points; covariances = covariances; sigma2 = sigma2 }

    let make_b_only ~points ~inv_b_chol_k ~model =
      make ~points ~covariances:(syrk ~trans:`T inv_b_chol_k) ~model

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
        points : Spec.Inducing.t;
        covariances : mat;
        sigma2 : float;
      }

      let calc model =
        let points = Common_model.get_inducing_points model in
        let inv_b_chol_k = solve_b_chol model ~k:(Common_model.get_km model) in
        let covariances = syrk ~trans:`T inv_b_chol_k in
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

    let calc_common ~kn_diag ~kmn ~inv_km_chol_kmn ~points ~model =
      let kernel = Common_model.get_kernel model in
      let covariances = Spec.Inputs.upper_no_diag kernel points in
      for i = 1 to Vec.dim kn_diag do
        covariances.{i, i} <- kn_diag.{i}
      done;
      ignore (syrk ~trans:`T ~alpha:(-1.) inv_km_chol_kmn ~beta:1. ~c:covariances);
      let inv_b_chol_kmn = solve_b_chol model ~k:kmn in
      ignore (syrk ~trans:`T ~alpha:1. inv_b_chol_kmn ~beta:1. ~c:covariances);
      make ~points ~covariances ~model

    let calc_model_inputs model =
      let kn_diag = model.Common_model.kn_diag in
      let kmn = model.Common_model.inputs.Inputs.kmn in
      let inv_km_chol_kmn = model.Common_model.inv_km_chol_kmn in
      let points = Common_model.get_input_points model in
      calc_common ~kn_diag ~kmn ~inv_km_chol_kmn ~points ~model

    let calc_induced model inputs =
      if Common_model.get_inducing model <> inputs.Inputs.inducing then
        failwith (
          "Make_common.FITC_covariances.calc_induced: \
          model and inputs disagree about inducing points");
      let kmn = inputs.Inputs.kmn in
      let inv_km_chol_kmn, kn_diag = Inputs.nystrom_chol_marginals inputs in
      let points = inputs.Inputs.points in
      calc_common ~kn_diag ~kmn ~inv_km_chol_kmn ~points ~model
  end

  module FIC_covariances = struct
    include Common_covariances

    let calc_model_inputs model =
      let points = Common_model.get_input_points model in
      let lam_diag = model.Common_model.lam_diag in
      let kmn = model.Common_model.inputs.Inputs.kmn in
      let inv_b_chol_kmn = solve_b_chol model ~k:kmn in
      let covariances = syrk ~trans:`T ~alpha:1. inv_b_chol_kmn in
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
      make_b_only ~points ~inv_b_chol_k:(solve_b_chol model ~k:kmt) ~model
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

module Make_traditional (Spec : Specs.Eval) = struct
  include Make_common (Spec)
  module Model = Common_model
end

module Make_variational (Spec : Specs.Eval) = struct
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

module Make_FITC (Spec : Specs.Eval) = struct
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

module Make_FIC (Spec : Specs.Eval) = struct
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

module Make_variational_FITC (Spec : Specs.Eval) = struct
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

module Make_variational_FIC (Spec : Specs.Eval) = struct
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

module Make (Spec : Specs.Eval) = struct
  module type Sig = Sigs.Eval with module Spec = Spec

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

module type Deriv_sig = functor (Spec : Specs.Eval_deriv) ->
  Sigs.Deriv
    with module Eval.Spec = Spec.Eval_spec
    with module Deriv.Spec = Spec.Deriv_spec

module Make_common_deriv (Spec : Specs.Eval_deriv) = struct
  open Spec

  module Eval_common = Make_common (Eval_spec)

  open Eval_common

  module Deriv_common = struct
    module Spec = Deriv_spec

    open Spec

    module Inducing = struct
      type t = {
        eval_inducing : Eval_common.Inducing.t;
        shared : Spec.Inducing.shared;
      }

      let calc kernel points =
        let km, shared = Spec.Inducing.calc_shared kernel points in
        let eval_inducing =
          Eval_common.Inducing.calc_internal kernel points km
        in
        {
          eval_inducing = eval_inducing;
          shared = shared;
        }

      let calc_eval inducing = inducing.eval_inducing

      let get_kernel inducing =
        Eval_common.Inducing.get_kernel inducing.eval_inducing
    end

    module Inputs = struct
      type t = {
        inducing : Inducing.t;
        eval_inputs : Eval_common.Inputs.t;
        shared_cross : Spec.Inputs.cross;
      }

      let calc inducing points =
        let eval_inducing = inducing.Inducing.eval_inducing in
        let kernel = eval_inducing.Eval_common.Inducing.kernel in
        let kmn, shared_cross =
          Spec.Inputs.calc_shared_cross kernel
            ~inducing:eval_inducing.Eval_common.Inducing.points ~inputs:points
        in
        let eval_inputs =
          Eval_common.Inputs.calc_internal eval_inducing points kmn
        in
        {
          inducing = inducing;
          eval_inputs = eval_inputs;
          shared_cross = shared_cross;
        }

      let calc_eval t = t.eval_inputs

      let get_kernel inputs = Inducing.get_kernel inputs.inducing
    end

    module Common_model = struct
      type hyper =
        [
        | `Hyper of Hyper.t * vec
        | `Sigma2
        ]

      type evidence = {
        hyper : hyper;
        evidence : float;
(*         dkm : mat; *)
      }

      type t = {
        inputs : Inputs.t;
        eval_model : Eval_common.Common_model.t;
        shared_diag : Spec.Inputs.diag;
        calc_evidence : Hyper.t -> evidence;
        calc_evidence_sigma2 : unit -> evidence;
      }

      module Eval_inducing = Eval_common.Inducing
      module Eval_inputs = Eval_common.Inputs
      module Eval_model = Eval_common.Common_model

      (* Checks whether a sparse row matrix is sane *)
      let check_sparse_sane mat = function
        | None -> ()
        | Some rows ->
            let m = Mat.dim1 mat in
            let n_rows = Array.length rows in
            if n_rows <> m then
              failwith (
                sprintf
                  "check_sparse_sane: number of rows in sparse matrix (%d) \
                  disagrees with size of row array (%d)" m n_rows);
            let rec loop ~i ~next =
              if i >= 0 then
                let rows_i = rows.(i) in
                if rows_i <= 0 then
                  failwith (
                    sprintf
                      "check_sparse_sane: sparse row %d contains \
                      illegal negative real row index %d" (i + 1) rows_i)
                else if rows_i >= next then
                  failwith (
                    sprintf
                      "check_sparse_sane: sparse row %d associated with \
                      real row index %d violates strict ordering (next: %d)"
                      (i + 1) rows_i next)
                else loop ~i:(i - 1) ~next:rows_i
            in
            loop ~i:(n_rows - 1) ~next:(Mat.dim2 mat + 1)

      (* Detriangularize sparse row matrix to make it symmetric *)
      let detri_sparse mat rows =
        let m = Mat.dim1 mat in
        let n = Mat.dim2 mat in
        if n < rows.(m - 1) then
          failwith "detri_sparse: sparse matrix cannot be square"
        else
          let rec loop r =
            if r > 1 then
              let r_1 = r - 1 in
              let real_r = rows.(r_1) in
              for r' = 1 to r_1 do
                mat.{r', real_r} <- mat.{r, rows.(r' - 1)}
              done;
              loop r_1
          in
          loop m

      (* Computes the symmetric additive decomposition of symmetric
         sparse matrices (in place) *)
      let symm_add_decomp_sparse mat rows =
        (* TODO: for not so sparse matrices we may want to multiply the
           non-shared elements by two instead.  Dependent algorithms
           need to be adapted as required.

           Cutoff:

             if m*m + n > real_m*real_m - m*m then
               multiply
             else
               devide
        *)
        let m = Mat.dim1 mat in
        let n = Mat.dim2 mat in
        let m_1 = m - 1 in
        for c = 1 to n do
          for i = 0 to m_1 do
            let r = rows.(i) in
            mat.{r, c} <- 0.5 *. mat.{r, c}
          done
        done

      let calc inputs =
        let kernel = Inputs.get_kernel inputs in
        let eval_inputs = inputs.Inputs.eval_inputs in
        let kn_diag, shared_diag =
          Spec.Inputs.calc_shared_diag kernel
            eval_inputs.Eval_common.Inputs.points
        in
        let eval_model =
          Eval_model.calc_internal inputs.Inputs.eval_inputs kn_diag
        in

        let km_chol = eval_inputs.Eval_inputs.inducing.Eval_inducing.km_chol in
        let inv_km_chol_kmn = eval_model.Eval_model.inv_km_chol_kmn in
        let inv_km_kmn = solve_triangular km_chol ~k:inv_km_chol_kmn in
        let inv_km = inv_chol km_chol in

        let b_chol = eval_model.Eval_model.b_chol in
        let kmn = eval_inputs.Eval_inputs.kmn in
        let inv_b_chol_kmn = solve_triangular ~trans:`T b_chol ~k:kmn in
        let inv_b_kmn = solve_triangular b_chol ~k:inv_b_chol_kmn in
        let inv_b_minus_inv_km = inv_chol b_chol in
        (* TODO: manipulate upper triangle only *)
        Mat.axpy ~alpha:(-1.) ~x:inv_km inv_b_minus_inv_km;
        Mat.detri inv_b_minus_inv_km;

        let m = Mat.dim1 inv_km_chol_kmn in
        let n = Mat.dim2 inv_km_chol_kmn in
        let inducing = inputs.Inputs.inducing in
        let inducing_shared = inducing.Inducing.shared in
        let inv_lam_sigma2_diag = eval_model.Eval_model.inv_lam_sigma2_diag in

        let update_prod_diag dst fact mat1 mat2 =
          for i = 1 to n do
            (* TODO: optimize dot and col *)
            let diag_i = dot ~x:(Mat.col mat1 i) (Mat.col mat2 i) in
            dst.{i} <- dst.{i} +. fact *. diag_i
          done
        in

        let calc_prod_trace mat1 mat2 =
          let rec loop trace i =
            if i = 0 then trace
            else
              (* TODO: optimize dot and col *)
              let diag_i = dot ~x:(Mat.col mat1 i) (Mat.col mat2 i) in
              loop (trace +. diag_i) (i - 1)
          in
          loop 0. (Mat.dim2 mat1)
        in

        let calc_evidence hyper =
          let dkm, dkm_rows = Spec.Inducing.calc_deriv inducing_shared hyper in
          check_sparse_sane dkm dkm_rows;
          let dkmn, dkmn_rows =
            Spec.Inputs.calc_deriv_cross inputs.Inputs.shared_cross hyper
          in
          check_sparse_sane dkmn dkmn_rows;
          let dlam_diag__ =
            match Spec.Inputs.calc_deriv_diag shared_diag hyper with
            | Some deriv_diag -> deriv_diag
            | None -> Vec.make0 n
          in
          let dkm_trace =
            match dkm_rows with
            | Some dkm_rows when Array.length dkm_rows <> m ->
                (* Only some inducing inputs depend on variable *)
                let _dkm_inv_km_kmn =
                  detri_sparse dkm dkm_rows;
                  symm_add_decomp_sparse dkm dkm_rows;
                  gemm dkm inv_km_kmn
                in
                (assert false (* XXX *))
            | _ ->
                (* All inducing inputs depend on variable *)
                let dkm_inv_km_kmn =
                  Mat.detri dkm;
                  symm dkm inv_km_kmn
                in
                update_prod_diag dlam_diag__ 1. inv_km_kmn dkm_inv_km_kmn;
                calc_prod_trace inv_b_minus_inv_km dkm
          in
          begin
            match dkmn_rows with
            | Some dkmn_rows when Array.length dkmn_rows <> m ->
                (assert false (* XXX *))
            | _ -> update_prod_diag dlam_diag__ (-2.) dkmn inv_km_kmn
          end;
          let dlam__trace =
            let rec loop trace i =
              if i = 0 then trace
              else
                let new_dlam_diag__ =
                  dlam_diag__.{i} *. inv_lam_sigma2_diag.{i}
                in
                dlam_diag__.{i} <- new_dlam_diag__;
                let new_trace = trace +. new_dlam_diag__ in
                loop new_trace (i - 1)
            in
            loop 0. n
          in
          let kmn_trace =
            match dkmn_rows with
            | Some dkmn_rows when Array.length dkmn_rows <> m ->
                (assert false (* XXX *))
            | _ ->
                let combined = Mat.create m n in
                for c = 1 to n do
                  let dlam__c = dlam_diag__.{c} in
                  for r = 1 to m do
                    combined.{r, c} <- 2. *. dkmn.{r, c} -. kmn.{r, c} *. dlam__c
                  done;
                done;
                calc_prod_trace inv_b_kmn combined
          in
          let trace_sum = kmn_trace +. dkm_trace +. dlam__trace in
          {
            hyper = `Hyper (hyper, dlam_diag__);
            evidence = 0.5 *. trace_sum;
          }
        in

        let calc_evidence_sigma2 () =
          let rec loop trace i =
            if i = 0 then
              {
                hyper = `Sigma2;
                evidence = 0.5 *. trace;
              }
            else
              let el =
                let inv_lam_sigma2_diag_i = inv_lam_sigma2_diag.{i} in
                (* TODO: optimize ssqr and col *)
                inv_lam_sigma2_diag_i *.
                  (1. -. Vec.ssqr (Mat.col inv_b_chol_kmn i))
              in
              loop (trace +. el) (i - 1)
          in
          loop 0. n
        in

        {
          inputs = inputs;
          eval_model = eval_model;
          shared_diag = shared_diag;
          calc_evidence = calc_evidence;
          calc_evidence_sigma2 = calc_evidence_sigma2;
        }

      let calc_eval t = t.eval_model
      let calc_evidence model hyper = model.calc_evidence hyper
      let calc_evidence_sigma2 model = model.calc_evidence_sigma2 ()
      let evidence e = e.evidence
    end

    module Variational_model = struct
      include Common_model

      let calc inputs =
        let common_model = calc inputs in
        let
          {
            Eval_model.
            inv_lam_sigma2_diag = inv_lam_sigma2_diag;
            lam_diag = lam_diag;
          } = common_model.Common_model.eval_model
        in
        let n = Vec.dim lam_diag in
        let update_trace_term traditional_evidence trace_term =
          {
            traditional_evidence with
            Common_model.evidence =
              traditional_evidence.Common_model.evidence -.
                0.5 *. trace_term;
          }
        in
        let calc_evidence hyper =
          let traditional_evidence =
            common_model.Common_model.calc_evidence hyper
          in
          match traditional_evidence.Common_model.hyper with
          | `Sigma2 -> assert false  (* impossible *)
          | `Hyper (_, dlam_diag__) ->
              let rec loop trace i =
                if i = 0 then update_trace_term traditional_evidence trace
                else
                  let el =
                    (1. -. lam_diag.{i} *. inv_lam_sigma2_diag.{i})
                      *. dlam_diag__.{i}
                  in
                  loop (trace +. el) (i - 1)
              in
              loop 0. n
        in
        let calc_evidence_sigma2 () =
          let traditional_evidence =
            common_model.Common_model.calc_evidence_sigma2 ()
          in
          let rec loop trace i =
            if i = 0 then update_trace_term traditional_evidence trace
            else
              let el =
                let inv_lam_sigma2_diag_i = inv_lam_sigma2_diag.{i} in
                (* NOTE: keep parenthesis for numeric stability *)
                inv_lam_sigma2_diag_i *.
                  (inv_lam_sigma2_diag_i *. lam_diag.{i})
              in
              loop (trace +. el) (i - 1)
          in
          loop 0. n
        in
        {
          common_model with
          calc_evidence = calc_evidence;
          calc_evidence_sigma2 = calc_evidence_sigma2;
        }
    end

    module Trained = struct
    end
  end
end

module Make_traditional_deriv (Spec : Specs.Eval_deriv) = struct
  include Make_common_deriv (Spec)

  module Eval = struct
    include Eval_common
    module Model = Common_model
  end

  module Deriv = struct
    include Deriv_common
    module Model = Common_model
  end
end

module Make_FITC_deriv (Spec : Specs.Eval_deriv) = struct
  include Make_common_deriv (Spec)

  module Eval = struct
    include Eval_common

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

  module Deriv = struct
    include Deriv_common
    module Model = Common_model
  end
end
