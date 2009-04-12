open Printf

open Lacaml.Impl.D

open Utils
open Interfaces
open Inducing_input_gpr
open Specs

(* TODO: dimension sanity checks *)

module type Sig = functor (Spec : Specs.Eval) ->
  Sigs.Eval with module Spec = Spec

module Make_common (Spec : Specs.Eval) = struct
  module Spec = Spec

  open Spec

  let jitter = !cholesky_jitter

  module Inducing = struct
    module Prepared = struct
      type t =
        {
          points : Spec.Inducing.t;
          upper : Spec.Inducing.Prepared.upper;
        }

      let calc points =
        {
          points = points;
          upper = Spec.Inducing.Prepared.calc_upper points;
        }
    end

    type t = {
      kernel : Kernel.t;
      prepared : Prepared.t;
      km : mat;
      km_chol : mat;
      log_det_km : float;
    }

    let calc_internal kernel prepared km =
      (* TODO: copy upper triangle only *)
      let km_chol = Mat.copy km in
      potrf ~jitter km_chol;
      let log_det_km = log_det km_chol in
      {
        kernel = kernel;
        prepared = prepared;
        km = km;
        km_chol = km_chol;
        log_det_km = log_det_km;
      }

    let calc kernel prepared =
      let km = Spec.Inducing.calc_upper kernel prepared.Prepared.upper in
      calc_internal kernel prepared km

    let get_kernel inducing = inducing.kernel
    let get_points inducing = inducing.prepared.Prepared.points
    let get_upper inducing = inducing.prepared.Prepared.upper
  end

  module Input = struct
    module Prepared = struct
      type t =
        {
          point : Spec.Input.t;
          inducing_points : Spec.Inducing.t;
          cross : Spec.Input.Prepared.cross;
        }

      let calc ind_prep input =
        let { Inducing.Prepared.points = points; upper = upper } = ind_prep in
        {
          point = input;
          inducing_points = points;
          cross = Spec.Input.Prepared.calc_cross upper input;
        }
    end

    type t = {
      inducing : Inducing.t;
      point : Spec.Input.t;
      k_m : vec;
    }

    let calc inducing prepared =
      if Inducing.get_points inducing != prepared.Prepared.inducing_points then
        failwith
          "Input.calc: inducing points incompatible with \
          those used for prepared";
      {
        inducing = inducing;
        point = prepared.Prepared.point;
        k_m = Spec.Input.eval inducing.Inducing.kernel prepared.Prepared.cross;
      }

    let get_kernel t = t.inducing.Inducing.kernel
  end

  module Inputs = struct
    module Prepared = struct
      type t =
        {
          points : Spec.Inputs.t;
          inducing_points : Spec.Inducing.t;
          cross : Spec.Inputs.Prepared.cross;
        }

      let calc ind_prep inputs =
        let { Inducing.Prepared.points = points; upper = upper } = ind_prep in
        {
          points = inputs;
          inducing_points = points;
          cross = Spec.Inputs.Prepared.calc_cross upper inputs;
        }
    end

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

    let calc inducing prepared =
      if Inducing.get_points inducing != prepared.Prepared.inducing_points
      then
        failwith
          "Inputs.calc: inducing points incompatible with \
          those used for prepared";
      let kernel = inducing.Inducing.kernel in
      let kmn = Inputs.calc_cross kernel prepared.Prepared.cross in
      calc_internal inducing prepared.Prepared.points kmn

    let get_kernel t = t.inducing.Inducing.kernel

    (* TODO: (re?)move *)
    let nystrom_chol_marginals inputs =
      let inducing = inputs.inducing in
      let kernel = get_kernel inputs in
      let kn_diag =
        match Inputs.calc_diag kernel inputs.points with
        | `Single kn_diag -> kn_diag
        | `Block _block_diag ->
            (* TODO *)
            (assert false (* XXX *))
      in
      let inv_km_chol_kmn =
        solve_triangular ~trans:`T inducing.Inducing.km_chol ~k:inputs.kmn
      in
      inv_km_chol_kmn, kn_diag
  end

  module Common_model = struct
    type t = {
      sigma2 : float;
      inputs : Inputs.t;
      kn_diag : vec;
      inv_km_chol_kmn : mat;
      lam_diag : vec;
      inv_lam_sigma2_diag : vec;
      b_chol : mat;
      log_evidence : float;
    }

    let check_sigma2 sigma2 =
      if sigma2 < 0. then failwith "Model.check_sigma2: sigma2 < 0"

    let calc_internal inputs sigma2 kn_diag =
      check_sigma2 sigma2;
      let inducing = inputs.Inputs.inducing in
      let kmn = inputs.Inputs.kmn in
      let inv_km_chol_kmn =
        solve_triangular ~trans:`T inducing.Inducing.km_chol ~k:kmn
      in
      let kmn_ = Mat.copy kmn in
      let n_inputs = Vec.dim kn_diag in
      let lam_diag = Vec.create n_inputs in
      let inv_lam_sigma2_diag = Vec.create n_inputs in
      let rec loop log_det_lam_sigma2 i =
        if i = 0 then log_det_lam_sigma2
        else
          let kn_diag_i = kn_diag.{i} in
          (* TODO: optimize sqr_nrm2 and col *)
          let qn_diag_i = Vec.sqr_nrm2 (Mat.col inv_km_chol_kmn i) in
          let lam_i = kn_diag_i -. qn_diag_i in
          lam_diag.{i} <- lam_i;
          let lam_sigma2_i = lam_i +. sigma2 in
          let inv_lam_sigma2_diag_i = 1. /. lam_sigma2_i in
          inv_lam_sigma2_diag.{i} <- inv_lam_sigma2_diag_i;
          (* TODO: optimize scal col *)
          scal (sqrt inv_lam_sigma2_diag_i) (Mat.col kmn_ i);
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
        sigma2 = sigma2;
        inputs = inputs;
        kn_diag = kn_diag;
        inv_km_chol_kmn = inv_km_chol_kmn;
        lam_diag = lam_diag;
        inv_lam_sigma2_diag = inv_lam_sigma2_diag;
        b_chol = b_chol;
        log_evidence = 0.5 *. (l1_2 -. float n_inputs *. log_2pi);
      }

    let calc inputs ~sigma2 =
      let kernel = Inputs.get_kernel inputs in
      let kn_diag =
        match Spec.Inputs.calc_diag kernel inputs.Inputs.points with
        | `Single kn_diag -> kn_diag
        | `Block _block_diag ->
            (* TODO *)
            (assert false (* XXX *))
      in
      calc_internal inputs sigma2 kn_diag

    let update_sigma2 model sigma2 =
      calc_internal model.inputs sigma2 model.kn_diag

    let calc_log_evidence model = model.log_evidence
    let get_inducing model = model.inputs.Inputs.inducing
    let get_inducing_points model = Inducing.get_points (get_inducing model)
    let get_upper model = Inducing.get_upper (get_inducing model)
    let get_input_points model = model.inputs.Inputs.points
    let get_kernel model = (get_inducing model).Inducing.kernel
    let get_sigma2 model = model.sigma2
    let get_kmn model = model.inputs.Inputs.kmn
    let get_lam_diag model = model.lam_diag
    let get_km model = (get_inducing model).Inducing.km
  end

  module Variational_model = struct
    include Common_model

    let calc_from_model model =
      let
        {
          lam_diag = lam_diag;
          inv_lam_sigma2_diag = x;
          log_evidence = log_evidence;
        } = model
      in
      { model with log_evidence = log_evidence -. 0.5 *. dot ~x lam_diag }

    let calc_internal inputs sigma2 kn_diag =
      calc_from_model (calc_internal inputs sigma2 kn_diag)

    let calc inputs ~sigma2 = calc_from_model (calc inputs ~sigma2)
  end

  module Trained = struct
    type t = {
      model : Common_model.t;
      targets : vec;
      y__ : vec;
      inv_b_chol_kmn_y__ : vec;
      log_evidence : float;
    }

    let calc_y__ssqr model targets =
      let n_targets = Vec.dim targets in
      let y__ = Vec.create n_targets in
      let inv_lam_sigma2_diag = model.Common_model.inv_lam_sigma2_diag in
      for i = 1 to n_targets do
        y__.{i} <- targets.{i} *. inv_lam_sigma2_diag.{i}
      done;
      y__, dot ~x:targets y__

    let make model targets y__ ssqr_y__ inv_b_chol_kmn_y__ =
      let fit_log_evidence =
        0.5 *. (Vec.sqr_nrm2 inv_b_chol_kmn_y__ -. ssqr_y__)
      in
      {
        model = model;
        targets = targets;
        y__ = y__;
        inv_b_chol_kmn_y__ = inv_b_chol_kmn_y__;
        log_evidence = model.Common_model.log_evidence +. fit_log_evidence;
      }

    let calc_internal model targets inv_b_chol_kmn =
      let y__, ssqr_y__ = calc_y__ssqr model targets in
      let inv_b_chol_kmn_y__ = gemv inv_b_chol_kmn y__ in
      make model targets y__ ssqr_y__ inv_b_chol_kmn_y__

    let calc model ~targets =
      let y__, ssqr_y__ = calc_y__ssqr model targets in
      let inv_b_chol_kmn_y__ = gemv (Common_model.get_kmn model) y__ in
      trsv ~trans:`T model.Common_model.b_chol inv_b_chol_kmn_y__;
      make model targets y__ ssqr_y__ inv_b_chol_kmn_y__

    let calc_log_evidence trained = trained.log_evidence

    let get_kmn trained = trained.model.Common_model.inputs.Inputs.kmn
  end

  module Weights = struct
    type t = { model : Common_model.t; coeffs : vec }

    let get_kernel weights = Common_model.get_kernel weights.model
    let get_inducing weights = Common_model.get_inducing weights.model
    let get_upper weights = Common_model.get_upper weights.model
    let get_coeffs weights = weights.coeffs

    let calc trained =
      let coeffs = copy trained.Trained.inv_b_chol_kmn_y__ in
      trsv trained.Trained.model.Common_model.b_chol coeffs;
      {
        model = trained.Trained.model;
        coeffs = coeffs;
      }

    let get_inducing_points weights =
      Common_model.get_inducing_points weights.model
  end

  module Mean = struct
    type t = { point : Spec.Input.t; value : float }

    let make ~point ~value = { point = point; value = value }

    let calc_input weights point =
      let upper = Weights.get_upper weights in
      let kernel = Weights.get_kernel weights in
      let coeffs = Weights.get_coeffs weights in
      let prepared = Spec.Input.Prepared.calc_cross upper point in
      let value = Spec.Input.weighted_eval kernel ~coeffs prepared in
      make ~point ~value

    let calc_induced weights input =
      if Weights.get_inducing weights <> input.Input.inducing then
        failwith
          "Mean.calc_induced: \
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
      let upper = Weights.get_upper weights in
      let prepared = Spec.Inputs.Prepared.calc_cross upper points in
      let values = Spec.Inputs.weighted_eval kernel ~coeffs prepared in
      make ~points ~values

    let calc_induced weights inputs =
      let { Inputs.points = points; kmn = kmn } = inputs in
      if Weights.get_inducing weights <> inputs.Inputs.inducing then
        failwith
          "Means.calc_induced: \
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
        (* TODO: optimize sqr_nrm2 and col *)
        variances.{i} <- variances.{i} +. Vec.sqr_nrm2 (Mat.col inv_b_chol_kmn i)
      done;
      make ~points:(Common_model.get_input_points model) ~variances ~model

    let calc_induced model inputs =
      if Common_model.get_inducing model <> inputs.Inputs.inducing then
        failwith
          "Variances.calc_induced: \
          model and inputs disagree about inducing points";
      let kmt = inputs.Inputs.kmn in
      let inv_km_chol_kmt, kt_diag = Inputs.nystrom_chol_marginals inputs in
      let variances = copy kt_diag in
      let inv_b_chol_kmt = solve_b_chol model ~k:kmt in
      let n = Mat.dim2 inv_b_chol_kmt in
      for i = 1 to n do
        let explained_variance =
          (* TODO: optimize sqr_nrm2 and col *)
          Vec.sqr_nrm2 (Mat.col inv_km_chol_kmt i) -.
            Vec.sqr_nrm2 (Mat.col inv_b_chol_kmt i)
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
          (* TODO: optimize sqr_nrm2 and col *)
          variances.{i} <- Vec.sqr_nrm2 (Mat.col inv_b_chol_km i)
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

    let calc_common ~kmn ~inv_km_chol_kmn ~points ~model =
      let kernel = Common_model.get_kernel model in
      let covariances = Spec.Inputs.calc_upper kernel points in
      ignore (syrk ~trans:`T ~alpha:(-1.) inv_km_chol_kmn ~beta:1. ~c:covariances);
      let inv_b_chol_kmn = solve_b_chol model ~k:kmn in
      ignore (syrk ~trans:`T ~alpha:1. inv_b_chol_kmn ~beta:1. ~c:covariances);
      make ~points ~covariances ~model

    let calc_model_inputs model =
      let kmn = model.Common_model.inputs.Inputs.kmn in
      let inv_km_chol_kmn = model.Common_model.inv_km_chol_kmn in
      let points = Common_model.get_input_points model in
      calc_common ~kmn ~inv_km_chol_kmn ~points ~model

    let calc_induced model inputs =
      if Common_model.get_inducing model <> inputs.Inputs.inducing then
        failwith (
          "FITC_covariances.calc_induced: \
          model and inputs disagree about inducing points");
      let kmn = inputs.Inputs.kmn in
      let inv_km_chol_kmn, _ = Inputs.nystrom_chol_marginals inputs in
      let points = inputs.Inputs.points in
      calc_common ~kmn ~inv_km_chol_kmn ~points ~model
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
          "FIC_covariances.calc_induced: \
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

module FIC_Loc = struct let loc = "FIC" end
module Variational_FITC_Loc = struct let loc = "Variational_FITC" end
module Variational_FIC_Loc = struct let loc = "Variational_FIC" end

let fitc_loc = "FITC"
let fic_loc = "FIC"
let variational_fitc_loc = "Variational_FITC"
let variational_fic_loc = "Variational_FIC"

module Make_FITC (Spec : Specs.Eval) = struct
  include Make_common (Spec)

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

module Make_FIC (Spec : Specs.Eval) = struct
  include Make_common (Spec)

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

module Make_variational_FITC (Spec : Specs.Eval) = struct
  include Make_common (Spec)

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

module Make_variational_FIC (Spec : Specs.Eval) = struct
  include Make_common (Spec)

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

module type Deriv_sig = functor (Spec : Specs.Deriv) ->
  Sigs.Deriv
    with module Eval.Spec = Spec.Eval
    with module Deriv.Spec = Spec

module Make_common_deriv (Spec : Specs.Deriv) = struct
  open Spec

  module Eval_common = Make_common (Spec.Eval)
  module Eval_inducing = Eval_common.Inducing
  module Eval_inputs = Eval_common.Inputs
  module Eval_model = Eval_common.Common_model
  module Eval_trained = Eval_common.Trained

  module Deriv_common = struct
    module Spec = Spec

    module Inducing = struct
      module Prepared = struct
        type t = {
          eval : Eval_inducing.Prepared.t;
          upper : Spec.Inducing.Prepared.upper;
        }

        let calc eval =
          {
            eval = eval;
            upper =
              Spec.Inducing.Prepared.calc_upper
                eval.Eval_inducing.Prepared.upper;
          }

        let get_points prepared = prepared.eval.Eval_inducing.Prepared.points
      end

      type t = {
        eval : Eval_inducing.t;
        shared : Spec.Inducing.shared;
      }

      let calc kernel { Prepared.eval = eval; upper = upper } =
        let km, shared = Spec.Inducing.calc_shared_upper kernel upper in
        {
          eval = Eval_inducing.calc_internal kernel eval km;
          shared = shared;
        }

      let calc_eval inducing = inducing.eval
      let get_kernel inducing = Eval_inducing.get_kernel inducing.eval
      let get_points inducing = Eval_inducing.get_points inducing.eval
    end

    module Inputs = struct
      module Prepared = struct
        type t = {
          inducing_prepared : Inducing.Prepared.t;
          eval : Eval_inputs.Prepared.t;
          cross : Spec.Inputs.Prepared.cross;
        }

        let calc inducing_prepared eval =
          if
            Inducing.Prepared.get_points inducing_prepared
              != eval.Eval_inputs.Prepared.inducing_points
          then
            failwith
              "Deriv.Inputs.Prepared.calc: prepared inducing points \
              incompatible with those used for prepared inputs";
          {
            inducing_prepared = inducing_prepared;
            eval = eval;
            cross =
              Spec.Inputs.Prepared.calc_cross
                inducing_prepared.Inducing.Prepared.upper
                eval.Eval_inputs.Prepared.cross;
          }

        let get_inducing_points prepared =
          Inducing.Prepared.get_points prepared.inducing_prepared
      end

      type t = {
        inducing : Inducing.t;
        eval : Eval_inputs.t;
        shared_cross : Spec.Inputs.cross;
      }

      let calc inducing prepared =
        if Inducing.get_points inducing != Prepared.get_inducing_points prepared
        then
          failwith
            "Deriv.Inputs.calc: inducing points \
            incompatible with those used for prepared inputs";
        let kernel = Inducing.get_kernel inducing in
        let kmn, shared_cross =
          Spec.Inputs.calc_shared_cross kernel prepared.Prepared.cross
        in
        let eval =
          Eval_inputs.calc_internal
            inducing.Inducing.eval
            prepared.Prepared.eval.Eval_inputs.Prepared.points
            kmn
        in
        { inducing = inducing; eval = eval; shared_cross = shared_cross; }

      let calc_eval t = t.eval

      let get_kernel inputs = Inducing.get_kernel inputs.inducing
    end

    type model_kind = Traditional | Variational

    module Common_model = struct
      type t = {
        inputs : Inputs.t;
        shared_diag : Spec.Inputs.diag;
        eval_model : Eval_model.t;
        inv_b_chol_kmn : mat;
        inv_b_kmn : mat;
        model_kind : model_kind;
      }

      let calc_internal inputs eval_model shared_diag model_kind =
        let eval_inputs = inputs.Inputs.eval in
        let b_chol = eval_model.Eval_model.b_chol in
        let kmn = eval_inputs.Eval_inputs.kmn in
        let inv_b_chol_kmn = solve_triangular ~trans:`T b_chol ~k:kmn in
        let inv_b_kmn = solve_triangular b_chol ~k:inv_b_chol_kmn in
        {
          inputs = inputs;
          shared_diag = shared_diag;
          eval_model = eval_model;
          inv_b_chol_kmn = inv_b_chol_kmn;
          inv_b_kmn = inv_b_kmn;
          model_kind = model_kind;
        }

      let calc_common model_kind inputs sigma2 =
        let kernel = Inputs.get_kernel inputs in
        let eval_inputs = inputs.Inputs.eval in
        let kn_diag, shared_diag =
          Spec.Inputs.calc_shared_diag kernel eval_inputs.Eval_inputs.points
        in
        let eval_model =
          match model_kind with
          | Traditional ->
              Eval_model.calc_internal inputs.Inputs.eval sigma2 kn_diag
          | Variational ->
              Eval_common.Variational_model.calc_internal
                inputs.Inputs.eval sigma2 kn_diag
        in
        calc_internal inputs eval_model shared_diag model_kind

      let calc_eval model = model.eval_model

      let calc inputs ~sigma2 = calc_common Traditional inputs sigma2

      let update_sigma2 ({ model_kind = model_kind } as model) sigma2 =
        let eval_model =
          match model_kind with
          | Traditional -> Eval_model.update_sigma2 model.eval_model sigma2
          | Variational ->
              Eval_common.Variational_model.update_sigma2
                model.eval_model sigma2
        in
        calc_internal model.inputs eval_model model.shared_diag model_kind

      (**)

      type log_evidence_sigma2 = {
        sigma2_model : t;
        log_evidence_sigma2 : float;
      }

      let calc_log_evidence_sigma2 ({ eval_model = eval_model } as model) =
        let inv_lam_sigma2_diag = eval_model.Eval_model.inv_lam_sigma2_diag in
        let inv_b_chol_kmn = model.inv_b_chol_kmn in
        let rec loop trace i =
          if i = 0 then
            let log_evidence_sigma2 = 0.5 *. trace in
            (
              log_evidence_sigma2,
              {
                sigma2_model = model;
                log_evidence_sigma2 = log_evidence_sigma2;
              }
            )
          else
            let el =
              let inv_lam_sigma2_diag_i = inv_lam_sigma2_diag.{i} in
              (* TODO: optimize sqr_nrm2 and col *)
              let sqr_nrm2 = Vec.sqr_nrm2 (Mat.col inv_b_chol_kmn i) in
              inv_lam_sigma2_diag_i *. (inv_lam_sigma2_diag_i *. sqr_nrm2 -. 1.)
            in
            loop (trace +. el) (i - 1)
        in
        loop 0. (Mat.dim2 inv_b_chol_kmn)

      (**)

      type hyper_t = {
        model : t;
        inv_km_kmn : mat;
        inv_km_minus_inv_b : mat;
      }

      let prepare_hyper model =
        let { inputs = inputs; eval_model = eval_model } = model in
        let eval_inputs = inputs.Inputs.eval in
        let km_chol = eval_inputs.Eval_inputs.inducing.Eval_inducing.km_chol in
        let inv_km_chol_kmn = eval_model.Eval_model.inv_km_chol_kmn in
        let inv_km_kmn = solve_triangular km_chol ~k:inv_km_chol_kmn in
        let b_chol = eval_model.Eval_model.b_chol in
        let inv_b = inv_chol b_chol in
        let inv_km_minus_inv_b = inv_chol km_chol in
        (* TODO: manipulate upper triangle only *)
        Mat.axpy ~alpha:(-1.) ~x:inv_b inv_km_minus_inv_b;
        Mat.detri inv_km_minus_inv_b;
        {
          model = model;
          inv_km_kmn = inv_km_kmn;
          inv_km_minus_inv_b = inv_km_minus_inv_b;
        }

      (**)

      type log_evidence = {
        hyper_model : hyper_t;
        var : Hyper.t;
        deriv_upper : symm_mat_deriv;
        dlam_diag__ : vec;
        deriv_cross : mat_deriv;
        log_evidence_hyper : float;
      }

      let calc_log_evidence hyper_model hyper =
        let
          {
            model =
              {
                inputs = inputs;
                shared_diag = shared_diag;
                eval_model = eval_model;
                inv_b_kmn = inv_b_kmn;
              };
            inv_km_kmn = inv_km_kmn;
            inv_km_minus_inv_b = inv_km_minus_inv_b;
          } = hyper_model
        in
        let kmn = inputs.Inputs.eval.Eval_inputs.kmn in
        let m = Mat.dim1 kmn in
        let n = Mat.dim2 kmn in
        let inv_lam_sigma2_diag = eval_model.Eval_model.inv_lam_sigma2_diag in
        let inducing_shared = inputs.Inputs.inducing.Inducing.shared in
        let dlam_diag__ =
          match Spec.Inputs.calc_deriv_diag shared_diag hyper with
          | `Vec deriv_diag -> deriv_diag
          | `Const c -> Vec.make n c
          | `Factor _ -> assert false (* XXX *)
        in
        let deriv_upper =
          Spec.Inducing.calc_deriv_upper inducing_shared hyper
        in
        let dkm_trace =
          match deriv_upper with
          | `Dense dkm ->
              (* All inducing inputs depend on variable *)
              Mat.detri dkm;
              let dkm_inv_km_kmn = symm dkm inv_km_kmn in
              update_prod_diag dlam_diag__ 1. inv_km_kmn dkm_inv_km_kmn;
              calc_prod_trace inv_km_minus_inv_b dkm
          | `Sparse_rows (dkm, dkm_rows) ->
              check_sparse_sane dkm dkm_rows ~real_m:m;
              (* Only some inducing inputs depend on variable *)
              let dkm_inv_km_kmn =
                detri_sparse dkm dkm_rows;
                symm_add_decomp_sparse dkm dkm_rows;
                gemm dkm inv_km_kmn
              in
              let n_rows = Array.length dkm_rows in
              let trace = ref 0. in
              for i = 1 to n_rows do
                let r = dkm_rows.(i - 1) in
                for c = 1 to n do
                  dlam_diag__.{c} <-
                    dlam_diag__.{c}
                      +. 2. *. inv_km_kmn.{r, c} *. dkm_inv_km_kmn.{i, c}
                done;
                for c = 1 to m do
                  trace := !trace +. inv_km_minus_inv_b.{r, c} *. dkm.{i, c}
                done
              done;
              !trace
          | `Diag_vec _vec ->
              (assert false (* XXX *))
          | `Diag_const _c ->
              (assert false (* XXX *))
          | `Const _c ->
              (assert false (* XXX *))
          | `Factor _c ->
              (assert false (* XXX *))
        in
        let deriv_cross =
          Spec.Inputs.calc_deriv_cross inputs.Inputs.shared_cross hyper
        in
        begin
          match deriv_cross with
          | `Dense dkmn -> update_prod_diag dlam_diag__ (-2.) dkmn inv_km_kmn
          | `Sparse_rows (dkmn, dkmn_rows) ->
              check_sparse_sane dkmn dkmn_rows ~real_m:m;
              let n_rows = Array.length dkmn_rows in
              let diag_el = ref 0. in
              for c = 1 to n do
                for i = 1 to n_rows do
                  let r = dkmn_rows.(i - 1) in
                  diag_el := !diag_el +. dkmn.{i, c} *. inv_km_kmn.{r, c}
                done;
                dlam_diag__.{c} <- dlam_diag__.{c} -. 2. *. !diag_el;
                diag_el := 0.
              done
          | `Const _c ->
              (assert false (* XXX *))
          | `Factor _c ->
              (assert false (* XXX *))
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
        let dkmn_trace =
          match deriv_cross with
          | `Dense dkmn ->
              let combined = Mat.create m n in
              for c = 1 to n do
                let dlam__c = dlam_diag__.{c} in
                let inv_lam_sigma2_diag_c = inv_lam_sigma2_diag.{c} in
                for r = 1 to m do
                  combined.{r, c} <-
                    inv_lam_sigma2_diag_c
                      *. (2. *. dkmn.{r, c} -. kmn.{r, c} *. dlam__c)
                done;
              done;
              calc_prod_trace inv_b_kmn combined
          | `Sparse_rows (dkmn, dkmn_rows) ->
              let combined =
                (* TODO: gift *)
                Mat.copy kmn
              in
              let n_rows = Array.length dkmn_rows in
              for c = 1 to n do
                (* TODO: optimize scal col *)
                scal (-. dlam_diag__.{c}) (Mat.col combined c);
                let inv_lam_sigma2_diag_c = inv_lam_sigma2_diag.{c} in
                for i = 1 to n_rows do
                  let r = dkmn_rows.(i - 1) in
                  combined.{r, c} <-
                    inv_lam_sigma2_diag_c *.
                      (2. *. dkmn.{i, c} +. combined.{r, c})
                done
              done;
              calc_prod_trace inv_b_kmn combined
          | `Const const ->
              let combined = Mat.create m n in
              let const2 = 2. *. const in
              for c = 1 to n do
                let dlam__c = dlam_diag__.{c} in
                let inv_lam_sigma2_diag_c = inv_lam_sigma2_diag.{c} in
                for r = 1 to m do
                  combined.{r, c} <-
                    inv_lam_sigma2_diag_c *. (const2 -. kmn.{r, c} *. dlam__c)
                done
              done;
              calc_prod_trace inv_b_kmn combined
          | `Factor _c ->
              (assert false (* XXX *))
        in
        let log_evidence_hyper =
          0.5 *. (dkm_trace -. dkmn_trace -. dlam__trace)
        in
        (
          log_evidence_hyper,
          {
            hyper_model = hyper_model;
            var = hyper;
            deriv_upper = deriv_upper;
            dlam_diag__ = dlam_diag__;
            deriv_cross = deriv_cross;
            log_evidence_hyper = log_evidence_hyper;
          }
        )
    end

    module Cm = Common_model

    module Variational_model = struct
      include Cm

      let calc inputs ~sigma2 = calc_common Variational inputs sigma2

      (**)

      let calc_log_evidence_sigma2 common_model =
        let
          {
            Eval_model.
            inv_lam_sigma2_diag = inv_lam_sigma2_diag;
            lam_diag = lam_diag;
          } = common_model.Cm.eval_model
        in
        let rec loop trace i =
          if i = 0 then
            let log_evidence, traditional =
              Cm.calc_log_evidence_sigma2 common_model
            in
            let variational_log_evidence = log_evidence -. 0.5 *. trace in
            (
              variational_log_evidence,
              {
                traditional with
                Cm.log_evidence_sigma2 = variational_log_evidence;
              }
            )
          else
            let el =
              let inv_lam_sigma2_diag_i = inv_lam_sigma2_diag.{i} in
              inv_lam_sigma2_diag_i *. (inv_lam_sigma2_diag_i *. lam_diag.{i})
            in
            loop (trace -. el) (i - 1)
        in
        loop 0. (Vec.dim lam_diag)

      (**)

      let calc_log_evidence common_hyper_model hyper =
        let
          {
            Eval_model.
            inv_lam_sigma2_diag = inv_lam_sigma2_diag;
            lam_diag = lam_diag;
          } = common_hyper_model.model.Cm.eval_model
        in
        let log_evidence, traditional =
          Cm.calc_log_evidence common_hyper_model hyper
        in
        let { dlam_diag__ = dlam_diag__ } = traditional in
        let rec loop trace i =
          if i = 0 then
            let variational_log_evidence = log_evidence -. 0.5 *. trace in
            (
              variational_log_evidence,
              {
                traditional with
                Cm.log_evidence_hyper = variational_log_evidence;
              }
            )
          else
            let el =
              dlam_diag__.{i} *. (1. -. lam_diag.{i} *. inv_lam_sigma2_diag.{i})
            in
            loop (trace +. el) (i - 1)
        in
        loop 0. (Vec.dim lam_diag)
    end

    module Trained = struct
      type t = {
        common_model : Cm.t;
        eval_trained : Eval_trained.t;
        inv_b_kmn_y__ : vec;
        knm_inv_b_kmn_y__ : vec;
        __knm_inv_b_kmn_y__ : vec;
      }

      let calc common_model ~targets =
        let eval_model = common_model.Cm.eval_model in
        let eval_trained =
          Eval_trained.calc_internal eval_model targets
            common_model.Cm.inv_b_chol_kmn
        in
        let y__ = eval_trained.Eval_trained.y__ in
        let inv_b_kmn_y__ = gemv common_model.Cm.inv_b_kmn y__ in
        let kmn = Eval_trained.get_kmn eval_trained in
        let knm_inv_b_kmn_y__ = gemv ~trans:`T kmn inv_b_kmn_y__ in
        let n = Vec.dim y__ in
        let __knm_inv_b_kmn_y__ = Vec.create n in
        let inv_lam_sigma2_diag =
          eval_trained.Eval_trained.model.Eval_model.inv_lam_sigma2_diag
        in
        for i = 1 to n do
          __knm_inv_b_kmn_y__.{i} <-
            inv_lam_sigma2_diag.{i} *. knm_inv_b_kmn_y__.{i}
        done;
        {
          common_model = common_model;
          eval_trained = eval_trained;
          inv_b_kmn_y__ = inv_b_kmn_y__;
          knm_inv_b_kmn_y__ = knm_inv_b_kmn_y__;
          __knm_inv_b_kmn_y__ = __knm_inv_b_kmn_y__;
        }

      let calc_eval trained = trained.eval_trained

      (**)

      let raise_incompatible_model loc =
        failwith (
          sprintf
            "Deriv.Trained.%s: model used for training log evidence \
            and model used for model log evidence are not the same" loc)

      let calc_log_evidence_sigma2 trained model_log_evidence_sigma2 =
        if
          trained.common_model != model_log_evidence_sigma2.Cm.sigma2_model
        then raise_incompatible_model "calc_log_evidence_sigma2";
        let eval_trained = trained.eval_trained in
        let y__ = eval_trained.Eval_trained.y__ in
        let __knm_inv_b_kmn_y__ = trained.__knm_inv_b_kmn_y__ in
        let rec loop sum i =
          if i = 0 then
            model_log_evidence_sigma2.Cm.log_evidence_sigma2 +. 0.5 *. sum
          else
            let new_sum =
              let z__i = __knm_inv_b_kmn_y__.{i} in
              let y__i = y__.{i} in
              sum +. y__i *. y__i +. z__i *. (z__i -. 2. *. y__i)
            in
            loop new_sum (i - 1)
        in
        loop 0. (Vec.dim __knm_inv_b_kmn_y__)

      (**)

      type hyper_t = {
        trained : t;
        hyper_model : Cm.hyper_t;
        dlam_factor : vec;
        dkmn_factor : vec;
      }

      let prepare_hyper trained hyper_model =
        let eval_trained = trained.eval_trained in
        let y__ = eval_trained.Eval_trained.y__ in
        let y = eval_trained.Eval_trained.targets in
        let knm_inv_b_kmn_y__ = trained.knm_inv_b_kmn_y__ in
        let __knm_inv_b_kmn_y__ = trained.__knm_inv_b_kmn_y__ in
        let dkm_multipliers = copy knm_inv_b_kmn_y__ in
        axpy ~alpha:(-1.) ~x:y__ dkm_multipliers;
        scal 2. dkm_multipliers;
        let n = Vec.dim y in
        let dlam_factor = Vec.create n in
        let dkmn_factor = Vec.create n in
        for i = 1 to n do
          let z_i = knm_inv_b_kmn_y__.{i} in
          let z__i = __knm_inv_b_kmn_y__.{i} in
          let y__i = y__.{i} in
          dlam_factor.{i} <- (2.*.z_i -. y.{i}) *. y__i -. z_i *. z__i;
          dkmn_factor.{i} <- 2. *. (z__i -. y__i);
        done;
        {
          trained = trained;
          hyper_model = hyper_model;
          dlam_factor = dlam_factor;
          dkmn_factor = dkmn_factor;
        }

      (**)

      let calc_log_evidence hyper_trained model_log_evidence =
        if hyper_trained.hyper_model != model_log_evidence.Cm.hyper_model then
          raise_incompatible_model "calc_log_evidence";
        let inv_b_kmn_y__ = hyper_trained.trained.inv_b_kmn_y__ in
        let deriv_upper = model_log_evidence.Cm.deriv_upper in
        let dkm_nll =
          match deriv_upper with
          | `Dense dkm -> dot ~x:(gemv dkm inv_b_kmn_y__) inv_b_kmn_y__
          | `Sparse_rows _ ->
              (assert false (* XXX *))
          | `Diag_vec _ ->
              (assert false (* XXX *))
          | `Diag_const _ ->
              (assert false (* XXX *))
          | `Const _ ->
              (assert false (* XXX *))
          | `Factor _c ->
              (assert false (* XXX *))
        in
        let deriv_cross = model_log_evidence.Cm.deriv_cross in
        let dkmn_factor = hyper_trained.dkmn_factor in
        let dkmn_nll =
          match deriv_cross with
          | `Dense dkmn ->
              dot ~x:(gemv ~trans:`T dkmn inv_b_kmn_y__) dkmn_factor
          | `Sparse_rows _ ->
              (assert false (* XXX *))
          | `Const _ ->
              (assert false (* XXX *))
          | `Factor _c ->
              (assert false (* XXX *))
        in
        let dlam_diag__ = model_log_evidence.Cm.dlam_diag__ in
        let dlam_nll = dot ~x:hyper_trained.dlam_factor dlam_diag__ in
        let nll = 0.5 *. (dkm_nll +. dkmn_nll +. dlam_nll) in
        model_log_evidence.Cm.log_evidence_hyper -. nll
    end
  end
end

module Make_FITC_deriv (Spec : Specs.Deriv) = struct
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

module Make_FIC_deriv (Spec : Specs.Deriv) = struct
  include Make_common_deriv (Spec)

  module Eval = struct
    include Eval_common

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

  module Deriv = struct
    include Deriv_common
    module Model = Common_model
  end
end

module Make_variational_FITC_deriv (Spec : Specs.Deriv) = struct
  include Make_common_deriv (Spec)

  module Eval = struct
    include Eval_common

    module Model = Variational_model

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
    module Model = Variational_model
  end
end

module Make_variational_FIC_deriv (Spec : Specs.Deriv) = struct
  include Make_common_deriv (Spec)

  module Eval = struct
    include Eval_common

    module Model = Variational_model

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

  module Deriv = struct
    include Deriv_common
    module Model = Variational_model
  end
end

module Make_deriv (Spec : Specs.Deriv) = struct
  module type Sig = Sigs.Deriv
    with module Eval.Spec = Spec.Eval
    with module Deriv.Spec = Spec

  module Common_deriv = Make_common_deriv (Spec)

  module FITC = struct
    include Common_deriv

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

  module FIC = struct
    include Common_deriv

    module Eval = struct
      include Eval_common

      module Model = Variational_model

      module Covariances = FIC_covariances

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
      module Model = Variational_model
    end
  end

  module Variational_FITC = struct
    include Common_deriv

    module Eval = struct
      include Eval_common

      module Model = Variational_model

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
      module Model = Variational_model
    end
  end

  module Variational_FIC = struct
    include Common_deriv

    module Eval = struct
      include Eval_common

      module Model = Variational_model

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

    module Deriv = struct
      include Deriv_common
      module Model = Variational_model
    end
  end
end
