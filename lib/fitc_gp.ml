open Format
open Bigarray

open Lacaml.Impl.D

open Utils
open Interfaces
open Specs

module type Sig = functor (Spec : Specs.Eval) ->
  Sigs.Eval with module Spec = Spec

(* Computations shared by FIC and FITC, and standard and variational
   version *)
module Make_common (Spec : Specs.Eval) = struct
  module Spec = Spec

  open Spec

  let jitter = !cholesky_jitter

  (* Evaluation of inducing points *)
  module Inducing = struct
    module Prepared = struct
      type t = {
        points : Spec.Inducing.t;
        upper : Spec.Inducing.Prepared.upper;
      }

      let calc points =
        {
          points = points;
          upper = Spec.Inducing.Prepared.calc_upper points;
        }

      let id x = x

      let check_n_inducing ~n_inducing inputs =
        let n_inputs = Spec.Inputs.get_n_inputs inputs in
        if n_inputs < 1 || n_inducing > n_inputs then
          failwith
            (sprintf
              "Gpr.Fitc_gp.Make_common.check_n_inducing: \
              violating 1 <= n_inducing (%d) <= n_inputs (%d)"
              n_inducing n_inputs)

      let choose kernel inputs indexes =
        let chosen_inputs = Spec.Inputs.choose_subset inputs indexes in
        calc (Spec.Inputs.create_inducing kernel chosen_inputs)

      let choose_n_first_inputs kernel ~n_inducing inputs =
        check_n_inducing ~n_inducing inputs;
        let indexes = Int_vec.create n_inducing in
        for i = 1 to n_inducing do indexes.{i} <- i done;
        choose kernel inputs indexes

      let choose_n_random_inputs
            ?(rnd_state = Random.get_state ()) kernel ~n_inducing inputs =
        check_n_inducing ~n_inducing inputs;
        let n_inputs = Spec.Inputs.get_n_inputs inputs in
        let indexes = Int_vec.create n_inputs in
        for i = 1 to n_inputs do indexes.{i} <- i done;
        for i = 1 to n_inducing do
          let rnd_index = Random.State.int rnd_state (n_inputs - i + 1) + 1 in
          let tmp = indexes.{rnd_index} in
          indexes.{rnd_index} <- indexes.{i};
          indexes.{i} <- tmp;
        done;
        let indexes = Array1.sub indexes 1 n_inducing in
        choose kernel inputs indexes

      let get_points inducing = inducing.points
    end

    type t = {
      kernel : Kernel.t;
      prepared : Prepared.t;
      km : mat;
      chol_km : mat;
      log_det_km : float;
    }

    let calc_internal kernel prepared km =
      let chol_km = lacpy ~uplo:`U km in
      potrf ~jitter chol_km;
      let log_det_km = log_det chol_km in
      {
        kernel = kernel;
        prepared = prepared;
        km = km;
        chol_km = chol_km;
        log_det_km = log_det_km;
      }

    let calc kernel prepared =
      let km = Spec.Inducing.calc_upper kernel prepared.Prepared.upper in
      calc_internal kernel prepared km

    let get_kernel inducing = inducing.kernel
    let get_points inducing = inducing.prepared.Prepared.points
    let get_prepared inducing = inducing.prepared
    let get_upper inducing = inducing.prepared.Prepared.upper
  end

  (* Evaluation of one input point *)
  module Input = struct
    module Prepared = struct
      type t = {
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

  (* Evaluation of input points *)
  module Inputs = struct
    module Prepared = struct
      type t = {
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

    let create_default_kernel inputs =
      let params = Spec.Inputs.create_default_kernel_params inputs in
      Kernel.create params

    let get_kernel t = t.inducing.Inducing.kernel
    let get_km t = t.inducing.Inducing.km
    let get_chol_km t = t.inducing.Inducing.chol_km
    let get_log_det_km t = t.inducing.Inducing.log_det_km
    let get_kmn t = t.kmn
    let get_points t = t.points

    let nystrom_marginals inputs =
      let kernel, inducing = get_kernel inputs, inputs.inducing in
      let kn_diag = Inputs.calc_diag kernel inputs.points in
      let v_mat = solve_tri ~trans:`T inducing.Inducing.chol_km inputs.kmn in
      v_mat, kn_diag
  end

  (* Model computations shared by standard and variational version *)
  module Common_model = struct
    type t = {
      sigma2 : float;
      inputs : Inputs.t;
      kn_diag : vec;
      v_mat : mat;
      r_vec : vec;
      s_vec : vec;
      is_vec : vec;
      chol_b : mat;
      l1 : float;
    }

    let check_sigma2 sigma2 =
      if sigma2 < 0. then failwith "Model.check_sigma2: sigma2 < 0"

    let calc_internal inputs sigma2 ~kn_diag ~v_mat ~r_vec =
      check_sigma2 sigma2;
      let n = Vec.dim kn_diag in
      let s_vec, is_vec, log_det_s_vec, chol_b, log_det_b =
        let s_vec = Vec.create n in
        let is_vec = Vec.create n in
        let rec loop log_det_s_vec i =
          if i = 0 then log_det_s_vec
          else
            let s_vec_i = r_vec.{i} +. sigma2 in
            s_vec.{i} <- s_vec_i;
            let is_vec_i = 1. /. s_vec_i in
            is_vec.{i} <- is_vec_i;
            loop (log_det_s_vec +. log s_vec_i) (i - 1)
        in
        let log_det_s_vec = loop 0. n in
        let km = Inputs.get_km inputs in
        let kmn_ = lacpy (Inputs.get_kmn inputs) in
        Mat.scal_cols kmn_ (Vec.sqrt is_vec);
        let chol_b = syrk kmn_ ~beta:1. ~c:(lacpy ~uplo:`U km) in
        potrf ~jitter chol_b;
        s_vec, is_vec, log_det_s_vec, chol_b, log_det chol_b
      in
      let l1 =
        let log_det_km = Inputs.get_log_det_km inputs in
        -0.5 *. (log_det_b -. log_det_km +. log_det_s_vec +. float n *. log_2pi)
      in
      {
        inputs = inputs;
        sigma2 = sigma2;
        chol_b = chol_b;
        kn_diag = kn_diag;
        v_mat = v_mat;
        r_vec = r_vec;
        s_vec = s_vec;
        is_vec = is_vec;
        l1 = l1;
      }

    let calc_r_vec ~kn_diag ~v_mat =
      Mat.syrk_diag ~trans:`T ~alpha:(-1.) v_mat ~beta:1. ~y:(copy kn_diag)

    let calc_with_kn_diag inputs sigma2 kn_diag =
      let chol_km = Inputs.get_chol_km inputs in
      let v_mat = solve_tri ~trans:`T chol_km inputs.Inputs.kmn in
      let r_vec = calc_r_vec ~kn_diag ~v_mat in
      calc_internal inputs sigma2 ~kn_diag ~v_mat ~r_vec

    let calc inputs ~sigma2 =
      let v_mat, kn_diag = Inputs.nystrom_marginals inputs in
      let r_vec = calc_r_vec ~kn_diag ~v_mat in
      calc_internal inputs sigma2 ~kn_diag ~v_mat ~r_vec

    let update_sigma2 model sigma2 =
      check_sigma2 sigma2;
      let { kn_diag = kn_diag; v_mat = v_mat; r_vec = r_vec } = model in
      calc_internal model.inputs sigma2 ~kn_diag ~v_mat ~r_vec

    let calc_log_evidence model = model.l1

    let get_chol_b model = model.chol_b
    let get_v_mat model = model.v_mat
    let get_inducing model = model.inputs.Inputs.inducing
    let get_inducing_points model = Inducing.get_points (get_inducing model)
    let get_upper model = Inducing.get_upper (get_inducing model)
    let get_inputs model = model.inputs
    let get_input_points model = model.inputs.Inputs.points
    let get_kernel model = (get_inducing model).Inducing.kernel
    let get_sigma2 model = model.sigma2
    let get_chol_km model = Inputs.get_chol_km model.inputs
    let get_r_vec model = model.r_vec
    let get_s_vec model = model.s_vec
    let get_is_vec model = model.is_vec
    let get_km model = (get_inducing model).Inducing.km
    let get_kmn model = model.inputs.Inputs.kmn
    let get_kn_diag model = model.kn_diag

    let calc_co_variance_coeffs model =
      let res = ichol (get_chol_km model) in
      Mat.axpy ~alpha:(-1.) ~x:(ichol model.chol_b) res;
      potrf ~jitter res;
      res
  end

  (* Model computation (variational version) *)
  module Variational_model = struct
    include Common_model

    let from_common ({ r_vec = r_vec; is_vec = is_vec; l1 = l1 } as model) =
      { model with l1 = l1 +. -0.5 *. dot ~x:is_vec r_vec }

    let calc_internal inputs sigma2 ~kn_diag ~v_mat ~r_vec =
      from_common (calc_internal inputs sigma2 ~kn_diag ~v_mat ~r_vec)

    let calc_with_kn_diag inputs sigma2 kn_diag =
      from_common (calc_with_kn_diag inputs sigma2 kn_diag)

    let update_sigma2 model sigma2 = from_common (update_sigma2 model sigma2)
    let calc inputs ~sigma2 = from_common (calc inputs ~sigma2)
  end

  (* Trained models *)
  module Trained = struct
    type t = {
      model : Common_model.t;
      y : vec;
      coeffs : vec;
      l : float;
    }

    let calc_internal model ~y ~t_vec ~l2 =
      {
        model = model;
        y = y;
        coeffs = t_vec;
        l = model.Common_model.l1 +. l2;
      }

    let calc model ~targets:y =
      let y__ = Vec.mul y model.Common_model.is_vec in
      let kmn_y__ = gemv (Common_model.get_kmn model) y__ in
      trsv ~trans:`T model.Common_model.chol_b kmn_y__;
      let l2 = -0.5 *. (dot ~x:y__ y -. Vec.sqr_nrm2 kmn_y__) in
      let t_vec = kmn_y__ in
      trsv model.Common_model.chol_b t_vec;
      calc_internal model ~y ~t_vec ~l2

    let calc_mean_coeffs trained = trained.coeffs
    let calc_log_evidence trained = trained.l

    let calc_means trained =
      gemv ~trans:`T (Common_model.get_kmn trained.model) trained.coeffs

    let calc_rmse trained =
      let means = calc_means trained in
      sqrt ((Vec.ssqr_diff trained.y means) /. float (Vec.dim trained.y))

    let get_kernel trained = Common_model.get_kernel trained.model
    let get_inducing trained = Common_model.get_inducing trained.model
    let get_upper trained = Common_model.get_upper trained.model
    let get_targets trained = trained.y
    let get_model trained = trained.model
  end

  module Mean_predictor = struct
    type t = {
      kernel : Spec.Kernel.t;
      sigma2 : float;
      inducing : Spec.Inducing.t;
      coeffs : vec;
    }

    let calc_trained trained =
      {
        kernel = Trained.get_kernel trained;
        sigma2 = Common_model.get_sigma2 (Trained.get_model trained);
        inducing = Inducing.get_points (Trained.get_inducing trained);
        coeffs = Trained.calc_mean_coeffs trained;
      }

    let calc kernel ~sigma2 inducing ~coeffs =
      {
        kernel = kernel;
        sigma2 = sigma2;
        inducing = inducing;
        coeffs = coeffs;
      }

    let get_kernel t = t.kernel
    let get_sigma2 t = t.sigma2
    let get_inducing t = t.inducing
    let get_coeffs t = t.coeffs
  end

  (* Prediction of mean for one input point *)
  module Mean = struct
    type t = { point : Spec.Input.t; value : float }

    let calc_induced mean_predictor input =
      if
        mean_predictor.Mean_predictor.inducing
        != Inducing.get_points input.Input.inducing
      then
        failwith
          "Mean.calc_induced: mean predictor and input disagree \
          about inducing points"
      else
        let value =
          dot ~x:input.Input.k_m mean_predictor.Mean_predictor.coeffs
        in
        { point = input.Input.point; value = value }

    let get mean = mean.value
  end

  (* Prediction of means for several input points *)
  module Means = struct
    type t = { points : Spec.Inputs.t; values : vec }

    let calc_induced mean_predictor inputs =
      if
        mean_predictor.Mean_predictor.inducing !=
        Inducing.get_points inputs.Inputs.inducing
      then
        failwith
          "Means.calc_induced: \
          trained and inputs disagree about inducing points"
      else
        let { Inputs.points = points; kmn = kmn } = inputs in
        let values = gemv ~trans:`T kmn mean_predictor.Mean_predictor.coeffs in
        { points = points; values = values }

    let get means = means.values

    module Inducing = struct
      type t = { points : Spec.Inducing.t; values : vec }

      let calc { Trained.coeffs = coeffs; model = model } =
        {
          points = Common_model.get_inducing_points model;
          values = symv (Common_model.get_km model) coeffs;
        }

      let get means = means.values
    end
  end

  module Co_variance_predictor = struct
    type t = {
      kernel : Spec.Kernel.t;
      inducing : Spec.Inducing.t;
      coeffs : mat;
    }

    let calc_model model =
      {
        kernel = Common_model.get_kernel model;
        inducing = Inducing.get_points (Common_model.get_inducing model);
        coeffs = Common_model.calc_co_variance_coeffs model;
      }

    let calc kernel inducing ~coeffs =
      {
        kernel = kernel;
        inducing = inducing;
        coeffs = coeffs;
      }
  end

  (* Prediction of variance for one input point *)
  module Variance = struct
    type t = { point : Spec.Input.t; variance : float; sigma2 : float }

    let calc_induced co_variance_predictor ~sigma2 input =
      if
        co_variance_predictor.Co_variance_predictor.inducing
        != Inducing.get_points input.Input.inducing
      then
        failwith
          "Variance.calc_induced: \
          co-variance predictor and input disagree about inducing points"
      else
        let { Input.point = point; k_m = k_m } = input in
        let kernel = co_variance_predictor.Co_variance_predictor.kernel in
        let prior_variance = Spec.Input.eval_one kernel point in
        let coeffs = co_variance_predictor.Co_variance_predictor.coeffs in
        let explained_variance =
          let tmp = copy k_m in
          trmv coeffs tmp;
          Vec.sqr_nrm2 tmp
        in
        let posterior_variance = prior_variance -. explained_variance in
        {
          point = point;
          variance = posterior_variance;
          sigma2 = sigma2;
        }

    let get ?predictive t =
      match predictive with
      | None | Some true -> t.variance +. t.sigma2
      | Some false -> t.variance
  end

  let solve_chol_b model k = solve_tri ~trans:`T model.Common_model.chol_b k

  (* Prediction of variance for several input points *)
  module Variances = struct
    type t = { points : Spec.Inputs.t; variances : vec; sigma2 : float }

    let make ~points ~variances ~sigma2 =
      { points = points; variances = variances; sigma2 = sigma2 }

    let calc_model_inputs model =
      let r_mat = solve_chol_b model (Common_model.get_kmn model) in
      let variances =
        Mat.syrk_diag ~trans:`T r_mat ~beta:1.
          ~y:(copy (Common_model.get_r_vec model))
      in
      let sigma2 = Common_model.get_sigma2 model in
      make ~points:(Common_model.get_input_points model) ~variances ~sigma2

    let calc_induced co_variance_predictor ~sigma2 inputs =
      if
        co_variance_predictor.Co_variance_predictor.inducing
        != Inducing.get_points inputs.Inputs.inducing
      then
        failwith
          "Variances.calc_induced: \
          co-variance predictor and inputs disagree about inducing points"
      else
        let { Inputs.points = points; kmn = kmt } = inputs in
        let kernel = co_variance_predictor.Co_variance_predictor.kernel in
        let prior_variances = Spec.Inputs.calc_diag kernel points in
        let coeffs = co_variance_predictor.Co_variance_predictor.coeffs in
        let posterior_variances =
          let tmp = lacpy kmt in
          trmm coeffs ~b:tmp;
          Mat.syrk_diag tmp ~alpha:(-1.) ~y:prior_variances
        in
        make ~points ~variances:posterior_variances ~sigma2

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

      let calc model =
        let km = lacpy ~uplo:`U (Common_model.get_km model) in
        Mat.detri km;
        let ichol_b_km = solve_chol_b model km in
        let variances = Mat.syrk_diag ~trans:`T ichol_b_km in
        let points = Common_model.get_inducing_points model in
        let sigma2 = Common_model.get_sigma2 model in
        { points = points; variances = variances; sigma2 = sigma2 }

      let get ?predictive { variances = variances; sigma2 = sigma2 } =
        get_common ?predictive ~variances ~sigma2
    end
  end

  (* Computations for predicting covariances shared by FIC and
     FITC, and standard and variational version *)
  module Common_covariances = struct
    type t = { points : Spec.Inputs.t; covariances : mat; sigma2 : float }

    let check_inducing ~loc co_variance_predictor inputs =
      if
        co_variance_predictor.Co_variance_predictor.inducing
        != Inducing.get_points inputs.Inputs.inducing
      then
        failwith (
          sprintf
            "%s_covariances.calc_induced: \
            co-variance predictor and inputs disagree about inducing points"
            loc)

    let make_with_prior co_variance_predictor sigma2 inputs prior_covariances =
      let coeffs = co_variance_predictor.Co_variance_predictor.coeffs in
      let { Inputs.points = points; kmn = kmt } = inputs in
      let covariances =
        let tmp = lacpy kmt in
        trmm coeffs ~b:tmp;
        syrk ~trans:`T ~alpha:(-1.) tmp ~c:prior_covariances
      in
      { points = points; covariances = covariances; sigma2 = sigma2 }

    let get_common ?predictive ~covariances ~sigma2 =
      match predictive with
      | None | Some true ->
          let res = lacpy ~uplo:`U covariances in
          for i = 1 to Mat.dim1 res do res.{i, i} <- res.{i, i} +. sigma2 done;
          res
      | Some false -> covariances

    let get ?predictive { covariances = covariances; sigma2 = sigma2 } =
      get_common ?predictive ~covariances ~sigma2

    let get_variances { points = points; covariances = covs; sigma2 = sigma2 } =
      {
        Variances.points = points;
        variances = Mat.copy_diag covs;
        sigma2 = sigma2;
      }

    module Inducing = struct
      type t = {
        points : Spec.Inducing.t;
        covariances : mat;
        sigma2 : float;
      }

      let calc model =
        let points = Common_model.get_inducing_points model in
        let km = lacpy ~uplo:`U (Common_model.get_km model) in
        Mat.detri km;
        let ichol_b_km = solve_chol_b model km in
        let covariances = syrk ~trans:`T ichol_b_km in
        let sigma2 = Common_model.get_sigma2 model in
        { points = points; covariances = covariances; sigma2 = sigma2 }

      let get ?predictive { covariances = covariances; sigma2 = sigma2 } =
        get_common ?predictive ~covariances ~sigma2

      let get_variances
            { points = points; covariances = covs; sigma2 = sigma2 } =
        {
          Variances.Inducing.
          points = points;
          variances = Mat.copy_diag covs;
          sigma2 = sigma2;
        }
    end
  end

  (* Predicting covariances with FITC (standard or variational) *)
  module FITC_covariances = struct
    include Common_covariances

    let calc_model_inputs model =
      let kmn = model.Common_model.inputs.Inputs.kmn in
      let v_mat = model.Common_model.v_mat in
      let points = Common_model.get_input_points model in
      let kernel = Common_model.get_kernel model in
      let covariances = Spec.Inputs.calc_upper kernel points in
      ignore (syrk ~trans:`T ~alpha:(-1.) v_mat ~beta:1. ~c:covariances);
      let r_mat = solve_chol_b model kmn in
      ignore (syrk ~trans:`T ~alpha:1. r_mat ~beta:1. ~c:covariances);
      let sigma2 = Common_model.get_sigma2 model in
      { points = points; covariances = covariances; sigma2 = sigma2 }

    let calc_induced co_variance_predictor ~sigma2 inputs =
      check_inducing ~loc:"FITC" co_variance_predictor inputs;
      let kernel = co_variance_predictor.Co_variance_predictor.kernel in
      let prior_covariances =
        Spec.Inputs.calc_upper kernel inputs.Inputs.points
      in
      make_with_prior co_variance_predictor sigma2 inputs prior_covariances
  end

  (* Predicting covariances with FIC (standard or variational) *)
  module FIC_covariances = struct
    include Common_covariances

    let calc_model_inputs model =
      let points = Common_model.get_input_points model in
      let r_vec = model.Common_model.r_vec in
      let kmn = model.Common_model.inputs.Inputs.kmn in
      let r_mat = solve_chol_b model kmn in
      let covariances = syrk ~trans:`T ~alpha:1. r_mat in
      for i = 1 to Vec.dim r_vec do
        covariances.{i, i} <- covariances.{i, i} +. r_vec.{i}
      done;
      let sigma2 = Common_model.get_sigma2 model in
      { points = points; covariances = covariances; sigma2 = sigma2 }

    let calc_induced co_variance_predictor ~sigma2 inputs =
      check_inducing ~loc:"FIC" co_variance_predictor inputs;
      let prior_covariances = 
        let v_mat, kn_diag = Inputs.nystrom_marginals inputs in
        let res = syrk ~trans:`T v_mat in
        for i = 1 to Vec.dim kn_diag do res.{i, i} <- kn_diag.{i} done;
        res
      in
      make_with_prior co_variance_predictor sigma2 inputs prior_covariances
  end

  (* Computations for sampling the marginal posterior GP distribution
     shared by standard and variational version *)
  module Common_sampler = struct
    type t = { mean : float; stddev : float }

    let calc ~loc ?predictive mean variance =
      if mean.Mean.point != variance.Variance.point then
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

  (* Computations for sampling the posterior GP distribution shared
     by FIC and FITC, and standard and variational version *)
  module Common_cov_sampler = struct
    type t = { means : vec; cov_chol : mat }

    let calc ~loc ?predictive means covariances =
      let module Covariances = Common_covariances in
      if means.Means.points != covariances.Covariances.points then
        failwith (
          loc ^
          ".Cov_sampler: means and covariances disagree about input points");
      let cov_chol = lacpy ~uplo:`U covariances.Covariances.covariances in
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
          samples.{row, col} <- samples.{row, col} +. means.{row}
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


(* Handling derivatives *)

module type Deriv_sig = functor (Spec : Specs.Deriv) ->
  Sigs.Deriv
    with module Eval.Spec = Spec.Eval
    with module Deriv.Spec = Spec

(* Computations shared by FIC and FITC, and standard and variational
   version for derivatives *)
module Make_common_deriv (Spec : Specs.Deriv) = struct
  open Spec

  (* Eval modules *)

  module Eval_common = Make_common (Spec.Eval)
  module Eval_inducing = Eval_common.Inducing
  module Eval_inputs = Eval_common.Inputs
  module Eval_model = Eval_common.Common_model
  module Eval_trained = Eval_common.Trained

  (* Kind of model *)
  type model_kind = Standard | Variational

  module Deriv_common = struct

    (* Derivative modules *)

    module Spec = Spec

    (* Derivative of inducing points *)
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

      type t = { eval : Eval_inducing.t; shared_upper : Spec.Inducing.upper }

      let calc kernel { Prepared.eval = eval; upper = upper } =
        let km, shared_upper = Spec.Inducing.calc_shared_upper kernel upper in
        {
          eval = Eval_inducing.calc_internal kernel eval km;
          shared_upper = shared_upper;
        }

      let calc_eval inducing = inducing.eval
      let get_kernel inducing = Eval_inducing.get_kernel inducing.eval
      let get_points inducing = Eval_inducing.get_points inducing.eval
    end

    (* Derivative of inputs *)
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
        { inducing = inducing; eval = eval; shared_cross = shared_cross }

      let calc_eval t = t.eval

      let get_kernel inputs = Inducing.get_kernel inputs.inducing
    end

    (* Derivative of hyper parameters *)
    module Shared = struct
      type shared = {
        km : mat;
        kmn : mat;
        kn_diag : vec;
        shared_upper : Spec.Inducing.upper;
        shared_cross : Spec.Inputs.cross;
        shared_diag : Spec.Inputs.diag;
      }

      type dfacts = { v_vec : vec; w_mat : mat; x_mat : mat }
      type hyper_t = { shared : shared; dfacts : dfacts }

      let calc_u_mat eval_model =
        let v_mat = Eval_model.get_v_mat eval_model in
        solve_tri (Eval_model.get_chol_km eval_model) v_mat

      let calc_dkn_diag_term ~v_vec ~kn_diag = function
        | `Vec dkn_diag -> dot ~x:v_vec dkn_diag
        | `Sparse_vec (svec, rows) ->
            check_sparse_vec_sane ~real_n:(Vec.dim v_vec) ~svec ~rows;
            let res_ref = ref 0. in
            for i = 1 to Vec.dim svec do
              res_ref := !res_ref +. v_vec.{rows.{i}} *. svec.{i}
            done;
            !res_ref
        | `Const 0. | `Factor 0. -> 0.
        | `Const c -> c *. Vec.sum v_vec
        | `Factor c -> c *. dot ~x:kn_diag v_vec

      let calc_dkm_term ~w_mat ~km = function
        | `Dense dkm -> Mat.symm2_trace w_mat dkm
        | `Sparse_rows (smat, rows) -> symm2_sparse_trace ~mat:w_mat ~smat ~rows
        | `Const 0. | `Factor 0. | `Diag_const 0. -> 0.
        | `Const c -> c *. sum_symm_mat w_mat
        | `Factor c -> c *. Mat.symm2_trace w_mat km
        | `Diag_vec ddkm ->
            let res_ref = ref 0. in
            for i = 1 to Mat.dim1 w_mat do
              res_ref := !res_ref +. ddkm.{i} *. w_mat.{i, i}
            done;
            !res_ref
        | `Diag_const c ->
            let res_ref = ref 0. in
            for i = 1 to Mat.dim1 w_mat do
              res_ref := !res_ref +. c *. w_mat.{i, i}
            done;
            !res_ref

      let calc_dkmn_term ~x_mat ~kmn = function
        | `Dense dkmn -> Mat.gemm_trace ~transa:`T x_mat dkmn
        | `Sparse_rows (sdkmn, rows) ->
            let real_m = Mat.dim1 x_mat in
            check_sparse_row_mat_sane ~real_m ~smat:sdkmn ~rows;
            let m = Int_vec.dim rows in
            let n = Mat.dim2 sdkmn in
            let res_ref = ref 0. in
            for r = 1 to m do
              let real_r = rows.{r} in
              for c = 1 to n do
                res_ref := !res_ref +. x_mat.{real_r, c} *. sdkmn.{r, c}
              done
            done;
            !res_ref
        | `Const 0. | `Factor 0. -> 0.
        | `Const c -> c *. sum_mat x_mat
        | `Factor c -> c *. Mat.gemm_trace ~transa:`T x_mat kmn
        | `Sparse_cols (sdkmn, cols) ->
            let real_n = Mat.dim2 x_mat in
            check_sparse_col_mat_sane ~real_n ~smat:sdkmn ~cols;
            let m = Mat.dim1 sdkmn in
            let n = Int_vec.dim cols in
            let res_ref = ref 0. in
            for c = 1 to n do
              let real_c = cols.{c} in
              for r = 1 to m do
                res_ref := !res_ref +. x_mat.{r, real_c} *. sdkmn.{r, c}
              done
            done;
            !res_ref

      let calc_log_evidence { shared = shared; dfacts = dfacts } hyper =
        let { v_vec = v_vec; w_mat = w_mat; x_mat = x_mat } = dfacts in
        let dkn_diag_term =
          let kn_diag = shared.kn_diag in
          let dkn_diag = Spec.Inputs.calc_deriv_diag shared.shared_diag hyper in
          calc_dkn_diag_term ~v_vec ~kn_diag dkn_diag
        in
        let dkm_term =
          let km = shared.km in
          let dkm = Spec.Inducing.calc_deriv_upper shared.shared_upper hyper in
          calc_dkm_term ~w_mat ~km dkm
        in
        let dkmn_term =
          let kmn = shared.kmn in
          let dkmn = Spec.Inputs.calc_deriv_cross shared.shared_cross hyper in
          calc_dkmn_term ~x_mat ~kmn dkmn
        in
        -0.5 *. (dkn_diag_term -. dkm_term) -. dkmn_term
    end

    (* Derivative of models *)
    module Common_model = struct

      (* Precomputations for all derivatives *)

      type t = {
        model_kind : model_kind;
        model_shared : Shared.shared;
        eval_model : Eval_model.t;
        inv_km : mat;
        r_diag : vec;
        s_mat : mat;
        t_mat : mat;
      }

      let calc_internal model_kind model_shared eval_model inv_km =
        let kmn = Eval_model.get_kmn eval_model in
        let r_mat = Eval_common.solve_chol_b eval_model kmn in
        let r_diag = Mat.syrk_diag ~trans:`T r_mat in
        let s_mat = r_mat in
        let chol_b = Eval_model.get_chol_b eval_model in
        trtrs chol_b s_mat;
        Mat.scal_cols s_mat eval_model.Eval_model.is_vec;
        let t_mat = lacpy ~uplo:`U inv_km in
        Mat.axpy ~alpha:(-1.) ~x:(ichol chol_b) t_mat;
        {
          model_kind = model_kind;
          model_shared = model_shared;
          eval_model = eval_model;
          inv_km = inv_km;
          r_diag = r_diag;
          s_mat = s_mat;
          t_mat = t_mat;
        }

      let calc_common model_kind inputs sigma2 =
        let kernel = Inputs.get_kernel inputs in
        let eval_inputs = inputs.Inputs.eval in
        let kn_diag, shared_diag =
          Spec.Inputs.calc_shared_diag kernel eval_inputs.Eval_inputs.points
        in
        let calc_with_kn_diag =
          match model_kind with
          | Standard -> Eval_model.calc_with_kn_diag
          | Variational -> Eval_common.Variational_model.calc_with_kn_diag
        in
        let eval_model = calc_with_kn_diag eval_inputs sigma2 kn_diag in
        let km = Eval_model.get_km eval_model in
        let kmn = Eval_model.get_kmn eval_model in
        let kn_diag = Eval_model.get_kn_diag eval_model in
        let chol_km = Eval_model.get_chol_km eval_model in
        let inv_km = ichol chol_km in
        let model_shared =
          {
            Shared.
            km = km;
            kmn = kmn;
            kn_diag = kn_diag;
            shared_upper = inputs.Inputs.inducing.Inducing.shared_upper;
            shared_cross = inputs.Inputs.shared_cross;
            shared_diag = shared_diag;
          }
        in
        calc_internal model_kind model_shared eval_model inv_km

      let calc_eval model = model.eval_model
      let calc inputs ~sigma2 = calc_common Standard inputs sigma2

      let update_sigma2 ({ model_kind = model_kind } as model) sigma2 =
        let update_sigma2 =
          match model_kind with
          | Standard -> Eval_model.update_sigma2
          | Variational -> Eval_common.Variational_model.update_sigma2
        in
        let eval_model = update_sigma2 model.eval_model sigma2 in
        calc_internal model_kind model.model_shared eval_model model.inv_km

      let calc_vc_vec model =
        let s_vec = Eval_model.get_s_vec model.eval_model in
        let vc_vec =
          match model.model_kind with
          | Standard -> copy s_vec
          | Variational ->
              let n = Vec.dim s_vec in
              let vc_vec = Vec.create n in
              let r_vec = Eval_model.get_r_vec model.eval_model in
              for i = 1 to n do
                let s_vec_i = s_vec.{i} in
                vc_vec.{i} <- s_vec_i +. (s_vec_i -. r_vec.{i})
              done;
              vc_vec
        in
        axpy ~alpha:(-1.) ~x:model.r_diag vc_vec;
        vc_vec

      let calc_v1_vec model =
        let z = calc_vc_vec model in
        let is_vec = Eval_model.get_is_vec model.eval_model in
        Vec.mul is_vec z ~z:(Vec.mul is_vec z ~z)

      (* Derivative of sigma2 *)

      let common_calc_log_evidence_sigma2 ({ eval_model = em } as model) v_vec =
        let sum_v1_vec = Vec.sum v_vec in
        let sum =
          match model.model_kind with
          | Standard -> sum_v1_vec
          | Variational -> sum_v1_vec -. Vec.sum em.Eval_model.is_vec
        in
        -0.5 *. sum

      let calc_log_evidence_sigma2 model =
        common_calc_log_evidence_sigma2 model (calc_v1_vec model)

      (* Prepare derivative of general hyper-parameters *)

      let prepare_hyper ({ eval_model = eval_model } as model) =
        let v_vec = calc_v1_vec model in
        let u_mat = Shared.calc_u_mat eval_model in
        let w_mat =
          let u1_mat = lacpy u_mat in
          Mat.scal_cols u1_mat (Vec.sqrt v_vec);
          syrk ~alpha:(-1.) u1_mat ~beta:1. ~c:(lacpy ~uplo:`U model.t_mat)
        in
        let x_mat = u_mat in
        Mat.scal_cols x_mat (Vec.neg v_vec);
        Mat.axpy ~x:model.s_mat x_mat;
        {
          Shared.
          shared = model.model_shared;
          dfacts = { Shared.v_vec = v_vec; w_mat = w_mat; x_mat = x_mat };
        }

      include Shared
    end

    module Cm = Common_model

    module Variational_model = struct
      include Cm

      let calc inputs ~sigma2 = calc_common Variational inputs sigma2
    end

    (* Derivative of trained models *)
    module Trained = struct
      type t = {
        common_model : Cm.t;
        eval_trained : Eval_trained.t;
        t_vec : vec;
        u_vec : vec;
        v_vec : vec;
      }

      let calc common_model ~targets:y =
        let t_vec = gemv common_model.Cm.s_mat y in
        let eval_model = common_model.Cm.eval_model in
        let kmn = Eval_model.get_kmn eval_model in
        let e_vec = gemv ~alpha:(-1.) ~beta:1. ~trans:`T kmn t_vec ~y:(copy y) in
        let sqr_e_vec = Vec.sqr e_vec in
        let is_vec = Eval_model.get_is_vec eval_model in
        let u_vec = Vec.mul is_vec e_vec ~z:e_vec in
        let eval_trained =
          let l2 = -0.5 *. dot ~x:u_vec y in
          Eval_trained.calc_internal eval_model ~y ~t_vec ~l2
        in
        let v_vec =
          let z = Cm.calc_vc_vec common_model in
          axpy ~alpha:(-1.) ~x:sqr_e_vec z;
          Vec.mul is_vec z ~z:(Vec.mul is_vec z ~z)
        in
        {
          common_model = common_model;
          eval_trained = eval_trained;
          t_vec = t_vec;
          u_vec = u_vec;
          v_vec = v_vec;
        }

      let calc_eval trained = trained.eval_trained

      (* Derivative of sigma2 *)

      let calc_log_evidence_sigma2 { common_model = cm; v_vec = v_vec } =
        Cm.common_calc_log_evidence_sigma2 cm v_vec

      (* Derivative of general hyper-parameters *)

      let prepare_hyper trained =
        let common_model = trained.common_model in
        let eval_model = common_model.Cm.eval_model in
        let u_mat = Shared.calc_u_mat eval_model in
        let t_vec = trained.t_vec in
        let u_vec = trained.u_vec in
        let w_mat =
          let w_mat =
            let t_mat = trained.common_model.Cm.t_mat in
            syr ~alpha:(-1.) t_vec (lacpy ~uplo:`U t_mat)
          in
          let u1_mat = lacpy u_mat in
          let sqrt_v1_vec = Cm.calc_v1_vec common_model in
          Mat.scal_cols u1_mat (Vec.sqrt ~y:sqrt_v1_vec sqrt_v1_vec);
          let w_mat = syrk ~alpha:(-1.) u1_mat ~beta:1. ~c:w_mat in
          let u2_mat = u1_mat in
          let u2_mat = lacpy u_mat ~b:u2_mat in
          Mat.scal_cols u2_mat u_vec;
          syrk u2_mat ~beta:1. ~c:w_mat
        in
        let v_vec = trained.v_vec in
        let x_mat =
          Mat.scal_cols u_mat v_vec;
          let x_mat = lacpy common_model.Cm.s_mat in
          Mat.axpy ~alpha:(-1.) ~x:u_mat x_mat;
          ger ~alpha:(-1.) t_vec u_vec x_mat
        in
        {
          Shared.
          shared = common_model.Cm.model_shared;
          dfacts = { Shared.v_vec = v_vec; w_mat = w_mat; x_mat = x_mat };
        }

      include Shared
    end

    module Test = struct
      let check_deriv_hyper kernel1 inputs_prepared hyper ~eps ~tol =
        let kernel2 =
          let hypers, values = Spec.Hyper.extract kernel1 in
          let rec loop i =
            if i < 0 then
              failwith
                "Gpr.Fitc_gp.Make_deriv.Test.check_deriv_hyper: \
                hyper variable is unknown"
            else
              if hypers.(i) = hyper then
                let i1 = i + 1 in
                let value = values.{i1} in
                let value_eps = value +. eps in
                let values_eps = copy values in
                values_eps.{i1} <- value_eps;
                Spec.Hyper.update kernel1 values_eps
              else loop (i - 1)
          in
          loop (Array.length hypers - 1)
        in
        let
          {
            Inputs.Prepared.
            inducing_prepared =
              ({ Inducing.Prepared.eval = eval_inducing_prepared }
              as inducing_prepared);
            eval = eval_inputs_prepared;
          } = inputs_prepared
        in
        let eval_inducing1 = Eval_inducing.calc kernel1 eval_inducing_prepared in
        let eval_cross1 = Eval_inputs.calc eval_inducing1 eval_inputs_prepared in
        let eval_inducing2 = Eval_inducing.calc kernel2 eval_inducing_prepared in
        let eval_cross2 = Eval_inputs.calc eval_inducing2 eval_inputs_prepared in
        let kmn1 = eval_cross1.Eval_inputs.kmn in
        let finite_diff_mat =
          let res = lacpy eval_cross2.Eval_inputs.kmn in
          Mat.axpy ~alpha:(-1.) ~x:kmn1 res;
          Mat.scal (1. /. eps) res;
          res
        in
        let inducing = Inducing.calc kernel1 inducing_prepared in
        let inputs = Inputs.calc inducing inputs_prepared in
        let check ~deriv ~r ~c =
          let finite = finite_diff_mat.{r, c} in
          if abs_float (finite -. deriv) > tol then
            failwith (
              sprintf
                "Gpr.Fitc_gp.Make_deriv.Test.check_deriv_hyper: \
                finite difference (%f) and derivative (%f) differ \
                by more than %f on dkmn.{%d, %d}" finite deriv tol r c)
        in
        match Spec.Inputs.calc_deriv_cross inputs.Inputs.shared_cross hyper with
        | `Dense deriv_mat ->
            let m = Mat.dim1 kmn1 in
            let n = Mat.dim2 kmn1 in
            for c = 1 to n do
              for r = 1 to m do check ~deriv:deriv_mat.{r, c} ~r ~c done
            done
        | `Sparse_rows (sdkmn, rows) ->
            let m = Int_vec.dim rows in
            let n = Mat.dim2 sdkmn in
            for r = 1 to m do
              let real_r = rows.{r} in
              for c = 1 to n do check ~deriv:sdkmn.{r, c} ~r:real_r ~c done
            done
        | `Const const ->
            let m = Mat.dim1 kmn1 in
            let n = Mat.dim2 kmn1 in
            for c = 1 to n do
              for r = 1 to m do check ~deriv:const ~r ~c done
            done
        | `Factor const ->
            let m = Mat.dim1 kmn1 in
            let n = Mat.dim2 kmn1 in
            for c = 1 to n do
              for r = 1 to m do check ~deriv:(const *. kmn1.{r, c}) ~r ~c done
            done
        | `Sparse_cols (sdkmn, cols) ->
            let m = Mat.dim1 sdkmn in
            let n = Int_vec.dim cols in
            for r = 1 to m do
              for c = 1 to n do
                let real_c = cols.{c} in
                check ~deriv:sdkmn.{r, c} ~r ~c:real_c done
            done
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
