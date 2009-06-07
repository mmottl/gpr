open Format

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

    let get_kernel t = t.inducing.Inducing.kernel
    let get_km t = t.inducing.Inducing.km
    let get_chol_km t = t.inducing.Inducing.chol_km
    let get_log_det_km t = t.inducing.Inducing.log_det_km
    let get_kmn t = t.kmn

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
        Mat.scal_cols (Vec.sqrt is_vec) kmn_;
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
    let get_input_points model = model.inputs.Inputs.points
    let get_kernel model = (get_inducing model).Inducing.kernel
    let get_sigma2 model = model.sigma2
    let get_kmn model = model.inputs.Inputs.kmn
    let get_chol_km model = Inputs.get_chol_km model.inputs
    let get_r_vec model = model.r_vec
    let get_s_vec model = model.s_vec
    let get_is_vec model = model.is_vec
    let get_km model = (get_inducing model).Inducing.km
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

    let get_coeffs trained = trained.coeffs
    let calc_log_evidence trained = trained.l

    let get_kernel trained = Common_model.get_kernel trained.model
    let get_inducing trained = Common_model.get_inducing trained.model
    let get_upper trained = Common_model.get_upper trained.model
  end

  (* Prediction of mean for one input point *)
  module Mean = struct
    type t = { point : Spec.Input.t; value : float }

    let make ~point ~value = { point = point; value = value }

    let calc_input trained point =
      let upper = Trained.get_upper trained in
      let prepared = Spec.Input.Prepared.calc_cross upper point in
      let kernel = Trained.get_kernel trained in
      let coeffs = Trained.get_coeffs trained in
      let value = Spec.Input.weighted_eval kernel ~coeffs prepared in
      make ~point ~value

    let calc_induced trained input =
      if Trained.get_inducing trained <> input.Input.inducing then
        failwith
          "Mean.calc_induced: trained and input disagree about inducing points";
      let value = dot ~x:input.Input.k_m trained.Trained.coeffs in
      make ~point:input.Input.point ~value

    let get mean = mean.value
  end

  (* Prediction of means for several input points *)
  module Means = struct
    type t = { points : Spec.Inputs.t; values : vec }

    let make ~points ~values = { points = points; values = values }

    let calc_model_inputs { Trained.coeffs = coeffs; model = model } =
      make
        ~points:(Common_model.get_input_points model)
        ~values:(gemv ~trans:`T (Common_model.get_kmn model) coeffs)

    let calc_inputs trained points =
      let upper = Trained.get_upper trained in
      let prepared = Spec.Inputs.Prepared.calc_cross upper points in
      let kernel = Trained.get_kernel trained in
      let coeffs = Trained.get_coeffs trained in
      let values = Spec.Inputs.weighted_eval kernel ~coeffs prepared in
      make ~points ~values

    let calc_induced trained inputs =
      let { Inputs.points = points; kmn = kmn } = inputs in
      if Trained.get_inducing trained <> inputs.Inputs.inducing then
        failwith
          "Means.calc_induced: \
          trained and inputs disagree about inducing points";
      make ~points ~values:(gemv ~trans:`T kmn trained.Trained.coeffs)

    let get means = means.values

    module Inducing = struct
      type t = { points : Spec.Inducing.t; values : vec }

      let make ~points ~values = { points = points; values = values }

      let calc { Trained.coeffs = coeffs; model = model } =
        make
          ~points:(Common_model.get_inducing_points model)
          ~values:(symv (Common_model.get_km model) coeffs)

      let get means = means.values
    end
  end

  (* Prediction of variance for one input point *)
  module Variance = struct
    type t = { point : Spec.Input.t; variance : float; sigma2 : float }

    let calc_induced model induced =
      let { Input.point = point; k_m = k_m } = induced in
      let kernel = Common_model.get_kernel model in
      let prior_variance = Spec.Input.eval_one kernel point in
      let ichol_km_k_m = copy k_m in
      let ichol_km_k_m_mat = Mat.from_col_vec ichol_km_k_m in
      let inducing = induced.Input.inducing in
      potrs ~factorize:false inducing.Inducing.chol_km ichol_km_k_m_mat;
      let km_arg = dot ~x:k_m ichol_km_k_m in
      let ichol_b_k_m = copy k_m ~y:ichol_km_k_m in
      let ichol_b_k_m_mat = ichol_km_k_m_mat in
      potrs ~factorize:false model.Common_model.chol_b ichol_b_k_m_mat;
      let b_arg = dot ~x:k_m ichol_b_k_m in
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

  let solve_chol_b model k = solve_tri ~trans:`T model.Common_model.chol_b k

  (* Prediction of variance for several input points *)
  module Variances = struct
    type t = { points : Spec.Inputs.t; variances : vec; sigma2 : float }

    let make ~points ~variances ~model =
      let sigma2 = Common_model.get_sigma2 model in
      { points = points; variances = variances; sigma2 = sigma2 }

    let calc_model_inputs model =
      let r_mat = solve_chol_b model (Common_model.get_kmn model) in
      let y = copy (Common_model.get_r_vec model) in
      let variances = Mat.syrk_diag ~trans:`T r_mat ~beta:1. ~y in
      make ~points:(Common_model.get_input_points model) ~variances ~model

    let calc_induced model inputs =
      if Common_model.get_inducing model <> inputs.Inputs.inducing then
        failwith
          "Variances.calc_induced: \
          model and inputs disagree about inducing points";
      let kmt = inputs.Inputs.kmn in
      let ichol_km_kmt, kt_diag = Inputs.nystrom_marginals inputs in
      let ichol_b_kmt = solve_chol_b model kmt in
      let variances =
        let y =
          Mat.syrk_diag ~trans:`T ichol_b_kmt
            ~beta:1. ~y:(copy kt_diag)
        in
        Mat.syrk_diag ~alpha:(-1.) ~trans:`T ichol_km_kmt ~beta:1. ~y
      in
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
        let km = lacpy ~uplo:`U (Common_model.get_km model) in
        Mat.detri km;
        let ichol_b_km = solve_chol_b model km in
        let variances = Mat.syrk_diag ~trans:`T ichol_b_km in
        make ~points:(Common_model.get_inducing_points model) ~variances ~model

      let get ?predictive { variances = variances; sigma2 = sigma2 } =
        get_common ?predictive ~variances ~sigma2
    end
  end

  (* Computations for predicting covariances shared by FIC and
     FITC, and standard and variational version *)
  module Common_covariances = struct
    type t = { points : Spec.Inputs.t; covariances : mat; sigma2 : float }

    let make ~points ~covariances ~model =
      let sigma2 = Common_model.get_sigma2 model in
      { points = points; covariances = covariances; sigma2 = sigma2 }

    let make_b_only ~points ~ichol_b_k ~model =
      make ~points ~covariances:(syrk ~trans:`T ichol_b_k) ~model

    let get_common ?predictive ~covariances ~sigma2 =
      match predictive with
      | None | Some true ->
          let res = lacpy ~uplo:`U covariances in
          for i = 1 to Mat.dim1 res do res.{i, i} <- res.{i, i} +. sigma2 done;
          res
      | Some false -> covariances

    let get ?predictive { covariances = covariances; sigma2 = sigma2 } =
      get_common ?predictive ~covariances ~sigma2

    let variances { points = points; covariances = covs; sigma2 = sigma2 } =
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

      let variances { points = points; covariances = covs; sigma2 = sigma2 } =
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

    let calc_common ~kmn ~v_mat ~points ~model =
      let kernel = Common_model.get_kernel model in
      let covariances = Spec.Inputs.calc_upper kernel points in
      ignore (syrk ~trans:`T ~alpha:(-1.) v_mat ~beta:1. ~c:covariances);
      let r_mat = solve_chol_b model kmn in
      ignore (syrk ~trans:`T ~alpha:1. r_mat ~beta:1. ~c:covariances);
      make ~points ~covariances ~model

    let calc_model_inputs model =
      let kmn = model.Common_model.inputs.Inputs.kmn in
      let v_mat = model.Common_model.v_mat in
      let points = Common_model.get_input_points model in
      calc_common ~kmn ~v_mat ~points ~model

    let calc_induced model inputs =
      if Common_model.get_inducing model <> inputs.Inputs.inducing then
        failwith (
          "FITC_covariances.calc_induced: \
          model and inputs disagree about inducing points");
      let kmn = inputs.Inputs.kmn in
      let v_mat, _ = Inputs.nystrom_marginals inputs in
      let points = inputs.Inputs.points in
      calc_common ~kmn ~v_mat ~points ~model
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
      make ~points ~covariances ~model

    let calc_induced model inputs =
      if Common_model.get_inducing model <> inputs.Inputs.inducing then
        failwith (
          "FIC_covariances.calc_induced: \
          model and inputs disagree about inducing points");
      let kmt = inputs.Inputs.kmn in
      let points = inputs.Inputs.points in
      make_b_only ~points ~ichol_b_k:(solve_chol_b model kmt) ~model
  end

  (* Computations for sampling the marginal posterior GP distribution
     shared by standard and variational version *)
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

  (* Computations for sampling the posterior GP distribution shared
     by FIC and FITC, and standard and variational version *)
  module Common_cov_sampler = struct
    type t = { means : vec; cov_chol : mat }

    let calc ~loc ?predictive means covariances =
      let module Covariances = Common_covariances in
      if means.Means.points <> covariances.Covariances.points then
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

      type t = { eval : Eval_inducing.t; shared : Spec.Inducing.shared }

      let calc kernel { Prepared.eval = eval; upper = upper } =
        let km, shared = Spec.Inducing.calc_shared_upper kernel upper in
        { eval = Eval_inducing.calc_internal kernel eval km; shared = shared }

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
        { inducing = inducing; eval = eval; shared_cross = shared_cross }

      let calc_eval t = t.eval

      let get_kernel inputs = Inducing.get_kernel inputs.inducing
    end

    type model_kind = Standard | Variational

    let calc_u_mat eval_model =
      let v_mat = Eval_model.get_v_mat eval_model in
      solve_tri (Eval_model.get_chol_km eval_model) v_mat

    let calc_dkn_term ~v_vec = function
      | `Vec dkn -> dot ~x:v_vec dkn
      | `Sparse_vec (_sdkn, _inds) -> (assert false (* XXX *))
      | `Const c -> c *. Vec.sum v_vec
      | `Factor _c -> (assert false (* XXX *))

    let calc_dkm_term ~w_mat = function
      | `Dense dkm -> Mat.symm2_trace w_mat dkm
      | `Sparse_rows (_sdkm, _inds) -> (assert false (* XXX *))
      | `Const _c -> (assert false (* XXX *))
      | `Factor _c -> (assert false (* XXX *))
      | `Diag_vec _ddkm -> (assert false (* XXX *))
      | `Diag_const _c -> (assert false (* XXX *))

    let calc_dkmn_term ~x_mat = function
      | `Dense dkm -> Mat.gemm_trace ~transa:`T x_mat dkm
      | `Sparse_rows (_sdkm, _inds) -> (assert false (* XXX *))
      | `Const _c -> (assert false (* XXX *))
      | `Factor _c -> (assert false (* XXX *))
      | `Sparse_cols (_sdkmn, _inds) -> (assert false (* XXX *))

    module Shared = struct
      type t = {
        inducing : Spec.Inducing.shared;
        cross : Spec.Inputs.cross;
        diag : Spec.Inputs.diag;
      }
    end

    module Dfacts = struct
      type t = { v_vec : vec; w_mat : mat; x_mat : mat }
    end

    let common_calc_log_evidence shared dfacts hyper =
      let { Dfacts.v_vec = v_vec; w_mat = w_mat; x_mat = x_mat } = dfacts in
      let dkn_term =
        let dkn = Spec.Inputs.calc_deriv_diag shared.Shared.diag hyper in
        calc_dkn_term ~v_vec dkn
      in
      let dkm_term =
        let dkm = Spec.Inducing.calc_deriv_upper shared.Shared.inducing hyper in
        calc_dkm_term ~w_mat dkm
      in
      let dkmn_term =
        let dkmn = Spec.Inputs.calc_deriv_cross shared.Shared.cross hyper in
        calc_dkmn_term ~x_mat dkmn
      in
      -0.5 *. (dkn_term -. dkm_term) -. dkmn_term

    module Common_model = struct

      (* Precomputations for all derivatives *)

      type t = {
        model_kind : model_kind;
        shared : Shared.t;
        eval_model : Eval_model.t;
        inv_km : mat;
        r_diag : vec;
        s_mat : mat;
        t_mat : mat;
      }

      let calc_internal model_kind shared eval_model inv_km =
        let kmn = Eval_model.get_kmn eval_model in
        let r_mat = Eval_common.solve_chol_b eval_model kmn in
        let r_diag = Mat.syrk_diag ~trans:`T r_mat in
        let s_mat = r_mat in
        let chol_b = Eval_model.get_chol_b eval_model in
        trtrs chol_b s_mat;
        Mat.scal_cols eval_model.Eval_model.is_vec s_mat;
        let t_mat = lacpy ~uplo:`U inv_km in
        Mat.axpy ~alpha:(-1.) ~x:(ichol chol_b) t_mat;
        {
          model_kind = model_kind;
          shared = shared;
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
        let eval_model = calc_with_kn_diag inputs.Inputs.eval sigma2 kn_diag in
        let chol_km = Eval_model.get_chol_km eval_model in
        let inv_km = ichol chol_km in
        let shared =
          {
            Shared.
            inducing = inputs.Inputs.inducing.Inducing.shared;
            cross = inputs.Inputs.shared_cross;
            diag = shared_diag;
          }
        in
        calc_internal model_kind shared eval_model inv_km

      let calc_eval model = model.eval_model
      let calc inputs ~sigma2 = calc_common Standard inputs sigma2

      let update_sigma2 ({ model_kind = model_kind } as model) sigma2 =
        let update_sigma2 =
          match model_kind with
          | Standard -> Eval_model.update_sigma2
          | Variational -> Eval_common.Variational_model.update_sigma2
        in
        let eval_model = update_sigma2 model.eval_model sigma2 in
        calc_internal model_kind model.shared eval_model model.inv_km

      let calc_common_v_vec model =
        let s_vec = Eval_model.get_s_vec model.eval_model in
        let common_v_vec =
          match model.model_kind with
          | Standard -> copy s_vec
          | Variational ->
              let n = Vec.dim s_vec in
              let common_v_vec = Vec.create n in
              let r_vec = Eval_model.get_r_vec model.eval_model in
              for i = 1 to n do
                let s_vec_i = s_vec.{i} in
                common_v_vec.{i} <- s_vec_i +. (s_vec_i -. r_vec.{i})
              done;
              common_v_vec
        in
        let r_diag = model.r_diag in
        axpy ~alpha:(-1.) ~x:r_diag common_v_vec;
        common_v_vec

      let calc_v1_vec model =
        let z = calc_common_v_vec model in
        let is_vec = Eval_model.get_is_vec model.eval_model in
        (* TODO: add more efficient in-place Vec multiplication to LACAML *)
        Vec.mul is_vec z ~z:(Vec.mul is_vec z ~z)

      (* Derivative of sigma2 *)

      let common_calc_log_evidence_sigma2 ({ eval_model = em } as model) v_vec =
        let sum_v1_vec = Vec.sum v_vec in
        match model.model_kind with
        | Standard -> sum_v1_vec
        | Variational -> sum_v1_vec -. Vec.sum em.Eval_model.is_vec

      let calc_log_evidence_sigma2 model =
        -0.5 *. common_calc_log_evidence_sigma2 model (calc_v1_vec model)

      (* Derivative of general hyper-parameters *)

      type hyper_t = { hyp_shared : Shared.t; dfacts : Dfacts.t }

      let prepare_hyper ({ eval_model = eval_model } as model) =
        let v_vec = calc_v1_vec model in
        let u_mat = calc_u_mat eval_model in
        let w_mat =
          let u1_mat = lacpy u_mat in
          Mat.scal_cols (Vec.sqrt v_vec) u1_mat;
          syrk ~alpha:(-1.) u1_mat ~beta:1. ~c:(lacpy ~uplo:`U model.t_mat)
        in
        let x_mat = u_mat in
        Mat.scal_cols v_vec x_mat;
        let m = Mat.dim1 x_mat in
        let n = Mat.dim2 x_mat in
        (* TODO: add inplace addition/subtraction of matrices to LACAML *)
        let s_mat = model.s_mat in
        for c = 1 to n do
          for r = 1 to m do x_mat.{r, c} <- s_mat.{r, c} -. x_mat.{r, c} done
        done;
        {
          hyp_shared = model.shared;
          dfacts = { Dfacts.  v_vec = v_vec; w_mat = w_mat; x_mat = x_mat };
        }

      let calc_log_evidence hyper_model hyper =
        common_calc_log_evidence hyper_model.hyp_shared hyper_model.dfacts hyper
    end

    module Cm = Common_model

    module Variational_model = struct
      include Cm

      let calc inputs ~sigma2 = calc_common Variational inputs sigma2
    end

    (**)

    module Trained = struct
      type t = {
        common_model : Cm.t;
        eval_trained : Eval_trained.t;
        t_vec : vec;
        u_vec : vec;
        v_vec : vec;
        y_minus_kmn_t : vec;
      }

      let calc common_model ~targets:y =
        let t_vec = gemv common_model.Cm.s_mat y in
        let eval_model = common_model.Cm.eval_model in
        let kmn = Eval_model.get_kmn eval_model in
        let y_minus_kmn_t =
          gemv ~alpha:(-1.) ~beta:1. ~trans:`T kmn t_vec ~y:(copy y)
        in
        let is_vec = Eval_model.get_is_vec eval_model in
        let u_vec = Vec.mul is_vec y_minus_kmn_t in
        let eval_trained =
          let l2 = -0.5 *. dot ~x:u_vec y in
          Eval_trained.calc_internal eval_model ~y ~t_vec ~l2
        in
        let v_vec =
          let z = Cm.calc_common_v_vec common_model in
          let sqrt_tmp = Vec.sqr y_minus_kmn_t in
          axpy ~alpha:(-1.) ~x:sqrt_tmp z;
          (* TODO: add more efficient in-place Vec multiplication to LACAML *)
          Vec.mul is_vec z ~z:(Vec.mul is_vec z ~z)
        in
        {
          common_model = common_model;
          eval_trained = eval_trained;
          t_vec = t_vec;
          u_vec = u_vec;
          v_vec = v_vec;
          y_minus_kmn_t = y_minus_kmn_t;
        }

      let calc_eval trained = trained.eval_trained

      (**)

      let calc_log_evidence_sigma2 { common_model = cm; v_vec = v_vec } =
        0.5 *. Cm.common_calc_log_evidence_sigma2 cm v_vec

      (**)

      type hyper_t = { shared : Shared.t; dfacts : Dfacts.t }

      let prepare_hyper trained =
        let common_model = trained.common_model in
        let eval_model = common_model.Cm.eval_model in
        let u_mat = calc_u_mat eval_model in
        let t_vec = trained.t_vec in
        let u_vec = trained.u_vec in
        let w_mat =
          let w_mat =
            let t_mat = trained.common_model.Cm.t_mat in
            syr ~alpha:(-1.) t_vec (lacpy ~uplo:`U t_mat)
          in
          let u1_mat = lacpy u_mat in
          let sqrt_v1_vec = Cm.calc_v1_vec common_model in
          (* TODO: add more efficient in-place Vec.sqrt to LACAML *)
          Mat.scal_cols (Vec.sqrt ~y:sqrt_v1_vec sqrt_v1_vec) u1_mat;
          let w_mat = syrk ~alpha:(-1.) u1_mat ~beta:1. ~c:w_mat in
          let u2_mat = u1_mat in
          let u2_mat = lacpy u_mat ~b:u2_mat in
          Mat.scal_cols u_vec u2_mat;
          syrk u2_mat ~beta:(1.) ~c:w_mat
        in
        let v_vec = trained.v_vec in
        let x_mat =
          Mat.scal_cols v_vec u_mat;
          let x_mat = lacpy common_model.Cm.s_mat in
          Mat.axpy ~alpha:(-1.) ~x:u_mat x_mat;
          ger ~alpha:(-1.) t_vec u_vec x_mat
        in
        {
          shared = common_model.Cm.shared;
          dfacts = { Dfacts.v_vec = v_vec; w_mat = w_mat; x_mat = x_mat };
        }


      (**)

      let calc_log_evidence hyp_trained hyper =
        common_calc_log_evidence hyp_trained.shared hyp_trained.dfacts hyper
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
