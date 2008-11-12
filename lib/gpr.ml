open Format

open Lacaml.Impl.D
open Lacaml.Io

open Interfaces

(* NOTE: for testing *)
let print_float name n = printf "%s: @[%f@]@.@." name n
let print_vec name vec = printf "%s: @[%a@]@.@." name pp_fvec vec
let print_mat name vec = printf "%s: @[%a@]@.@." name pp_fmat vec

module type FIC = sig
  module Common : sig
    type kernel

    type t = {
      kernel : kernel;
      sigma2 : float;
      inducing_points : mat;
    }
  end

  module Trained : sig
    type t

    val train : ?jitter : float -> Common.t -> inputs : mat -> targets : vec -> t
    val neg_log_likelihood : t -> float
  end

  module Mean_predictor : sig
    type t = private {
      common : Common.t;
      coeffs : vec;
    }

    val of_trained : Trained.t -> t
    val mean : t -> input : vec -> float
    val means : t -> inputs : mat -> dst : vec -> unit
  end

  module Full_predictor : sig
    type t
    type fmt = [ `Inv | `Chols | `Packed_inv | `Packed_chols ]
    type kind = [ `Fic | `Fitc ]

    val of_trained : kind -> ?fmt : fmt -> Trained.t -> t
    val mean_predictor : t -> Mean_predictor.t
    val mean_variance : ?with_noise : bool -> t -> input : vec -> float * float
    val means_variances : ?with_noise : bool -> t -> inputs : mat -> vec * vec
    val means_covariances : ?with_noise : bool -> t -> inputs : mat -> vec * mat
  end
end

module type Kernel = sig
  type t

  include Kernel_fun
    with type input = vec
    with type inputs = mat
end

(* Assumes Cholesky factorized matrix *)
let log_det mat =
  let n = Mat.dim1 mat in
  if Mat.dim2 mat <> n then failwith "log_det: not a square matrix";
  let rec loop acc i =
    if i = 0 then acc
    else loop (acc +. log mat.{i, i}) (i - 1)
  in
  2. *. loop 0. n

let pi = 4. *. atan 1.
let log_2pi = log (2. *. pi)

module Make_FIC (Kernel : Kernel) : FIC with type Common.kernel = Kernel.t =
struct
  module Common = struct
    type kernel = Kernel.t

    type t = {
      kernel : kernel;
      sigma2 : float;
      inducing_points : mat;
    }
  end

  module Trained = struct
    type t = {
      common : Common.t;
      km_chol : mat;
      kmn_y_ : vec;
      b_chol : mat;
      neg_log_likelihood : float;
    }

    let calc_kernel_mats ~inputs ~inducing_points =
      let n_inducing_points = Mat.dim2 inducing_points in
      let n_inputs = Mat.dim2 inputs in
      let km = Mat.create n_inducing_points n_inducing_points in
      let kmn = Mat.create n_inducing_points n_inputs in
      Kernel.upper inducing_points ~dst:km;
      Kernel.cross inducing_points inputs ~dst:kmn;
      km, kmn

    let common_train ~inputs ~targets ~sigma2 ~km_chol ~kmn =
      let n_inputs = Mat.dim2 inputs in
      let kn = Vec.create n_inputs in
      Kernel.diag_vec inputs ~dst:kn;
      let inv_lam_sigma2 = kn in
      let sqrt_inv_lam_sigma2 = Vec.create n_inputs in
      let y_ = Vec.create n_inputs in
      (* TODO: copying may not be necessary *)
      let kmn_ = Mat.copy kmn in
      let km_2_kmn = Mat.copy kmn in
      (* TODO: most expensive *)
      trtrs ~trans:`T km_chol km_2_kmn;
      (* TODO: loop somewhat expensive *)
      let rec loop_lambda sum_log_lam_sigma2_diag i =
        if i = 0 then sum_log_lam_sigma2_diag
        else
          let kn_i = kn.{i} in
          (* TODO: optimize ssqr and col *)
          let qn_i = Vec.ssqr (Mat.col km_2_kmn i) in
          let lam_sigma2 = kn_i -. qn_i +. sigma2 in
          let alpha = 1. /. lam_sigma2 in
          inv_lam_sigma2.{i} <- alpha;
          let sqrt_alpha = sqrt alpha in
          sqrt_inv_lam_sigma2.{i} <- sqrt_alpha;
          y_.{i} <- sqrt_alpha *. targets.{i};
          (* TODO: optimize col *)
          scal sqrt_alpha (Mat.col kmn_ i);
          loop_lambda (sum_log_lam_sigma2_diag +. log lam_sigma2) (i - 1)
      in
      loop_lambda 0. n_inputs, y_, kmn_, inv_lam_sigma2

    let train ?(jitter = 10e-9) common ~inputs ~targets =
      let { Common.inducing_points = inducing_points } = common in
      let km, kmn = calc_kernel_mats ~inputs ~inducing_points in

      (* TODO: copy symmetric *)
      let km_chol = Mat.copy km in
      potrf ~up:true ~jitter km_chol;

      let sigma2 = common.Common.sigma2 in

      let log_det_lam_sigma2, y_, kmn_, inv_lam_sigma2 =
        common_train ~inputs ~targets ~sigma2 ~km_chol ~kmn
      in

      (* TODO: copying may not be necessary; copy symmetric *)
      let b_chol = syrk ~beta:1. ~c:(Mat.copy km) kmn_ in
      potrf ~jitter b_chol;

      let ssqr_y_ = Vec.ssqr y_ in

      let kmn_y_ = gemv kmn_ y_ in

      (* TODO: copying may not be necessary *)
      let b_kmn_y_sol_ = copy kmn_y_ in
      trsv ~trans:`T b_chol b_kmn_y_sol_;

      let log_det_b = log_det b_chol in

      let log_det_km = log_det km_chol in
      let l1_2 = log_det_b -. log_det_km +. log_det_lam_sigma2 in
      let l2_2 = ssqr_y_ -. Vec.ssqr b_kmn_y_sol_ in

      let n_inputs = Mat.dim2 inputs in
      let f_inputs = float n_inputs in
      let neg_log_likelihood = (l1_2 +. l2_2 +. f_inputs *. log_2pi) /. 2. in

      {
        common = common;
        km_chol = km_chol;
        kmn_y_ = kmn_y_;
        b_chol = b_chol;
        neg_log_likelihood = neg_log_likelihood;
      }

      let neg_log_likelihood t = t.neg_log_likelihood
  end

  module Mean_predictor = struct
    type t = {
      common : Common.t;
      coeffs : vec;
    }

    let of_trained trained =
      let
        {
          Trained.
          common = common;
          kmn_y_ = kmn_y_;
          b_chol = b_chol;
        } = trained
      in
      let coeffs = copy kmn_y_ in
      let coeffs_mat = Mat.from_col_vec coeffs in
      potrs ~factorize:false b_chol coeffs_mat;
      {
        common = common;
        coeffs = coeffs;
      }

    let mean t ~input =
      Kernel.weighted_eval
        input ~coeffs:t.coeffs t.common.Common.inducing_points

    let means t ~inputs ~dst =
      Kernel.weighted_evals
        ~coeffs:t.coeffs t.common.Common.inducing_points inputs ~dst
  end

  module Full_predictor = struct
    type fmt = [ `Inv | `Chols | `Packed_inv | `Packed_chols ]
    type kind = [ `Fic | `Fitc ]

    type 'storage_type chols = { km_chol : 'storage_type; b_chol : 'storage_type }

    module Fic_invs = struct
      type 'storage_type t = {
        inv_km : 'storage_type;
        inv : 'storage_type;
      }
    end

    type covariance_factor =
      | Cf_chols of kind * mat chols
      | Cf_packed_chols of kind * vec chols
      | Fic_invs of mat Fic_invs.t
      | Fic_packed_invs of vec Fic_invs.t
      | Fitc_inv of mat
      | Fitc_packed_inv of vec

    type t =
      {
        mean_predictor : Mean_predictor.t;
        covariance_factor : covariance_factor;
      }

    let inv_copy_chol_mat mat =
      (* TODO: copy triangle *)
      let inv = Mat.copy mat in
      potri ~factorize:false inv;
      inv

    let calc_inv ~km_chol ~b_chol =
      let inv_km = inv_copy_chol_mat km_chol in
      let inv_b = inv_copy_chol_mat b_chol in
      let res = inv_b in
      Mat.axpy ~x:inv_km res;
      res, inv_km

    let of_trained kind ?(fmt =`Inv) trained =
      let { Trained.common = common; km_chol = km_chol; b_chol = b_chol } =
        trained
      in
      let mean_predictor = Mean_predictor.of_trained trained in
      let covariance_factor =
        match kind, fmt with
        | `Fic,`Inv ->
            let inv, inv_km = calc_inv ~km_chol ~b_chol in
            Fic_invs { Fic_invs.inv_km = inv_km; inv = inv }
        | `Fic, `Packed_inv ->
            let inv, inv_km = calc_inv ~km_chol ~b_chol in
            let packed_inv = Mat.packed inv in
            let packed_inv_km = Mat.packed inv_km in
            Fic_packed_invs {
              Fic_invs.
              inv_km = packed_inv_km;
              inv = packed_inv;
            }
        | `Fitc,`Inv -> Fitc_inv (fst (calc_inv ~km_chol ~b_chol))
        | `Fitc, `Packed_inv ->
            let inv = fst (calc_inv ~km_chol ~b_chol) in
            Fitc_packed_inv (Mat.packed inv)
        | _, `Chols -> Cf_chols (kind, { km_chol = km_chol; b_chol = b_chol })
        | _, `Packed_chols ->
            Cf_packed_chols (
              kind,
              {
                km_chol = Mat.packed km_chol;
                b_chol = Mat.packed b_chol;
              })
      in
      {
        mean_predictor = mean_predictor;
        covariance_factor = covariance_factor;
      }

    let mean_predictor t = t.mean_predictor

    let solve_chol_vec chol vec =
      let inv_vec = Mat.of_col_vecs [| vec |] in
      potrs ~factorize:false chol inv_vec;
      Mat.col inv_vec 1

    let cev_chols ~ks ~km_chol ~b_chol =
      let inv_km_ks = solve_chol_vec km_chol ks in
      let inv_b_ks = solve_chol_vec b_chol ks in
      let res = inv_km_ks in
      axpy ~x:inv_b_ks res;
      res

    let calc_explained_variance t ks =
      let right_term =
        match t.covariance_factor with
        | Fic_invs { Fic_invs.inv = inv } | Fitc_inv inv -> symv inv ks
        | Fic_packed_invs { Fic_invs.inv = pinv } | Fitc_packed_inv pinv ->
            (* TODO: interface pptrf, pptrs, pptri *)
            symv (Mat.unpacked pinv) ks
        | Cf_chols (_, { km_chol = km_chol; b_chol = b_chol }) ->
            cev_chols ~ks ~km_chol ~b_chol
        | Cf_packed_chols (_, pchols) ->
            (* TODO: interface pptrf, pptrs, pptri *)
            let km_chol = Mat.unpacked pchols.km_chol in
            let b_chol = Mat.unpacked pchols.b_chol in
            cev_chols ~ks ~km_chol ~b_chol
      in
      dot ks right_term

    let mean_variance ?(with_noise = true) t ~input =
      let { Common.inducing_points = inducing_points } as common =
        t.mean_predictor.Mean_predictor.common
      in
      let ks = Vec.create (Mat.dim2 inducing_points) in
      Kernel.evals input inducing_points ~dst:ks;
      let mean = dot ks t.mean_predictor.Mean_predictor.coeffs in
      let prior_variance = Kernel.eval input input in
      let explained_variance = calc_explained_variance t ks in
      let posterior_variance = prior_variance -. explained_variance in
      let variance =
        if with_noise then posterior_variance +. common.Common.sigma2
        else posterior_variance
      in
      mean, variance

    let calc_kmt ~inducing_points ~inputs =
      let n_inputs = Mat.dim2 inputs in
      let n_inducing_points = Mat.dim2 inducing_points in
      let ks = Mat.create n_inducing_points n_inputs in
      Kernel.cross inducing_points inputs ~dst:ks;
      ks

    let solve_chol_mat chol mat =
      let sol = Mat.copy mat in
      potrs ~factorize:false chol sol;
      sol

    let sub_diag_transa_prod mat1 mat2 ~dst =
      (* TODO: optimize away col, dot *)
      for i = 1 to Vec.dim dst do
        dst.{i} <- dst.{i} -. dot (Mat.col mat1 i) (Mat.col mat2 i)
      done

    let cevs_chols ~prior ~kmt ~km_chol ~b_chol =
      let res = solve_chol_mat km_chol kmt in
      let b_sol = solve_chol_mat b_chol kmt in
      Mat.axpy ~x:b_sol res;
      for i = 1 to Vec.dim prior do
        (* TODO: optimze ssqr + col *)
        prior.{i} <- prior.{i} -. Vec.ssqr (Mat.col res i)
      done

    let calc_explained_variances t ~prior ~kmt =
      match t.covariance_factor with
      | Fic_invs { Fic_invs.inv = inv } | Fitc_inv inv ->
          let sol = symm inv kmt in
          sub_diag_transa_prod sol kmt ~dst:prior
      | Fic_packed_invs { Fic_invs.inv = pinv } | Fitc_packed_inv pinv ->
          (* TODO: interface pptrf, pptrs, pptri *)
          let sol = symm ~side:`R kmt (Mat.unpacked pinv) in
          sub_diag_transa_prod kmt sol ~dst:prior
      | Cf_chols (_, { km_chol = km_chol; b_chol = b_chol }) ->
          cevs_chols ~prior ~kmt ~km_chol ~b_chol
      | Cf_packed_chols (_, pchols) ->
          (* TODO: interface pptrf, pptrs, pptri *)
          let km_chol = Mat.unpacked pchols.km_chol in
          let b_chol = Mat.unpacked pchols.b_chol in
          cevs_chols ~prior ~kmt ~km_chol ~b_chol

    (* TODO: pass in preallocated stuff *)
    let means_variances ?(with_noise = true) t ~inputs =
      let { Common.inducing_points = inducing_points } as common =
        t.mean_predictor.Mean_predictor.common
      in
      let kmt = calc_kmt ~inducing_points ~inputs in
      let means = gemv ~trans:`T kmt t.mean_predictor.Mean_predictor.coeffs in
      let n_inputs = Mat.dim2 inputs in
      let variances = Vec.create n_inputs in
      Kernel.diag_vec inputs ~dst:variances;
      calc_explained_variances t ~prior:variances ~kmt;
      if with_noise then begin
        let sigma2 = common.Common.sigma2 in
        for i = 1 to n_inputs do
          variances.{i} <- variances.{i} +. sigma2
        done;
      end;
      means, variances

    let cecs_chols ~prior ~kmt ~km_chol ~b_chol =
      let res = solve_chol_mat km_chol kmt in
      let b_sol = solve_chol_mat b_chol kmt in
      Mat.axpy ~x:b_sol res;
      syrk ~alpha:(-1.) res ~beta:1. ~c:prior

    let calc_explained_covariances t ~prior ~kmt =
      match t.covariance_factor with
      | Fic_invs { Fic_invs.inv = inv } | Fitc_inv inv ->
          gemm kmt (symm inv kmt) ~beta:1. ~c:prior
      | Fic_packed_invs { Fic_invs.inv = pinv } | Fitc_packed_inv pinv ->
          (* TODO: interface pptrf, pptrs, pptri *)
          gemm kmt (symm (Mat.unpacked pinv) kmt) ~beta:1. ~c:prior
      | Cf_chols (_, { km_chol = km_chol; b_chol = b_chol }) ->
          cecs_chols ~prior ~kmt ~km_chol ~b_chol
      | Cf_packed_chols (_, pchols) ->
          (* TODO: interface pptrf, pptrs, pptri *)
          let km_chol = Mat.unpacked pchols.km_chol in
          let b_chol = Mat.unpacked pchols.b_chol in
          cecs_chols ~prior ~kmt ~km_chol ~b_chol

    let cov_maybe_add_noise ?(with_noise = true) common covariances =
      if with_noise then
        let sigma2 = common.Common.sigma2 in
        for i = 1 to Mat.dim1 covariances do
          covariances.{i, i} <- covariances.{i, i} +. sigma2
        done

    let means_covariances ?with_noise t ~inputs =
      let { Common.inducing_points = inducing_points } as common =
        t.mean_predictor.Mean_predictor.common
      in
      let kmt = calc_kmt ~inducing_points ~inputs in
      let means = gemv kmt t.mean_predictor.Mean_predictor.coeffs in
      let n_inputs = Mat.dim2 inputs in
      let prior_covariances = Mat.create n_inputs n_inputs in
      let kind =
        match t.covariance_factor with
        | Fic_invs invs ->
            ignore (gemm ~transa:`T kmt (symm invs.Fic_invs.inv_km kmt)
              ~c:prior_covariances);
            `Fic
        | Fic_packed_invs pinvs ->
            (* TODO: interface pptrf, pptrs, pptri *)
            ignore (
              gemm ~transa:`T kmt (symm (Mat.unpacked pinvs.Fic_invs.inv_km) kmt)
                ~c:prior_covariances);
            `Fic
        | Fitc_inv _ | Fitc_packed_inv _ ->
            Kernel.upper inputs ~dst:prior_covariances;
            `Fitc
        | Cf_chols (kind, { km_chol = km_chol }) ->
            ignore (syrk (solve_chol_mat km_chol kmt) ~c:prior_covariances);
            kind
        | Cf_packed_chols (kind, { km_chol = km_chol }) ->
            (* TODO: interface pptrf, pptrs, pptri *)
            ignore (syrk (solve_chol_mat (Mat.unpacked km_chol) kmt)
              ~c:prior_covariances);
            kind
      in
      if kind = `Fic then Kernel.diag_mat inputs ~dst:prior_covariances;
      let covariances =
        calc_explained_covariances t ~prior:prior_covariances ~kmt
      in
      cov_maybe_add_noise ?with_noise common covariances;
      means, covariances
  end
end
