open Lacaml.Impl.D

open Interfaces
open Utils
open Fitc_impl

module type Spec = sig
  include Fitc_impl.Spec

  val update_sigma2 : Eval_spec.kernel -> float -> Eval_spec.kernel
end

module type Sig = functor (SPGP_spec : Spec) ->
  Inducing_input_gpr.Sigs.Eval with module Spec = SPGP_spec.Eval_spec

let sigma2_deriv ~inv_lam_sigma2_diag ~kmn ~b_chol ~inv_b_chol_kmn_y__ ~y__ =
  let inv_b_chol_kmn = Mat.copy kmn in
  trtrs ~trans:`T b_chol inv_b_chol_kmn;
  let inv_b_kmn_y__ = copy inv_b_chol_kmn_y__ in
  trsv b_chol inv_b_kmn_y__;
  let knm_inv_b_kmn_y__ = gemv ~trans:`T kmn inv_b_kmn_y__ in
  let rec loop res i =
    if i = 0 then 0.5 *. res
    else
      let new_res =
        (* TODO: optimize ssqr and col *)
        let ssqr = Vec.ssqr (Mat.col inv_b_chol_kmn i) in
        let inv_lam_sigma2_i = inv_lam_sigma2_diag.{i} in
        let l1 = inv_lam_sigma2_i *. (1. -. ssqr *. inv_lam_sigma2_i) in
        let t1 = inv_lam_sigma2_i *. knm_inv_b_kmn_y__.{i} in
        let y__i = y__.{i} in
        let l2 = y__i *. y__i +. t1 *. (t1 -. 2. *. y__i) in
        res +. l1 -. l2
      in
      loop new_res (i - 1)
  in
  loop 0. (Mat.dim2 kmn)

(*
let () =
  Lacaml.Io.pp_float_el_default := (fun ppf n -> Format.fprintf ppf "%.9f" n);
  let m = 10 in
  let n = 20 in
  let kmn = Mat.init_cols m n (fun r c -> float (c * m + r)) in
  let b_chol = Mat.hilbert m in
  for i = 1 to m do b_chol.{i, i} <- b_chol.{i, i} +. 100. done;
  potrf b_chol;
  let sigma2 = 1.3 in
  let inv_lam_sigma2_diag = Vec.init n (fun i -> 1. /. (float i +. sigma2)) in
  let y__ = Vec.init n float in
  let inv_b_chol_kmn_y__ = gemv kmn y__ in
  trsv ~trans:`T b_chol inv_b_chol_kmn_y__;
  let res =
    sigma2_deriv ~inv_lam_sigma2_diag ~kmn ~b_chol ~inv_b_chol_kmn_y__ ~y__
  in
  print_float "res" res;
*)
