open Format

open Lacaml.Impl.D
open Lacaml.Io

(* NOTE: for testing *)
let print_float name n = printf "%s: @[%f@]@.@." name n
let print_vec name vec = printf "%s: @[%a@]@.@." name pp_fvec vec
let print_mat name vec = printf "%s: @[%a@]@.@." name pp_fmat vec

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

let inv_copy_chol_mat mat =
  (* TODO: copy triangle *)
  let inv = Mat.copy mat in
  potri ~factorize:false inv;
  inv

let inv_copy_chol_vec chol vec =
  let inv = copy vec in
  trsv chol inv;
  inv

let sub_diag_transa_prod mat1 mat2 ~dst =
  (* TODO: optimize away col, dot *)
  for i = 1 to Vec.dim dst do
    dst.{i} <- dst.{i} -. dot ~x:(Mat.col mat1 i) (Mat.col mat2 i)
  done

let default_rng = Gsl_rng.make (Gsl_rng.default ())

(*
let get_unit_gaussian_vec ~dst =
  let sum = ref 0. in
  let n = Vec.dim dst in
  for i = 1 to n do
    let x = Gsl_randist.gaussian rng ~sigma:1. in
    sum := !sum +. x;
    dst.{i} <- x
  done;
  let ssqr = ref 0. in
  let f_n = float n in
  let mean = !sum /. f_n in
  for i = 1 to n do
    let x = dst.{i} in
    let diff = x -. mean in
    dst.{i} <- diff;
    ssqr := !ssqr +. diff *. diff
  done;
  let norm = sqrt !ssqr in
  scal (1. /. norm) dst

let draw_samples_chol ?(jitter = 10E-9) ?dst ~means_dst ~covs_chol n =
  let dst =
    match dst with
    | None -> Vec.create n
    | Some dst -> Lacaml.Utils.check_vec "draw_samples_chol" "dst" dst n; dst
  in
  get_unit_gaussian_vec ~dst;
  trmv ~trans:`T covs_chol dst

let draw_samples ?(jitter = 10E-9) ~means ~covs ~n =
  let covs_chol = Mat.copy covs in
  potrf ~jitter covs_chol;
  let u =
    (assert false (* XXX *))
  in
  trmv ~trans:`T covs_chol vec
*)
