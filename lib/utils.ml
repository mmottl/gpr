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
