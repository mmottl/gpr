open Format

open Lacaml.Impl.D
open Lacaml.Io

(* NOTE: for testing *)
let print_float name n = printf "%s: @[%.9f@]@.@." name n
let print_vec name vec = printf "%s: @[%a@]@.@." name pp_fvec vec
let print_mat name vec = printf "%s: @[%a@]@.@." name pp_fmat vec

let timing name f =
  let t1 = Unix.times () in
  let res = f () in
  let t2 = Unix.times () in
  printf "%s %.2f@." name (t2.Unix.tms_utime -. t1.Unix.tms_utime);
  res

let gen_write pp file obj =
  let oc = open_out (Filename.concat "data" file) in
  fprintf (formatter_of_out_channel oc) "%a@." pp obj;
  close_out oc

let write_float file = gen_write pp_print_float file
let write_vec file = gen_write pp_fvec file
let write_mat file = gen_write pp_fmat file

(* Assumes Cholesky factorized matrix *)
let log_det mat =
  let n = Mat.dim1 mat in
  if Mat.dim2 mat <> n then failwith "log_det: not a square matrix";
  let rec loop acc i =
    if i = 0 then 2. *. acc
    else loop (acc +. log mat.{i, i}) (i - 1)
  in
  loop 0. n

let pi = 4. *. atan 1.
let log_2pi = log (2. *. pi)

let default_rng = Gsl_rng.make (Gsl_rng.default ())

let cholesky_jitter = ref 10e-9
