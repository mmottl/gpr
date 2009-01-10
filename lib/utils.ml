open Format

open Lacaml.Impl.D
open Lacaml.Io

(* NOTE: for testing *)
let print_float name n = printf "%s: @[%.9f@]@.@." name n
let print_vec name vec = printf "%s: @[%a@]@.@." name pp_fvec vec
let print_mat name mat = printf "%s: @[%a@]@.@." name pp_fmat mat

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

let cholesky_jitter = ref 10e-6

let solve_triangular ?trans chol ~k =
  (* TODO: special case for copying symmetric matrices? *)
  let inv_chol_k = Mat.copy k in
  trtrs ?trans chol inv_chol_k;
  inv_chol_k

let inv_chol chol =
  (* TODO: copy upper triangle only *)
  let inv = Mat.copy chol in
  potri ~factorize:false inv;
  inv


(* Sparse row matrices *)

(* Checks whether a sparse row matrix is sane *)
let check_sparse_sane mat rows =
  let m = Mat.dim1 mat in
  let n_rows = Array.length rows in
  if n_rows <> m then
    failwith (
      sprintf
        "Deriv.Model.check_sparse_sane: number of rows in \
        sparse matrix (%d) disagrees with size of row array (%d)"
        m n_rows);
  let rec loop ~i ~next =
    if i >= 0 then
      let rows_i = rows.(i) in
      if rows_i <= 0 then
        failwith (
          sprintf
            "Deriv.Model.check_sparse_sane: sparse row %d contains \
            illegal negative real row index %d" (i + 1) rows_i)
      else if rows_i >= next then
        failwith (
          sprintf
            "Deriv.Model.check_sparse_sane: sparse row %d \
            associated with real row index %d violates strict \
            ordering (next: %d)"
            (i + 1) rows_i next)
      else loop ~i:(i - 1) ~next:rows_i
  in
  loop ~i:(n_rows - 1) ~next:(Mat.dim2 mat + 1)

(* Detriangularize sparse row matrix to make it symmetric *)
let detri_sparse mat rows =
  let m = Mat.dim1 mat in
  let n = Mat.dim2 mat in
  if n < rows.(m - 1) then
    failwith "Deriv.Model.detri_sparse: sparse matrix cannot be square"
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

let update_prod_diag dst fact mat1 mat2 =
  for i = 1 to (Mat.dim2 mat1) do
    (* TODO: optimize dot and col *)
    let diag_i = dot ~x:(Mat.col mat1 i) (Mat.col mat2 i) in
    dst.{i} <- dst.{i} +. fact *. diag_i
  done

let calc_prod_trace mat1 mat2 =
  let rec loop trace i =
    if i = 0 then trace
    else
      (* TODO: optimize dot and col *)
      let diag_i = dot ~x:(Mat.col mat1 i) (Mat.col mat2 i) in
      loop (trace +. diag_i) (i - 1)
  in
  loop 0. (Mat.dim2 mat1)
