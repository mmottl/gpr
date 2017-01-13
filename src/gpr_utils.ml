(* File: utils.ml

   OCaml-GPR - Gaussian Processes for OCaml

     Copyright (C) 2009-  Markus Mottl
     email: markus.mottl@gmail.com
     WWW:   http://www.ocaml.info

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*)

open Core.Std
open Bigarray
open Lacaml.D

(* Global definitions *)

module Int_vec = struct
  type t = (int, int_elt, fortran_layout) Array1.t
  let create n : t = Array1.create int fortran_layout n
  let dim (t : t) = Array1.dim t
  let sub (t : t) n = Array1.sub t n
end

let debug = ref false
let cholesky_jitter = ref 1e-6

type fast_float_ref = { mutable x : float }

let pi = 4. *. atan 1.
let log_2pi = log (pi +. pi)

let default_rng = Gsl.Rng.make (Gsl.Rng.default ())


(* Testing and I/O functionality *)

let print_int name n = Format.printf "%s: @[%d@]@.@." name n
let print_float name n = Format.printf "%s: @[%.9f@]@.@." name n
let print_vec name vec = Format.printf "%s: @[%a@]@.@." name pp_vec vec
let print_mat name mat = Format.printf "%s: @[%a@]@.@." name pp_mat mat

let timing name f =
  let t1 = Unix.times () in
  let res = f () in
  let t2 = Unix.times () in
  Format.printf "%s %.2f@." name (t2.Unix.tms_utime -. t1.Unix.tms_utime);
  res


(* General matrix functions *)

(* Choose columns of a matrix *)
let choose_cols mat indexes =
  let m = Mat.dim1 mat in
  let n = Mat.dim2 mat in
  let k = Int_vec.dim indexes in
  let res = Mat.create m k in
  for c = 1 to k do
    let real_c = indexes.{c} in
    if real_c < 1 || real_c > n then
      failwithf
        "Gpr.Gpr_utils.choose_cols: violating 1 <= index (%d) <= dim (%d)"
        real_c n ()
    else for r = 1 to m do res.{r, c} <- mat.{r, real_c} done
  done;
  res

(* Compute the sum of all elements in a matrix *)
let sum_mat mat = Vec.sum (Mat.as_vec mat)

(* Compute the sum of all elements in a symmetric matrix *)
let sum_symm_mat mat =
  let diag_ref = ref 0. in
  let rest_ref = ref 0. in
  let n = Mat.dim1 mat in
  for c = 1 to n do
    for r = 1 to c - 1 do rest_ref := !rest_ref +. mat.{r, c} done;
    diag_ref := !diag_ref +. mat.{c, c}
  done;
  let rest = !rest_ref in
  rest +. !diag_ref +. rest

(* Computes logarithm of determinant; assumes Cholesky factorized matrix *)
let log_det mat =
  let n = Mat.dim1 mat in
  if Mat.dim2 mat <> n then failwith "log_det: not a square matrix";
  let rec loop acc i =
    if i = 0 then acc +. acc
    else loop (acc +. log mat.{i, i}) (i - 1)
  in
  loop 0. n

(* Solve triangular system *)
let solve_tri ?trans chol mat =
  let ichol_mat = lacpy mat in
  trtrs ?trans chol ichol_mat;
  ichol_mat

(* Compute the inverse of a matrix using the cholesky factor *)
let ichol chol =
  let inv = lacpy ~uplo:`U chol in
  potri ~factorize:false inv;
  inv


(* Sparse matrices and vectors *)

(* Checks whether a sparse row matrix is sane *)
let check_sparse_row_mat_sane ~real_m ~smat ~rows =
  if !debug then begin
    if real_m < 0 then
      failwith "Gpr.Gpr_utils.check_sparse_row_mat_sane: real_m < 0";
    let m = Mat.dim1 smat in
    let n_rows = Int_vec.dim rows in
    if n_rows <> m then
      failwithf
        "Gpr.Gpr_utils.check_sparse_row_mat_sane: number of rows in \
        sparse matrix (%d) disagrees with size of row array (%d)"
        m n_rows ();
    let rec loop ~i ~limit =
      if i > 0 then
        let rows_i = rows.{i} in
        if rows_i <= 0 then
          failwithf
            "Gpr.Gpr_utils.check_sparse_row_mat_sane: sparse row %d contains \
            illegal negative real row index %d" i rows_i ()
        else if rows_i > limit then
          failwithf
            "Gpr.Gpr_utils.check_sparse_row_mat_sane: sparse row %d \
            associated with real row index %d violates consistency \
            (current row limit: %d)"
            i rows_i limit ()
        else loop ~i:(i - 1) ~limit:rows_i
    in
    loop ~i:n_rows ~limit:real_m
  end

(* Checks whether a sparse column matrix is sane *)
let check_sparse_col_mat_sane ~real_n ~smat ~cols =
  if !debug then begin
    if real_n < 0 then
      failwith "Gpr.Gpr_utils.check_sparse_col_mat_sane: real_n < 0";
    let n = Mat.dim2 smat in
    let n_cols = Int_vec.dim cols in
    if n_cols <> n then
      failwithf
        "Gpr.Gpr_utils.check_sparse_col_mat_sane: number of cols in \
        sparse matrix (%d) disagrees with size of col array (%d)"
        n n_cols ();
    let rec loop ~i ~limit =
      if i > 0 then
        let cols_i = cols.{i} in
        if cols_i <= 0 then
          failwithf
            "Gpr.Gpr_utils.check_sparse_col_mat_sane: sparse col %d contains \
            illegal negative real col index %d" i cols_i ()
        else if cols_i > limit then
          failwithf
            "Gpr.Gpr_utils.check_sparse_col_mat_sane: sparse col %d \
            associated with real col index %d violates consistency \
            (current col limit: %d)"
            i cols_i limit ()
        else loop ~i:(i - 1) ~limit:cols_i
    in
    loop ~i:n_cols ~limit:real_n
  end

(* Checks whether a parse vector is sane *)
let check_sparse_vec_sane ~real_n ~svec ~rows =
  if !debug then
    let k = Vec.dim svec in
    if Int_vec.dim rows <> k then
      failwith
        "Gpr.Gpr_utils.check_sparse_vec_sane: \
        size of sparse vector disagrees with indexes";
    let rec loop ~last i =
      if i > 0 then
        let ind = rows.{i} in
        if ind >= last || ind <= 0 then
          failwith "Gpr.Gpr_utils.check_sparse_vec_sane: rows inconsistent"
        else loop ~last:ind (i - 1)
    in
    loop ~last:real_n (Int_vec.dim rows)

(* Computes the trace of the product of a symmetric and sparse
   symmetric matrix *)
let symm2_sparse_trace ~mat ~smat ~rows =
  let m = Int_vec.dim rows in
  let n = Mat.dim2 smat in
  let full_ref = ref 0. in
  let half_ref = ref 0. in
  let rows_ix_ref = ref 1 in
  for sparse_r = 1 to m do
    let c = rows.{sparse_r} in
    for r = 1 to n do
      let mat_el = if r > c then mat.{c, r} else mat.{r, c} in
      let rows_ix = !rows_ix_ref in
      if
        rows_ix > m ||
        let rows_el = rows.{rows_ix} in
        r < rows_el || c < rows_el
      then full_ref := !full_ref +. mat_el *. smat.{sparse_r, r}
      else begin
        half_ref := !half_ref +. mat_el *. smat.{rows_ix, c};
        incr rows_ix_ref
      end
    done;
    rows_ix_ref := 1
  done;
  let full = !full_ref in
  full +. !half_ref +. full
