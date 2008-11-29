open Lacaml.Impl.D
open Lacaml.Io

module type From_all_vec = sig
  type kernel

  val eval_one : kernel -> vec -> float
  val eval : kernel -> vec -> vec -> float
  val eval_vec_col : kernel -> vec -> mat -> int -> float
  val eval_mat_cols : kernel -> mat -> int -> mat -> int -> float
  val eval_mat_col : kernel -> mat -> int -> float
end

module Make_from_all_vec (Spec : From_all_vec) = struct
  include Spec

  module Inducing = struct
    type t = mat

    let size points = Mat.dim2 points

    let upper kernel mat =
      let n = Mat.dim2 mat in
      let dst = Mat.create n n in
      for i = 1 to n do
        dst.{i, i} <- eval_mat_col kernel mat i;
        for j = i + 1 to n do
          dst.{i, j} <- eval_mat_cols kernel mat i mat j;
        done
      done;
      dst
  end

  module Input = struct
    type t = vec

    let eval_one = eval_one

    let eval kernel ~inducing ~input =
      let n = Mat.dim2 inducing in
      let dst = Vec.create n in
      for col = 1 to n do
        dst.{col} <- eval_vec_col kernel input inducing col
      done;
      dst

    let weighted_eval kernel ~coeffs ~inducing ~input =
      let n = Vec.dim input in
      let res = ref 0. in
      for col = 1 to n do
        res := !res +. coeffs.{col} *. eval_vec_col kernel input inducing col
      done;
      !res
  end

  module Inputs = struct
    include Inducing

    let weighted_eval kernel ~coeffs ~inducing ~inputs =
      let n_inducing_inputs = Mat.dim2 inducing in
      let n_inputs = Mat.dim2 inputs in
      let dst = Vec.create n_inputs in
      for i = 1 to n_inputs do
        dst.{i} <- coeffs.{1} *. eval_mat_cols kernel inducing 1 inputs i;
        for j = 2 to n_inducing_inputs do
          dst.{i} <-
            dst.{i} +.  coeffs.{j} *. eval_mat_cols kernel inputs i inducing j
        done
      done;
      dst

    let upper_no_diag kernel mat =
      let n = Mat.dim2 mat in
      let dst = Mat.create n n in
      for i = 1 to n do
        for j = i + 1 to n do
          dst.{i, j} <- eval_mat_cols kernel mat i mat j;
        done
      done;
      dst

    let cross kernel ~inducing ~inputs =
      let n1 = Mat.dim2 inducing in
      let n2 = Mat.dim2 inputs in
      let dst = Mat.create n1 n2 in
      for i = 1 to n1 do
        for j = 1 to n2 do
          dst.{i, j} <- eval_mat_cols kernel inducing i inputs j
        done
      done;
      dst

    let diag kernel mat =
      let n = Mat.dim2 mat in
      let dst = Vec.create n in
      for i = 1 to n do
        dst.{i} <- eval_mat_col kernel mat i
      done;
      dst
  end

(* TODO: this is surprisingly faster; maybe implement weighted, etc.,
   dot product operations on matrix columns, etc. in C.

  let evals vec1 mat ~dst =
    let n = Mat.dim2 mat in
    for i = 1 to n do
      let vec2 = Mat.col mat i in
      dst.{i} <- eval vec1 vec2
    done

  let weighted_eval ~coeffs vec1 mat =
    let n = Mat.dim2 mat in
    let res = ref 0. in
    for i = 1 to n do
      let vec2 = Mat.col mat i in
      res := !res +. coeffs.{i} *. eval vec1 vec2
    done;
    !res

  let weighted_evals ~coeffs mat1 mat2 ~dst =
    let n1 = Mat.dim2 mat1 in
    let n2 = Mat.dim2 mat2 in
    for i = 1 to n1 do
      let vec1 = Mat.col mat1 i in
      dst.{i} <- 0.;
      for j = 1 to n2 do
        let vec2 = Mat.col mat2 j in
        dst.{i} <- dst.{i} +. coeffs.{i} *. eval vec1 vec2
      done
    done

  let upper mat ~dst =
    let n = Mat.dim2 mat in
    for i = 1 to n do
      let vec1 = Mat.col mat i in
      for j = i to n do
        let vec2 = Mat.col mat j in
        dst.{j, i} <- eval vec1 vec2
      done
    done

  let cross mat1 mat2 ~dst =
    let n1 = Mat.dim2 mat1 in
    let n2 = Mat.dim2 mat2 in
    for i = 1 to n1 do
      let vec1 = Mat.col mat1 i in
      for j = 1 to n2 do
        let vec2 = Mat.col mat2 j in
        dst.{i, j} <- eval vec1 vec2
      done
    done

  let diag mat ~dst =
    let n = Mat.dim2 mat in
    for i = 1 to n do
      let vec = Mat.col mat i in
      dst.{i} <- eval vec vec
    done
*)
end

module Gauss_all_vec_spec = struct
  type kernel = float * float

  let eval_rbf2 (a, b) r = exp (a +. b *. r)

  let eval_one k vec = eval_rbf2 k (Vec.ssqr vec)

  let eval k vec1 vec2 = eval_rbf2 k (Vec.ssqr_diff vec1 vec2)

  let eval_mat_col = fun (a, _b) _mat _col -> exp a

  let eval_vec_col k vec mat col =
    let d = Vec.dim vec in
    let r2 = ref 0. in
    for i = 1 to d do
      let diff = vec.{i} -. mat.{i, col} in
      r2 := !r2 +. diff *. diff
    done;
    eval_rbf2 k !r2

  let eval_mat_cols k mat1 col1 mat2 col2 =
    let d = Mat.dim1 mat1 in
    let r2 = ref 0. in
    for i = 1 to d do
      let diff = mat1.{i, col1} -. mat2.{i, col2} in
      r2 := !r2 +. diff *. diff
    done;
    eval_rbf2 k !r2
end

module Gauss_all_vec = Make_from_all_vec (Gauss_all_vec_spec)

module Wiener_all_vec_spec = struct
  type kernel = float

  let eval_one k vec = exp k *. sqrt (Vec.ssqr vec)

  let eval k vec1 vec2 = exp k *. sqrt (min (Vec.ssqr vec1) (Vec.ssqr vec2))

  let eval_mat_col k mat col = eval_one k (Mat.col mat col)

  let eval_vec_col k vec mat col = eval k vec (Mat.col mat col)

  let eval_mat_cols k mat1 col1 mat2 col2 =
    eval k (Mat.col mat1 col1) (Mat.col mat2 col2)
end

module Wiener_all_vec = Make_from_all_vec (Wiener_all_vec_spec)
