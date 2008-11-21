open Lacaml.Impl.D
open Lacaml.Io

module type From_vec = sig
  type t
  type input = vec
  type inputs = mat

  val get_n_inputs : inputs -> int
  val eval_one : t -> input -> float
  val eval : t -> input -> input -> float
  val eval_vec_col : t -> input -> inputs -> int -> float
  val eval_mat_cols : t -> inputs -> int -> inputs -> int -> float
  val eval_mat_col : t -> inputs -> int -> float
end

module Make_from_vec (From_vec : From_vec) = struct
  include From_vec

  let evals t vec mat =
    let n = Mat.dim2 mat in
    let dst = Vec.create n in
    for col = 1 to n do
      dst.{col} <- eval_vec_col t vec mat col
    done;
    dst

  let weighted_eval t ~weights vec mat =
    let n = Vec.dim vec in
    let res = ref 0. in
    for col = 1 to n do
      res := !res +. weights.{col} *. eval_vec_col t vec mat col
    done;
    !res

  let weighted_evals t ~weights inducing_inputs inputs =
    let n_inducing_inputs = Mat.dim2 inducing_inputs in
    let n_inputs = Mat.dim2 inputs in
    let dst = Vec.create n_inputs in
    for i = 1 to n_inputs do
      dst.{i} <- weights.{1} *. eval_mat_cols t inducing_inputs 1 inputs i;
      for j = 2 to n_inducing_inputs do
        dst.{i} <-
          dst.{i} +. weights.{j} *. eval_mat_cols t inputs i inducing_inputs j
      done
    done;
    dst

  let upper t mat =
    let n = Mat.dim2 mat in
    let dst = Mat.create n n in
    for i = 1 to n do
      dst.{i, i} <- eval_mat_col t mat i;
      for j = i + 1 to n do
        dst.{i, j} <- eval_mat_cols t mat i mat j;
      done
    done;
    dst

  let upper_no_diag t mat =
    let n = Mat.dim2 mat in
    let dst = Mat.create n n in
    for i = 1 to n do
      for j = i + 1 to n do
        dst.{i, j} <- eval_mat_cols t mat i mat j;
      done
    done;
    dst

  let cross t mat1 mat2 =
    let n1 = Mat.dim2 mat1 in
    let n2 = Mat.dim2 mat2 in
    let dst = Mat.create n1 n2 in
    for i = 1 to n1 do
      for j = 1 to n2 do
        dst.{i, j} <- eval_mat_cols t mat1 i mat2 j
      done
    done;
    dst

  let diag_mat t mat ~dst =
    let n = Mat.dim2 mat in
    for i = 1 to n do
      dst.{i, i} <- eval_mat_col t mat i
    done

  let diag_vec t mat =
    let n = Mat.dim2 mat in
    let dst = Vec.create n in
    for i = 1 to n do
      dst.{i} <- eval_mat_col t mat i
    done;
    dst

(* TODO: this is surprisingly faster; maybe implement weighted, etc.,
   dot product operations on matrix columns, etc. in C.

  let evals vec1 mat ~dst =
    let n = Mat.dim2 mat in
    for i = 1 to n do
      let vec2 = Mat.col mat i in
      dst.{i} <- eval vec1 vec2
    done

  let weighted_eval ~weights vec1 mat =
    let n = Mat.dim2 mat in
    let res = ref 0. in
    for i = 1 to n do
      let vec2 = Mat.col mat i in
      res := !res +. weights.{i} *. eval vec1 vec2
    done;
    !res

  let weighted_evals ~weights mat1 mat2 ~dst =
    let n1 = Mat.dim2 mat1 in
    let n2 = Mat.dim2 mat2 in
    for i = 1 to n1 do
      let vec1 = Mat.col mat1 i in
      dst.{i} <- 0.;
      for j = 1 to n2 do
        let vec2 = Mat.col mat2 j in
        dst.{i} <- dst.{i} +. weights.{i} *. eval vec1 vec2
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

  let diag_mat mat ~dst =
    let n = Mat.dim2 mat in
    for i = 1 to n do
      let vec = Mat.col mat i in
      dst.{i, i} <- eval vec vec
    done

  let diag_vec mat ~dst =
    let n = Mat.dim2 mat in
    for i = 1 to n do
      let vec = Mat.col mat i in
      dst.{i} <- eval vec vec
    done
*)
end

module Gauss_vec = struct
  type t = float * float
  type input = vec
  type inputs = mat

  let get_n_inputs inputs = Mat.dim2 inputs

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

module Gauss = Make_from_vec (Gauss_vec)
