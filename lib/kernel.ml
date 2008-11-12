module type From_vec = sig
  type input = vec
  type inputs = mat

  val eval : input -> input -> float
  val eval_vec_col : input -> inputs -> int -> float
  val eval_mat_cols : inputs -> int -> inputs -> int -> float
  val eval_mat_col : inputs -> int -> float
end

module Make_from_vec (From_vec : From_vec) = struct
  include From_vec

  let evals vec mat ~dst =
    let n = Mat.dim2 mat in
    for col = 1 to n do
      dst.{col} <- eval_vec_col vec mat col
    done

  let weighted_eval ~coeffs vec mat =
    let n = Vec.dim vec in
    let res = ref 0. in
    for col = 1 to n do
      res := !res +. coeffs.{col} *. eval_vec_col vec mat col
    done;
    !res

  let weighted_evals ~coeffs inducing_inputs inputs ~dst =
    let n_inducing_inputs = Mat.dim2 inducing_inputs in
    let n_inputs = Mat.dim2 inputs in
    for i = 1 to n_inputs do
      dst.{i} <- coeffs.{1} *. eval_mat_cols inducing_inputs 1 inputs i;
      for j = 2 to n_inducing_inputs do
        dst.{i} <-
          dst.{i} +. coeffs.{j} *. eval_mat_cols inputs i inducing_inputs j
      done
    done

  let upper mat ~dst =
    let n = Mat.dim2 mat in
    for i = 1 to n do
      dst.{i, i} <- eval_mat_col mat i;
      for j = i + 1 to n do
        dst.{i, j} <- eval_mat_cols mat i mat j;
      done
    done

  let cross mat1 mat2 ~dst =
    let n1 = Mat.dim2 mat1 in
    let n2 = Mat.dim2 mat2 in
    for i = 1 to n1 do
      for j = 1 to n2 do
        dst.{i, j} <- eval_mat_cols mat1 i mat2 j
      done
    done

  let diag_mat mat ~dst =
    let n = Mat.dim2 mat in
    for i = 1 to n do
      dst.{i, i} <- eval_mat_col mat i
    done

  let diag_vec mat ~dst =
    let n = Mat.dim2 mat in
    for i = 1 to n do
      dst.{i} <- eval_mat_col mat i
    done

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
  type input = vec
  type inputs = mat

  let eval_rbf2 r = exp (-0.5 *. r)

  let eval vec1 vec2 = eval_rbf2 (Vec.ssqr_diff vec1 vec2)

  let res0 = eval_rbf2 0.

  let eval_mat_col = fun _mat _col -> res0

  let eval_vec_col vec mat col =
    let d = Vec.dim vec in
    let r2 = ref 0. in
    for i = 1 to d do
      let diff = vec.{i} -. mat.{i, col} in
      r2 := !r2 +. diff *. diff
    done;
    eval_rbf2 !r2

  let eval_mat_cols mat1 col1 mat2 col2 =
    let d = Mat.dim1 mat1 in
    let r2 = ref 0. in
    for i = 1 to d do
      let diff = mat1.{i, col1} -. mat2.{i, col2} in
      r2 := !r2 +. diff *. diff
    done;
    eval_rbf2 !r2
end

module Gauss_mat = Make_from_vec (Gauss_vec)

module Kernel = struct
  type t = unit

  include Gauss_mat
end
