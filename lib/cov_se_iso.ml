open Printf
open Bigarray
open Lacaml.Impl.D

open Interfaces
open Utils

module Params = struct type t = { log_ell : float; log_sf2 : float } end

type inducing_hyper = { ind : int; dim : int }

module Eval = struct
  module Kernel = struct
    type params = Params.t

    type t = {
      params : params;
      inv_ell2 : float;
      inv_ell2_05 : float;
      log_sf2 : float;
      sf2 : float;
    }

    let create params =
      let log_sf2 = params.Params.log_sf2 in
      let inv_ell2 = exp (-2. *. params.Params.log_ell) in
      {
        params = params;
        inv_ell2 = inv_ell2;
        inv_ell2_05 = -0.5 *. inv_ell2;
        log_sf2 = log_sf2;
        sf2 = exp log_sf2;
      }

    let get_params k = k.params
  end

  open Kernel

  module Inducing = struct
    type t = mat

    module Prepared = struct
      type upper = { sq_diff_mat : mat; ssqr_inducing : vec; inducing : t }

      let calc_upper inducing =
        let sq_diff_mat = syrk ~trans:`T inducing in
        let m = Mat.dim2 sq_diff_mat in
        let ssqr_inducing = Vec.create m in
        for i = 1 to m do ssqr_inducing.{i} <- sq_diff_mat.{i, i} done;
        for c = 2 to m do
          let ssqr_inducing_c = ssqr_inducing.{c} in
          for r = 1 to c - 1 do
            sq_diff_mat.{r, c} <-
              ssqr_inducing.{r} -. 2. *. sq_diff_mat.{r, c} +. ssqr_inducing_c
          done
        done;
        for i = 1 to m do sq_diff_mat.{i, i} <- 0. done;
        {
          sq_diff_mat = sq_diff_mat;
          ssqr_inducing = ssqr_inducing;
          inducing = inducing;
        }
    end

    let calc_upper k { Prepared.sq_diff_mat = sq_diff_mat } =
      let m = Mat.dim2 sq_diff_mat in
      let res = Mat.create m m in
      let { inv_ell2_05 = inv_ell2_05; log_sf2 = log_sf2; sf2 = sf2 } = k in
      for c = 2 to m do
        for r = 1 to c - 1 do
          res.{r, c} <- exp (log_sf2 +. inv_ell2_05 *. sq_diff_mat.{r, c});
        done
      done;
      for i = 1 to m do res.{i, i} <- sf2 done;
      res
  end

  module Input = struct
    type t = vec

    module Prepared = struct
      type cross = { sq_diff_vec : vec; inducing : Inducing.t }

      let calc_cross inducing_prepared input =
        let
          {
            Inducing.Prepared.ssqr_inducing = ssqr_inducing;
            inducing = inducing
          } = inducing_prepared
        in
        let sq_diff_vec = gemv ~trans:`T inducing input in
        let ssqr_input = Vec.sqr_nrm2 input in
        for i = 1 to Vec.dim sq_diff_vec do
          sq_diff_vec.{i} <-
            ssqr_inducing.{i} -. 2. *. sq_diff_vec.{i} +. ssqr_input
        done;
        { sq_diff_vec = sq_diff_vec; inducing = inducing }
    end

    let eval k prepared =
      let sq_diff_vec = prepared.Prepared.sq_diff_vec in
      let m = Vec.dim sq_diff_vec in
      let res = Vec.create m in
      let { Kernel.inv_ell2_05 = inv_ell2_05; log_sf2 = log_sf2 } = k in
      for i = 1 to m do
        res.{i} <- exp (log_sf2 +. inv_ell2_05 *. sq_diff_vec.{i})
      done;
      res

    let weighted_eval k ~coeffs prepared =
      let sq_diff_vec = prepared.Prepared.sq_diff_vec in
      let m = Vec.dim sq_diff_vec in
      if Vec.dim coeffs <> m then
        failwith "Gpr.Cov_se_iso.Eval.Input.weighted_eval: dim(coeffs) <> m";
      let { Kernel.inv_ell2_05 = inv_ell2_05; log_sf2 = log_sf2 } = k in
      let rec loop acc i =
        if i = 0 then acc
        else
          let new_acc =
            coeffs.{i} *. exp (log_sf2 +. inv_ell2_05 *. sq_diff_vec.{i})
          in
          loop new_acc (i - 1)
      in
      loop 0. m

    let eval_one k _input = k.Kernel.sf2
  end

  module Inputs = struct
    type t = mat

    let get_n_inputs = Mat.dim2
    let choose_subset inputs indexes = choose_cols inputs indexes
    let create_inducing _kernel inputs = inputs

    let create_default_kernel_params _inputs =
      { Params.log_ell = 0.; log_sf2 = 0. }

    module Prepared = struct
      type cross = {
        sq_diff_mat : mat;
        inducing : Inducing.t;
        inputs : t;
      }

      let calc_cross inducing_prepared inputs =
        let
          {
            Inducing.Prepared.
            ssqr_inducing = ssqr_inducing;
            inducing = inducing;
          } = inducing_prepared
        in
        let m = Mat.dim2 inducing in
        let n = Mat.dim2 inputs in
        let sq_diff_mat = gemm ~transa:`T inducing inputs in
        let ssqr_inputs = Vec.create n in
        for i = 1 to n do
          (* TODO: optimize sqr_nrm2 and col *)
          ssqr_inputs.{i} <- Vec.sqr_nrm2 (Mat.col inputs i)
        done;
        for c = 1 to n do
          let ssqr_inputs_c = ssqr_inputs.{c} in
          for r = 1 to m do
            sq_diff_mat.{r, c} <-
              ssqr_inducing.{r} -. 2. *. sq_diff_mat.{r, c} +. ssqr_inputs_c
          done
        done;
        { sq_diff_mat = sq_diff_mat; inducing = inducing; inputs = inputs }
    end

    let calc_upper k inputs =
      Inducing.calc_upper k (Inducing.Prepared.calc_upper inputs)

    let calc_diag k inputs =
      (* TODO: return sparse data? *)
      Vec.make (Mat.dim2 inputs) k.Kernel.sf2

    let calc_cross k cross =
      let { Kernel.inv_ell2_05 = inv_ell2_05; log_sf2 = log_sf2 } = k in
      let sq_diff_mat = cross.Prepared.sq_diff_mat in
      let m = Mat.dim1 sq_diff_mat in
      let n = Mat.dim2 sq_diff_mat in
      let res = Mat.create m n in
      for c = 1 to n do
        for r = 1 to m do
          res.{r, c} <- exp (log_sf2 +. inv_ell2_05 *. sq_diff_mat.{r, c})
        done
      done;
      res

    let weighted_eval k ~coeffs cross =
      let sq_diff_mat = cross.Prepared.sq_diff_mat in
      let m = Mat.dim1 sq_diff_mat in
      let n = Mat.dim2 sq_diff_mat in
      if Vec.dim coeffs <> m then
        failwith "Gpr.Cov_se_iso.Eval.Inputs.weighted_eval: dim(coeffs) <> m";
      let { Kernel.inv_ell2_05 = inv_ell2_05; log_sf2 = log_sf2 } = k in
      let rec loop c acc r =
        if r = 0 then acc
        else
          let el =
            coeffs.{r} *. exp (log_sf2 +. inv_ell2_05 *. sq_diff_mat.{r, c})
          in
          loop c (acc +. el) (r - 1)
      in
      Vec.init n (fun c -> loop c 0. m)
  end
end

module Deriv = struct
  module Eval = Eval

  type gen_deriv = [ `Log_ell | `Log_sf2 ]

  module Hyper = struct
    type t = [ gen_deriv | `Inducing_hyper of inducing_hyper ]

    let extract { Eval.Kernel.params = params } =
      let values = Vec.create 2 in
      values.{1} <- params.Params.log_ell;
      values.{2} <- params.Params.log_sf2;
      [| `Log_ell; `Log_sf2 |], values

    let update _kernel (values : vec) =
      Eval.Kernel.create { Params.log_ell = values.{1}; log_sf2 = values.{2} }
  end

  type deriv_common = {
    kernel : Eval.Kernel.t;
    sq_diff_mat : mat;
    eval_mat : mat;
  }

  let calc_gen_deriv
        ({ sq_diff_mat = sq_diff_mat; eval_mat = eval_mat } as sh) =
    function
    | `Log_sf2 -> `Factor 1.
    | `Log_ell ->
        let m = Mat.dim1 sq_diff_mat in
        let n = Mat.dim2 sq_diff_mat in
        let res = Mat.create m n in
        let { Eval.Kernel.inv_ell2 = inv_ell2 } = sh.kernel in
        (* TODO: Lacaml-version of element-wise multiplication? *)
        for c = 1 to n do
          for r = 1 to m do
            res.{r, c} <- eval_mat.{r, c} *. sq_diff_mat.{r, c} *. inv_ell2
          done
        done;
        `Dense res

  module Inducing = struct
    module Prepared = struct
      type upper = Eval.Inducing.Prepared.upper

      let calc_upper upper = upper
    end

    type upper = Eval.Inducing.t * deriv_common

    let calc_shared_upper kernel prepared_upper =
      let module EI = Eval.Inducing in
      let module EIP = EI.Prepared in
      let upper = EI.calc_upper kernel prepared_upper in
      let shared =
        (
          prepared_upper.EIP.inducing,
          {
            kernel = kernel;
            sq_diff_mat = prepared_upper.EIP.sq_diff_mat;
            eval_mat = upper;
          }
        )
      in
      upper, shared

    let calc_deriv_upper (inducing, common) = function
      | #gen_deriv as hyper -> calc_gen_deriv common hyper
      | `Inducing_hyper inducing_hyper ->
          let { ind = ind; dim = dim } = inducing_hyper in
          let eval_mat = common.eval_mat in
          let m = Mat.dim2 eval_mat in
          let res = Mat.create 1 m in
          let indx_d = inducing.{dim, ind} in
          let inv_ell2 = common.kernel.Eval.Kernel.inv_ell2 in
          for i = 1 to ind - 1 do
            let ind_d = inducing.{dim, i} in
            res.{1, i} <- inv_ell2 *. (ind_d -. indx_d) *. eval_mat.{i, ind}
          done;
          res.{1, ind} <- 0.;
          for i = ind + 1 to m do
            let ind_d = inducing.{dim, i} in
            res.{1, i} <- inv_ell2 *. (ind_d -. indx_d) *. eval_mat.{ind, i}
          done;
          let rows = Sparse_indices.create 1 in
          rows.{1} <- ind;
          `Sparse_rows (res, rows)
  end

  module Inputs = struct
    module Prepared = struct
      type cross = Eval.Inputs.Prepared.cross

      let calc_cross _upper cross = cross
    end

    type diag = Eval.Kernel.t

    type cross = Eval.Inducing.t * Eval.Inputs.t * deriv_common

    let calc_shared_diag k diag_eval_inputs =
      Eval.Inputs.calc_diag k diag_eval_inputs, k

    let calc_shared_cross kernel prepared_cross =
      let module EI = Eval.Inputs in
      let module EIP = EI.Prepared in
      let cross = EI.calc_cross kernel prepared_cross in
      let shared =
        (
          prepared_cross.EIP.inducing,
          prepared_cross.EIP.inputs,
          {
            kernel = kernel;
            sq_diff_mat = prepared_cross.EIP.sq_diff_mat;
            eval_mat = cross;
          }
        )
      in
      cross, shared

    let calc_deriv_diag _diag = function
      | `Log_sf2 -> `Factor 1.
      | `Log_ell | `Inducing_hyper _ -> `Const 0.

    let calc_deriv_cross (inducing, inputs, common) = function
      | #gen_deriv as hyper -> calc_gen_deriv common hyper
      | `Inducing_hyper inducing_hyper ->
          let { ind = ind; dim = dim } = inducing_hyper in
          let eval_mat = common.eval_mat in
          let n = Mat.dim2 eval_mat in
          let res = Mat.create 1 n in
          let indx_d = inducing.{dim, ind} in
          let inv_ell2 = common.kernel.Eval.Kernel.inv_ell2 in
          for c = 1 to n do
            let inp_d = inputs.{dim, c} in
            res.{1, c} <- inv_ell2 *. (inp_d -. indx_d) *. eval_mat.{ind, c}
          done;
          let rows = Sparse_indices.create 1 in
          rows.{1} <- ind;
          `Sparse_rows (res, rows)
  end
end

module SPGP = struct
  module Eval = Eval
  module Deriv = Deriv

  module Inducing_hypers = struct
    let extract inducing =
      let d = Mat.dim1 inducing in
      let n = Mat.dim2 inducing in
      let all = d * n in
      let hypers = Array.create all (`Inducing_hyper { ind = 1; dim = 1 }) in
      for ind = 1 to n do
        let indd = (ind - 1) * d in
        for dim = 1 to d do
          let inducing_hyper = { ind = ind; dim = dim } in
          hypers.(indd + dim - 1) <- (`Inducing_hyper inducing_hyper)
        done
      done;
      hypers, Mat.as_vec inducing

    let update inducing values =
      let gen = genarray_of_array1 values in
      reshape_2 gen (Mat.dim1 inducing) (Mat.dim2 inducing)
  end
end
