open Printf
open Lacaml.Impl.D
open Lacaml.Io

open Utils

module Params = struct type t = { log_ells : vec } end

module Eval = struct
  module Kernel = struct
    type params = Params.t
    type t = { params : params; consts : vec }

    let create params =
      let log_ells = params.Params.log_ells in
      let d = Vec.dim log_ells in
      let consts = Vec.create d in
      for i = 1 to d do consts.{i} <- exp (-. log_ells.{i}) done;
      { params = params; consts = consts }

    let get_params k = k.params
  end

  module Inducing = struct
    type t = mat

    let get_n_points = Mat.dim2
    let calc_upper _k inducing = syrk ~trans:`T inducing
  end

  module Input = struct
    type t = vec

    let calc_ard_input { Kernel.consts = consts } input =
      let d = Vec.dim input in
      let ard_input = Vec.create d in
      for i = 1 to d do ard_input.{i} <- consts.{i} *. input.{i} done;
      ard_input

    let eval k inducing input =
      gemv ~trans:`T inducing (calc_ard_input k input)

    let weighted_eval k inducing ~coeffs input =
      dot ~x:coeffs (eval k inducing input)

    let eval_one { Kernel.consts = consts } input =
      let rec loop res i =
        if i = 0 then res
        else
          let x = consts.{i} *. input.{i} in
          loop (res +. x *. x) (i - 1)
      in
      loop 0. (Vec.dim input)
  end

  module Inputs = struct
    type t = mat

    let get_n_points = Mat.dim2
    let choose_subset inputs indexes = choose_cols inputs indexes

    let create_inducing { Kernel.consts = consts } inputs =
      let m = Mat.dim1 inputs in
      let n = Mat.dim2 inputs in
      let res = Mat.create m n in
      (* TODO: implement Mat.scal_rows in Lacaml *)
      for r = 1 to m do
        let const = consts.{r} in
        for c = 1 to n do res.{r, c} <- const *. inputs.{r, c} done;
      done;
      res

    let create_default_kernel_params ~n_inducing:_ inputs =
      { Params.log_ells = Vec.make (Mat.dim1 inputs) 0. }

    let calc_ard_inputs { Kernel.consts = consts } inputs =
      let d = Mat.dim1 inputs in
      let n = Mat.dim2 inputs in
      let ard_inputs = Mat.create d n in
      for c = 1 to n do
        for r = 1 to d do ard_inputs.{r, c} <- consts.{r} *. inputs.{r, c} done
      done;
      ard_inputs

    let calc_upper k inputs = syrk ~trans:`T (calc_ard_inputs k inputs)
    let calc_diag k inputs = Mat.syrk_diag ~trans:`T (calc_ard_inputs k inputs)

    let calc_cross k inducing inputs =
      gemm ~transa:`T inducing (calc_ard_inputs k inputs)

    let weighted_eval k inducing ~coeffs inputs =
      gemv ~trans:`T (calc_cross k inducing inputs) coeffs
  end
end

module Deriv = struct
  module Eval = Eval

  module Hyper = struct
    type t = [ `Log_ell of int ]

    let get_all { Eval.Kernel.params = params } _inducing =
      Array.init (Vec.dim params.Params.log_ells) (fun d ->
        `Log_ell (d + 1))

    let get_value { Eval.Kernel.params = params } _inducing = function
      | `Log_ell i -> params.Params.log_ells.{i}

    let set_values kernel inducing hypers values =
      let { Eval.Kernel.params = params } = kernel in
      let log_ells_lazy = lazy (copy params.Params.log_ells) in
      for i = 1 to Array.length hypers do
        match hypers.(i - 1) with
        | `Log_ell d -> (Lazy.force log_ells_lazy).{d} <- values.{i}
      done;
      let new_kernel =
        if Lazy.lazy_is_val log_ells_lazy then
          Eval.Kernel.create { Params.log_ells = Lazy.force log_ells_lazy }
        else kernel
      in
      new_kernel, inducing
  end

  module Inducing = struct
    type upper = Eval.Inducing.t

    let calc_shared_upper k eval_inducing =
      let upper = Eval.Inducing.calc_upper k eval_inducing in
      upper, eval_inducing

    let calc_deriv_upper inducing (`Log_ell d) =
      let m = Mat.dim2 inducing in
      let res = Mat.create m m in
      for c = 1 to m do
        for r = 1 to c do
          let prod = inducing.{d, r} *. inducing.{d, c} in
          res.{r, c} <- -. prod -. prod
        done
      done;
      `Dense res
  end

  module Inputs = struct
    type diag = Eval.Kernel.t * Eval.Inputs.t
    type cross = Eval.Kernel.t * Eval.Inducing.t * Eval.Inputs.t

    let calc_shared_diag k eval_inputs =
      Eval.Inputs.calc_diag k eval_inputs, (k, eval_inputs)

    let calc_shared_cross k eval_inducing eval_inputs =
      (
        Eval.Inputs.calc_cross k eval_inputs eval_inputs,
        (k, eval_inducing, eval_inputs)
      )

    let calc_const k d = let kd = k.Eval.Kernel.consts.{d} in -2. *. kd *. kd

    let calc_deriv_diag (k, inputs) (`Log_ell d) =
      let n = Mat.dim2 inputs in
      let res = Vec.create n in
      let const = calc_const k d in
      for i = 1 to n do
        let el = inputs.{d, i} in
        res.{i} <- (const *. el) *. el
      done;
      `Vec res

    let calc_deriv_cross (k, inducing, inputs) (`Log_ell d) =
      let m = Mat.dim2 inducing in
      let n = Mat.dim2 inputs in
      let res = Mat.create m n in
      let const = calc_const k d in
      for c = 1 to n do
        for r = 1 to m do
          res.{r, c} <- const *. inducing.{d, r} *. inputs.{d, c}
        done
      done;
      `Dense res
  end
end
