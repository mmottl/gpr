open Bigarray
open Printf
open Lacaml.Impl.D
open Lacaml.Io

module Params = struct type t = { log_theta : float } end

module Eval = struct
  module Kernel = struct
    type params = Params.t
    type t = { params : params; const : float }

    let create params =
      { params = params; const = exp (-2. *. params.Params.log_theta) }

    let get_params k = k.params
  end

  module Inducing = struct
    type t = int

    let calc_upper k m = Mat.make m m k.Kernel.const
    let get_n_points m = m
  end

  module Input = struct
    type t = unit

    let eval k m () = Vec.make m k.Kernel.const
    let weighted_eval k _ ~coeffs () = k.Kernel.const *. Vec.sum coeffs
    let eval_one k () = k.Kernel.const
  end

  module Inputs = struct
    type t = int

    let get_n_points n = n
    let choose_subset _inputs indexes = Array1.dim indexes
    let create_inducing _kernel n = n

    let create_default_kernel_params ~n_inducing:_ _inputs =
      { Params.log_theta = 0. }

    let calc_upper = Inducing.calc_upper
    let calc_diag k n = Vec.make n k.Kernel.const
    let calc_cross k m n = Mat.make m n k.Kernel.const

    let weighted_eval k _ ~coeffs _ =
      let res = copy coeffs in
      scal k.Kernel.const res;
      res
  end
end

module Deriv = struct
  module Eval = Eval

  module Hyper = struct
    type t = [ `Log_theta ]

    let extract { Eval.Kernel.params = params } =
      let values = Vec.create 1 in
      values.{1} <- params.Params.log_theta;
      [| `Log_theta |], values

    let update _kernel (values : vec) =
      Eval.Kernel.create { Params.log_theta = values.{1} }
  end

  let calc_const_deriv k = -2. *. k.Eval.Kernel.const

  module Inducing = struct
    type upper = { m : int; deriv_const : float }

    let calc_shared_upper k m =
      Eval.Inducing.calc_upper k m, { m = m; deriv_const = calc_const_deriv k }

    let calc_deriv_upper shared `Log_theta = `Const shared.deriv_const
  end

  module Inputs = struct
    type diag = { diag_eval_inputs : Eval.Inputs.t; diag_const_deriv : float }
    type cross = { cross_const_deriv : float }

    let calc_shared_diag k diag_eval_inputs =
      (
        Eval.Inputs.calc_diag k diag_eval_inputs,
        {
          diag_eval_inputs = diag_eval_inputs;
          diag_const_deriv = calc_const_deriv k;
        }
      )

    let calc_shared_cross k eval_inducing eval_inputs =
      (
        Eval.Inputs.calc_cross k eval_inducing eval_inputs,
        { cross_const_deriv = calc_const_deriv k }
      )

    let calc_deriv_diag diag `Log_theta = `Const diag.diag_const_deriv
    let calc_deriv_cross cross `Log_theta = `Const cross.cross_const_deriv
  end
end
