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

    module Prepared = struct
      type upper = t

      let calc_upper points = points
    end

    (* TODO: build upper triangle only *)
    let calc_upper k m = Mat.make m m k.Kernel.const
  end

  module Input = struct
    type t = unit

    module Prepared = struct
      type cross = int

      let calc_cross m _input = m
    end

    let eval k m = Vec.make m k.Kernel.const

    let weighted_eval k ~coeffs m =
      if Vec.dim coeffs <> m then
        failwith "Gpr.Cov_const.Eval.Input.weighted_eval: dim(coeffs) <> m";
      k.Kernel.const *. Vec.sum coeffs

    let eval_one k () = k.Kernel.const
  end

  module Inputs = struct
    type t = int

    let get_n_inputs n = n
    let choose_subset _inputs indexes = Array1.dim indexes
    let create_inducing _kernel n = n
    let create_default_kernel_params _inputs = { Params.log_theta = 0. }

    module Prepared = struct
      type cross = { m : int; n : int }

      let calc_cross m n = { m = m; n = n }
    end

    let calc_upper = Inducing.calc_upper

    let calc_diag k n = Vec.make n k.Kernel.const
    let calc_cross k { Prepared.m = m; n = n } = Mat.make m n k.Kernel.const

    let weighted_eval k ~coeffs { Prepared.m = m } =
      if Vec.dim coeffs <> m then
        failwith "Gpr.Cov_const.Eval.Inputs.weighted_eval: dim(coeffs) <> m";
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
    module Prepared = struct
      type upper = Eval.Inducing.Prepared.upper

      let calc_upper upper = upper
    end

    type upper = { m : Prepared.upper; deriv_const : float }

    let calc_shared_upper k m =
      Eval.Inducing.calc_upper k m, { m = m; deriv_const = calc_const_deriv k }

    let calc_deriv_upper shared `Log_theta = `Const shared.deriv_const
  end

  module Inputs = struct
    module Prepared = struct
      type cross = Eval.Inputs.Prepared.cross

      let calc_cross m cross =
        if m <> cross.Eval.Inputs.Prepared.m then
          failwith
            "Gpr.Cov_const.Deriv.Inputs.Prepared.calc_cross: \
            dimension mismatch";
        { cross with Eval.Inputs.Prepared.m = m }
    end

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

    let calc_shared_cross k prepared_cross =
      (
        Eval.Inputs.calc_cross k prepared_cross,
        { cross_const_deriv = calc_const_deriv k }
      )

    let calc_deriv_diag diag `Log_theta = `Const diag.diag_const_deriv
    let calc_deriv_cross cross `Log_theta = `Const cross.cross_const_deriv
  end
end
