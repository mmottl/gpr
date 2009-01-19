open Lacaml.Impl.D
open Lacaml.Io

module Eval = struct
  module Kernel = struct
    type params = < log_theta : float >
    type t = float

    let create params = exp (-2. *. params#log_theta)
  end

  module Inducing = struct
    type t = < m : int >

    module Prepared = struct
      type upper = t

      let calc_upper points = points
    end

    let calc_upper k upper =
      let m = upper#m in
      (* TODO: build upper triangle only *)
      Mat.make m m k
  end

  module Input = struct
    type t = unit

    module Prepared = struct
      type cross = int

      let calc_cross upper _input = upper#m
    end

    let eval k m = Vec.make m k

    let weighted_eval k ~coeffs m =
      if Vec.dim coeffs <> m then
        failwith
          "Gpr.Cov_const.Deriv.Eval.Input.weighted_eval: dim(coeffs) <> m";
      k *. Vec.sum coeffs

    let eval_one k () = k
  end

  module Inputs = struct
    type t = < n : int >

    module Prepared = struct
      type cross = { m : int; n : int }

      let calc_cross upper inputs = { m = upper#m; n = inputs#n }
    end

    let calc_upper k inputs =
      let m = inputs#n in
      Inducing.calc_upper k (object method m = m end)

    let calc_diag k inputs = Vec.make inputs#n k
    let calc_cross k { Prepared.m = m; n = n } = Mat.make m n k

    let weighted_eval k ~coeffs { Prepared.m = m } =
      if Vec.dim coeffs <> m then
        failwith
          "Gpr.Cov_const.Deriv.Eval.Inputs.weighted_eval: dim(coeffs) <> m";
      let res = copy coeffs in
      scal k res;
      res
  end
end

module Hyper = struct type t = [ `Log_theta ] end

let calc_const_deriv const = -2. *. const

module Inducing = struct
  module Prepared = struct
    type upper = Eval.Inducing.Prepared.upper

    let calc_upper upper = upper
  end

  type shared = { m : Prepared.upper; deriv_const : float }

  let calc_shared_upper k m =
    Eval.Inducing.calc_upper k m, { m = m; deriv_const = calc_const_deriv k }

  let calc_deriv_upper shared `Log_theta = `Const shared.deriv_const
end

module Inputs = struct
  module Prepared = struct
    type cross = Eval.Inputs.Prepared.cross

    let calc_cross upper cross =
      let m = upper#m in
      if m <> cross.Eval.Inputs.Prepared.m then
        failwith
          "Gpr.Cov_const.Deriv.Inputs.Prepared.calc_cross: dimension mismatch";
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
