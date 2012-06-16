(* File: cov_const.ml

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
open Lacaml.D

module Params = struct type t = { log_theta : float } end

module Eval = struct
  module Kernel = struct
    type params = Params.t
    type t = { params : params; const : float }

    let create params = { params; const = exp (-2. *. params.Params.log_theta) }

    let get_params k = k.params
  end

  module Inducing = struct
    type t = int

    let calc_upper k m = Mat.make m m k.Kernel.const
    let get_n_points m = m
  end

  module Input = struct
    type t = unit

    let eval k () m = Vec.make m k.Kernel.const
    let weighted_eval k () _ ~coeffs = k.Kernel.const *. Vec.sum coeffs
    let eval_one k () = k.Kernel.const
  end

  module Inputs = struct
    type t = int

    let create = Array.length
    let get_n_points n = n
    let choose_subset _inputs indexes = Bigarray.Array1.dim indexes
    let create_inducing _kernel n = n

    let create_default_kernel_params _inputs ~n_inducing:_ =
      { Params.log_theta = 0. }

    let calc_upper = Inducing.calc_upper
    let calc_diag k n = Vec.make n k.Kernel.const
    let calc_cross k ~inputs:n ~inducing:m = Mat.make n m k.Kernel.const

    let weighted_eval k ~inputs:_ ~inducing:_ ~coeffs =
      let res = copy coeffs in
      scal k.Kernel.const res;
      res
  end
end

module Deriv = struct
  module Eval = Eval

  module Hyper = struct
    type t = [ `Log_theta ]

    let get_all _kernel _inducing _inputs = [| `Log_theta |]

    let get_value { Eval.Kernel.params } _inducing _inputs = function
      | `Log_theta -> params.Params.log_theta

    let set_values kernel inducing inputs hypers values =
      let { Eval.Kernel.params } = kernel in
      let log_theta_ref = ref params.Params.log_theta in
      let kernel_changed_ref = ref false in
      for i = 1 to Array.length hypers do
        match hypers.(i - 1) with
        | `Log_theta -> log_theta_ref := values.{i}; kernel_changed_ref := true
      done;
      let new_kernel =
        if !kernel_changed_ref then
          Eval.Kernel.create { Params.log_theta = !log_theta_ref }
        else kernel
      in
      new_kernel, inducing, inputs
  end

  let calc_const_deriv k = -2. *. k.Eval.Kernel.const

  module Inducing = struct
    type upper = { m : int; deriv_const : float }

    let calc_shared_upper k m =
      Eval.Inducing.calc_upper k m, { m; deriv_const = calc_const_deriv k }

    let calc_deriv_upper shared `Log_theta = `Const shared.deriv_const
  end

  module Inputs = struct
    type diag = { diag_eval_inputs : Eval.Inputs.t; diag_const_deriv : float }
    type cross = { cross_const_deriv : float }

    let calc_shared_diag k diag_eval_inputs =
      (
        Eval.Inputs.calc_diag k diag_eval_inputs,
        { diag_eval_inputs; diag_const_deriv = calc_const_deriv k }
      )

    let calc_shared_cross k ~inputs ~inducing =
      (
        Eval.Inputs.calc_cross k ~inputs ~inducing,
        { cross_const_deriv = calc_const_deriv k }
      )

    let calc_deriv_diag diag `Log_theta = `Const diag.diag_const_deriv
    let calc_deriv_cross cross `Log_theta = `Const cross.cross_const_deriv
  end
end
