(* File: cov_lin_ard.ml

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
      { params; consts }

    let get_params k = k.params
  end

  module Inducing = struct
    type t = mat

    let get_n_points = Mat.dim2
    let calc_upper _k inducing = syrk ~trans:`T inducing
  end

  module Input = struct
    type t = vec

    let calc_ard_input { Kernel.consts } input =
      let d = Vec.dim input in
      let ard_input = Vec.create d in
      for i = 1 to d do ard_input.{i} <- consts.{i} *. input.{i} done;
      ard_input

    let eval k input inducing =
      gemv ~trans:`T inducing (calc_ard_input k input)

    let weighted_eval k input inducing ~coeffs =
      dot ~x:coeffs (eval k input inducing)

    let eval_one { Kernel.consts } input =
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

    let create = Mat.of_col_vecs
    let get_n_points = Mat.dim2
    let choose_subset inputs indexes = Gpr_utils.choose_cols inputs indexes

    let calc_ard_inputs { Kernel.consts } inputs =
      let res = lacpy inputs in
      Mat.scal_rows consts res;
      res

    let create_inducing = calc_ard_inputs

    let create_default_kernel_params inputs ~n_inducing:_ =
      { Params.log_ells = Vec.make (Mat.dim1 inputs) 0. }

    let calc_upper k inputs = syrk ~trans:`T (calc_ard_inputs k inputs)
    let calc_diag k inputs = Mat.syrk_diag ~trans:`T (calc_ard_inputs k inputs)

    let calc_cross k ~inputs ~inducing =
      gemm ~transa:`T (calc_ard_inputs k inputs) inducing

    let weighted_eval k ~inputs ~inducing ~coeffs =
      gemv (calc_cross k ~inputs ~inducing) coeffs
  end
end

module Deriv = struct
  module Eval = Eval

  module Hyper = struct
    type t = [ `Log_ell of int ]

    let get_all { Eval.Kernel.params } _inducing _inputs =
      Array.init (Vec.dim params.Params.log_ells) ~f:(fun d ->
        `Log_ell (d + 1))

    let get_value { Eval.Kernel.params } _inducing _inputs = function
      | `Log_ell i -> params.Params.log_ells.{i}

    let set_values k inducing inputs hypers values =
      let { Eval.Kernel.params } = k in
      let log_ells_lazy = lazy (copy params.Params.log_ells) in
      for i = 1 to Array.length hypers do
        match hypers.(i - 1) with
        | `Log_ell d -> (Lazy.force log_ells_lazy).{d} <- values.{i}
      done;
      let new_kernel =
        if Lazy.is_val log_ells_lazy then
          Eval.Kernel.create { Params.log_ells = Lazy.force log_ells_lazy }
        else k
      in
      new_kernel, inducing, inputs
  end

  module Inducing = struct
    type upper = Eval.Inducing.t

    let calc_shared_upper k eval_inducing =
      let upper = Eval.Inducing.calc_upper k eval_inducing in
      upper, eval_inducing

    let calc_deriv_upper _inducing = function
      | `Log_ell _ -> `Const 0.
  end

  module Inputs = struct
    type diag = Eval.Kernel.t * Eval.Inputs.t
    type cross = Eval.Kernel.t * Eval.Inputs.t* Eval.Inducing.t

    let calc_shared_diag k eval_inputs =
      Eval.Inputs.calc_diag k eval_inputs, (k, eval_inputs)

    let calc_shared_cross k ~inputs ~inducing =
      (
        Eval.Inputs.calc_cross k ~inputs ~inducing,
        (k, inputs, inducing)
      )

    let calc_deriv_diag (k, inputs) (`Log_ell d) =
      let n = Mat.dim2 inputs in
      let res = Vec.create n in
      let const = -2. *. k.Eval.Kernel.consts.{d} in
      for i = 1 to n do
        let el = inputs.{d, i} in
        res.{i} <- (const *. el) *. el
      done;
      `Vec res

    let calc_deriv_cross (k, inputs, inducing) (`Log_ell d) =
      let m = Mat.dim2 inducing in
      let n = Mat.dim2 inputs in
      let res = Mat.create n m in
      let const = -. k.Eval.Kernel.consts.{d} in
      for c = 1 to m do
        for r = 1 to n do
          res.{r, c} <- const *. inducing.{d, c} *. inputs.{d, r}
        done
      done;
      `Dense res
  end
end
