(* File: cov_se_iso.ml

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

open Interfaces

open Core.Std
open Lacaml.D

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

    let create ({ Params.log_sf2; log_ell } as params) =
      let inv_ell2 = exp (-2. *. log_ell) in
      let inv_ell2_05 = -0.5 *. inv_ell2 in
      { params; inv_ell2; inv_ell2_05; log_sf2; sf2 = exp log_sf2 }

    let get_params k = k.params
  end

  open Kernel

  module Inducing = struct
    type t = mat

    let get_n_points = Mat.dim2

    let calc_sqr_diff_mat inducing =
      let d = Mat.dim1 inducing in
      let m = Mat.dim2 inducing in
      let res = Mat.create m m in
      let ssqr_diff_ref = ref 0. in
      for c = 1 to m do
        for r = 1 to c - 1 do
          for i = 1 to d do
            let diff = inducing.{i, c} -. inducing.{i, r} in
            ssqr_diff_ref := !ssqr_diff_ref +. diff *. diff
          done;
          res.{r, c} <- !ssqr_diff_ref;
          ssqr_diff_ref := 0.
        done;
        res.{c, c} <- 0.
      done;
      res

    let calc_upper_with_sqr_diff_mat k sqr_diff_mat =
      let m = Mat.dim2 sqr_diff_mat in
      let res = Mat.create m m in
      let { inv_ell2_05; log_sf2; sf2 } = k in
      for c = 1 to m do
        for r = 1 to c - 1 do
          res.{r, c} <- exp (log_sf2 +. inv_ell2_05 *. sqr_diff_mat.{r, c});
        done;
        res.{c, c} <- sf2;
      done;
      res

    let calc_upper k inducing =
      calc_upper_with_sqr_diff_mat k (calc_sqr_diff_mat inducing)
  end

  module Input = struct
    type t = vec

    let eval { Kernel.inv_ell2_05; log_sf2 } input inducing =
      let d = Mat.dim1 inducing in
      let m = Mat.dim2 inducing in
      let res = Vec.create m in
      let ssqr_diff_ref = ref 0. in
      for c = 1 to m do
        for i = 1 to d do
          let diff = input.{i} -. inducing.{i, c} in
          ssqr_diff_ref := !ssqr_diff_ref +. diff *. diff
        done;
        res.{c} <- exp (log_sf2 +. inv_ell2_05 *. !ssqr_diff_ref);
        ssqr_diff_ref := 0.;
      done;
      res

    let weighted_eval k input inducing ~coeffs =
      dot ~x:coeffs (eval k input inducing)

    let eval_one k _input = k.Kernel.sf2
  end

  module Inputs = struct
    type t = mat

    let create = Mat.of_col_vecs
    let get_n_points = Mat.dim2
    let choose_subset inputs indexes = Gpr_utils.choose_cols inputs indexes
    let create_inducing _kernel inputs = inputs

    let create_default_kernel_params _inputs ~n_inducing:_ =
      { Params.log_ell = 0.; log_sf2 = 0. }

    let calc_upper k inputs = Inducing.calc_upper k inputs
    let calc_diag k inputs = Vec.make (Mat.dim2 inputs) k.Kernel.sf2

    let calc_sqr_diff_mat ~inputs ~inducing =
      let d = Mat.dim1 inducing in
      let m = Mat.dim2 inducing in
      let n = Mat.dim2 inputs in
      let res = Mat.create n m in
      let ssqr_diff_ref = ref 0. in
      for c = 1 to m do
        for r = 1 to n do
          for i = 1 to d do
            let diff = inputs.{i, r} -. inducing.{i, c} in
            ssqr_diff_ref := !ssqr_diff_ref +. diff *. diff
          done;
          res.{r, c} <- !ssqr_diff_ref;
          ssqr_diff_ref := 0.
        done
      done;
      res

    let calc_cross_with_sqr_diff_mat k sqr_diff_mat =
      let { Kernel.inv_ell2_05; log_sf2 } = k in
      let n = Mat.dim1 sqr_diff_mat in
      let m = Mat.dim2 sqr_diff_mat in
      let res = Mat.create n m in
      for c = 1 to m do
        for r = 1 to n do
          res.{r, c} <- exp (log_sf2 +. inv_ell2_05 *. sqr_diff_mat.{r, c})
        done
      done;
      res

    let calc_cross k ~inputs ~inducing =
      calc_cross_with_sqr_diff_mat k (calc_sqr_diff_mat ~inputs ~inducing)

    let weighted_eval k ~inputs ~inducing ~coeffs =
      let sqr_diff_mat = calc_sqr_diff_mat ~inputs ~inducing in
      let n = Mat.dim1 sqr_diff_mat in
      let m = Mat.dim2 sqr_diff_mat in
      if Vec.dim coeffs <> m then
        failwith "Gpr.Cov_se_iso.Eval.Inputs.weighted_eval: dim(coeffs) <> m";
      let { Kernel.inv_ell2_05; log_sf2 } = k in
      let rec loop r acc c =
        if c = 0 then acc
        else
          let el =
            coeffs.{c} *. exp (log_sf2 +. inv_ell2_05 *. sqr_diff_mat.{r, c})
          in
          loop r (acc +. el) (c - 1)
      in
      Vec.init n (fun r -> loop r 0. m)
  end
end

module Deriv = struct
  module Eval = Eval

  type gen_deriv = [ `Log_ell | `Log_sf2 ]

  module Hyper = struct
    type t = [ gen_deriv | `Inducing_hyper of inducing_hyper ]

    let get_all _kernel inducing _inputs =
      let d = Mat.dim1 inducing in
      let m = Mat.dim2 inducing in
      let n_inducing_hypers = d * m in
      let n_all_hypers = 2 + n_inducing_hypers in
      let hypers = Array.create ~len:n_all_hypers `Log_ell in
      hypers.(1) <- `Log_sf2 ;
      for ind = 1 to m do
        let indd = (ind - 1) * d in
        for dim = 1 to d do
          let inducing_hyper = { ind; dim } in
          hypers.(1 + indd + dim) <- `Inducing_hyper inducing_hyper
        done
      done;
      hypers

    let get_value { Eval.Kernel.params } inducing _inputs = function
      | `Log_ell -> params.Params.log_ell
      | `Log_sf2 -> params.Params.log_sf2
      | `Inducing_hyper { ind; dim } -> inducing.{dim, ind}

    let set_values { Eval.Kernel.params } inducing inputs hypers values =
      let { Params.log_ell; log_sf2 } = params in
      let log_ell_ref = ref log_ell in
      let log_sf2_ref = ref log_sf2 in
      let inducing_lazy = lazy (lacpy inducing) in
      for i = 1 to Array.length hypers do
        match hypers.(i - 1) with
        | `Log_ell -> log_ell_ref := values.{i}
        | `Log_sf2 -> log_sf2_ref := values.{i}
        | `Inducing_hyper { ind; dim } ->
            (Lazy.force inducing_lazy).{dim, ind} <- values.{i}
      done;
      let new_kernel =
        let log_ell = !log_ell_ref in
        Eval.Kernel.create { Params.log_ell; log_sf2 = !log_sf2_ref }
      in
      let lift lazy_value value =
        if Lazy.is_val lazy_value then Lazy.force lazy_value
        else value
      in
      let new_inducing = lift inducing_lazy inducing in
      new_kernel, new_inducing, inputs
  end

  type deriv_common = {
    kernel : Eval.Kernel.t;
    sqr_diff_mat : mat;
    eval_mat : mat;
  }

  module Inducing = struct
    type upper = Eval.Inducing.t * deriv_common

    let calc_shared_upper kernel eval_inducing =
      let module EI = Eval.Inducing in
      let sqr_diff_mat = EI.calc_sqr_diff_mat eval_inducing in
      let eval_mat = EI.calc_upper_with_sqr_diff_mat kernel sqr_diff_mat in
      eval_mat, (eval_inducing, { kernel; sqr_diff_mat; eval_mat })

    let calc_deriv_upper (inducing, common) = function
      | `Log_sf2 -> `Factor 1.
      | `Log_ell ->
          let { sqr_diff_mat; eval_mat; kernel } = common in
          let m = Mat.dim1 sqr_diff_mat in
          let res = Mat.create m m in
          let { Eval.Kernel.inv_ell2 } = kernel in
          for c = 1 to m do
            for r = 1 to c - 1 do
              res.{r, c} <- eval_mat.{r, c} *. sqr_diff_mat.{r, c} *. inv_ell2
            done;
            res.{c, c} <- 0.;
          done;
          `Dense res
        | `Inducing_hyper { ind; dim } ->
            let eval_mat = common.eval_mat in
            let m = Mat.dim2 eval_mat in
            let res = Mat.create 1 m in
            let inducing_dim = inducing.{dim, ind} in
            let inv_ell2 = common.kernel.Eval.Kernel.inv_ell2 in
            for i = 1 to ind - 1 do
              let ind_d = inducing.{dim, i} in
              res.{1, i} <-
                inv_ell2 *. (ind_d -. inducing_dim) *. eval_mat.{i, ind}
            done;
            res.{1, ind} <- 0.;
            for i = ind + 1 to m do
              let ind_d = inducing.{dim, i} in
              res.{1, i} <-
                inv_ell2 *. (ind_d -. inducing_dim) *. eval_mat.{ind, i}
            done;
            let rows = Sparse_indices.create 1 in
            rows.{1} <- ind;
            `Sparse_rows (res, rows)
  end

  module Inputs = struct
    type diag = Eval.Kernel.t

    type cross = Eval.Inputs.t * Eval.Inducing.t * deriv_common

    let calc_shared_diag k diag_eval_inputs =
      Eval.Inputs.calc_diag k diag_eval_inputs, k

    let calc_shared_cross kernel ~inputs ~inducing =
      let module EI = Eval.Inputs in
      let sqr_diff_mat = EI.calc_sqr_diff_mat ~inputs ~inducing in
      let eval_mat = EI.calc_cross_with_sqr_diff_mat kernel sqr_diff_mat in
      let shared = inputs, inducing, { kernel; sqr_diff_mat; eval_mat } in
      eval_mat, shared

    let calc_deriv_diag _diag = function
      | `Log_sf2 -> `Factor 1.
      | `Log_ell | `Inducing_hyper _ -> `Const 0.

    let calc_deriv_cross (inputs, inducing, common) = function
      | `Log_sf2 -> `Factor 1.
      | `Log_ell ->
          let { sqr_diff_mat; eval_mat; kernel } = common in
          let n = Mat.dim1 sqr_diff_mat in
          let m = Mat.dim2 sqr_diff_mat in
          let res = Mat.create n m in
          let { Eval.Kernel.inv_ell2 } = kernel in
          for c = 1 to m do
            for r = 1 to n do
              res.{r, c} <- eval_mat.{r, c} *. sqr_diff_mat.{r, c} *. inv_ell2
            done
          done;
          `Dense res
      | `Inducing_hyper { ind; dim } ->
          let eval_mat = common.eval_mat in
          let n = Mat.dim1 eval_mat in
          let res = Mat.create n 1 in
          let indx_d = inducing.{dim, ind} in
          let inv_ell2 = common.kernel.Eval.Kernel.inv_ell2 in
          for r = 1 to n do
            let inp_d = inputs.{dim, r} in
            res.{r, 1} <- inv_ell2 *. (inp_d -. indx_d) *. eval_mat.{r, ind}
          done;
          let cols = Sparse_indices.create 1 in
          cols.{1} <- ind;
          `Sparse_cols (res, cols)
  end
end
