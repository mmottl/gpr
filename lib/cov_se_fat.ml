(* File: cov_se_fat.ml

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
open Gpr_utils

open Core.Std
open Lacaml.D

let option_map ~f = function None -> None | Some v -> Some (f v)
let option_iter ~f = function None -> () | Some v -> f v

module Params = struct
  type params = {
    d : int;
    log_sf2 : float;
    tproj : mat option;
    log_hetero_skedasticity : vec option;
    log_multiscales_m05 : mat option;
  }

  type t = params

  let create (params : params) =
    let check v_dim name v =
      let n = v_dim v in
      if n <> params.d then
        failwithf
          "Cov_se_fat.Params.create: %s projection (%d) disagrees \
          with target dimension d (%d)" name n params.d ()
    in
    option_iter params.tproj ~f:(check Mat.dim2 "tproj");
    params
end

module Eval = struct
  module Kernel = struct
    type params = Params.t

    type t = {
      params : params;
      sf2 : float;
      hetero_skedasticity : vec option;
      multiscales : mat option;
    }

    let create params =
      let hetero_skedasticity =
        option_map params.Params.log_hetero_skedasticity ~f:(Vec.map exp)
      in
      let multiscales =
        let f v = exp v +. 0.5 in
        option_map params.Params.log_multiscales_m05 ~f:(Mat.map f)
      in
      {
        params;
        sf2 = exp params.Params.log_sf2;
        hetero_skedasticity;
        multiscales;
      }

    let get_params k = k.params
  end

  let calc_res_el ~log_sf2 tmp =
    let x = tmp.x in
    tmp.x <- 0.;
    exp (log_sf2 -. 0.5 *. x)

  let calc_upper_vanilla k mat =
    let { Kernel.sf2; params = { Params.d; log_sf2 } } = k in
    let n = Mat.dim2 mat in
    let res = Mat.create n n in
    let tmp = { x = 0. } in
    for c = 1 to n do
      for r = 1 to c - 1 do
        for i = 1 to d do
          let diff = mat.{i, r} -. mat.{i, c} in
          tmp.x <- tmp.x +. diff *. diff
        done;
        res.{r, c} <- calc_res_el ~log_sf2 tmp;
      done;
      res.{c, c} <- sf2;
    done;
    res

  let update_tmp_sum ~tmp ~diff ~scale =
    tmp.x <- tmp.x +. diff *. (diff /. scale) +. log scale

  module Inducing = struct
    type t = mat

    let get_n_points = Mat.dim2

    let calc_upper k inducing =
      let m = Mat.dim2 inducing in
      let res =
        match k.Kernel.multiscales with
        | None -> calc_upper_vanilla k inducing
        | Some multiscales ->
            let { Kernel.params = { Params.d; log_sf2 } } = k in
            let res = Mat.create m m in
            let tmp = { x = 0. } in
            for c = 1 to m do
              for r = 1 to c - 1 do
                for i = 1 to d do
                  let diff = inducing.{i, r} -. inducing.{i, c} in
                  let scale = multiscales.{i, r} +. multiscales.{i, c} -. 1. in
                  update_tmp_sum ~tmp ~diff ~scale
                done;
                res.{r, c} <- calc_res_el ~log_sf2 tmp
              done;
              for i = 1 to d do
                let multiscale = multiscales.{i, c} in
                tmp.x <- tmp.x +. log (multiscale +. multiscale -. 1.)
              done;
              res.{c, c} <- calc_res_el ~log_sf2 tmp;
            done;
            res
      in
      match k.Kernel.hetero_skedasticity with
      | None -> res
      | Some hetero_skedasticity ->
          for i = 1 to m do
            res.{i, i} <- res.{i, i} +. hetero_skedasticity.{i}
          done;
          res
  end

  module Input = struct
    type t = vec

    let eval k input inducing =
      let
        { Kernel.multiscales; params = { Params.d; log_sf2; tproj } } = k
      in
      let projection =
        match tproj with
        | None -> input
        | Some tproj -> gemv ~trans:`T tproj input
      in
      let m = Mat.dim2 inducing in
      let res = Vec.create m in
      let tmp = { x = 0. } in
      begin match multiscales with
      | None ->
          for c = 1 to m do
            for i = 1 to d do
              let diff = projection.{i} -. inducing.{i, c} in
              tmp.x <- tmp.x +. diff *. diff
            done;
            res.{c} <- calc_res_el ~log_sf2 tmp;
          done;
      | Some multiscales ->
          for c = 1 to m do
            for i = 1 to d do
              let diff = projection.{i} -. inducing.{i, c} in
              let scale = multiscales.{i, c} in
              update_tmp_sum ~tmp ~diff ~scale
            done;
            res.{c} <- calc_res_el ~log_sf2 tmp;
          done;
      end;
      res

    let weighted_eval k input inducing ~coeffs =
      dot ~x:(eval k input inducing) coeffs

    let eval_one k _input = k.Kernel.sf2
  end

  module Inputs = struct
    type t = mat

    let create = Mat.of_col_vecs
    let get_n_points = Mat.dim2
    let choose_subset = choose_cols

    let create_default_kernel_params inputs ~n_inducing =
      let big_dim = Mat.dim1 inputs in
      let n_inputs = Mat.dim2 inputs in
      let d = min big_dim 10 in
      let tproj = Mat.create big_dim d in
      let factor = float n_inputs /. float big_dim in
      for r = 1 to big_dim do
        let sum_ref = ref 0. in
        for c = 1 to n_inputs do sum_ref := !sum_ref +. inputs.{r, c} done;
        let mean_factor = factor /. !sum_ref in
        for c = 1 to d do
          tproj.{r, c} <-  mean_factor *. (Random.float 2. -. 1.)
        done;
      done;
      {
        Params.
        d;
        log_sf2 = Random.float 2. -. 1.;
        tproj = Some tproj;
        log_hetero_skedasticity = Some (Vec.make n_inducing ~-.5.);
        log_multiscales_m05 = Some (Mat.make0 d n_inducing);
      }

    let project k inputs =
      match k.Kernel.params.Params.tproj with
      | None -> inputs
      | Some tproj -> gemm ~transa:`T tproj inputs

    let create_inducing = project
    let calc_upper k inputs = calc_upper_vanilla k (project k inputs)
    let calc_diag k inputs = Vec.make (Mat.dim2 inputs) k.Kernel.sf2

    let calc_cross_with_projections k ~projections ~inducing =
      let { Kernel.multiscales; params = { Params.d; log_sf2 } } = k in
      let m = Mat.dim2 inducing in
      let n = Mat.dim2 projections in
      let res = Mat.create n m in
      let tmp = { x = 0. } in
      begin match multiscales with
      | None ->
          for c = 1 to m do
            for r = 1 to n do
              for i = 1 to d do
                let diff = projections.{i, r} -. inducing.{i, c} in
                tmp.x <- tmp.x +. diff *. diff
              done;
              res.{r, c} <- calc_res_el ~log_sf2 tmp;
            done;
          done;
      | Some multiscales ->
          for c = 1 to m do
            for r = 1 to n do
              for i = 1 to d do
                let diff = projections.{i, r} -. inducing.{i, c} in
                let scale = multiscales.{i, c} in
                update_tmp_sum ~tmp ~diff ~scale;
              done;
              res.{r, c} <- calc_res_el ~log_sf2 tmp;
            done;
          done;
      end;
      res

    let calc_cross k ~inputs ~inducing =
      let projections = project k inputs in
      calc_cross_with_projections k ~projections ~inducing

    let weighted_eval k ~inputs ~inducing ~coeffs =
      gemv (calc_cross k ~inputs ~inducing) coeffs
  end
end

module Proj_hyper = struct type t = { big_dim : int; small_dim : int } end
module Dim_hyper = struct type t = int end
module Inducing_hyper = struct type t = { ind : int; dim : int } end

module Hyper_repr = struct
  type t = [
    | `Log_sf2
    | `Proj of Proj_hyper.t
    | `Log_hetero_skedasticity of Dim_hyper.t
    | `Log_multiscale_m05 of Inducing_hyper.t
    | `Inducing_hyper of Inducing_hyper.t
  ]
end

module Deriv = struct
  module Eval = Eval

  module Hyper = struct
    type t = Hyper_repr.t

    let get_all { Eval.Kernel.params } inducing _inputs =
      let
        { Params.d; tproj; log_hetero_skedasticity; log_multiscales_m05 } =
          params
      in
      let m = Mat.dim2 inducing in
      let n_mandatory_hypers = 1 + d * m in
      let n_hypers_ref = ref n_mandatory_hypers in
      let update_count_mat maybe_mat =
        option_iter maybe_mat ~f:(fun mat ->
          n_hypers_ref := !n_hypers_ref + Mat.dim1 mat * Mat.dim2 mat)
      in
      let update_count_vec maybe_vec =
        option_iter maybe_vec ~f:(fun vec ->
          n_hypers_ref := !n_hypers_ref + Vec.dim vec)
      in
      update_count_mat tproj;
      update_count_vec log_hetero_skedasticity;
      update_count_mat log_multiscales_m05;
      let n_hypers = !n_hypers_ref in
      let hypers = Array.create ~len:n_hypers `Log_sf2 in
      for ind = 1 to m do
        let indd = (ind - 1) * d in
        for dim = 1 to d do
          let inducing_hyper = { Inducing_hyper.ind; dim } in
          hypers.(indd + dim) <- `Inducing_hyper inducing_hyper
        done
      done;
      let pos_ref = ref n_mandatory_hypers in
      option_iter tproj ~f:(fun tproj ->
        let dim = Mat.dim1 tproj in
        for big_dim = 1 to dim do
          for small_dim = 1 to d do
            let pos = !pos_ref in
            pos_ref := pos + 1;
            hypers.(pos) <- `Proj { Proj_hyper.big_dim; small_dim };
          done;
        done);
      option_iter log_hetero_skedasticity ~f:(fun log_hetero_skedasticity ->
        let m = Vec.dim log_hetero_skedasticity in
        for i = 1 to m do
          let pos = !pos_ref in
          pos_ref := pos + 1;
          hypers.(pos) <- `Log_hetero_skedasticity i;
        done);
      option_iter log_multiscales_m05 ~f:(fun log_multiscales_m05 ->
        for ind = 1 to Mat.dim2 log_multiscales_m05 do
          for dim = 1 to d do
            let pos = !pos_ref in
            pos_ref := pos + 1;
            hypers.(pos) <- `Log_multiscale_m05 { Inducing_hyper.ind; dim };
          done;
        done);
      hypers

    let option_get_value name = function
      | None ->
          failwithf "Deriv.Hyper.option_get_value: %s not supported" name ()
      | Some v -> v

    let get_value { Eval.Kernel.params } inducing _inputs = function
      | `Log_sf2 -> params.Params.log_sf2
      | `Proj { Proj_hyper.big_dim; small_dim } ->
          (option_get_value "tproj" params.Params.tproj).{big_dim, small_dim}
      | `Log_hetero_skedasticity dim ->
          (option_get_value "log_hetero_skedasticity"
            params.Params.log_hetero_skedasticity).{dim}
      | `Log_multiscale_m05 { Inducing_hyper.ind; dim } ->
          (option_get_value
            "log_multiscales_m05" params.Params.log_multiscales_m05).{dim, ind}
      | `Inducing_hyper { Inducing_hyper.ind; dim } -> inducing.{dim, ind}

    let set_values { Eval.Kernel.params } inducing inputs hypers values =
      let log_sf2_ref = ref params.Params.log_sf2 in
      let lazy_opt name f opt_v = lazy (f (option_get_value name opt_v)) in
      let tproj_lazy = lazy_opt "tproj" lacpy params.Params.tproj in
      let log_hetero_skedasticity_lazy =
        lazy_opt "log_hetero_skedasticity"
          copy params.Params.log_hetero_skedasticity
      in
      let log_multiscales_m05_lazy =
        lazy_opt "log_multiscales_m05" lacpy params.Params.log_multiscales_m05
      in
      let inducing_lazy = lazy (lacpy inducing) in
      for i = 1 to Array.length hypers do
        match hypers.(i - 1) with
        | `Log_sf2 -> log_sf2_ref := values.{i}
        | `Proj { Proj_hyper.big_dim; small_dim } ->
            (Lazy.force tproj_lazy).{big_dim, small_dim} <- values.{i}
        | `Log_hetero_skedasticity dim ->
            (Lazy.force log_hetero_skedasticity_lazy).{dim} <- values.{i}
        | `Log_multiscale_m05 { Inducing_hyper.ind; dim } ->
            (Lazy.force log_multiscales_m05_lazy).{dim, ind} <- values.{i}
        | `Inducing_hyper { Inducing_hyper.ind; dim } ->
            (Lazy.force inducing_lazy).{dim, ind} <- values.{i}
      done;
      let lift_opt lazy_value value =
        if Lazy.is_val lazy_value then Some (Lazy.force lazy_value)
        else value
      in
      let lift lazy_value value =
        if Lazy.is_val lazy_value then Lazy.force lazy_value
        else value
      in
      let new_kernel =
        Eval.Kernel.create
          {
            Params.
            d = params.Params.d;
            log_sf2 = !log_sf2_ref;
            tproj = lift_opt tproj_lazy params.Params.tproj;
            log_hetero_skedasticity =
              lift_opt log_hetero_skedasticity_lazy
                params.Params.log_hetero_skedasticity;
            log_multiscales_m05 =
              lift_opt
                log_multiscales_m05_lazy params.Params.log_multiscales_m05;
          }
      in
      let new_inducing = lift inducing_lazy inducing in
      new_kernel, new_inducing, inputs
  end

  type deriv_common = { kernel : Eval.Kernel.t; eval_mat : mat }

  module Inducing = struct
    type upper = Eval.Inducing.t * deriv_common

    let calc_shared_upper kernel inducing =
      let eval_mat = Eval.Inducing.calc_upper kernel inducing in
      eval_mat, (inducing, { kernel; eval_mat })

    let calc_deriv_upper (inducing, { kernel; eval_mat }) hyper =
      match hyper with
      | `Log_sf2 ->
          begin
            match kernel.Eval.Kernel.hetero_skedasticity with
            | None -> `Factor 1.
            | Some hetero_skedasticity ->
                let res = lacpy eval_mat in
                for i = 1 to Mat.dim1 res do
                  res.{i, i} <- res.{i, i} -. hetero_skedasticity.{i}
                done;
                `Dense res
          end
      | `Proj _ -> `Const 0.
      | `Log_hetero_skedasticity dim ->
          begin
            match kernel.Eval.Kernel.hetero_skedasticity with
            | None ->
                failwith (
                    "Cov_se_fat.Deriv.Inducing.calc_deriv_upper: \
                    heteroskedastic modeling disabled, \
                    cannot calculate derivative")
            | Some hetero_skedasticity ->
                let deriv = Vec.make0 (Vec.dim hetero_skedasticity) in
                deriv.{dim} <- hetero_skedasticity.{dim};
                (* TODO: sparse diagonal derivatives? *)
                `Diag_vec deriv
          end
      | `Log_multiscale_m05 { Inducing_hyper.ind; dim } ->
          begin match kernel.Eval.Kernel.multiscales with
          | None ->
              failwith (
                  "Cov_se_fat.Deriv.Inducing.calc_deriv_upper: \
                  multiscale modeling disabled, cannot calculate derivative")
          | Some multiscales ->
              let m = Mat.dim2 eval_mat in
              let res = Mat.create 1 m in
              let inducing_dim = inducing.{dim, ind} in
              let multiscale = multiscales.{dim, ind} in
              let multiscale_const = multiscale -. 1. in
              let h = 0.5 in
              let multiscale_h = h -. multiscale in
              let multiscale_factor = h *. multiscale_h in
              for i = 1 to ind - 1 do
                let diff = inducing.{dim, i} -. inducing_dim in
                let iscale = 1. /. (multiscales.{dim, i} +. multiscale_const) in
                let sdiff = diff *. iscale in
                let sdiff2 = sdiff *. sdiff in
                let inner = (iscale -. sdiff2) *. multiscale_factor in
                res.{1, i} <- inner *. eval_mat.{i, ind}
              done;
              begin match kernel.Eval.Kernel.hetero_skedasticity with
              | None ->
                  res.{1, ind} <-
                    multiscale_h /. (multiscale +. multiscale_const)
                      *. eval_mat.{ind, ind};
              | Some hetero_skedasticity ->
                  res.{1, ind} <-
                    multiscale_h /. (multiscale +. multiscale_const)
                      *. (eval_mat.{ind, ind} -. hetero_skedasticity.{ind});
              end;
              for i = ind + 1 to m do
                let diff = inducing.{dim, i} -. inducing_dim in
                let iscale = 1. /. (multiscales.{dim, i} +. multiscale_const) in
                let sdiff = diff *. iscale in
                let sdiff2 = sdiff *. sdiff in
                let inner = (iscale -. sdiff2) *. multiscale_factor in
                res.{1, i} <- inner *. eval_mat.{ind, i}
              done;
              let rows = Sparse_indices.create 1 in
              rows.{1} <- ind;
              `Sparse_rows (res, rows)
          end
      | `Inducing_hyper { Inducing_hyper.ind; dim } ->
          let m = Mat.dim2 eval_mat in
          let res = Mat.create 1 m in
          let inducing_dim = inducing.{dim, ind} in
          begin match kernel.Eval.Kernel.multiscales with
          | None ->
              for i = 1 to ind - 1 do
                let diff = inducing.{dim, i} -. inducing_dim in
                res.{1, i} <- diff *. eval_mat.{i, ind}
              done;
              res.{1, ind} <- 0.;
              for i = ind + 1 to m do
                let diff = inducing.{dim, i} -. inducing_dim in
                res.{1, i} <- diff *. eval_mat.{ind, i}
              done
          | Some multiscales ->
              let multiscale_const = multiscales.{dim, ind} -. 1. in
              for i = 1 to ind - 1 do
                let diff = inducing.{dim, i} -. inducing_dim in
                let scale = multiscales.{dim, i} +. multiscale_const in
                res.{1, i} <- diff /. scale *. eval_mat.{i, ind}
              done;
              res.{1, ind} <- 0.;
              for i = ind + 1 to m do
                let diff = inducing.{dim, i} -. inducing_dim in
                let scale = multiscales.{dim, i} +. multiscale_const in
                res.{1, i} <- diff /. scale *. eval_mat.{ind, i}
              done;
          end;
          let rows = Sparse_indices.create 1 in
          rows.{1} <- ind;
          `Sparse_rows (res, rows)
  end

  module Inputs = struct
    (* Diag *)

    type diag = Eval.Kernel.t

    let calc_shared_diag k diag_eval_inputs =
      Eval.Inputs.calc_diag k diag_eval_inputs, k

    let calc_deriv_diag _diag = function
      | `Log_sf2 -> `Factor 1.
      | `Proj _ | `Log_hetero_skedasticity _ | `Log_multiscale_m05 _
      | `Inducing_hyper _ -> `Const 0.

    (* Cross *)

    module Cross = struct
      type t = {
        common : deriv_common;
        inputs : Eval.Inputs.t;
        inducing : Eval.Inducing.t;
        projections : Eval.Inducing.t;
      }
    end

    type cross = Cross.t

    let calc_shared_cross kernel ~inputs ~inducing =
      let projections = Eval.Inputs.project kernel inputs in
      let eval_mat =
        Eval.Inputs.calc_cross_with_projections kernel ~projections ~inducing
      in
      let shared =
        { Cross.common = { kernel; eval_mat }; inputs; inducing; projections }
      in
      eval_mat, shared

    let check_tproj_available = function
      | None ->
          failwith
            "Cov_se_fat.Deriv.Inputs.calc_deriv_cross: \
            tproj disabled, cannot calculate derivative"
      | Some _ -> ()

    let calc_deriv_cross cross hyper =
      let
        { Cross.common = { kernel; eval_mat }; inputs; inducing; projections } =
          cross
      in
      match hyper with
      | `Log_sf2 -> `Factor 1.
      | `Proj { Proj_hyper.big_dim; small_dim } ->
          check_tproj_available kernel.Eval.Kernel.params.Params.tproj;
          let m = Mat.dim2 inducing in
          let n = Mat.dim2 inputs in
          let res = Mat.create n m in
          begin match kernel.Eval.Kernel.multiscales with
          | None ->
              for c = 1 to m do
                let ind_el = inducing.{small_dim, c} in
                for r = 1 to n do
                  let alpha = inputs.{big_dim, r} in
                  let proj = projections.{small_dim, r} in
                  res.{r, c} <- alpha *. (ind_el -. proj) *. eval_mat.{r, c}
                done
              done;
          | Some multiscales ->
              for c = 1 to m do
                let ind_el = inducing.{small_dim, c} in
                let multiscale = multiscales.{small_dim, c} in
                for r = 1 to n do
                  let alpha = inputs.{big_dim, r} in
                  let proj = projections.{small_dim, r} in
                  res.{r, c} <-
                    alpha *. ((ind_el -. proj) /. multiscale) *. eval_mat.{r, c}
                done
              done;
          end;
          `Dense res
      | `Log_hetero_skedasticity _ -> `Const 0.
      | `Log_multiscale_m05 { Inducing_hyper.ind; dim } ->
          begin match kernel.Eval.Kernel.multiscales with
          | None ->
              failwith (
                  "Cov_se_fat.Deriv.Inputs.calc_deriv_cross: \
                  multiscale modeling disabled, cannot calculate derivative")
          | Some multiscales ->
            let n = Mat.dim1 eval_mat in
            let res = Mat.create n 1 in
            let inducing_dim = inducing.{dim, ind} in
            let multiscale = multiscales.{dim, ind} in
            let h = 0.5 in
            let multiscale_h = h -. multiscale in
            let multiscale_factor = h *. multiscale_h in
            for r = 1 to n do
              let diff = projections.{dim, r} -. inducing_dim in
              let iscale = 1. /. multiscales.{dim, ind} in
              let sdiff = diff *. iscale in
              let sdiff2 = sdiff *. sdiff in
              let inner = (iscale -. sdiff2) *. multiscale_factor in
              res.{r, 1} <- inner *. eval_mat.{r, ind}
            done;
            let cols = Sparse_indices.create 1 in
            cols.{1} <- ind;
            `Sparse_cols (res, cols)
          end
      | `Inducing_hyper { Inducing_hyper.ind; dim } ->
          let n = Mat.dim1 eval_mat in
          let res = Mat.create n 1 in
          let inducing_dim = inducing.{dim, ind} in
          begin match kernel.Eval.Kernel.multiscales with
          | None ->
              for r = 1 to n do
                let diff = projections.{dim, r} -. inducing_dim in
                res.{r, 1} <- diff *. eval_mat.{r, ind}
              done;
          | Some multiscales ->
              let multiscale_factor = 1. /. multiscales.{dim, ind} in
              for r = 1 to n do
                let diff = projections.{dim, r} -. inducing_dim in
                res.{r, 1} <- multiscale_factor *. diff *. eval_mat.{r, ind}
              done;
          end;
          let cols = Sparse_indices.create 1 in
          cols.{1} <- ind;
          `Sparse_cols (res, cols)
  end
end
