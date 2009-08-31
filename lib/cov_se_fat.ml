open Printf
open Bigarray
open Lacaml.Impl.D

open Interfaces

let option_map ~f = function None -> None | Some v -> Some (f v)
let option_iter ~f = function None -> () | Some v -> f v

module Params = struct
  type params = {
    d : int;
    log_sf2 : float;
    tproj : mat option;
    log_hetero_skedasticity : vec option;
    log_multiscales : mat option;
  }

  type t = params

  let create params =
    let check v_dim name v =
      let n = v_dim v in
      if n <> params.d then
        failwith (
          sprintf
            "Cov_se_fat.Params.create: %s projection (%d) disagrees \
            with target dimension d (%d)" name n params.d)
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
        option_map params.Params.log_multiscales ~f:(Mat.map exp)
      in
      {
        params = params;
        sf2 = exp params.Params.log_sf2;
        hetero_skedasticity = hetero_skedasticity;
        multiscales = multiscales;
      }

    let get_params k = k.params
    let get_d k = k.params.Params.d
  end

  open Kernel

  let calc_upper_vanilla k inputs =
    let m = Mat.dim2 inputs in
    let res = Mat.create m m in
    let { Kernel.params = { Params.d = d; log_sf2 = log_sf2 } } = k in
    let res_ref = ref 0. in
    for c = 2 to m do
      for r = 1 to c - 1 do
        for i = 1 to d do
          let diff = inputs.{i, r} -. inputs.{i, c} in
          res_ref := !res_ref +. diff *. diff
        done;
        res.{r, c} <- exp (log_sf2 -. 0.5 *. !res_ref);
        res_ref := 0.;
      done;
    done;
    res

  module Inducing = struct
    type t = mat

    let get_n_points = Mat.dim2

    let calc_upper k inducing =
      let m = Mat.dim2 inducing in
      let
        {
          Kernel.
          sf2 = sf2;
          hetero_skedasticity = hetero_skedasticity;
          multiscales = multiscales;
          params = { Params.d = d; log_sf2 = log_sf2 }
        } = k
      in
      let res_ref = ref 0. in
      let res =
        match multiscales with
        | None -> calc_upper_vanilla k inducing
        | Some multiscales ->
            (* TODO: save squared differences for later *)
            let res = Mat.create m m in
            for c = 2 to m do
              for r = 1 to c - 1 do
                for i = 1 to d do
                  let diff = inducing.{i, r} -. inducing.{i, c} in
                  let scale = multiscales.{i, r} +. multiscales.{i, c} -. 1. in
                  res_ref := !res_ref +. diff *. (diff /. scale)
                done;
                res.{r, c} <- exp (log_sf2 -. 0.5 *. !res_ref);
                res_ref := 0.;
              done;
            done;
            res
      in
      begin match hetero_skedasticity with
      | None -> for i = 1 to m do res.{i, i} <- sf2 done
      | Some hetero_skedasticity ->
          for i = 1 to m do res.{i, i} <- sf2 +. hetero_skedasticity.{i} done
      end;
      res
  end

  module Input = struct
    type t = vec

    let eval k inducing input =
      let
        {
          Kernel.
          multiscales = multiscales;
          params = { Params.d = d; log_sf2 = log_sf2; tproj = tproj }
        } = k
      in
      let projection =
        match tproj with
        | None -> input
        | Some tproj -> gemv ~trans:`T tproj input
      in
      let m = Mat.dim2 inducing in
      let res = Vec.create m in
      let res_ref = ref 0. in
      begin match multiscales with
      | None ->
          for c = 1 to m do
            for i = 1 to d do
              let diff = projection.{i} -. inducing.{i, c} in
              res_ref := !res_ref +. diff *. diff
            done;
            res.{c} <- exp (log_sf2 -. 0.5 *. !res_ref);
            res_ref := 0.;
          done;
      | Some multiscales ->
          for c = 1 to m do
            for i = 1 to d do
              let diff = projection.{i} -. inducing.{i, c} in
              let scale = multiscales.{i, c} in
              res_ref := !res_ref +. diff *. (diff /. scale)
            done;
            res.{c} <- exp (log_sf2 -. 0.5 *. !res_ref);
            res_ref := 0.;
          done;
      end;
      res

    let weighted_eval k inducing ~coeffs input =
      dot ~x:(eval k inducing input) coeffs

    let eval_one k _input = k.Kernel.sf2
  end

  module Inputs = struct
    type t = mat

    let get_n_points = Mat.dim2

    let choose_subset = Utils.choose_cols

    let create_default_kernel_params ~n_inducing inputs =
      let big_dim = Mat.dim1 inputs in
      let small_dim = min big_dim 10 in
      let _n = Mat.dim2 inputs in
      {
        Params.
        d = small_dim;
        log_sf2 = Random.float 1.;
        tproj = None;
(*         tproj = Some (Mat.make big_dim small_dim (1. /. float n)); *)
(*         log_hetero_skedasticity = None; *)
        log_hetero_skedasticity = Some (Vec.random n_inducing);
(*         log_multiscales = None; *)
        log_multiscales = Some (Mat.random small_dim n_inducing);
      }

    let project { Kernel.params = { Params.tproj = tproj } } inputs =
      match tproj with
      | None -> inputs
      | Some tproj -> gemm ~transa:`T tproj inputs

    let create_inducing = project

    let calc_upper k inputs =
      let res = calc_upper_vanilla k (project k inputs) in
      let sf2 = k.Kernel.sf2 in
      for i = 1 to Mat.dim2 res do res.{i, i} <- sf2 done;
      res

    let calc_diag k inputs = Vec.make (Mat.dim2 inputs) k.Kernel.sf2

    let calc_cross_with_projections k ~inducing ~projections =
      let
        {
          Kernel.
          multiscales = multiscales;
          params = { Params.d = d; log_sf2 = log_sf2 }
        } = k
      in
      let n = Mat.dim2 projections in
      let m = Mat.dim2 inducing in
      let res = Mat.create m n in
      let res_ref = ref 0. in
      begin match multiscales with
      | None ->
          for c = 1 to n do
            for r = 1 to m do
              for i = 1 to d do
                let diff = projections.{i, c} -. inducing.{i, r} in
                res_ref := !res_ref +. diff *. diff
              done;
              res.{r, c} <- exp (log_sf2 -. 0.5 *. !res_ref);
              res_ref := 0.;
            done;
          done;
      | Some multiscales ->
          for c = 1 to n do
            for r = 1 to m do
              for i = 1 to d do
                let diff = projections.{i, c} -. inducing.{i, r} in
                let scale = multiscales.{i, r} in
                res_ref := !res_ref +. diff *. (diff /. scale)
              done;
              res.{r, c} <- exp (log_sf2 -. 0.5 *. !res_ref);
              res_ref := 0.;
            done;
          done;
      end;
      res

    let calc_cross k inducing inputs =
      let projections = project k inputs in
      calc_cross_with_projections k ~inducing ~projections

    let weighted_eval k inducing ~coeffs inputs =
      gemv (calc_cross k inducing inputs) coeffs
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
    | `Log_multiscales of Inducing_hyper.t
    | `Inducing_hyper of Inducing_hyper.t
  ]
end

module Deriv = struct
  module Eval = Eval

  module Hyper = struct
    type t = Hyper_repr.t

    let get_all { Eval.Kernel.params = params } inducing =
      let
        {
          Params.
          d = d;
          tproj = tproj;
          log_hetero_skedasticity = log_hetero_skedasticity;
          log_multiscales = log_multiscales;
        } = params
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
      update_count_mat log_multiscales;
      let n_hypers = !n_hypers_ref in
      let hypers = Array.create n_hypers `Log_sf2 in
      for ind = 1 to m do
        let indd = (ind - 1) * d in
        for dim = 1 to d do
          let inducing_hyper = { Inducing_hyper.ind = ind; dim = dim } in
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
            hypers.(pos) <-
              `Proj { Proj_hyper.big_dim = big_dim; small_dim = small_dim };
          done;
        done);
      option_iter log_hetero_skedasticity ~f:(fun log_hetero_skedasticity ->
        let m = Vec.dim log_hetero_skedasticity in
        for i = 1 to m do
          let pos = !pos_ref in
          pos_ref := pos + 1;
          hypers.(pos) <- `Log_hetero_skedasticity i;
        done);
      option_iter log_multiscales ~f:(fun log_multiscales ->
        for ind = 1 to Mat.dim2 log_multiscales do
          for dim = 1 to d do
            let pos = !pos_ref in
            pos_ref := pos + 1;
            hypers.(pos) <-
              `Log_multiscales { Inducing_hyper.ind = ind; dim = dim };
          done;
        done);
      hypers

    let option_get_value name = function
      | None ->
          failwith (
            sprintf "Deriv.Hyper.option_get_value: %s not supported" name)
      | Some v -> v

    let get_value { Eval.Kernel.params = params } inducing = function
      | `Log_sf2 -> params.Params.log_sf2
      | `Proj { Proj_hyper.big_dim = big_dim; small_dim = small_dim } ->
          (option_get_value "tproj" params.Params.tproj).{big_dim, small_dim}
      | `Log_hetero_skedasticity dim ->
          (option_get_value "log_hetero_skedasticity"
            params.Params.log_hetero_skedasticity).{dim}
      | `Log_multiscales { Inducing_hyper.ind = ind; dim = dim } ->
          (option_get_value
            "log_multiscales" params.Params.log_multiscales).{dim, ind}
      | `Inducing_hyper { Inducing_hyper.ind = ind; dim = dim } ->
          inducing.{dim, ind}

    let set_values kernel inducing hypers values =
      let { Eval.Kernel.params = params } = kernel in
      let log_sf2_ref = ref params.Params.log_sf2 in
      let lazy_opt name f opt_v = lazy (f (option_get_value name opt_v)) in
      let tproj_lazy = lazy_opt "tproj" lacpy params.Params.tproj in
      let log_hetero_skedasticity_lazy =
        lazy_opt "log_hetero_skedasticity"
          copy params.Params.log_hetero_skedasticity
      in
      let log_multiscales_lazy =
        lazy_opt "log_multiscales" lacpy params.Params.log_multiscales
      in
      let inducing_lazy = lazy (lacpy inducing) in
      for i = 1 to Array.length hypers do
        match hypers.(i - 1) with
        | `Log_sf2 -> log_sf2_ref := values.{i}
        | `Proj { Proj_hyper.big_dim = big_dim; small_dim = small_dim } ->
            (Lazy.force tproj_lazy).{big_dim, small_dim} <- values.{i}
        | `Log_hetero_skedasticity dim ->
            (Lazy.force log_hetero_skedasticity_lazy).{dim} <- values.{i}
        | `Log_multiscales { Inducing_hyper.ind = ind; dim = dim } ->
            (Lazy.force log_multiscales_lazy).{dim, ind} <- values.{i}
        | `Inducing_hyper { Inducing_hyper.ind = ind; dim = dim } ->
            (Lazy.force inducing_lazy).{dim, ind} <- values.{i}
      done;
      let lift_opt lazy_value value =
        if Lazy.lazy_is_val lazy_value then Some (Lazy.force lazy_value)
        else value
      in
      let lift lazy_value value =
        if Lazy.lazy_is_val lazy_value then Lazy.force lazy_value
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
            log_multiscales =
              lift_opt log_multiscales_lazy params.Params.log_multiscales;
          }
      in
      let new_inducing = lift inducing_lazy inducing in
      new_kernel, new_inducing
  end

  type deriv_common = { kernel : Eval.Kernel.t; eval_mat : mat }

  module Inducing = struct
    type upper = Eval.Inducing.t * deriv_common

    let calc_shared_upper kernel inducing =
      let upper = Eval.Inducing.calc_upper kernel inducing in
      let shared = inducing, { kernel = kernel; eval_mat = upper } in
      upper, shared

    let calc_deriv_upper (inducing, common) = function
      | `Log_sf2 ->
          begin
            match common.kernel.Eval.Kernel.hetero_skedasticity with
            | None -> `Factor 1.
            | Some hetero_skedasticity ->
                let res = lacpy common.eval_mat in
                for i = 1 to Mat.dim1 res do
                  res.{i, i} <- res.{i, i} -. hetero_skedasticity.{i}
                done;
                `Dense res
          end
      | `Proj _ -> `Const 0.
      | `Log_hetero_skedasticity dim ->
          begin
            match common.kernel.Eval.Kernel.hetero_skedasticity with
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
      | `Log_multiscales { Inducing_hyper.ind = ind; dim = dim } ->
          begin match common.kernel.Eval.Kernel.multiscales with
          | None ->
              failwith (
                  "Cov_se_fat.Deriv.Inducing.calc_deriv_upper: \
                  multiscale modeling disabled, cannot calculate derivative")
          | Some multiscales ->
              let eval_mat = common.eval_mat in
              let m = Mat.dim2 eval_mat in
              let res = Mat.create 1 m in
              let inducing_dim = inducing.{dim, ind} in
              let multiscale = multiscales.{dim, ind} in
              let multiscale_2 = 0.5 *. multiscale in
              let multiscale_const = multiscale -. 1. in
              for i = 1 to ind - 1 do
                let diff = inducing.{dim, i} -. inducing_dim in
                let scale = multiscales.{dim, i} +. multiscale_const in
                let scaled_diff = diff /. scale in
                let scaled_diff_multi = scaled_diff *. multiscale_2 in
                let scaled_diff_eval = scaled_diff *. eval_mat.{i, ind} in
                res.{1, i} <- scaled_diff_multi *. scaled_diff_eval
              done;
              res.{1, ind} <- 0.;
              for i = ind + 1 to m do
                let diff = inducing.{dim, i} -. inducing_dim in
                let scale = multiscales.{dim, i} +. multiscale_const in
                let scaled_diff = diff /. scale in
                let scaled_diff_multi = scaled_diff *. multiscale_2 in
                let scaled_diff_eval = scaled_diff *. eval_mat.{ind, i} in
                res.{1, i} <- scaled_diff_multi *. scaled_diff_eval
              done;
              let rows = Sparse_indices.create 1 in
              rows.{1} <- ind;
              `Sparse_rows (res, rows)
          end
      | `Inducing_hyper { Inducing_hyper.ind = ind; dim = dim } ->
          let eval_mat = common.eval_mat in
          let m = Mat.dim2 eval_mat in
          let res = Mat.create 1 m in
          let inducing_dim = inducing.{dim, ind} in
          begin match common.kernel.Eval.Kernel.multiscales with
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
      | `Proj _ | `Log_hetero_skedasticity _ | `Log_multiscales _
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

    let calc_shared_cross k inducing inputs =
      let projections = Eval.Inputs.project k inputs in
      let eval_mat =
        Eval.Inputs.calc_cross_with_projections k ~inducing ~projections
      in
      let shared =
        {
          Cross.
          common = { kernel = k; eval_mat = eval_mat };
          inducing = inducing;
          inputs = inputs;
          projections = projections;
        }
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
        {
          Cross.
          common = { kernel = kernel; eval_mat = eval_mat };
          inducing = inducing;
          inputs = inputs;
          projections = projections;
        } = cross
      in
      match hyper with
      | `Log_sf2 -> `Factor 1.
      | `Proj { Proj_hyper.big_dim = big_dim; small_dim = small_dim } ->
          check_tproj_available kernel.Eval.Kernel.params.Params.tproj;
          let m = Mat.dim2 inducing in
          let n = Mat.dim2 inputs in
          let res = Mat.create m n in
          for c = 1 to n do
            let proj = projections.{small_dim, big_dim} in
            for r = 1 to m do
              let alpha = inputs.{big_dim, c} in
              let ind_el = inducing.{small_dim, r} in
              res.{r, c} <- (alpha *. (ind_el -. proj)) *. eval_mat.{r, c}
            done
          done;
          `Dense res
      | `Log_hetero_skedasticity _ -> `Const 0.
      | `Log_multiscales { Inducing_hyper.ind = ind; dim = dim } ->
          begin match kernel.Eval.Kernel.multiscales with
          | None ->
              failwith (
                  "Cov_se_fat.Deriv.Inputs.calc_deriv_cross: \
                  multiscale modeling disabled, cannot calculate derivative")
          | Some multiscales ->
              let n = Mat.dim2 eval_mat in
              let res = Mat.create 1 n in
              let inducing_dim = inducing.{dim, ind} in
              let multiscale_factor = 0.5 /. multiscales.{dim, ind} in
              for c = 1 to n do
                let diff = projections.{dim, c} -. inducing_dim in
                let diff_multi = diff *. multiscale_factor in
                let diff_eval = diff *. eval_mat.{ind, c} in
                res.{1, c} <- diff_multi *. diff_eval
              done;
              let rows = Sparse_indices.create 1 in
              rows.{1} <- ind;
              `Sparse_rows (res, rows)
          end
      | `Inducing_hyper { Inducing_hyper.ind = ind; dim = dim } ->
          let n = Mat.dim2 eval_mat in
          let res = Mat.create 1 n in
          let inducing_dim = inducing.{dim, ind} in
          begin match kernel.Eval.Kernel.multiscales with
          | None ->
              for c = 1 to n do
                let diff = projections.{dim, c} -. inducing_dim in
                res.{1, c} <- diff *. eval_mat.{ind, c}
              done;
          | Some multiscales ->
              let multiscale_factor = 1. /. multiscales.{dim, ind} in
              for c = 1 to n do
                let diff = projections.{dim, c} -. inducing_dim in
                res.{1, c} <- multiscale_factor *. diff *. eval_mat.{ind, c}
              done;
          end;
          let rows = Sparse_indices.create 1 in
          rows.{1} <- ind;
          `Sparse_rows (res, rows)
  end
end
