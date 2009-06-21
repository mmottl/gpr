open Printf
open Lacaml.Impl.D
open Lacaml.Io

open Utils

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
    type t = mat

    module Prepared = struct
      type upper = { upper : mat; inducing : t }

      let calc_km points =
        let m = Mat.dim2 points in
        (* TODO: make upper triangle only *)
        syrk ~trans:`T points ~beta:1. ~c:(Mat.make m m 1.)

      let calc_upper points = { upper = calc_km points; inducing = points }
    end

    let calc_upper_mat k upper =
      let res = lacpy ~uplo:`U upper in
      (* TODO: scale upper triangle only *)
      Mat.scal k.Kernel.const res;
      res

    let calc_upper k prepared_upper =
      calc_upper_mat k prepared_upper.Prepared.upper
  end

  module Input = struct
    type t = vec

    module Prepared = struct
      type cross = t

      let calc_cross { Inducing.Prepared.inducing = inducing } input =
        gemv ~trans:`T inducing
          input ~beta:1. ~y:(Vec.make (Mat.dim2 inducing) 1.)
    end

    let eval k cross =
      let res = copy cross in
      scal k.Kernel.const res;
      res

    let weighted_eval k ~coeffs cross =
      if Vec.dim coeffs <> Vec.dim cross then
        failwith
          "Gpr.Cov_lin_one.Eval.Input.weighted_eval: dim(coeffs) <> m";
      k.Kernel.const *. dot ~x:coeffs cross

    let eval_one k input = k.Kernel.const *. (Vec.sqr_nrm2 input +. 1.)
  end

  module Inputs = struct
    type t = mat

    let get_n_inputs = Mat.dim2
    let choose_subset inputs indexes = choose_cols inputs indexes
    let create_inducing _kernel inputs = inputs
    let create_default_kernel_params _inputs = { Params.log_theta = 0. }

    module Prepared = struct
      type cross = t

      let calc_cross { Inducing.Prepared.inducing = inducing } inputs =
        let m = Mat.dim2 inducing in
        let n = Mat.dim2 inputs in
        gemm ~transa:`T inducing inputs ~beta:1. ~c:(Mat.make m n 1.)
    end

    let calc_upper k inputs =
      Inducing.calc_upper_mat k (Inducing.Prepared.calc_km inputs)

    let calc_diag k inputs =
      let n = Mat.dim2 inputs in
      let res = Mat.syrk_diag ~trans:`T inputs ~beta:1. ~y:(Vec.make n 1.) in
      scal k.Kernel.const res;
      res

    let calc_cross k cross =
      let res = lacpy cross in
      Mat.scal k.Kernel.const res;
      res

    let weighted_eval k ~coeffs cross =
      if Vec.dim coeffs <> Mat.dim1 cross then
        failwith
          "Gpr.Cov_lin_one.Eval.Inputs.weighted_eval: dim(coeffs) <> m";
      gemv ~alpha:k.Kernel.const ~trans:`T cross coeffs
  end
end

module Deriv = struct
  module Eval = Eval

  module Hyper = struct
    type t = [ `Log_theta ]

    let n_hypers = 1
    let get_n_hypers _kernel = n_hypers

    let of_index _kernel ~index =
      match index with
      | 1 -> `Log_theta
      | _ ->
          failwith (
            sprintf
              "Gpr.Cov_lin_one.Deriv.Hyper.of_index: index (%d) > n_hypers (%d)"
              index n_hypers)
  end

  let calc_deriv_common () `Log_theta = `Factor (-2.)

  module Inducing = struct
    module Prepared = struct
      type upper = Eval.Inducing.Prepared.upper

      let calc_upper upper = upper
    end

    type upper = unit

    let calc_shared_upper k prepared_upper =
      Eval.Inducing.calc_upper k prepared_upper, ()

    let calc_deriv_upper = calc_deriv_common
  end

  module Inputs = struct
    module Prepared = struct
      type cross = Eval.Inputs.Prepared.cross

      let calc_cross _upper cross = cross
    end

    type diag = unit
    type cross = unit

    let calc_shared_diag k diag_eval_inputs =
      Eval.Inputs.calc_diag k diag_eval_inputs, ()

    let calc_shared_cross k cross_eval_inputs =
      Eval.Inputs.calc_cross k cross_eval_inputs, ()

    let calc_deriv_diag = calc_deriv_common
    let calc_deriv_cross = calc_deriv_common
  end
end
