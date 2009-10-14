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

    let get_n_points = Mat.dim2

    let calc_upper { Kernel.const = alpha } inducing =
      let m = Mat.dim2 inducing in
      syrk ~alpha ~trans:`T inducing ~beta:1. ~c:(Mat.make m m alpha)
  end

  module Input = struct
    type t = vec

    let eval { Kernel.const = alpha } inducing input =
      gemv ~alpha ~trans:`T inducing input
        ~beta:1. ~y:(Vec.make (Mat.dim2 inducing) alpha)

    let weighted_eval k inducing ~coeffs input =
      dot ~x:coeffs (eval k inducing input)

    let eval_one k input = k.Kernel.const *. (Vec.sqr_nrm2 input +. 1.)
  end

  module Inputs = struct
    type t = mat

    let get_n_points = Mat.dim2
    let choose_subset inputs indexes = choose_cols inputs indexes
    let create_inducing _kernel inputs = inputs

    let create_default_kernel_params ~n_inducing:_ _inputs =
      { Params.log_theta = 0. }

    let calc_upper = Inducing.calc_upper

    let calc_diag { Kernel.const = alpha } inputs =
      Mat.syrk_diag ~alpha ~trans:`T inputs ~beta:1.
        ~y:(Vec.make (Mat.dim2 inputs) alpha)

    let calc_cross { Kernel.const = alpha } inducing inputs =
      let m = Mat.dim2 inducing in
      let n = Mat.dim2 inputs in
      gemm ~alpha ~transa:`T inputs inducing ~beta:1. ~c:(Mat.make n m alpha)

    let weighted_eval k inducing ~coeffs inputs =
      gemv (calc_cross k inducing inputs) coeffs
  end
end

module Deriv = struct
  module Eval = Eval

  module Hyper = struct
    type t = [ `Log_theta ]

    let get_all _kernel _inducing = [| `Log_theta |]

    let get_value { Eval.Kernel.params = params } _inducing = function
      | `Log_theta -> params.Params.log_theta

    let set_values kernel inducing hypers values =
      let { Eval.Kernel.params = params } = kernel in
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
      new_kernel, inducing
  end

  let calc_deriv_common () `Log_theta = `Factor ~-.2.

  module Inducing = struct
    type upper = unit

    let calc_shared_upper k inducing = Eval.Inducing.calc_upper k inducing, ()
    let calc_deriv_upper = calc_deriv_common
  end

  module Inputs = struct
    type diag = unit
    type cross = unit

    let calc_shared_diag k diag_eval_inputs =
      Eval.Inputs.calc_diag k diag_eval_inputs, ()

    let calc_shared_cross k inducing cross_eval_inputs =
      Eval.Inputs.calc_cross k inducing cross_eval_inputs, ()

    let calc_deriv_diag = calc_deriv_common
    let calc_deriv_cross = calc_deriv_common
  end
end
