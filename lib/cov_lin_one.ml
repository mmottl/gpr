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
      gemm ~alpha ~transa:`T inducing inputs ~beta:1. ~c:(Mat.make m n alpha)

    let weighted_eval k inducing ~coeffs inputs =
      gemv ~trans:`T (calc_cross k inducing inputs) coeffs
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

  let calc_deriv_common () `Log_theta = `Factor (-2.)

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
