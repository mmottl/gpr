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
    type t = mat

    module Prepared = struct
      type upper = {
        upper : mat;
        inducing : t;
      }

      let calc_km points =
        let m = Mat.dim2 points in
        (* TODO: make upper triangle only *)
        syrk ~trans:`T points ~beta:1. ~c:(Mat.make m m 1.)

      let calc_upper points = { upper = calc_km points; inducing = points }
    end

    let calc_upper_mat k upper =
      (* TODO: copy and scale upper triangle only *)
      let res = Mat.copy upper in
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
      let res = Vec.create n in
      let const = k.Kernel.const in
      for i = 1 to n do
        (* TODO: optimize sqr_nrm2 and col *)
        res.{i} <- const *. (Vec.sqr_nrm2 (Mat.col inputs i) +. 1.);
      done;
      res

    let calc_cross k cross =
      let res = Mat.copy cross in
      Mat.scal k.Kernel.const res;
      res

    let weighted_eval k ~coeffs cross =
      if Vec.dim coeffs <> Mat.dim1 cross then
        failwith
          "Gpr.Cov_lin_one.Eval.Inputs.weighted_eval: dim(coeffs) <> m";
      gemv ~alpha:k.Kernel.const ~trans:`T cross coeffs
  end
end

module Hyper = struct type t = [ `Log_theta ] end

let calc_deriv_mat mat =
  (* TODO: copy and scale upper only *)
  (* TODO: even better: introduce passing through matrix and factor *)
  let res = Mat.copy mat in
  Mat.scal (-2.) res;
  `Dense res

module Inducing = struct
  module Prepared = struct
    type upper = Eval.Inducing.Prepared.upper

    let calc_upper upper = upper
  end

  type shared = Eval.Inducing.t

  let calc_shared_upper k prepared_upper =
    let upper = Eval.Inducing.calc_upper k prepared_upper in
    upper, upper

  let calc_deriv_upper upper `Log_theta = calc_deriv_mat upper
end

module Inputs = struct
  module Prepared = struct
    type cross = Eval.Inputs.Prepared.cross

    let calc_cross _upper cross = cross
  end

  type diag = vec
  type cross = Eval.Inputs.t

  let calc_shared_diag k diag_eval_inputs =
    let diag = Eval.Inputs.calc_diag k diag_eval_inputs in
    diag, diag

  let calc_shared_cross k cross_eval_inputs =
    let cross = Eval.Inputs.calc_cross k cross_eval_inputs in
    cross, cross

  let calc_deriv_diag diag `Log_theta =
    let res = copy diag in
    scal (-2.) res;
    `Vec res

  let calc_deriv_cross cross `Log_theta = calc_deriv_mat cross
end
