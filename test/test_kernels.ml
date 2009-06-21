open Lacaml.Impl.D

open Gpr

module Const = struct
  include Fitc.Make_deriv (Cov_const.Deriv)

  let params = { Cov_const.Params.log_theta = log 1. }
  let kernel = Cov_const.Eval.Kernel.create params
end

module Lin_one = struct
  include Fitc.Make_deriv (Cov_lin_one.Deriv)

  let params = { Cov_lin_one.Params.log_theta = log 1. }
  let kernel = Cov_lin_one.Eval.Kernel.create params
end

module Lin_ard = struct
  include Fitc.Make_deriv (Cov_lin_ard.Deriv)

  (* Any other value with ARD is useless for graphical demonstration *)
  let params = { Cov_lin_ard.Params.log_ells = Vec.make 1 0. }
  let kernel = Cov_lin_ard.Eval.Kernel.create params
end

module SE_iso = struct
  include Fitc.Make_deriv (Cov_se_iso.Deriv)

  let params = { Cov_se_iso.Params.log_ell = log 1.; log_sf2 = log 1. }
  let kernel = Cov_se_iso.Eval.Kernel.create params
end
