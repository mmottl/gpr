open Lacaml.Impl.D

include Interfaces.Inducing_input_gpr.Specs.Deriv
  with type Eval.Kernel.params = < log_theta : float >
  with type Eval.Inducing.t = mat
  with type Eval.Input.t = vec
  with type Eval.Inputs.t = mat
  with type Hyper.t = [ `Log_theta ]
