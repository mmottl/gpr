include Interfaces.Inducing_input_gpr.Specs.Deriv
  with type Eval.Kernel.params = < log_theta : float >
  with type Eval.Inducing.t = < m : int >
  with type Eval.Input.t = unit
  with type Eval.Inputs.t = < n : int >
  with type Hyper.t = [ `Log_theta ]
