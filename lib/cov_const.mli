module Params : sig type t = { log_theta : float } end

include Interfaces.Inducing_input_gpr.Specs.Deriv
  with type Eval.Kernel.params = Params.t
  with type Eval.Inducing.t = int
  with type Eval.Input.t = unit
  with type Eval.Inputs.t = int
  with type Hyper.t = [ `Log_theta ]
