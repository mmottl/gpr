open Interfaces.Specs

module Params : sig type t = { log_theta : float } end

module Eval :
  Eval
    with type Kernel.params = Params.t
    with type Inducing.t = int
    with type Input.t = unit
    with type Inputs.t = int

module Deriv :
  Deriv
    with module Eval = Eval
    with type Hyper.t = [ `Log_theta ]
