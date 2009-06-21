open Lacaml.Impl.D

open Interfaces.Specs

module Params : sig type t = { log_theta : float } end

module Eval :
  Eval
    with type Kernel.params = Params.t
    with type Inducing.t = mat
    with type Input.t = vec
    with type Inputs.t = mat

module Deriv :
  Deriv
    with module Eval = Eval
    with type Hyper.t = [ `Log_theta ]
