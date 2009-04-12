open Lacaml.Impl.D

module Params : sig type t = { log_ell : float; log_sf2 : float } end

type inducing_hyper = { ind : int; dim : int }

include Interfaces.Inducing_input_gpr.Specs.Deriv
  with type Eval.Kernel.params = Params.t
  with type Eval.Inducing.t = mat
  with type Eval.Input.t = vec
  with type Eval.Inputs.t = mat
  with type Hyper.t = [ `Log_ell | `Log_sf2 | `Inducing_hyper of inducing_hyper ]
