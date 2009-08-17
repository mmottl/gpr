open Lacaml.Impl.D

open Interfaces.Specs

module Params : sig type t = { log_ell : float; log_sf2 : float } end

type inducing_hyper = { ind : int; dim : int }

module Eval :
  Eval
    with type Kernel.params = Params.t
    with type Inducing.t = mat
    with type Input.t = vec
    with type Inputs.t = mat

module Deriv :
  Deriv
    with module Eval = Eval
    with type Hyper.t =
      [ `Log_ell | `Log_sf2 | `Inducing_hyper of inducing_hyper ]
