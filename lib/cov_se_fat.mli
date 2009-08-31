open Lacaml.Impl.D

open Interfaces.Specs

module Params : sig
  type params = {
    d : int;
    log_sf2 : float;
    tproj : mat option;
    log_hetero_skedasticity : vec option;
    log_multiscales : mat option;
  }

  type t = private params

  val create : params -> t
end

module Eval :
  Eval
    with type Kernel.params = Params.t
    with type Inducing.t = mat
    with type Input.t = vec
    with type Inputs.t = mat


(* Derivatives *)

(* module Proj_hyper : sig type t = private { big_dim : int; small_dim : int } end *)
module Proj_hyper : sig type t = { big_dim : int; small_dim : int } end
module Dim_hyper : sig type t = private int end
module Inducing_hyper : sig type t = private { ind : int; dim : int } end

module Hyper_repr : sig
  type t =
    [
    | `Log_sf2
    | `Proj of Proj_hyper.t
    | `Log_hetero_skedasticity of Dim_hyper.t
    | `Inducing_hyper of Inducing_hyper.t
    | `Log_multiscales of Inducing_hyper.t
    ]
end

module Deriv :
  Deriv
    with module Eval = Eval
    with type Hyper.t = Hyper_repr.t
