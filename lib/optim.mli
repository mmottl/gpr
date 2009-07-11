open Lacaml.Impl.D

open Interfaces

(** Hyper parameter optimization with the GNU Scientific Library *)
module Gsl : sig

  exception Optim_exception of exn

  (** Ordinary hyper parameter optimization *)
  module Make (Spec : Sigs.Deriv) : sig
    open Spec

    val train :
      ?step : float ->
      ?tol : float ->
      ?epsabs : float ->
      ?report_trained_model : (iter : int -> Eval.Trained.t -> unit) ->
      ?report_gradient_norm : (iter : int -> float -> unit) ->
      ?kernel : Eval.Spec.Kernel.t ->
      ?sigma2 : float ->
      ?inducing : Eval.Spec.Inducing.t ->
      ?n_rand_inducing : int ->
      inputs : Eval.Spec.Inputs.t ->
      targets : vec ->
      unit ->
      Eval.Trained.t
  end


  (** SPGP *)
  module Make_SPGP
    (Deriv : Sigs.Deriv)
    (Spec : Specs.SPGP
      with module Eval = Deriv.Eval.Spec
      with module Deriv = Deriv.Deriv.Spec)
    :
      Sigs.SPGP
        with module Eval = Deriv.Eval
        with module Deriv = Deriv
        with module SPGP.Spec = Spec
end
