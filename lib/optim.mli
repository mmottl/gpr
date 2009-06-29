open Lacaml.Impl.D

open Interfaces

(** Hyper parameter optimization with the GNU Scientific Library *)
module Gsl : sig

  (** Ordinary hyper parameter optimization *)
  module Make (Spec : Sigs.Deriv) : sig
    open Spec

    val train :
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
