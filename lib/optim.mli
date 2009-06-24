open Lacaml.Impl.D

open Interfaces

(* Ordinary hyper parameter optimization *)

module Make (Spec : Sigs.Deriv) : sig
  open Spec

  module Solution : sig
    type t = {
      sigma2 : float;
      kernel : Eval.Spec.Kernel.t;
      coeffs : vec;
      log_evidence : float;
    }
  end

  val solve :
    ?kernel : Eval.Spec.Kernel.t ->
    ?sigma2 : float ->
    ?inducing : Eval.Spec.Inducing.t ->
    ?n_rand_inducing : int ->
    inputs : Eval.Spec.Inputs.t ->
    targets : vec ->
    Solution.t
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
