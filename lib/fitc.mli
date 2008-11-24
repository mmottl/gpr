open Interfaces

module type Spec = sig
  module Eval_spec : Inducing_input_gpr.Specs.Eval

  val get_sigma2 : Eval_spec.kernel -> float
  val jitter : float
end

module type Sig = functor (Spec : Spec) ->
  Inducing_input_gpr.Sigs.Eval with module Spec = Spec.Eval_spec

module Make_FITC : Sig
module Make_FIC : Sig
module Make_variational_FITC : Sig
module Make_variational_FIC : Sig

module Make (Spec : Spec) : sig
  module type Sig =
    Inducing_input_gpr.Sigs.Eval with module Spec = Spec.Eval_spec

  module FITC : Sig

  module FIC :
    Sig
      with module Inducing = FITC.Inducing
      with module Input = FITC.Input
      with module Inputs = FITC.Inputs
      with module Model = FITC.Model
      with module Trained = FITC.Trained
      with module Weights = FITC.Weights
      with module Mean = FITC.Mean
      with module Means = FITC.Means
      with module Variance = FITC.Variance
      with module Variances = FITC.Variances
      with module Sampler = FITC.Sampler

  module Variational_FITC :
    Sig
      with module Inducing = FITC.Inducing
      with module Input = FITC.Input
      with module Inputs = FITC.Inputs

  module Variational_FIC :
    Sig
      with module Inducing = Variational_FITC.Inducing
      with module Input = Variational_FITC.Input
      with module Inputs = Variational_FITC.Inputs
      with module Model = Variational_FITC.Model
      with module Trained = Variational_FITC.Trained
      with module Weights = Variational_FITC.Weights
      with module Mean = Variational_FITC.Mean
      with module Means = Variational_FITC.Means
      with module Variance = Variational_FITC.Variance
      with module Variances = Variational_FITC.Variances
      with module Sampler = Variational_FITC.Sampler
end
