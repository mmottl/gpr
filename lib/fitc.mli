open Interfaces

module type Spec = sig
  module Kernel : Kernel

  val get_sigma2 : Kernel.t -> float
  val jitter : float
end

module type Sig = functor (Spec : Spec) ->
  Inducing_input_gpr with module Kernel = Spec.Kernel

module Make_FITC : Sig
module Make_FIC : Sig
module Make_variational_FITC : Sig
module Make_variational_FIC : Sig

module Make (Spec : Spec) : sig
  module type Sig = Inducing_input_gpr with module Kernel = Spec.Kernel

  module FITC : Sig

  module FIC :
    Sig
      with module Inducing = FITC.Inducing
      with module Induced = FITC.Induced
      with module Induceds = FITC.Induceds
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
      with module Induced = FITC.Induced
      with module Induceds = FITC.Induceds

  module Variational_FIC :
    Sig
      with module Inducing = Variational_FITC.Inducing
      with module Induced = Variational_FITC.Induced
      with module Induceds = Variational_FITC.Induceds
      with module Model = Variational_FITC.Model
      with module Trained = Variational_FITC.Trained
      with module Weights = Variational_FITC.Weights
      with module Mean = Variational_FITC.Mean
      with module Means = Variational_FITC.Means
      with module Variance = Variational_FITC.Variance
      with module Variances = Variational_FITC.Variances
      with module Sampler = Variational_FITC.Sampler
end
