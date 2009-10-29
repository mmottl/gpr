(* File: fitc_gp.mli

   OCaml-GPR - Gaussian Processes for OCaml

     Copyright (C) 2009-  Markus Mottl
     email: markus.mottl@gmail.com
     WWW:   http://www.ocaml.info

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*)

open Interfaces

(** Evaluation *)

module type Sig = functor (Spec : Specs.Eval) ->
  Sigs.Eval with module Spec = Spec

module Make_FITC : Sig
module Make_FIC : Sig
module Make_variational_FITC : Sig
module Make_variational_FIC : Sig

module Make (Spec : Specs.Eval) : sig
  module type Sig = Sigs.Eval with module Spec = Spec

  module FITC : Sig

  module FIC :
    Sig
      with module Inducing = FITC.Inducing
      with module Input = FITC.Input
      with module Inputs = FITC.Inputs
      with module Model = FITC.Model
      with module Trained = FITC.Trained
      with module Mean_predictor = FITC.Mean_predictor
      with module Co_variance_predictor = FITC.Co_variance_predictor
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
      with module Mean_predictor = Variational_FITC.Mean_predictor
      with module Co_variance_predictor = Variational_FITC.Co_variance_predictor
      with module Mean = Variational_FITC.Mean
      with module Means = Variational_FITC.Means
      with module Variance = Variational_FITC.Variance
      with module Variances = Variational_FITC.Variances
      with module Sampler = Variational_FITC.Sampler
end


(** Derivatives *)

module type Deriv_sig = functor (Spec : Specs.Deriv) ->
  Sigs.Deriv
    with module Eval.Spec = Spec.Eval
    with module Deriv.Spec = Spec

module Make_FITC_deriv : Deriv_sig
module Make_FIC_deriv : Deriv_sig
module Make_variational_FITC_deriv : Deriv_sig
module Make_variational_FIC_deriv : Deriv_sig

module Make_deriv (Spec : Specs.Deriv) : sig
  module type Sig = Sigs.Deriv
    with module Eval.Spec = Spec.Eval
    with module Deriv.Spec = Spec

  module FITC : Sig

  module FIC :
    Sig
      with module Eval.Inducing = FITC.Eval.Inducing
      with module Eval.Input = FITC.Eval.Input
      with module Eval.Inputs = FITC.Eval.Inputs
      with module Eval.Model = FITC.Eval.Model
      with module Eval.Trained = FITC.Eval.Trained
      with module Eval.Mean_predictor = FITC.Eval.Mean_predictor
      with module Eval.Co_variance_predictor = FITC.Eval.Co_variance_predictor
      with module Eval.Mean = FITC.Eval.Mean
      with module Eval.Means = FITC.Eval.Means
      with module Eval.Variance = FITC.Eval.Variance
      with module Eval.Variances = FITC.Eval.Variances
      with module Eval.Sampler = FITC.Eval.Sampler

      with module Deriv.Inducing = FITC.Deriv.Inducing
      with module Deriv.Inputs = FITC.Deriv.Inputs
      with module Deriv.Model = FITC.Deriv.Model
      with module Deriv.Trained = FITC.Deriv.Trained

  module Variational_FITC :
    Sig
      with module Eval.Inducing = FITC.Eval.Inducing
      with module Eval.Input = FITC.Eval.Input
      with module Eval.Inputs = FITC.Eval.Inputs
      with module Deriv.Inducing = FITC.Deriv.Inducing
      with module Deriv.Inputs = FITC.Deriv.Inputs

  module Variational_FIC :
    Sig
      with module Eval.Inducing = Variational_FITC.Eval.Inducing
      with module Eval.Input = Variational_FITC.Eval.Input
      with module Eval.Inputs = Variational_FITC.Eval.Inputs
      with module Eval.Model = Variational_FITC.Eval.Model
      with module Eval.Trained = Variational_FITC.Eval.Trained
      with module Eval.Mean_predictor = Variational_FITC.Eval.Mean_predictor
      with module Eval.Co_variance_predictor =
        Variational_FITC.Eval.Co_variance_predictor
      with module Eval.Mean = Variational_FITC.Eval.Mean
      with module Eval.Means = Variational_FITC.Eval.Means
      with module Eval.Variance = Variational_FITC.Eval.Variance
      with module Eval.Variances = Variational_FITC.Eval.Variances
      with module Eval.Sampler = Variational_FITC.Eval.Sampler

      with module Deriv.Inducing = Variational_FITC.Deriv.Inducing
      with module Deriv.Inputs = Variational_FITC.Deriv.Inputs
      with module Deriv.Model = Variational_FITC.Deriv.Model
      with module Deriv.Trained = Variational_FITC.Deriv.Trained
end
