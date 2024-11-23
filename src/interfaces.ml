(* OCaml-GPR - Gaussian Processes for OCaml

   Copyright © 2009- Markus Mottl <markus.mottl@gmail.com>

   This library is free software; you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the Free
   Software Foundation; either version 2.1 of the License, or (at your option)
   any later version.

   This library is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
   details.

   You should have received a copy of the GNU Lesser General Public License
   along with this library; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA *)

open Core
open Lacaml.D
open Utils

(** {2 Representations of (sparse) derivative matrices} *)

module Sparse_indices = Int_vec
(** Representation of indices into sparse matrices *)

type common_mat_deriv =
  [ `Dense of mat
  | `Sparse_rows of mat * Sparse_indices.t
  | `Const of float
  | `Factor of float ]
(** Derivative representations for both symmetric and unsymmetric matrices.

    - Dense: matrix is dense.
    - Sparse_rows: matrix is zero everywhere except for rows whose index is
      stored in the sparse index argument. The rows in the matrix correspond to
      the given indices.
    - Const: matrix is constant everywhere.
    - Factor: matrix is the non-derived matrix times the given factor (useful
      with exponential functions). *)

type mat_deriv = [ common_mat_deriv | `Sparse_cols of mat * Sparse_indices.t ]
(** Only general matrices support sparse column representations.

    - Sparse_cols: matrix is zero everywhere except for columns whose index is
      stored in the sparse index argument. The columns in the matrix correspond
      to the given indices. *)

type symm_mat_deriv =
  [ common_mat_deriv | `Diag_vec of vec | `Diag_const of float ]
(** Only symmetric (square) matrices support diagonal vectors and diagonal
    constants as derivatives.

    - Diag_vec: matrix is zero everywhere except for the diagonal whose values
      are given in the argument.
    - Diag_const: matrix is zero everywhere except for the diagonal whose values
      are set to the given constant.

    Note that sparse rows do not need to compute or store all elements for
    symmetric matrices. Entries that have already appeared in previous rows by
    symmetry can be left uninitialized. *)

type diag_deriv =
  [ `Vec of vec
  | `Sparse_vec of vec * Sparse_indices.t
  | `Const of float
  | `Factor of float ]
(** Derivatives of diagonal matrices.

    - Vec: the derivatives of the diagonal given in a dense vector.
    - Sparse_vec: matrix is zero everywhere except at those indices along the
      diagonal that are mentioned in the sparse indices argument. The element
      associated with such an index is stored in the vector argument.
    - Const: the derivative of the diagonal matrix is a constant.
    - Factor: the derivative of the diagonal is the the non-derived diagonal
      matrix times the given factor (useful with exponential functions). *)

(** Specifications of covariance functions (= kernels) and their derivatives *)
module Specs = struct
  (** Signature of kernels and their parameters *)
  module type Kernel = sig
    type t
    (** Type of kernel *)

    type params
    (** Type of kernel parameters *)

    val create : params -> t
    (** [create params]

        @return kernel given parameters [params]. *)

    val get_params : t -> params
    (** [get_params kernel]

        @return parameters used to parameterize [kernel]. *)
  end

  (** Evaluation of covariance functions *)
  module type Eval = sig
    module Kernel : Kernel
    (** Kernel used for evaluation *)

    (** Signature for evaluating inducing inputs *)
    module Inducing : sig
      type t

      val get_n_points : t -> int
      (** [get_n_points inducing]

          @return number of inducing points. *)

      val calc_upper : Kernel.t -> t -> mat
      (** [calc_upper kernel inducing]

          @return
            upper triangle of covariance matrix of [inducing] inputs given
            [kernel]. *)
    end

    (** Signature for evaluating single inputs *)
    module Input : sig
      type t
      (** Type of input point *)

      val eval : Kernel.t -> t -> Inducing.t -> vec
      (** [eval kernel input inducing]

          @return
            (row) vector of covariance evaluations between [input] and
            [inducing] inputs given [kernel]. *)

      val weighted_eval : Kernel.t -> t -> Inducing.t -> coeffs:vec -> float
      (** [weighted_eval kernel input inducing ~coeffs]

          @return
            [coeff]-weighted sum of covariances between [input] and [inducing]
            inputs given [kernel]. *)

      val eval_one : Kernel.t -> t -> float
      (** [eval_one kernel point]

          @return variance of [point] given [kernel]. *)
    end

    (** Signature for evaluating multiple inputs *)
    module Inputs : sig
      type t
      (** Type of input points *)

      val create : Input.t array -> t
      (** [create inputs]

          @return inputs given an array of single [inputs]. *)

      val get_n_points : t -> int
      (** [get_n_points inputs]

          @return number of input points. *)

      val choose_subset : t -> Int_vec.t -> t
      (** [choose_subset inputs indexes]

          @return subset of input points from [inputs] having [indexes]. *)

      val create_inducing : Kernel.t -> t -> Inducing.t
      (** [create_inducing kernel inputs]

          @return inducing points made from [inputs] and given [kernel]. *)

      val create_default_kernel_params : t -> n_inducing:int -> Kernel.params
      (** [create_default_kernel_params inputs ~n_inducing]

          @return
            default kernel parameters to be used with [n_inducing] inducing
            points and [inputs]. *)

      val calc_upper : Kernel.t -> t -> mat
      (** [calc_upper kernel inputs]

          @return
            upper triangle of covariance matrix of [inputs] given [kernel]. *)

      val calc_diag : Kernel.t -> t -> vec
      (** [calc_diag kernel inputs]

          @return diagonal of covariance matrix of [inputs] given [kernel]. *)

      val calc_cross : Kernel.t -> inputs:t -> inducing:Inducing.t -> mat
      (** [calc_cross kernel ~inputs ~inducing]

          @return
            cross-covariance matrix of [inputs] (indexing rows) and [inducing]
            points (indexing columns). *)

      val weighted_eval :
        Kernel.t -> inputs:t -> inducing:Inducing.t -> coeffs:vec -> vec
      (** [weighted_eval kernel ~inputs ~inducing ~coeffs]

          @return
            vector of [coeff]-weighted sums of covariances between [inputs] and
            [inducing] inputs given [kernel]. *)
    end
  end

  (** Derivatives of covariance functions *)
  module type Deriv = sig
    module Eval : Eval
    (** Derivatives always require evaluation functions *)

    (** Hyper parameters that have derivatives *)
    module Hyper : sig
      type t
      (** Type of hyper parameter *)

      val get_all : Eval.Kernel.t -> Eval.Inducing.t -> Eval.Inputs.t -> t array
      (** [get_all kernel inducing inputs]

          @return
            array of all hyper parameters of [kernel] and/or ([inducing])
            [inputs] for which derivatives can be computed. *)

      val get_value :
        Eval.Kernel.t -> Eval.Inducing.t -> Eval.Inputs.t -> t -> float
      (** [get_value kernel inducing inputs hyper]

          @return
            value of hyper parameter [hyper] of [kernel] and/or ([inducing])
            [inputs]. *)

      val set_values :
        Eval.Kernel.t ->
        Eval.Inducing.t ->
        Eval.Inputs.t ->
        t array ->
        vec ->
        Eval.Kernel.t * Eval.Inducing.t * Eval.Inputs.t
      (** [set_values kernel inducing inputs hypers values]

          @return
            triple of [(kernel, inducing, inputs)] in which [hypers] have been
            substituted with [values] position-wise. *)
    end

    (** Derivatives of the covariance matrix of inducing inputs *)
    module Inducing : sig
      type upper
      (** Representation of precomputed data for calculating the upper triangle
          of the derivative of the covariance matrix of inducing inputs. *)

      val calc_shared_upper : Eval.Kernel.t -> Eval.Inducing.t -> mat * upper
      (** [calc_shared_upper kernel inducing]

          @return
            the pair [(eval, upper)], where [eval] is the upper triangle of the
            covariance matrix of inducing inputs for [kernel], and [upper] is
            the precomputed data needed for taking derivatives. *)

      val calc_deriv_upper : upper -> Hyper.t -> symm_mat_deriv
      (** [calc_deriv_upper upper hyper]

          @return
            the derivative of the (symmetric) covariance matrix of inducing
            inputs given precomputed data [upper] and the [hyper]-variable. *)
    end

    (** Derivatives of the (cross-) covariance matrix of inputs. *)
    module Inputs : sig
      type diag
      (** Representation of precomputed data for calculating the derivative of
          the diagonal of the covariance matrix of inputs. *)

      type cross
      (** Representation of precomputed data for calculating the derivative of
          the cross-covariance matrix between inputs and inducing inputs. *)

      val calc_shared_diag : Eval.Kernel.t -> Eval.Inputs.t -> vec * diag
      (** [calc_shared_diag kernel inputs]

          @return
            the pair [(eval, diag)], where [eval] is the diagonal of the
            covariance matrix of [inputs] for [kernel], and [diag] is the
            precomputed data needed for taking derivatives. *)

      val calc_shared_cross :
        Eval.Kernel.t ->
        inputs:Eval.Inputs.t ->
        inducing:Eval.Inducing.t ->
        mat * cross
      (** [calc_shared_cross kernel ~inputs ~inducing]

          @return
            the pair [(eval,
          cross)], where [eval] is the
            cross-covariance matrix of inputs and inducing inputs for [kernel],
            and [diag] is the precomputed data needed for taking derivatives. *)

      val calc_deriv_diag : diag -> Hyper.t -> diag_deriv
      (** [calc_deriv_diag diag hyper]

          @return
            the derivative of the diagonal of the covariance matrix of inputs
            given precomputed data [diag] and the [hyper]-variable. *)

      val calc_deriv_cross : cross -> Hyper.t -> mat_deriv
      (** [calc_deriv_cross cross hyper]

          @return
            the derivative of the cross-covariance matrix of the inputs and
            inducing inputs given precomputed data [cross] and the
            [hyper]-variable. *)
    end
  end

  (** Derivatives of inputs for global optimization. *)
  module type Optimizer = sig
    module Eval : Eval
    (** Derivatives always require evaluation functions *)

    (** Input parameters that have derivatives *)
    module Var : sig
      type t
      (** Type of input parameter *)
    end

    module Input : sig
      val get_vars : Eval.Input.t -> Var.t array
      (** [get_vars input]

          @return
            array of all input parameters for which derivatives can be computed
            given [input]. *)

      val get_value : Eval.Input.t -> Var.t -> float
      (** [get_value input var]

          @return value of input parameter [var] for [input]. *)

      val set_values : Eval.Input.t -> Var.t array -> vec -> Eval.Input.t
      (** [set_values input vars values]

          @return
            input in which [vars] have been substituted with [values]
            position-wise. *)
    end

    module Inputs : sig
      val get_vars : Eval.Inputs.t -> Var.t array
      (** [get_vars inputs]

          @return
            array of all input parameters for which derivatives can be computed
            given [inputs]. *)

      val get_value : Eval.Inputs.t -> Var.t -> float
      (** [get_value inputs var]

          @return value of input parameter [var] for [inputs]. *)

      val set_values : Eval.Inputs.t -> Var.t array -> vec -> Eval.Inputs.t
      (** [set_values inputs vars values]

          @return
            inputs in which [vars] have been substituted with [values]
            position-wise. *)
    end
  end
end

(** Signatures for learning sparse Gaussian processes with inducing inputs *)
module Sigs = struct
  (** Modules for learning without derivatives of covariance functions. *)
  module type Eval = sig
    module Spec : Specs.Eval
    (** Specification of covariance function *)

    (** Evaluating inducing inputs *)
    module Inducing : sig
      type t
      (** Type of inducing inputs *)

      val choose_n_first_inputs :
        Spec.Kernel.t -> Spec.Inputs.t -> n_inducing:int -> Spec.Inducing.t
      (** [choose_n_first_inputs kernel inputs ~n_inducing]

          @return
            the first [n_inducing] inputs in [inputs] as inducing points given
            [kernel]. *)

      val choose_n_random_inputs :
        ?rnd_state:Random.State.t ->
        Spec.Kernel.t ->
        Spec.Inputs.t ->
        n_inducing:int ->
        Spec.Inducing.t
      (** [choose_n_random_inputs ?rnd_state kernel inputs ~n_inducing]

          @return
            [n_inducing] random inputs in [inputs] as inducing points given
            [kernel] and (optional) random state [rnd_state].

          @param rnd_state default = default used by the Random module *)

      val calc : Spec.Kernel.t -> Spec.Inducing.t -> t
      (** [calc kernel inducing_points]

          @return
            inducing inputs (= precomputed data) prepared using
            [inducing_points] and [kernel]. *)

      val get_points : t -> Spec.Inducing.t
      (** [get_points kernel inducing]

          @return
            inducing points associated with the prepared [inducing] inputs. *)
    end

    (** Evaluating single inputs *)
    module Input : sig
      type t
      (** Type of single input *)

      val calc : Inducing.t -> Spec.Input.t -> t
      (** [calc inducing point]

          @return
            input (= precomputed data) prepared using [inducing] inputs and
            input [point]. *)
    end

    (** Evaluating (multiple) inputs *)
    module Inputs : sig
      type t
      (** Type of (multiple) inputs *)

      val create_default_kernel :
        Spec.Inputs.t -> n_inducing:int -> Spec.Kernel.t
      (** [create_default_kernel points]

          @return
            a default kernel given input [points] and [n_inducing] inducing
            inputs. *)

      val calc : Spec.Inputs.t -> Inducing.t -> t
      (** [create points inducing]

          @return
            inputs (= precomputed data) prepared using [inducing] inputs and
            input [points]. *)

      val get_points : t -> Spec.Inputs.t
      (** [get_points kernel inputs]

          @return points associated with the prepared [inputs]. *)
    end

    (** (Untrained) model - does not require targets *)
    module Model : sig
      type t
      (** Type of models *)

      type co_variance_coeffs
      (** Type of covariance coefficients *)

      val calc : Inputs.t -> sigma2:float -> t
      (** [calc inputs ~sigma2]

          @return
            model given [inputs] and noise level [sigma2] (= variance, i.e.
            squared standard deviation). *)

      val update_sigma2 : t -> float -> t
      (** [update_sigma2 model sigma2]

          @return model by updating [model] with new noise level [sigma2]. *)

      val calc_log_evidence : t -> float
      (** [calc_log_evidence model]

          @return
            the contribution to the log evidence (= log marginal likelihood) of
            [model]. *)

      val calc_co_variance_coeffs : t -> co_variance_coeffs
      (** [calc_co_variance_coeffs model]

          @return
            the coefficients required for computing posterior (co-)variances for
            [model]. *)

      val get_kernel : t -> Spec.Kernel.t
      (** [get_kernel model]

          @return the kernel associated with [model]. *)

      val get_sigma2 : t -> float
      (** [get_sigma2 model]

          @return the noise level associated with [model]. *)

      val get_inputs : t -> Inputs.t
      (** [get_inputs model]

          @return the inputs associated with [model]. *)

      val get_inducing : t -> Inducing.t
      (** [get_inputs model]

          @return the inducing inputs associated with [model]. *)
    end

    (** Trained model - requires targets *)
    module Trained : sig
      type t
      (** Type of trained models *)

      val calc : Model.t -> targets:vec -> t
      (** [calc model ~targets]

          @return trained model given [model] and [targets]. *)

      val calc_mean_coeffs : t -> vec
      (** [calc_mean_coeffs trained]

          @return the vector of coefficients for computing posterior means. *)

      val calc_log_evidence : t -> float
      (** [calc_log_evidence trained]

          @return
            the log evidence for the trained model (includes contribution to log
            evidence by underlying model). *)

      val get_model : t -> Model.t
      (** [get_model trained]

          @return the model associated with the [trained] model. *)

      val get_targets : t -> vec
      (** [get_targets trained]

          @return targets used for training [trained]. *)
    end

    (** Statistics derived from trained models *)
    module Stats : sig
      type t = {
        n_samples : int;  (** Number of samples used for training *)
        target_variance : float;  (** Variance of targets *)
        sse : float;  (** Sum of squared errors *)
        mse : float;  (** Mean sum of squared errors *)
        rmse : float;  (** Root mean sum of squared errors *)
        smse : float;  (** Standardized mean squared error *)
        msll : float;  (** Mean standardized log loss *)
        mad : float;  (** Mean absolute deviation *)
        maxad : float;  (** Maximum absolute deviation *)
      }
      (** Type of full statistics *)

      val calc_n_samples : Trained.t -> int
      (** [calc_n_samples trained]

          @return number of samples used for training [trained]. *)

      val calc_target_variance : Trained.t -> float
      (** [calc_target_variance trained]

          @return variance of targets used for training [trained]. *)

      val calc_sse : Trained.t -> float
      (** [calc_sse trained]

          @return the sum of squared errors of the [trained] model. *)

      val calc_mse : Trained.t -> float
      (** [calc_mse trained]

          @return the mean sum of squared errors of the [trained] model. *)

      val calc_rmse : Trained.t -> float
      (** [calc_sse trained]

          @return
            the root of the mean sum of squared errors of the [trained] model. *)

      val calc_smse : Trained.t -> float
      (** [calc_smse trained]

          @return
            the standardized mean squared error of the [trained] model. This is
            equivalent to the mean squared error divided by the target variance. *)

      val calc_msll : Trained.t -> float
      (** [calc_msll trained]

          @return
            the mean standardized log loss. This is equivalent to subtracting
            the log evidence of the trained model from the log evidence of a
            normal distribution fit to the targets, and dividing the result by
            the number of samples. *)

      val calc_mad : Trained.t -> float
      (** [calc_mad trained]

          @return the mean absolute deviation of the [trained] model. *)

      val calc_maxad : Trained.t -> float
      (** [calc_mad trained]

          @return the maximum absolute deviation of the [trained] model. *)

      val calc : Trained.t -> t
      (** [calc trained]

          @return
            the full set of statistics associated with the [trained] model. *)
    end

    (** Module for making mean predictions *)
    module Mean_predictor : sig
      type t
      (** Type of mean predictors *)

      val calc : Spec.Inducing.t -> coeffs:vec -> t
      (** [calc inducing_points ~coeffs]

          @return
            a mean predictor given [inducing_points] and coefficients [coeffs]. *)

      val calc_trained : Trained.t -> t
      (** [calc_trained trained]

          @return a mean predictor given the [trained] model. *)

      val get_inducing : t -> Spec.Inducing.t
      (** [get_inducing mean_predictor]

          @return inducing points associated with [mean_predictor]. *)

      val get_coeffs : t -> vec
      (** [get_coeffs mean_predictor]

          @return coefficients associated with [mean_predictor]. *)
    end

    (** Posterior mean for a single input *)
    module Mean : sig
      type t
      (** Type of mean *)

      val calc : Mean_predictor.t -> Input.t -> t
      (** [calc mean_predictor input]

          @return mean for [input] given [mean_predictor]. *)

      val get : t -> float
      (** [get mean]

          @return the mean as a float. *)
    end

    (** Posterior means for (multiple) inputs *)
    module Means : sig
      type t
      (** Type of means *)

      val calc : Mean_predictor.t -> Inputs.t -> t
      (** [calc mean_predictor inputs]

          @return means for [inputs] given [mean_predictor]. *)

      val get : t -> vec
      (** [get means]

          @return the means as a vector. *)
    end

    (** Module for making (co-)variance predictions *)
    module Co_variance_predictor : sig
      type t
      (** Type of (co-)variance predictor *)

      val calc :
        Spec.Kernel.t -> Spec.Inducing.t -> Model.co_variance_coeffs -> t
      (** [calc kernel inducing_points co_variance_coeffs]

          @return
            (co-)variance predictor given [kernel], [inducing_points], and the
            (co-)variance coefficients [co_variance_coeffs]. *)

      val calc_model : Model.t -> t
      (** [calc_model model]

          @return (co-)variance predictor given the (untrained) [model]. *)
    end

    (** Posterior variance for a single input *)
    module Variance : sig
      type t
      (** Type of variance *)

      val calc : Co_variance_predictor.t -> sigma2:float -> Input.t -> t
      (** [calc co_variance_predictor ~sigma2 input]

          @return
            variance for [input] given [mean_predictor] and noise level
            [sigma2]. *)

      val get : ?predictive:bool -> t -> float
      (** [get ?predictive variance]

          @return
            the [variance] as a float. If [predictive] is [true], then the noise
            level will be added.

          @param predictive default = [true] *)
    end

    (** Posterior variances for (multiple) inputs *)
    module Variances : sig
      type t
      (** Type of variances *)

      val calc_model_inputs : Model.t -> t
      (** [calc_model_inputs model]

          @return variances for all inputs used in [model]. *)

      val calc : Co_variance_predictor.t -> sigma2:float -> Inputs.t -> t
      (** [calc co_variance_predictor ~sigma2 inputs]

          @return
            variances for [inputs] given [co_variance_predictor] and noise level
            [sigma2]. *)

      val get : ?predictive:bool -> t -> vec
      (** [get ?predictive variances]

          @return
            the [variances] as a vector. If [predictive] is [true], then the
            noise level will be added.

          @param predictive default = [true] *)
    end

    (** Posterior covariances *)
    module Covariances : sig
      type t
      (** Type of covariances *)

      val calc_model_inputs : Model.t -> t
      (** [calc_model_inputs model]

          @return
            covariances for all inputs used in [model]. This may be extremely
            expensive (O(N^2)) for large numbers of model inputs. *)

      val calc : Co_variance_predictor.t -> sigma2:float -> Inputs.t -> t
      (** [calc co_variance_predictor ~sigma2 inputs]

          @return
            posterior covariances for [inputs] given [co_variance_predictor] and
            noise level [sigma2]. This may be extremely expensive (O(N^2)) for
            large numbers of inputs. *)

      val get : ?predictive:bool -> t -> mat
      (** [get ?predictive covariances]

          @return
            the [covariances] as a matrix. If [predictive] is [true], then the
            noise level will be added (to the diagonal only).

          @param predictive default = [true] *)

      val get_variances : t -> Variances.t
      (** [get_variances covariances]

          @return the variances in [covariances]. *)
    end

    (** Module for sampling single points from the posterior distribution *)
    module Sampler : sig
      type t
      (** Type of sampler *)

      val calc : ?predictive:bool -> Mean.t -> Variance.t -> t
      (** [calc ?predictive mean variance]

          @return
            sampler given [mean] and [variance]. If [predictive] is true, the
            samples will be noisy. *)

      val sample : ?rng:Gsl.Rng.t -> t -> float
      (** [sample ?rng sampler]

          @return
            a sample from the posterior distribution given [sampler] and GSL
            random number generator [rng].

          @param rng default = GSL default *)

      val samples : ?rng:Gsl.Rng.t -> t -> n:int -> vec
      (** [samples ?rng sampler ~n]

          @return [n] samples from the posterior distribution given [sampler].

          @param rng default = GSL default *)
    end

    (** Module for sampling (multiple) points from the posterior distribution
        accounting for their covariance *)
    module Cov_sampler : sig
      type t
      (** Type of covariance sampler *)

      val calc : ?predictive:bool -> Means.t -> Covariances.t -> t
      (** [calc ?predictive mean variance]

          @return
            sampler given [means] and [covariances]. If [predictive] is true,
            the samples will be noisy. *)

      val sample : ?rng:Gsl.Rng.t -> t -> vec
      (** [sample ?rng sampler]

          @return
            a sample vector from the posterior distribution given [sampler] and
            GSL random number generator [rng].

          @param rng default = GSL default *)

      val samples : ?rng:Gsl.Rng.t -> t -> n:int -> mat
      (** [samples ?rng sampler ~n]

          @return
            matrix of [n] sample vectors (stored row-wise) from the posterior
            distribution given [sampler].

          @param rng default = GSL default *)
    end
  end

  (** Modules for learning with derivatives of the log evidence (evidence
      maximization framework) *)
  module type Deriv = sig
    module Eval : Eval
    (** Sub-modules for learning without derivatives. *)

    (** Sub-modules for learning with derivatives. *)
    module Deriv : sig
      module Spec : Specs.Deriv with module Eval = Eval.Spec
      (** Specification of covariance function derivatives *)

      (** Module for inducing inputs with derivatives *)
      module Inducing : sig
        type t
        (** Type of inducing inputs with derivatives *)

        val calc : Eval.Spec.Kernel.t -> Eval.Spec.Inducing.t -> t
        (** [calc kernel inducing_points]

            @return
              inducing inputs with derivative information given [kernel] and
              [inducing_points]. *)

        val calc_eval : t -> Eval.Inducing.t
        (** [calc_eval inducing]

            @return inducing inputs without derivative information. *)
      end

      (** Module for inputs with derivatives *)
      module Inputs : sig
        type t
        (** Type of inputs with derivatives *)

        val calc : Inducing.t -> Eval.Spec.Inputs.t -> t
        (** [calc inducing points]

            @return
              inputs with derivative information given [inducing] inputs and
              input [points]. *)

        val calc_eval : t -> Eval.Inputs.t
        (** [calc_eval inputs]

            @return inputs without derivative information. *)
      end

      (** (Untrained) model with derivative information *)
      module Model : sig
        type t
        (** Type of models with derivatives *)

        type hyper_t
        (** Type of models for general hyper parameters *)

        val calc : Inputs.t -> sigma2:float -> t
        (** [calc inputs ~sigma2]

            @return
              model with derivative information given [inputs] and noise level
              [sigma2]. *)

        val update_sigma2 : t -> float -> t
        (** [update_sigma2 model sigma2]

            @return
              model with derivative information by updating [model] with new
              noise level [sigma2]. *)

        val calc_eval : t -> Eval.Model.t
        (** [calc_eval model]

            @return model without derivative information given [model]. *)

        val calc_log_evidence_sigma2 : t -> float
        (** [calc_log_evidence_sigma2 model]

            @return
              the derivative of the log evidence of [model] with respect to the
              noise level (sigma2). *)

        val prepare_hyper : t -> hyper_t
        (** [prepare_hyper model]

            @return
              the model prepared for calculating derivatives for arbitrary hyper
              parameters. *)

        val calc_log_evidence : hyper_t -> Spec.Hyper.t -> float
        (** [calc_log_evidence hyper_t hyper]

            @return
              the derivative of the log evidence given prepared model [hyper_t]
              with respect to the [hyper] variable. *)
      end

      (** Trained model with derivative information *)
      module Trained : sig
        type t
        (** Type of trained models with derivatives *)

        type hyper_t
        (** Type of trained models for general hyper parameters *)

        val calc : Model.t -> targets:vec -> t
        (** [calc model ~targets]

            @return
              trained model with derivative information given the untrained
              [model] and [targets]. *)

        val calc_eval : t -> Eval.Trained.t
        (** [calc_eval trained]

            @return
              trained model without derivative information given [trained]. *)

        val calc_log_evidence_sigma2 : t -> float
        (** [calc_log_evidence_sigma2 trained]

            @return
              the derivative of the log evidence for the [trained] model with
              respect to the noise level (sigma2). This includes the
              contribution to the derivative by [model]. *)

        val prepare_hyper : t -> hyper_t
        (** [prepare_hyper trained]

            @return
              the trained model prepared for calculating derivatives for
              arbitrary hyper parameters. *)

        val calc_log_evidence : hyper_t -> Spec.Hyper.t -> float
        (** [calc_log_evidence hyper_t hyper]

            @return
              the derivative of the log evidence given prepared, trained model
              [hyper_t] with respect to the [hyper] variable. *)
      end

      (** Module for testing derivative code *)
      module Test : sig
        val check_deriv_hyper :
          ?eps:float ->
          ?tol:float ->
          Eval.Spec.Kernel.t ->
          Eval.Spec.Inducing.t ->
          Eval.Spec.Inputs.t ->
          Spec.Hyper.t ->
          unit
        (** [check_deriv_hyper ?eps ?tol kernel inducing_points points hyper]
            will raise [Failure] if the derivative code provided in the
            specification of the covariance function given parameter [hyper],
            the [kernel], [inducing_points] and input [points] exceeds the
            tolerance [tol] when compared to finite differences using epsilon
            [eps]. The failure exception will contain details on which
            derivative matrix was incorrect and indicate the matrix element.

            @param eps default = [1e-8]
            @param tol default = [1e-2] *)

        val self_test :
          ?eps:float ->
          ?tol:float ->
          Eval.Spec.Kernel.t ->
          Eval.Spec.Inducing.t ->
          Eval.Spec.Inputs.t ->
          sigma2:float ->
          targets:vec ->
          [ `Sigma2 | `Hyper of Spec.Hyper.t ] ->
          unit
        (** [self_test ?eps ?tol kernel inducing_points points ~sigma2 ~targets
            hyper]
            will raise [Failure] if the internal derivative code for the log
            evidence given parameter [hyper], the [kernel], [inducing_points],
            input [points], noise level [sigma2] and [targets] exceeds the
            tolerance [tol] when compared to finite differences using epsilon
            [eps].

            @param eps default = [1e-8]
            @param tol default = [1e-2] *)
      end

      (** Optimization module for evidence maximization *)
      module Optim : sig
        (** Optimization with the GNU Scientific library (GSL) *)
        module Gsl : sig
          exception Optim_exception of exn
          (** [Optim_exception exn] is raised when an internal exception occurs,
              e.g. because GSL fails, or because a callback raised it. *)

          val train :
            ?step:float ->
            ?tol:float ->
            ?epsabs:float ->
            ?report_trained_model:(iter:int -> Eval.Trained.t -> unit) ->
            ?report_gradient_norm:(iter:int -> float -> unit) ->
            ?kernel:Eval.Spec.Kernel.t ->
            ?sigma2:float ->
            ?inducing:Eval.Spec.Inducing.t ->
            ?n_rand_inducing:int ->
            ?learn_sigma2:bool ->
            ?hypers:Spec.Hyper.t array ->
            inputs:Eval.Spec.Inputs.t ->
            targets:vec ->
            unit ->
            Eval.Trained.t
          (** [train ?step ?tol ?epsabs ?report_trained_model
              ?report_gradient_norm ?kernel ?sigma2 ?inducing ?n_rand_inducing
              ?learn_sigma2 ?hypers ~inputs ~targets ()]
              takes the optional initial optimizer step size [step], the
              optimizer line search tolerance [tol], the minimum gradient norm
              [epsabs] to achieve by the optimizer, callbacks for reporting
              intermediate results [report_trained_model] and
              [report_gradient_norm], an optional [kernel], noise level
              [sigma2], inducing inputs [inducing], number of randomly chosen
              inducing inputs [n_rand_inducing], a flag for whether the noise
              level should be learnt [learn_sigma2], an array of optional hyper
              parameters [hypers] which should be optimized, and the [inputs]
              and [targets].

              @return
                the trained model obtained by evidence maximization (= type II
                maximum likelihood).

              @param step default = [1e-1]
              @param tol default = [1e-1]
              @param epsabs default = [1e-1]
              @param report_trained_model default = ignored
              @param report_gradient_norm default = ignored
              @param kernel default = default kernel computed from specification
              @param sigma2 default = target variance
              @param inducing default = randomly selected subset of inputs
              @param n_rand_inducing default = 10% of inputs, at most 1000
              @param learn_sigma2 default = [true]
              @param hypers default = all hyper parameters *)
        end

        module SGD : sig
          type t

          val create :
            ?tau:float ->
            ?eta0:float ->
            ?step:int ->
            ?kernel:Eval.Spec.Kernel.t ->
            ?sigma2:float ->
            ?inducing:Eval.Spec.Inducing.t ->
            ?n_rand_inducing:int ->
            ?learn_sigma2:bool ->
            ?hypers:Spec.Hyper.t array ->
            inputs:Eval.Spec.Inputs.t ->
            targets:vec ->
            unit ->
            t

          val step : t -> t
          val gradient_norm : t -> float
          val get_trained : t -> Eval.Trained.t
          val get_eta : t -> float
          val get_step : t -> int

          val test :
            ?epsabs:float -> ?max_iter:int -> ?report:(t -> unit) -> t -> t
        end

        module SMD : sig
          type t

          val create :
            ?eps:float ->
            ?lambda:float ->
            ?mu:float ->
            ?eta0:vec ->
            ?nu0:vec ->
            ?kernel:Eval.Spec.Kernel.t ->
            ?sigma2:float ->
            ?inducing:Eval.Spec.Inducing.t ->
            ?n_rand_inducing:int ->
            ?learn_sigma2:bool ->
            ?hypers:Spec.Hyper.t array ->
            inputs:Eval.Spec.Inputs.t ->
            targets:vec ->
            unit ->
            t

          val step : t -> t
          val gradient_norm : t -> float
          val get_trained : t -> Eval.Trained.t
          val get_eta : t -> vec
          val get_nu : t -> vec

          val test :
            ?epsabs:float -> ?max_iter:int -> ?report:(t -> unit) -> t -> t
        end
      end

      (* (** Online learning *) module Online : sig type t

         val sgd : ?capacity : int -> ?eta0 : float -> ?tau : float ->
         Spec.Eval.Kernel.t -> t

         val smd : ?capacity : int -> ?eta0 : vec -> ?mu : float -> ?lam : float
         -> Spec.Eval.Kernel.t -> t

         val train : t -> Spec.Eval.Input.t -> target : float -> t

         val calc_mean_predictor : t -> Eval.Mean_predictor.t val
         calc_co_variance_predictor : t -> Eval.Co_variance_predictor.t end *)
    end
  end

  (** Modules for global optimization *)
  module type Optimizer = sig
    module Eval : Eval
    (** Sub-modules for learning without derivatives. *)

    (** Sub-modules for global optimization. *)
    module Optimizer : sig
      module Spec : Specs.Optimizer with module Eval = Eval.Spec

      type t

      val create : ?max_memory:int -> Spec.Eval.Kernel.t -> t
      val learn : t -> (Spec.Eval.Input.t * float) array -> t
      val calc_mpi_criterion : t -> Spec.Eval.Input.t -> float
      val calc_mpi_deriv : t -> Spec.Eval.Input.t
    end
  end
end
