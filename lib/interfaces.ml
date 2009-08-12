open Bigarray
open Lacaml.Impl.D

open Utils

module Sparse_indices = Int_vec

type common_mat_deriv = [
  | `Dense of mat
  | `Sparse_rows of mat * Sparse_indices.t
  | `Const of float
  | `Factor of float
]

type mat_deriv = [
  | common_mat_deriv
  | `Sparse_cols of mat * Sparse_indices.t
]

type symm_mat_deriv = [
  | common_mat_deriv
  | `Diag_vec of vec
  | `Diag_const of float
]

type diag_deriv = [
  | `Vec of vec
  | `Sparse_vec of vec * Sparse_indices.t
  | `Const of float
  | `Factor of float
]

module Specs = struct
  module type Kernel = sig
    type t
    type params

    val create : params -> t
    val get_params : t -> params
  end

  module type Eval = sig
    module Kernel : Kernel

    module Inducing : sig
      type t

      val calc_upper : Kernel.t -> t -> mat
      val get_n_points : t -> int
    end

    module Input : sig
      type t

      val eval : Kernel.t -> Inducing.t -> t -> vec
      val weighted_eval : Kernel.t -> Inducing.t -> coeffs : vec -> t -> float
      val eval_one : Kernel.t -> t -> float
    end

    module Inputs : sig
      type t

      val get_n_points : t -> int
      val choose_subset : t -> Int_vec.t -> t
      val create_inducing : Kernel.t -> t -> Inducing.t
      val create_default_kernel_params : n_inducing : int -> t -> Kernel.params

      val calc_upper : Kernel.t -> t -> mat
      val calc_diag : Kernel.t -> t -> vec
      val calc_cross : Kernel.t -> Inducing.t -> t -> mat

      val weighted_eval : Kernel.t -> Inducing.t -> coeffs : vec -> t -> vec
    end
  end

  module type Deriv = sig
    module Eval : Eval

    module Hyper : sig
      type t

      val extract : Eval.Kernel.t -> t array * vec
      val update : Eval.Kernel.t -> vec -> Eval.Kernel.t
    end

    module Inducing : sig
      type upper

      val calc_shared_upper : Eval.Kernel.t -> Eval.Inducing.t -> mat * upper
      val calc_deriv_upper : upper -> Hyper.t -> symm_mat_deriv
    end

    module Inputs : sig
      type diag
      type cross

      val calc_shared_diag : Eval.Kernel.t -> Eval.Inputs.t -> vec * diag

      val calc_shared_cross :
        Eval.Kernel.t -> Eval.Inducing.t -> Eval.Inputs.t -> mat * cross

      val calc_deriv_diag : diag -> Hyper.t -> diag_deriv
      val calc_deriv_cross : cross -> Hyper.t -> mat_deriv
    end
  end

  module type SPGP = sig
    module Eval : Eval
    module Deriv : Deriv with module Eval = Eval

    module Inducing_hypers : sig
      val extract : Deriv.Eval.Inducing.t -> Deriv.Hyper.t array * vec
      val update : Deriv.Eval.Inducing.t -> vec -> Deriv.Eval.Inducing.t
    end
  end
end

module Sigs = struct
  module type Eval = sig
    module Spec : Specs.Eval

    module Inducing : sig
      type t

      val choose_n_first_inputs :
        Spec.Kernel.t -> n_inducing : int -> Spec.Inputs.t -> Spec.Inducing.t

      val choose_n_random_inputs :
        ?rnd_state : Random.State.t ->
        Spec.Kernel.t ->
        n_inducing : int ->
        Spec.Inputs.t ->
        Spec.Inducing.t

      val calc : Spec.Kernel.t -> Spec.Inducing.t -> t

      val get_points : t -> Spec.Inducing.t
    end

    module Input : sig
      type t

      val calc : Inducing.t -> Spec.Input.t -> t
    end

    module Inputs : sig
      type t

      val create_default_kernel :
        n_inducing : int -> Spec.Inputs.t -> Spec.Kernel.t

      val calc : Inducing.t -> Spec.Inputs.t -> t
      val get_points : t -> Spec.Inputs.t
    end

    module Model : sig
      type t

      val calc : Inputs.t -> sigma2 : float -> t
      val update_sigma2 : t -> float -> t
      val calc_log_evidence : t -> float
      val calc_co_variance_coeffs : t -> mat

      val get_kernel : t -> Spec.Kernel.t
      val get_sigma2 : t -> float
      val get_inputs : t -> Inputs.t
      val get_inducing : t -> Inducing.t
    end

    module Trained : sig
      type t

      val calc : Model.t -> targets : vec -> t
      val calc_mean_coeffs : t -> vec
      val calc_log_evidence : t -> float
      val calc_rmse : t -> float

      val get_model : t -> Model.t
      val get_targets : t -> vec
    end

    module Mean_predictor : sig
      type t

      val calc :
        Spec.Kernel.t -> sigma2 : float -> Spec.Inducing.t -> coeffs : vec -> t

      val calc_trained : Trained.t -> t

      val get_kernel : t -> Spec.Kernel.t
      val get_sigma2 : t -> float
      val get_inducing : t -> Spec.Inducing.t
      val get_coeffs : t -> vec
    end

    module Mean : sig
      type t

      val calc_induced : Mean_predictor.t -> Input.t -> t
      val get : t -> float
    end

    module Means : sig
      type t

      val calc_induced : Mean_predictor.t -> Inputs.t -> t
      val get : t -> vec

      module Inducing : sig
        type t

        val calc : Trained.t -> t
        val get : t -> vec
      end
    end

    module Co_variance_predictor : sig
      type t

      val calc : Spec.Kernel.t -> Spec.Inducing.t -> coeffs : mat -> t
      val calc_model : Model.t -> t
    end

    module Variance : sig
      type t

      val calc_induced :
        Co_variance_predictor.t -> sigma2 : float -> Input.t -> t

      val get : ?predictive : bool -> t -> float
    end

    module Variances : sig
      type t

      val calc_model_inputs : Model.t -> t

      val calc_induced :
        Co_variance_predictor.t -> sigma2 : float -> Inputs.t -> t

      val get : ?predictive : bool -> t -> vec

      module Inducing : sig
        type t

        val calc : Model.t -> t
        val get : ?predictive : bool -> t -> vec
      end
    end

    module Covariances : sig
      type t

      val calc_model_inputs : Model.t -> t

      val calc_induced :
        Co_variance_predictor.t -> sigma2 : float -> Inputs.t -> t

      val get : ?predictive : bool -> t -> mat
      val get_variances : t -> Variances.t

      module Inducing : sig
        type t

        val calc : Model.t -> t
        val get : ?predictive : bool -> t -> mat
        val get_variances : t -> Variances.Inducing.t
      end
    end

    module Sampler : sig
      type t

      val calc : ?predictive : bool -> Mean.t -> Variance.t -> t
      val sample : ?rng : Gsl_rng.t -> t -> float
      val samples : ?rng : Gsl_rng.t -> t -> n : int -> vec
    end

    module Cov_sampler : sig
      type t

      val calc : ?predictive : bool -> Means.t -> Covariances.t -> t
      val sample : ?rng : Gsl_rng.t -> t -> vec
      val samples : ?rng : Gsl_rng.t -> t -> n : int -> mat
    end
  end

  module type Deriv = sig
    module Eval : Eval

    module Deriv : sig
      module Spec : Specs.Deriv with module Eval = Eval.Spec

      module Inducing : sig
        type t

        val calc : Eval.Spec.Kernel.t -> Eval.Spec.Inducing.t -> t
        val calc_eval : t -> Eval.Inducing.t
      end

      module Inputs : sig
        type t

        val calc : Inducing.t -> Eval.Spec.Inputs.t -> t
        val calc_eval : t -> Eval.Inputs.t
      end

      module Model : sig
        type t
        type hyper_t

        val calc : Inputs.t -> sigma2 : float -> t
        val update_sigma2 : t -> float -> t
        val calc_eval : t -> Eval.Model.t

        val calc_log_evidence_sigma2 : t -> float

        val prepare_hyper : t -> hyper_t
        val calc_log_evidence : hyper_t -> Spec.Hyper.t -> float
      end

      module Trained : sig
        type t
        type hyper_t

        val calc : Model.t -> targets : vec -> t
        val calc_eval : t -> Eval.Trained.t

        val calc_log_evidence_sigma2 : t -> float

        val prepare_hyper : t -> hyper_t
        val calc_log_evidence : hyper_t -> Spec.Hyper.t -> float
      end

      module Test : sig
        val check_deriv_hyper :
          Eval.Spec.Kernel.t ->
          Eval.Spec.Inducing.t ->
          Eval.Spec.Inputs.t ->
          Spec.Hyper.t ->
          eps : float ->
          tol : float ->
          unit
      end
    end
  end

  module type SPGP = sig
    module Eval : Eval
    module Deriv : Deriv with module Eval = Eval

    module SPGP : sig
      module Spec : Specs.SPGP
        with module Eval = Deriv.Eval.Spec
        with module Deriv = Deriv.Deriv.Spec

      module Gsl : sig
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
    end
  end
end
