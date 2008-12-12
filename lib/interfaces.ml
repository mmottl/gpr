open Lacaml.Impl.D

(* TODO: consider loghyper variables *)

module Inducing_input_gpr = struct
  module Specs = struct
    module type Eval = sig
      module Kernel : sig
        type t

        val get_sigma2 : t -> float
      end

      module Inducing : sig
        type t

        val upper : Kernel.t -> t -> mat
      end

      module Input : sig
        type t

        val eval_one : Kernel.t -> t -> float

        val eval : Kernel.t -> inducing : Inducing.t -> input : t -> vec

        val weighted_eval :
          Kernel.t -> coeffs : vec -> inducing : Inducing.t -> input : t -> float
      end

      module Inputs : sig
        type t

        val weighted_eval :
          Kernel.t -> coeffs : vec -> inducing : Inducing.t -> inputs : t -> vec

        val upper : Kernel.t -> t -> mat
        val upper_no_diag : Kernel.t -> t -> mat
        val diag : Kernel.t -> t -> vec
        val cross : Kernel.t -> inducing : Inducing.t -> inputs : t -> mat
      end
    end

    module type Deriv = sig
      module Kernel : sig
        type t
      end

      module Hyper : sig
        type t
      end

      module Inducing : sig
        type t
        type shared

        val calc_shared : Kernel.t -> t -> mat * shared
        val calc_deriv : shared -> Hyper.t -> mat -> int array -> int
      end

      module Inputs : sig
        type t
        type upper
        type diag
        type cross

        val calc_shared_upper : Kernel.t -> t -> mat * upper

        val calc_shared_diag : Kernel.t -> t -> vec * diag

        val calc_shared_cross :
          Kernel.t -> inducing : Inducing.t -> inputs : t -> mat * cross

        val calc_deriv_upper : upper -> Hyper.t -> mat option
        val calc_deriv_diag : diag -> Hyper.t -> vec option
        val calc_deriv_cross : cross -> Hyper.t -> mat * int array * int
      end
    end
  end

  module Sigs = struct
    module type Eval = sig
      module Spec : Specs.Eval

      open Spec

      module Inducing : sig
        type t

        val calc : Kernel.t -> Spec.Inducing.t -> t
      end

      module Input : sig
        type t

        val calc : Inducing.t -> Spec.Input.t -> t
      end

      module Inputs : sig
        type t

        val calc : Inducing.t -> Spec.Inputs.t -> t
      end

      module Model : sig
        type t

        val calc : Inputs.t -> t
        val calc_evidence : t -> float
      end

      module Trained : sig
        type t

        val calc : Model.t -> targets : vec -> t
        val calc_evidence : t -> float
      end

      module Weights : sig
        type t

        val calc : Trained.t -> t

        val get_kernel : t -> Kernel.t
        val get_inducing_points : t -> Spec.Inducing.t
        val get_coeffs : t -> vec
      end

      module Mean : sig
        type t

        val calc_input : Weights.t -> Spec.Input.t -> t
        val calc_induced : Weights.t -> Input.t -> t

        val get : t -> float
      end

      module Means : sig
        type t

        val calc_model_inputs : Weights.t -> t
        val calc_inputs : Weights.t -> Spec.Inputs.t -> t
        val calc_induced : Weights.t -> Inputs.t -> t

        val get : t -> vec

        module Inducing : sig
          type t

          val calc : Weights.t -> t
          val get : t -> vec
        end
      end

      module Variance : sig
        type t

        (* NOTE: there is deliberately no [calc_input] function. *)
        val calc_induced : Model.t -> Input.t -> t
        val get : ?predictive : bool -> t -> float
      end

      module Variances : sig
        type t

        val calc_model_inputs : Model.t -> t
        (* NOTE: there is deliberately no [calc_inputs] function. *)
        val calc_induced : Model.t -> Inputs.t -> t
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
        val calc_induced : Model.t -> Inputs.t -> t
        val get : ?predictive : bool -> t -> mat
        val variances : t -> Variances.t

        module Inducing : sig
          type t

          val calc : Model.t -> t
          val get : ?predictive : bool -> t -> mat
          val variances : t -> Variances.Inducing.t
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
        module Spec :
          Specs.Deriv
            with type Kernel.t = Eval.Spec.Kernel.t
            with type Inputs.t = Eval.Spec.Inputs.t

        open Spec

        module Inducing : sig
          type t

          val calc : Kernel.t -> Inputs.t -> t
          val calc_eval : t -> Eval.Inducing.t
        end

        module Inputs : sig
          type t

          val calc : Inducing.t -> Inputs.t -> t
          val calc_eval : t -> Eval.Inputs.t
        end

        module Model : sig
          type t

          val calc : Inputs.t -> t
          val calc_eval : t -> Eval.Model.t
          val calc_evidence : t -> Hyper.t -> float
          val calc_evidence_sigma2 : t -> float
        end

        module Trained : sig
          type t

          val calc : Model.t -> targets : vec -> t
          val calc_eval : t -> Eval.Trained.t
          val calc_evidence : t -> Hyper.t -> float
          val calc_evidence_sigma2 : t -> float
        end
      end
    end
  end
end
