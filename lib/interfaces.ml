open Lacaml.Impl.D

(* TODO: consider loghyper variables *)

module Inducing_input_gpr = struct
  module Specs = struct
    module type Eval = sig
      type kernel

      module Inducing : sig
        type t

        val size : t -> int
        val upper : kernel -> t -> mat
      end

      module Input : sig
        type t

        val eval_one : kernel -> t -> float

        val eval : kernel -> inducing : Inducing.t -> input : t -> vec

        val weighted_eval :
          kernel -> coeffs : vec -> inducing : Inducing.t -> input : t -> float
      end

      module Inputs : sig
        type t

        val size : t -> int

        val weighted_eval :
          kernel -> coeffs : vec -> inducing : Inducing.t -> inputs : t -> vec

        val upper : kernel -> t -> mat
        val upper_no_diag : kernel -> t -> mat
        val diag : kernel -> t -> vec
        val cross : kernel -> inducing : Inducing.t -> inputs : t -> mat
      end
    end

    module type Deriv = sig
      module Kernel : sig
        type t

        val get_n_hypers : t -> int
      end

      module Inducing : sig
        type t
        type upper

        val shared_upper : Kernel.t -> t -> upper
      end

      module Inputs : sig
        type t
        type upper
        type diag
        type cross

        val get_n_hypers : t -> int array

        val shared_upper : Kernel.t -> t -> upper

        val shared_cross :
          Kernel.t -> inducing : Inducing.t -> inputs : t -> cross

        val shared_diag : Kernel.t -> t -> diag

        val eval_upper : upper -> mat
        val eval_cross : upper -> mat
        val eval_diag : upper -> vec

        val deriv_upper : upper -> int -> mat
        val deriv_cross : cross -> int -> mat
        val deriv_diag : diag -> int -> vec
      end
    end
  end

  module Sigs = struct
    module type Eval = sig
      module Spec : Specs.Eval
      open Spec

      module Inducing : sig
        type t

        val calc : kernel -> Spec.Inducing.t -> t
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
        val neg_log_marginal_likelihood : t -> float
      end

      module Trained : sig
        type t

        val calc : Model.t -> targets : vec -> t
        val neg_log_marginal_likelihood : t -> float
      end

      module Weights : sig
        type t

        val calc : Trained.t -> t

        val get_kernel : t -> kernel
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
  end

(*
    module type Deriv = sig
      module Eval : Eval

      module Deriv : sig
        module Kernel :
          Deriv_kernel
            with type t = Eval.Kernel.t
            with type inputs = Eval.Kernel.inputs

        module Inducing : sig
          type t

          val calc : Kernel.t -> Kernel.inputs -> t
          val to_eval : t -> Eval.Inducing.t
        end

        module Induceds : sig
          type t

          val calc : Inducing.t -> Kernel.inputs -> t
          val to_eval : t -> Eval.Induceds.t
        end

        module Model : sig
          type t

          val calc : Induceds.t -> t
          val neg_log_marginal_likelihood : t -> vec
          val to_eval : t -> Eval.Model.t
        end

        module Trained : sig
          type t

          val calc : Model.t -> targets : vec -> t
          val neg_log_marginal_likelihood : t -> vec
          val to_eval : t -> Eval.Trained.t
        end
      end
    end
*)
end
