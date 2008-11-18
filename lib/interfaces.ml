open Lacaml.Impl.D

(* TODO: consider loghyper variables *)
(* TODO: store Kmn *)

module type Kernel = sig
  type t
  type input
  type inputs

  val get_n_inputs : inputs -> int

  val eval_one : t -> input -> float
  val eval : t -> input -> input -> float
  val evals : t -> input -> inputs -> vec

  val weighted_eval : t -> weights : vec -> input -> inputs -> float
  val weighted_evals : t -> weights : vec -> inputs -> inputs -> vec

  val upper : t -> inputs -> mat
  val upper_no_diag : t -> inputs -> mat
  val cross : t -> inputs -> inputs -> mat
  val diag_mat : t -> inputs -> dst : mat -> unit
  val diag_vec : t -> inputs -> vec
end

module type Vec_kernel =
  Kernel
    with type input = vec
    with type inputs = mat

module type Inducing_input_gpr = sig
  module Spec : sig
    type kernel
    type input
    type inputs

    type t = {
      kernel : kernel;
      sigma2 : float;
      inducing_inputs : inputs;
    }
  end

  open Spec

  module Trained : sig
    type t

    val train : Spec.t -> inputs : inputs -> targets : vec -> t

    val neg_log_likelihood : t -> float
  end

  module Mean_predictor : sig
    type t

    val of_trained : Trained.t -> t

    val mean : t -> input -> float
    val means : t -> inputs -> vec
    val inducing_means : t -> Trained.t -> vec

    val spec : t -> Spec.t
    val weights : t -> vec
  end

  module Full_predictor : sig
    type t

    val of_trained : Trained.t -> t
    val mean_predictor : t -> Mean_predictor.t
    val mean_variance : ?with_noise : bool -> t -> input -> float * float
    val means_variances : ?with_noise : bool -> t -> inputs -> vec * vec
    val means_covariances : ?with_noise : bool -> t -> inputs -> vec * mat
  end
end
