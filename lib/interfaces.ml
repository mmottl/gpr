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
  module Kernel : Kernel

  open Kernel

  module Inducing : sig
    type t

    val calc : Kernel.t -> inputs -> t
  end

  module Induced : sig
    type t

    val calc : Inducing.t -> input -> t
  end

  module Induceds : sig
    type t

    val calc : Inducing.t -> inputs -> t
  end

  module Model : sig
    type t

    val calc : Induceds.t -> t
    val neg_log_likelihood : t -> float
  end

  module Trained : sig
    type t

    val calc : Model.t -> targets : vec -> t
    val neg_log_likelihood : t -> float
  end

  module Weights : sig
    type t

    val calc : Trained.t -> t
    val copy : t -> vec
  end

  module Mean : sig
    type t

    val calc : Weights.t -> input -> t
    val calc_induced : Weights.t -> Induced.t -> t
    val copy : t -> float
  end

  module Means : sig
    type t

    val calc_inducing : Weights.t -> Model.t -> t
    val calc_inputs : Weights.t -> Model.t -> t
    val calc : Weights.t -> inputs -> t
    val calc_induceds : Weights.t -> Induceds.t -> t
    val copy : t -> vec
  end

  module Variance : sig
    type t

    val calc_induced : Model.t -> Induced.t -> t
    val copy : ?predictive : bool -> t -> float
  end

  module Variances : sig
    type t

    val calc_inducing : Model.t -> t
    val calc_inputs : Model.t -> t
    val calc_induceds : Model.t -> Induceds.t -> t
    val copy : ?predictive : bool -> t -> vec
  end

  module Covariances : sig
    type t

    val calc_inducing : Model.t -> t
    val calc_inputs : Model.t -> t
    val calc_induceds : Model.t -> Induceds.t -> t
    val variances : t -> Variances.t
    val copy : ?predictive : bool -> t -> mat
  end

  module Sampler : sig
    type t

    val calc : ?predictive : bool -> Mean.t -> Variance.t -> t
    val sample : ?rng : Gsl_rng.t -> t -> float
  end

  module Samplers : sig
    type t

    val calc : ?predictive : bool -> Means.t -> Covariances.t -> t
    val sample : ?rng : Gsl_rng.t -> t -> vec
  end
end
