open Lacaml.Impl.D

module type Kernel = sig
  type t
  type input
  type inputs

  val get_n_inputs : inputs -> int

  val eval_one : t -> input -> float
  val eval : t -> input -> input -> float
  val evals : t -> input -> inputs -> dst : vec -> unit

  val weighted_eval : t -> coeffs : vec -> input -> inputs -> float

  val weighted_evals :
    t -> coeffs : vec -> inputs -> inputs -> dst : vec -> unit

  val upper : t -> inputs -> dst : mat -> unit
  val cross : t -> inputs -> inputs -> dst : mat -> unit
  val diag_mat : t -> inputs -> dst : mat -> unit
  val diag_vec : t -> inputs -> dst : vec -> unit
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
    val means : t -> inputs -> means : vec -> unit
    val inducing_means : t -> Trained.t -> vec

    val spec : t -> Spec.t
    val coeffs : t -> vec
  end

  module Full_predictor : sig
    type t

    val of_trained : Trained.t -> t

    val mean_predictor : t -> Mean_predictor.t

    val mean_variance : ?with_noise : bool -> t -> input -> float * float

    val means_variances :
      ?with_noise : bool ->
      t ->
      inputs ->
      means : vec ->
      variances : vec
      -> unit

    val means_covariances :
      ?with_noise : bool ->
      t ->
      inputs ->
      means : vec ->
      covariances : mat
      -> unit

    val means_covariances_packed :
      ?with_noise : bool ->
      t ->
      inputs ->
      means : vec ->
      covariances : vec
      -> unit

    val inducing_means_variances :
      ?with_noise : bool -> t -> Trained.t -> vec * vec

    val inducing_means_covariances :
      ?with_noise : bool -> t -> Trained.t -> vec * mat
  end
end

(*
module type Deriv_kernel_fun = sig
  include Kernel_fun

  val n_hyper : int

  val upper_deriv_precomp : hyper : int -> x : inputs -> dst : mat -> unit
  val cross_deriv_precomp : hyper : int -> x : inputs -> y : inputs -> dst : mat -> unit
  val deriv_diag_mat_precomp : hyper : int -> x : inputs -> dst : mat -> unit
  val deriv_diag_vec_precomp : hyper : int -> x : inputs -> dst : vec -> unit

  val upper_deriv : hyper : int -> x : inputs -> dst : mat -> unit
  val cross_deriv : hyper : int -> x : inputs -> y : inputs -> dst : mat -> unit
  val deriv_diag_mat : hyper : int -> x : inputs -> dst : mat -> unit
  val deriv_diag_vec : hyper : int -> x : inputs -> dst : vec -> unit
end

module type Pseudo_input_kernel = sig
  include Kernel_fun

  val get_pseudo_init : n : int -> inputs -> inputs

  val upper_deriv_precomp : hyper : int -> x : inputs -> dst : mat -> unit
  val cross_deriv_precomp : hyper : int -> x : inputs -> y : inputs -> dst : mat -> unit
  val deriv_diag_mat_precomp : hyper : int -> x : inputs -> dst : mat -> unit
  val deriv_diag_vec_precomp : hyper : int -> x : inputs -> dst : vec -> unit

  val upper_deriv : hyper : int -> x : inputs -> dst : mat -> unit
  val cross_deriv : hyper : int -> x : inputs -> y : inputs -> dst : mat -> unit
  val deriv_diag_mat : hyper : int -> x : inputs -> dst : mat -> unit
  val deriv_diag_vec : hyper : int -> x : inputs -> dst : vec -> unit
end
*)
