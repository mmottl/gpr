module type Kernel = sig
  type t
  type input
  type inputs

  val eval : input -> input -> float
  val evals : input -> inputs -> dst : vec -> unit

  val weighted_eval : coeffs : vec -> input -> inputs -> float
  val weighted_evals : coeffs : vec -> inputs -> inputs -> dst : vec -> unit

  val upper : inputs -> dst : mat -> unit
  val cross : inputs -> inputs -> dst : mat -> unit
  val diag_mat : inputs -> dst : mat -> unit
  val diag_vec : inputs -> dst : vec -> unit
end

module type Inducing_input_gp : sig
  type kernel
  type input
  type inputs

  module Spec : sig
    type t = private {
      kernel : Kernel.t;
      sigma2 : float;
      inducing_inputs : inputs;
    }

    val make :
      ?jitter : float ->
      kernel ->
      sigma2 : float ->
      inducing_inputs : inputs
      -> t
  end

  module Trained : sig
    type t

    val train : Spec.t -> inputs : inputs -> targets : vec -> t

    val neg_log_likelihood : t -> float

    val mean_squared_error : t -> float

    val inducing_means : t -> vec
    val inducing_means_variances : t -> vec * vec
    val inducing_means_covariances : t -> vec * mat
    val inducing_means_covariances_packed : t -> vec * vec
    (* TODO: etc.; can be done efficiently *)
  end

  module Mean_predictor : sig
    type t

    val of_trained : Trained.t -> t

    val mean : t -> input -> float
    val means : t -> inputs -> means : vec -> unit

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
