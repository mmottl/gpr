open Lacaml.Impl.D

open Interfaces
open Utils

module type Inducing_input_gpr = sig
  module Kernel : Kernel

  open Kernel

  module Inducing : sig
    type t

    val calc : Kernel.t -> inputs -> t
  end

  module Prepared : sig
    type t

    val calc : Inducing.t -> input -> t
  end

  module Prepareds : sig
    type t

    val calc : Inducing.t -> inputs -> t
  end

  module Model : sig
    type t

    val calc : Prepareds.t -> t
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
    val calc_prepared : Weights.t -> Prepared.t -> t
    val copy : t -> float
  end

  module Means : sig
    type t

    val calc_inducing : Weights.t -> Model.t -> t
    val calc_inputs : Weights.t -> Model.t -> t
    val calc : Weights.t -> inputs -> t
    val calc_prepareds : Weights.t -> Prepareds.t -> t
    val copy : t -> vec
  end

  module Variance : sig
    type t

    val calc_prepared : Model.t -> Prepared.t -> t
    val copy : ?predictive : bool -> t -> float
  end

  module Sampler : sig
    type t

    val calc : ?predictive : bool -> Mean.t -> Variance.t -> t
    val sample : ?rng : Gsl_rng.t -> t -> float
  end

  module Variances : sig
    type t

    val calc_inducing : Model.t -> t
    val calc_inputs : Model.t -> t
    val calc_prepareds : Model.t -> Prepareds.t -> t
    val copy : ?predictive : bool -> t -> vec
  end

  module Covariances : sig
    type t

    val calc_inducing : Model.t -> t
    val calc_inputs : Model.t -> t
    val calc_prepareds : Model.t -> Prepareds.t -> t
    val variances : t -> Variances.t
    val copy : ?predictive : bool -> t -> mat
  end

  module Samplers : sig
    type t

    val calc : ?predictive : bool -> Means.t -> Covariances.t -> t
    val sample : ?rng : Gsl_rng.t -> t -> vec
  end
end

(* TODO: dimension sanity checks *)

module type FITC_spec = sig
  module Kernel : Kernel

  val get_sigma2 : Kernel.t -> float
  val jitter : float
end

module Make_FITC_common (Loc : sig val loc : string end) (Spec : FITC_spec) =
struct
  include Spec

  open Kernel

  module Inducing = struct
    type t =
      {
        kernel : Kernel.t;
        inputs : inputs;
        km : mat;
        km_chol : mat;
        log_det_km : float;
      }

    let calc kernel inputs =
      let km = Kernel.upper kernel inputs in
      (* TODO: copy upper triangle only *)
      let km_chol = Mat.copy km in
      potrf ~jitter km_chol;
      let log_det_km = log_det km_chol in
      {
        kernel = kernel;
        inputs = inputs;
        km = km;
        km_chol = km_chol;
        log_det_km = log_det_km;
      }
  end

  module Prepared = struct
    type t =
      {
        inducing : Inducing.t;
        input : input;
        k_m : vec;
      }

    let calc inducing input =
      let kernel = inducing.Inducing.kernel in
      {
        inducing = inducing;
        input = input;
        k_m = Kernel.evals kernel input inducing.Inducing.inputs;
      }

    let get_kernel t = t.inducing.Inducing.kernel
  end

  module Prepareds = struct
    type t =
      {
        inducing : Inducing.t;
        inputs : inputs;
        kmn : mat;
      }

    let calc inducing inputs =
      let kernel = inducing.Inducing.kernel in
      let kmn = Kernel.cross kernel inducing.Inducing.inputs inputs in
      {
        inducing = inducing;
        inputs = inputs;
        kmn = kmn;
      }

    let get_kernel t = t.inducing.Inducing.kernel
  end

  let calc_basis_2_k basis_chol ~k =
    (* TODO: consider symmetric matrices *)
    let b_2_k = Mat.copy k in
    trtrs ~trans:`T basis_chol b_2_k;
    b_2_k

  (* Compute square root of Nystrom approximation, and diagonal of
     marginal variances *)
  let nystrom_2_marginals prepareds =
    let inducing = prepareds.Prepareds.inducing in
    let kernel = Prepareds.get_kernel prepareds in
    let kn_diag = Kernel.diag_vec kernel prepareds.Prepareds.inputs in
    let kmn = prepareds.Prepareds.kmn in
    let km_2_kmn = calc_basis_2_k inducing.Inducing.km_chol ~k:kmn in
    kn_diag, km_2_kmn, kmn

  module Model = struct
    type t =
      {
        prepareds : Prepareds.t;
        kn_diag : vec;
        km_2_kmn : mat;
        lam_diag : vec;
        inv_lam_sigma2_diag : vec;
        b_chol : mat;
        neg_log_likelihood : float;
      }

    let calc prepareds =
      let inducing = prepareds.Prepareds.inducing in
      let kernel = Prepareds.get_kernel prepareds in
      let sigma2 = get_sigma2 kernel in
      let kn_diag, km_2_kmn, kmn = nystrom_2_marginals prepareds in
      let kmn_= Mat.copy kmn in
      let n_inputs = Vec.dim kn_diag in
      let lam_diag = Vec.create n_inputs in
      let inv_lam_sigma2_diag = Vec.create n_inputs in
      let rec loop log_det_lam_sigma2 i =
        if i = 0 then log_det_lam_sigma2
        else
          let kn_diag_i = kn_diag.{i} in
          (* TODO: optimize ssqr and col *)
          let qn_diag_i = Vec.ssqr (Mat.col km_2_kmn i) in
          let lam_i = kn_diag_i -. qn_diag_i in
          lam_diag.{i} <- lam_i;
          let lam_sigma2_i = lam_i +. sigma2 in
          let inv_lam_sigma2_i = 1. /. lam_sigma2_i in
          inv_lam_sigma2_diag.{i} <- inv_lam_sigma2_i;
          (* TODO: optimize scal col *)
          scal (sqrt inv_lam_sigma2_i) (Mat.col kmn_ i);
          loop (log_det_lam_sigma2 +. log lam_sigma2_i) (i - 1)
      in
      let log_det_lam_sigma2 = loop 0. n_inputs in
      (* TODO: copy upper triangle only *)
      let b_chol = syrk kmn_ ~beta:1. ~c:(Mat.copy inducing.Inducing.km) in
      potrf ~jitter b_chol;
      let log_det_km = inducing.Inducing.log_det_km in
      let log_det_b = log_det b_chol in
      let l1_2 = log_det_b -. log_det_km +. log_det_lam_sigma2 in
      let neg_log_likelihood = (l1_2 +. float n_inputs *. log_2pi) /. 2. in
      {
        prepareds = prepareds;
        kn_diag = kn_diag;
        km_2_kmn = km_2_kmn;
        lam_diag = lam_diag;
        inv_lam_sigma2_diag = inv_lam_sigma2_diag;
        b_chol = b_chol;
        neg_log_likelihood = neg_log_likelihood;
      }

    let neg_log_likelihood model = model.neg_log_likelihood

    let get_inducing model = model.prepareds.Prepareds.inducing
    let get_inducing_inputs model = (get_inducing model).Inducing.inputs
    let get_inputs model = model.prepareds.Prepareds.inputs
    let get_kernel model = (get_inducing model).Inducing.kernel
    let get_sigma2 model = get_sigma2 (get_kernel model)
    let get_kmn model = model.prepareds.Prepareds.kmn
    let get_lam_diag model = model.lam_diag
    let get_km model = (get_inducing model).Inducing.km
  end

  module Trained = struct
    type t =
      {
        kernel : Kernel.t;
        inducing_inputs : inputs;
        b_chol : mat;
        inv_b_chol_kmn_y__ : vec;
        neg_log_likelihood : float;
      }

    let calc model ~targets =
      let n_targets = Vec.dim targets in
      let y__ = Vec.create n_targets in
      let inv_lam_sigma2_diag = model.Model.inv_lam_sigma2_diag in
      for i = 1 to n_targets do
        y__.{i} <- targets.{i} *. inv_lam_sigma2_diag.{i}
      done;
      let inv_b_chol_kmn_y__ = gemv (Model.get_kmn model) y__ in
      let ssqr_y__ = Vec.ssqr y__ in
      let b_chol = model.Model.b_chol in
      trsv ~trans:`T b_chol inv_b_chol_kmn_y__;
      let fit_neg_log_likelihood =
        (ssqr_y__ -. Vec.ssqr inv_b_chol_kmn_y__) /. 2.
      in
      {
        kernel = Model.get_kernel model;
        inducing_inputs = (Model.get_inducing model).Inducing.inputs;
        b_chol = b_chol;
        inv_b_chol_kmn_y__ = inv_b_chol_kmn_y__;
        neg_log_likelihood =
          model.Model.neg_log_likelihood +. fit_neg_log_likelihood;
      }

    let neg_log_likelihood trained = trained.neg_log_likelihood
  end

  module Weights = struct
    type t =
      {
        kernel : Kernel.t;
        inducing_inputs : inputs;
        weights : vec;
      }

    let calc trained =
      let weights = copy trained.Trained.inv_b_chol_kmn_y__ in
      trsv trained.Trained.b_chol weights;
      {
        kernel = trained.Trained.kernel;
        inducing_inputs = trained.Trained.inducing_inputs;
        weights = weights;
      }

    let copy t = copy t.weights
  end

  module Mean = struct
    type t = { input : input; mean : float }

    let make ~input ~mean = { input = input; mean = mean }


    let calc w input =
      let module W = Weights in
      let mean =
        Kernel.weighted_eval
          w.W.kernel ~weights:w.W.weights input w.W.inducing_inputs
      in
      make ~input ~mean

    let calc_prepared w { Prepared.input = input; k_m = k_m } =
      make ~input ~mean:(dot ~x:k_m w.Weights.weights)

    let copy m = m.mean
  end

  module Means = struct
    type t = { inputs : inputs; means : vec }

    let make ~inputs ~means = { inputs = inputs; means = means }

    let calc_inducing { Weights.weights = weights } model =
      make
        ~inputs:(Model.get_inputs model)
        ~means:(symv (Model.get_km model) weights)

    let calc_inputs { Weights.weights = weights } model =
      make
        ~inputs:(Model.get_inputs model)
        ~means:(gemv ~trans:`T (Model.get_kmn model) weights)

    let calc w inputs =
      let module W = Weights in
      let means =
        Kernel.weighted_evals
          w.W.kernel ~weights:w.W.weights inputs w.W.inducing_inputs
      in
      make ~inputs ~means

    let calc_prepareds w { Prepareds.inputs = inputs; kmn = kmn } =
      make ~inputs ~means:(gemv ~trans:`T kmn w.Weights.weights)

    let copy m = copy m.means
  end

  module Variance = struct
    type t = { input : input; variance : float; sigma2 : float }

    let calc_prepared model prepared =
      let { Prepared.input = input; k_m = k_m } = prepared in
      let kernel = Model.get_kernel model in
      let prior_variance = Kernel.eval_one kernel input in
      let inv_km_chol_k_m = copy k_m in
      let inv_km_chol_k_m_mat = Mat.from_col_vec inv_km_chol_k_m in
      let inducing = prepared.Prepared.inducing in
      potrs ~factorize:false inducing.Inducing.km_chol inv_km_chol_k_m_mat;
      let km_arg = dot ~x:k_m inv_km_chol_k_m in
      let inv_b_chol_k_m = copy k_m ~y:inv_km_chol_k_m in
      let inv_b_chol_k_m_mat = inv_km_chol_k_m_mat in
      potrs ~factorize:false model.Model.b_chol inv_b_chol_k_m_mat;
      let b_arg = dot ~x:k_m inv_b_chol_k_m in
      let explained_variance = km_arg -. b_arg in
      let variance = prior_variance -. explained_variance in
      { input = input; variance = variance; sigma2 = Model.get_sigma2 model }

    let copy ?predictive t =
      match predictive with
      | None | Some true -> t.variance +. t.sigma2
      | Some false -> t.variance
  end

  let calc_b_2_k model ~k = calc_basis_2_k model.Model.b_chol ~k

  module Variances = struct
    type t = { inputs : inputs; variances : vec; sigma2 : float }

    let make ~inputs ~variances ~model =
      let sigma2 = Model.get_sigma2 model in
      { inputs = inputs; variances = variances; sigma2 = sigma2 }

    let calc_inducing model =
      let b_2_km = calc_b_2_k model ~k:(Model.get_km model) in
      let m = Mat.dim2 b_2_km in
      let variances = Vec.create m in
      for i = 1 to m do
        (* TODO: optimize ssqr and col *)
        variances.{i} <- Vec.ssqr (Mat.col b_2_km i)
      done;
      make ~inputs:(Model.get_inducing_inputs model) ~variances ~model

    let calc_inputs model =
      let variances = copy (Model.get_lam_diag model) in
      let b_2_kmn = calc_b_2_k model ~k:(Model.get_kmn model) in
      let n = Mat.dim2 b_2_kmn in
      for i = 1 to n do
        (* TODO: optimize ssqr and col *)
        variances.{i} <- variances.{i} +. Vec.ssqr (Mat.col b_2_kmn i)
      done;
      make ~inputs:(Model.get_inputs model) ~variances ~model

    let calc_prepareds model prepareds =
      let kt_diag, km_2_kmt, kmt = nystrom_2_marginals prepareds in
      let variances = copy kt_diag in
      let b_2_kmt = calc_b_2_k model ~k:kmt in
      let n = Mat.dim2 b_2_kmt in
      for i = 1 to n do
        let explained_variance =
          (* TODO: optimize ssqr and col *)
          Vec.ssqr (Mat.col km_2_kmt i) -. Vec.ssqr (Mat.col b_2_kmt i)
        in
        variances.{i} <- variances.{i} -. explained_variance
      done;
      make ~inputs:prepareds.Prepareds.inputs ~variances ~model

    let copy ?predictive { variances = variances; sigma2 = sigma2 } =
      match predictive with
      | None | Some true ->
          let predictive_variances = Vec.make (Vec.dim variances) sigma2 in
          axpy ~x:variances predictive_variances;
          predictive_variances
      | Some false -> copy variances
  end

  module Sampler = struct
    type t = { mean : float; stddev : float }

    let calc ?predictive mean variance =
      if mean.Mean.input != variance.Variance.input then
        failwith (Loc.loc ^ ".Sampler: mean and variance disagree about input");
      let used_variance =
        match predictive with
        | None | Some true ->
            variance.Variance.variance +. variance.Variance.sigma2
        | Some false -> variance.Variance.variance
      in
      { mean = mean.Mean.mean; stddev = sqrt used_variance }

    let sample ?(rng = default_rng) sampler =
      let noise = Gsl_randist.gaussian rng ~sigma:sampler.stddev in
      sampler.mean +. noise
  end
end

module Make_FITC (Spec : FITC_spec) : Inducing_input_gpr = struct
  module Loc = struct let loc = "FITC" end

  include Make_FITC_common (Loc) (Spec)

  open Kernel

  module Covariances = struct
    type t = { inputs : inputs; covariances : mat; sigma2 : float }

    let make ~inputs ~covariances ~model =
      let sigma2 = Model.get_sigma2 model in
      { inputs = inputs; covariances = covariances; sigma2 = sigma2 }

    let calc_inducing model =
      let b_2_km = calc_b_2_k model ~k:(Model.get_km model) in
      let inputs = Model.get_inducing_inputs model in
      make ~inputs ~covariances:(syrk b_2_km) ~model

    let calc_common ~kn_diag ~kmn ~km_2_kmn ~inputs ~model =
      let kernel = Model.get_kernel model in
      let covariances = Kernel.upper_no_diag kernel inputs in
      for i = 1 to Vec.dim kn_diag do
        covariances.{i, i} <- kn_diag.{i}
      done;
      ignore (syrk ~alpha:(-1.) km_2_kmn ~beta:1. ~c:covariances);
      let b_2_kmn = calc_b_2_k model ~k:kmn in
      ignore (syrk ~alpha:(1.) b_2_kmn ~beta:1. ~c:covariances);
      make ~inputs ~covariances ~model

    let calc_inputs model =
      let kn_diag = model.Model.kn_diag in
      let kmn = model.Model.prepareds.Prepareds.kmn in
      let km_2_kmn = model.Model.km_2_kmn in
      let inputs = Model.get_inputs model in
      calc_common ~kn_diag ~kmn ~km_2_kmn ~inputs ~model

    let calc_prepareds model prepareds =
      let inputs = prepareds.Prepareds.inputs in
      let kn_diag, km_2_kmn, kmn = nystrom_2_marginals prepareds in
      calc_common ~kn_diag ~kmn ~km_2_kmn ~inputs ~model

    let variances { inputs = inputs; covariances = covs; sigma2 = sigma2 } =
      { Variances.inputs = inputs; variances = Mat.diag covs; sigma2 = sigma2 }

    let copy ?predictive { covariances = covariances; sigma2 = sigma2 } =
      (* TODO: copy upper triangle only *)
      let res = Mat.copy covariances in
      match predictive with
      | None | Some true ->
          for i = 1 to Mat.dim1 res do res.{i, i} <- res.{i, i} +. sigma2 done;
          res
      | Some false -> res
  end

  module Samplers = struct
    type t = { means : vec; cov_chol : mat }

    let calc ?predictive means covariances =
      if means.Means.inputs != covariances.Covariances.inputs then
        failwith
          (Loc.loc ^ ".Samplers: means and covariances disagree about inputs");
      (* TODO: copy upper triangle only *)
      let cov_chol = Mat.copy covariances.Covariances.covariances in
      begin
        match predictive with
        | None | Some true ->
            let sigma2 = covariances.Covariances.sigma2 in
            for i = 1 to Mat.dim1 cov_chol do
              cov_chol.{i, i} <- cov_chol.{i, i} +. sigma2
            done
        | Some false -> ()
      end;
      potrf ~jitter cov_chol;
      { means = means.Means.means; cov_chol = cov_chol }

    let sample ?(rng = default_rng) samplers =
      let n = Vec.dim samplers.means in
      let noises = Vec.init n (fun _ -> Gsl_randist.gaussian rng ~sigma:1.) in
      trmv ~trans:`T samplers.cov_chol noises;
      axpy ~x:samplers.means noises;
      noises
  end
end

(*
      let kn = Kernel.upper
      let b_2_kmn = calc_b_2_k model ~k:(Model.get_kmn model) in
      let covariances = syrk b_2_kmn in
      let n = Vec.dim lam_diag in
      for i = 1 to n do
        covariances.{i, i} <- lam_diag.{i} +. covariances.{i, i}
      done;
      covariances, Model.get_sigma2 model
*)

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
(*
let trace_prod =
  (assert false (* XXX *))

module Squared_exponential_ard : Kernel = struct
  type hyper =
    {
      alpha2 : float;
      lambdas : vec;
    }

  let calc_self_ard_r ~lambdas ~x ~col =
    let m = Mat.dim1 x in
    let rec loop ~res ~der =
      for row = 1 to m do
        res := !res
      done
    in
    loop ~res:0. ~der:0.

  let calc_self hyper ~x ~col =
    let m = Mat.dim1 x in
    for row = 1 to m do
    done

  let calc_diag_mat ~x ~dst =
    let n = Mat.dim1 x in
    for j = 1 to n do
      dst.{j, j} <- calc_self ~x ~col
    done

  let calc_upper_gram_mat ~x ~dst =
    calc_diag_mat ~x ~dst;
    let m = Mat.dim1 x in
    let n = Mat.dim2 x in
    for col1 = 1 to n do
      for col2 = col1 + 1 to n do
        calc_slice ~x ~col1 ~col2 ~dst
      done
    done
end

type ('data_spec, 'data) kernel =
  {
    get_num_hyps : 'data_spec -> int;
    eval : 'data ->
  }

let add_kernels k1 k2 =
  {
    get_num_hyps = fun data_spec ->
      k1.get_num_hyps data_spec + k2.get_num_hyps data_spec;
    eval = fun hyper ->
  }


module type Kernel = sig
  type t
end

module Squared_Hyp = struct
  type t
end

type k =
  | SqConst
  | EuclidianProd
  | Prod of k * k
  | Sum of k * k
  | Exp of any

type any =
  | Const of float

type kernel =
  | `SqConst of float
  | `Sum of kernel * kernel
  | `Prod of kernel * kernel
  | `InnerProd
  | `WeightedInnerProd

type env =
  | `SqConst of float
  | `Sum of env * env
  | `Prod of env * env
  | `InnerProd
  | `WeightedInnerProd

module type KERNEL = sig
  type hyper
  type data_set
  type data_point
  type kernel

  val prod : kernel -> kernel -> kernel
  val sum : kernel -> kernel -> kernel

  val create_hyper :
end



val prod :
  ('hyper1, 'data1) kernel -> ('hyper2, 'data2) kernel
  -> ('hyper1 * 'hyper2, 'data1 * 'data2) kernel

val sum :
  ('hyper1, 'data1) kernel -> ('hyper2, 'data2) kernel
  -> ('hyper1 * 'hyper2, 'data1 * 'data2) kernel

val sq_const : float -> (unit, 'a) kernel

val sq_var : unit -> (float, 'a) kernel

val inner_prod : 'ip_space -> (unit, 'ip_space) kernel

val weighted_inner_prod : unit kernel

val create :
  ?upper : (hyper : 'hyper -> data : 'data -> dst : mat)
  ?single : (hyper : 'hyper -> data : 'data -> dst : mat)


let _ =
  prod (sq_const 1.3)
    (inner_prod


(* Autom. generate OCaml-code from kernel spec *)
*)
