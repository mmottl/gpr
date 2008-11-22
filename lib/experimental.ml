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
