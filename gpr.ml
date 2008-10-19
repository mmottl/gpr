(*
open Lacaml.Impl.D

module type Kernel_spec = sig
  type data

  val get_n_samples : data -> int
end

module type Kernel_fun = sig
  type data
  type hyper

  val upper : hyper : hyper -> x : data -> dst : mat -> unit
  val cross : hyper : hyper -> x : data -> y : data -> dst : mat -> unit
  val diag_mat : hyper : hyper -> x : data -> dst : mat -> unit
  val diag_vec : hyper : hyper -> x : data -> dst : vec -> unit
end

module type Deriv_kernel_fun = sig
  include Kernel_fun

  val upper_deriv_precomp : hyper : hyper -> x : data -> dst : mat -> unit
  val cross_deriv_precomp : hyper : hyper -> x : data -> y : data -> dst : mat -> unit
  val deriv_diag_mat_precomp : hyper : hyper -> x : data -> dst : mat -> unit
  val deriv_diag_vec_precomp : hyper : hyper -> x : data -> dst : vec -> unit

  val upper_deriv : hyper : hyper -> x : data -> dst : mat -> unit
  val cross_deriv : hyper : hyper -> x : data -> y : data -> dst : mat -> unit
  val deriv_diag_mat : hyper : hyper -> x : data -> dst : mat -> unit
  val deriv_diag_vec : hyper : hyper -> x : data -> dst : vec -> unit
end

module type Pseudo_input_kernel = sig
  include Kernel_fun
end

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

let main () =
  ()

let () = main ()
*)

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

val prod :
  ('env1, 'data1) kernel -> ('env2, 'data2) kernel
  -> ('env1 * 'env2, 'data1 * 'data2) kernel

val sum :
  ('env1, 'data1) kernel -> ('env2, 'data2) kernel
  -> ('env1 * 'env2, 'data1 * 'data2) kernel

val sq_const : float -> (unit, 'a) kernel

val sq_var : unit -> (float, 'a) kernel

val inner_prod : 'ip_space -> (unit, 'ip_space) kernel

val weighted_inner_prod : unit kernel

let _ =
  prod (sq_const 1.3)
    (inner_prod 


(* Autom. generate OCaml-code from kernel spec *)
