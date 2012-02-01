(* File: test_derivatives.ml

   OCaml-GPR - Gaussian Processes for OCaml

     Copyright (C) 2009-  Markus Mottl
     email: markus.mottl@gmail.com
     WWW:   http://www.ocaml.info

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*)

open Core.Std
open Lacaml.D

open Gpr

module GP = Fitc_gp.Make_deriv (Cov_se_fat.Deriv)
module FITC = GP.FITC

let main () =
  Random.self_init ();
  let n_inducing = 5 in
  let training_inputs = Mat.random 3 10 in
  let training_targets = Vec.random 10 in
  let params =
    Cov_se_fat.Eval.Inputs.create_default_kernel_params
      ~n_inducing training_inputs
  in
  let kernel = Cov_se_fat.Eval.Kernel.create params in
  let inducing_points =
    FITC.Eval.Inducing.choose_n_random_inputs kernel ~n_inducing training_inputs
  in
  let all_hypers =
    FITC.Deriv.Spec.Hyper.get_all kernel inducing_points training_inputs
  in
  let sigma2 = 1. in
  FITC.Deriv.Test.self_test
    kernel inducing_points training_inputs
    ~sigma2 ~targets:training_targets `Sigma2;
  Array.iter all_hypers ~f:(fun hyper ->
    let module Csf = Cov_se_fat in
      let hyper_str =
        match hyper with
        | `Log_sf2 -> "Log_sf2"
        | `Inducing_hyper { Csf.Inducing_hyper.ind; dim } ->
            sprintf "Inducing_hyper { ind = %d; dim = %d }" ind dim
        | `Proj { Csf.Proj_hyper.big_dim; small_dim } ->
            sprintf "Proj { big_dim = %d; small_dim = %d }" big_dim small_dim
        | `Log_hetero_skedasticity dim ->
            sprintf "Log_hetero_skedasticity %d" dim
        | `Log_multiscale_m05 { Csf.Inducing_hyper.ind; dim } ->
            sprintf "Log_multiscale_m05 { ind = %d; dim = %d }" ind dim
      in
      printf "-------- testing finite differences for hyper: %s\n%!" hyper_str;
      FITC.Deriv.Test.check_deriv_hyper
        kernel inducing_points training_inputs hyper;
      FITC.Deriv.Test.self_test
        kernel inducing_points training_inputs
        ~sigma2 ~targets:training_targets (`Hyper hyper)
    )

let () = main ()
