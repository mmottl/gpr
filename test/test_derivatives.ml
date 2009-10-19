open Format
open Lacaml.Impl.D

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
  let all_hypers = FITC.Deriv.Spec.Hyper.get_all kernel inducing_points in
  let sigma2 = 1. in
  FITC.Deriv.Test.self_test
    kernel inducing_points training_inputs
    ~sigma2 ~targets:training_targets `Sigma2;
  Array.iter (fun hyper ->
    let module Csf = Cov_se_fat in
      let hyper_str =
        match hyper with
        | `Log_sf2 -> "Log_sf2"
        | `Inducing_hyper { Csf.Inducing_hyper.ind = ind; dim = dim } ->
            sprintf "Inducing_hyper { ind = %d; dim = %d }" ind dim
        | `Proj { Csf.Proj_hyper.big_dim = big_dim; small_dim = small_dim } ->
            sprintf "Proj { big_dim = %d; small_dim = %d }" big_dim small_dim
        | `Log_hetero_skedasticity dim ->
            sprintf "Log_hetero_skedasticity %d" dim
        | `Log_multiscale_m05 { Csf.Inducing_hyper.ind = ind; dim = dim } ->
            sprintf "Log_multiscale_m05 { ind = %d; dim = %d }" ind dim
      in
      printf "-------- testing finite differences for hyper: %s\n%!" hyper_str;
      FITC.Deriv.Test.check_deriv_hyper
        kernel inducing_points training_inputs hyper;
      FITC.Deriv.Test.self_test
        kernel inducing_points training_inputs
        ~sigma2 ~targets:training_targets (`Hyper hyper)
    )
    all_hypers

let () = main ()
