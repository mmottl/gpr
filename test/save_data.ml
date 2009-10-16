open Format

open Lacaml.Impl.D
open Lacaml.Io

open Gpr
open Utils

open Gen_data

module K = Fitc_gp.Make_deriv (Cov_se_iso.Deriv)
module FITC_all = K.FITC
module FITC = K.FITC.Eval
module FIC = K.FIC.Eval

let main () =
  Random.self_init ();
  begin try Unix.mkdir "data" 0o755 with _ -> () end;
  write_mat "inputs" training_inputs;
  write_vec "targets" training_targets;

  let params =
    Cov_se_iso.Eval.Inputs.create_default_kernel_params
      ~n_inducing training_inputs
  in
  let kernel = Cov_se_iso.Eval.Kernel.create params in
  let trained =
    let report_trained_model ~iter trained =
      let log_evidence = FITC.Trained.calc_log_evidence trained in
      let rmse = FITC.Trained.calc_rmse trained in
      printf "iter %4d:  log evidence: %.5f  rmse: %.5f\n%!"
        iter log_evidence rmse
    in
    let report_gradient_norm ~iter norm =
      Printf.printf "iter %4d:  |gradient| = %.5f\n%!" iter norm
    in
    let inducing_points =
      FITC.Inducing.choose_n_random_inputs kernel ~n_inducing training_inputs
    in
    let all_hypers =
      FITC_all.Deriv.Spec.Hyper.get_all kernel inducing_points
    in
    FITC_all.Deriv.Test.self_test
      kernel inducing_points training_inputs
      ~sigma2:noise_sigma2 ~targets:training_targets `Sigma2;
    Array.iter (fun hyper ->
        let hyper_str =
          match hyper with
          | `Log_sf2 -> "Log_sf2"
          | `Log_ell -> "Log_ell"
          | `Inducing_hyper _ -> "Inducing_hyper"
(*
          | `Proj _ -> "Proj"
          | `Log_hetero_skedasticity _ -> "Log_hetero_skedasticity"
          | `Log_multiscale_m05 _ -> "Log_multiscale_m05"
*)
        in
        printf "-------- testing finite difference for hyper: %s\n%!" hyper_str;
        FITC_all.Deriv.Test.check_deriv_hyper
          kernel inducing_points training_inputs hyper;
        FITC_all.Deriv.Test.self_test
          kernel inducing_points training_inputs
          ~sigma2:noise_sigma2 ~targets:training_targets (`Hyper hyper)
      )
      all_hypers;
    FITC_all.Deriv.Optim.Gsl.train
      ~report_trained_model ~report_gradient_norm
      ~kernel ~n_rand_inducing:n_inducing
      ~tol:0.1 ~step:0.1 ~epsabs:3.
      ~inputs:training_inputs ~targets:training_targets ()
(*
    let inputs =
      let inducing_points =
        FITC.Inducing.choose_n_random_inputs kernel ~n_inducing training_inputs
      in
      let inducing = FITC.Inducing.calc kernel inducing_points in
      FITC.Inputs.calc inducing training_inputs
    in
    let model = FITC.Model.calc inputs ~sigma2:noise_sigma2 in
    FITC.Trained.calc model ~targets:training_targets
*)
  in

  let model = FITC.Trained.get_model trained in
  printf "model log evidence: %.9f@." (FITC.Model.calc_log_evidence model);
  printf "log evidence: %.9f@." (FITC.Trained.calc_log_evidence trained);

  let inputs = FITC.Model.get_inputs model in

  let sigma2 = FITC.Model.get_sigma2 (FITC.Trained.get_model trained) in
  write_float "sigma2" sigma2;
  write_float "noise_sigma2" noise_sigma2;

  let inducing = FITC.Model.get_inducing model in
  let inducing_points = FITC.Inducing.get_points inducing in
  let inducing_inputs = FITC.Inputs.calc inducing inducing_points in

  write_mat "inducing_points" inducing_points;

(*
  write_float "log_sf2"
    (params :> Cov_se_fat.Params.params).Cov_se_fat.Params.log_sf2;
*)

  write_float "log_sf2" params.Cov_se_iso.Params.log_sf2;
  write_float "log_ell" params.Cov_se_iso.Params.log_ell;

  let mean_predictor = FITC.Mean_predictor.calc_trained trained in
  let co_variance_predictor = FITC.Co_variance_predictor.calc_model model in

  let means = FITC.Means.calc mean_predictor inputs in
  write_vec "means" (FITC.Means.get means);

  let inducing_means =
    FITC.Means.get (FITC.Means.calc mean_predictor inducing_inputs)
  in
  write_vec "inducing_means" inducing_means;

  let mean, variance =
    let input = Mat.col inducing_points n_inducing in
    write_vec "one_inducing" input;
    let induced = FITC.Input.calc inducing input in
    let mean = FITC.Mean.get (FITC.Mean.calc mean_predictor induced) in
    let variance =
      FITC.Variance.get ~predictive:false
        (FITC.Variance.calc co_variance_predictor ~sigma2 induced)
    in
    mean, variance
  in
  write_float "one_mean" mean;
  write_float "one_variance" variance;

  let variances =
    FITC.Variances.get ~predictive:false
      (FITC.Variances.calc_model_inputs model)
  in
  write_vec "variances" variances;

  let inducing_variances =
    FITC.Variances.get ~predictive:false
      (FITC.Variances.calc co_variance_predictor ~sigma2 inducing_inputs)
  in
  write_vec "inducing_variances" inducing_variances;

  let fitc_covariances = FITC.Covariances.calc_model_inputs model in
  let fitc_sampler =
    FITC.Cov_sampler.calc ~predictive:false means fitc_covariances
  in
  let fitc_samples = FITC.Cov_sampler.samples fitc_sampler ~n:3 in
  write_vec "sample1" (Mat.col fitc_samples 1);
  write_vec "sample2" (Mat.col fitc_samples 2);
  write_vec "sample3" (Mat.col fitc_samples 3);

  let fic_covariances = FIC.Covariances.calc_model_inputs model in
  let fic_sampler =
    FIC.Cov_sampler.calc ~predictive:false means fic_covariances
  in
  let fic_samples = FIC.Cov_sampler.samples fic_sampler ~n:3 in
  write_vec "fic_sample1" (Mat.col fic_samples 1);
  write_vec "fic_sample2" (Mat.col fic_samples 2);
  write_vec "fic_sample3" (Mat.col fic_samples 3)

let () = main ()
