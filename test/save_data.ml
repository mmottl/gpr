(* File: save_data.ml

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

open Gen_data

module GP = Fitc_gp.Make_deriv (Cov_se_iso.Deriv)
module FITC = GP.FITC.Eval
module FIC = GP.FIC.Eval
module SMD = GP.FITC.Deriv.Optim.SMD

let main () =
  Random.self_init ();

  begin try Unix.mkdir "test/data" ~perm:0o755 with _ -> () end;

  write_mat "inputs" training_inputs;
  write_vec "targets" training_targets;

  let trained =
    let params =
      Cov_se_iso.Eval.Inputs.create_default_kernel_params
        ~n_inducing training_inputs
    in
    let kernel = Cov_se_iso.Eval.Kernel.create params in
    let smd =
      SMD.create
        ~kernel ~n_rand_inducing:n_inducing
        ~inputs:training_inputs ~targets:training_targets ()
    in
    let smd =
      let step = ref 1 in
      let report smd =
        let trained = SMD.get_trained smd in
        let le = FITC.Trained.calc_log_evidence trained in
        let rmse = FITC.Stats.calc_rmse trained in
        printf "iter %4d:  log evidence: %.5f  rmse: %.5f  |gradient|: %f\n%!"
          !step le rmse (SMD.gradient_norm smd);
        incr step;
      in
      SMD.test ~epsabs:3. ~max_iter:1000 ~report smd
    in
    SMD.get_trained smd
(*
    GP.FITC.Deriv.Optim.Gsl.train
      ~report_trained_model:(fun ~iter trained ->
        let le = FITC.Trained.calc_log_evidence trained in
        let rmse = FITC.Stats.calc_rmse trained in
        printf "iter %4d:  log evidence: %.5f  rmse: %.5f\n%!" iter le rmse)
      ~report_gradient_norm:(fun ~iter norm ->
        printf "iter %4d:  |gradient| = %.5f\n%!" iter norm)
      ~kernel ~n_rand_inducing:n_inducing
      ~tol:0.1 ~step:0.1 ~epsabs:3.
      ~inputs:training_inputs ~targets:training_targets ()
*)
  in
  let params =
    let model = FITC.Trained.get_model trained in
    let kernel = FITC.Model.get_kernel model in
    Cov_se_iso.Eval.Kernel.get_params kernel
  in

  let model = FITC.Trained.get_model trained in
  printf "model log evidence: %.9f\n%!" (FITC.Model.calc_log_evidence model);
  printf "log evidence: %.9f\n%!" (FITC.Trained.calc_log_evidence trained);

  let inputs = FITC.Model.get_inputs model in

  let sigma2 = FITC.Model.get_sigma2 (FITC.Trained.get_model trained) in
  write_float "sigma2" sigma2;
  write_float "noise_sigma2" noise_sigma2;

  let inducing = FITC.Model.get_inducing model in
  let inducing_points = FITC.Inducing.get_points inducing in
  let inducing_inputs = FITC.Inputs.calc inducing_points inducing in

  write_mat "inducing_points" inducing_points;

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
