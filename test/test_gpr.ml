open Format

open Lacaml.Impl.D
open Lacaml.Io

open Gpr
open Kernel_impl
open Utils
open Fitc

let n_inputs = 10
let n_inducing_inputs = 5
let noise_sigma = 1.5
let noise_sigma2 = noise_sigma *. noise_sigma

module Gauss = struct
  include Fitc.Make (Gauss_all_vec)

  let k =
    {
      Gauss_all_vec_spec.Kernel.
      a = -0.5;
      b = -1.0;
      sigma2 = noise_sigma2;
    }
end

module Wiener = struct
  include Fitc.Make (Wiener_all_vec)

  let k =
    {
      Wiener_all_vec_spec.Kernel.
      a = 0.;
      sigma2 = noise_sigma2;
    }
end

module All = Gauss

open All

let f ?(with_noise = false) x =
  let v = sin (3. *. x) /. x +. (x -. 3.) /. (x *. x +. 1.) in
  if with_noise then v +. Gsl_randist.gaussian default_rng ~sigma:noise_sigma
  else v

let get_data ?with_noise n =
  let inputs = Mat.create 1 n in
  let targets = Vec.create n in
  for i = 1 to n do
    let x = Random.float 10. -. 5. in
    inputs.{1, i} <- x;
    targets.{i} <- f ?with_noise x;
  done;
  inputs, targets

let get_training () =
  let training_inputs, training_targets = get_data ~with_noise:true n_inputs in
  write_float "sigma2" noise_sigma2;
  write_mat "inputs" training_inputs;
  write_vec "targets" training_targets;
  let inducing_inputs = Mat.create 1 n_inducing_inputs in
  for i = 1 to n_inducing_inputs do
    inducing_inputs.{1, i} <-
      -5. +. float i *. 10. /. float (n_inducing_inputs + 1)
  done;
  write_mat "inducing_inputs" inducing_inputs;
  training_inputs, training_targets, inducing_inputs

let main () =
  let training_inputs, training_targets, inducing_inputs = get_training () in
  let inducing = FITC.Inducing.calc k inducing_inputs in
  let inputs = FITC.Inputs.calc inducing training_inputs in

  let model = FITC.Model.calc inputs in
  printf "model evidence: %.9f@." (FITC.Model.calc_evidence model);

  let trained = FITC.Trained.calc model ~targets:training_targets in
  printf "evidence: %.9f@." (FITC.Trained.calc_evidence trained);

  let weights = FITC.Weights.calc trained in

  let means = FITC.Means.calc_model_inputs weights in
  let inducing_means =
    FITC.Means.Inducing.get (FITC.Means.Inducing.calc weights)
  in
  write_vec "inducing_means" inducing_means;
  let means_vec = FITC.Means.get means in
  write_vec "means" means_vec;

  let variances =
    FITC.Variances.get ~predictive:false
      (FITC.Variances.calc_model_inputs model)
  in
  write_vec "variances" variances;

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

module Gauss_deriv = struct
  include Make_FITC_deriv (Gauss_deriv_all_vec)

  let k =
    {
      Gauss_all_vec_spec.Kernel.
      a = -0.5;
      b = -1.0;
      sigma2 = noise_sigma2;
    }
end

open Gauss_deriv

let training_inputs, training_targets, inducing_inputs = get_training ()

let with_k k =
  let inducing = Deriv.Inducing.calc k inducing_inputs in
  let inputs = Deriv.Inputs.calc inducing training_inputs in
  let model = Deriv.Model.calc inputs in
  let trained = Deriv.Trained.calc model ~targets:training_targets in
  let eval_trained = Deriv.Trained.calc_eval trained in
  let eval_trained_evidence = Eval.Trained.calc_evidence eval_trained in
  printf "eval trained evidence: %.9f@." eval_trained_evidence;
  let deriv_evidence =
    Deriv.Trained.calc_evidence trained (Deriv.Model.calc_evidence_sigma2 model)
  in
  printf "eval deriv evidence: %.9f@." deriv_evidence;
  eval_trained_evidence

let main () =
  let e1 = with_k k in
  let epsilon = 10e-10 in
  let new_k =
    {
      k with
      Gauss_all_vec_spec.Kernel.
      sigma2 = k.Gauss_all_vec_spec.Kernel.sigma2 +. epsilon;
    }
  in
  let e2 = with_k new_k in
  printf "finite: %.9f\n" ((e2 -. e1) /. epsilon)

let () = main ()
