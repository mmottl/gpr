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

  let k = { Gauss_all_vec_spec.Kernel.a = -0.5; b = -1.0 }
end

module Wiener = struct
  include Fitc.Make (Wiener_all_vec)

  let k = { Wiener_all_vec_spec.Kernel.a = 0. }
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

let training_inputs, training_targets, inducing_inputs = get_training ()

let main () =
  let inducing = FITC.Inducing.calc k inducing_inputs in
  let inputs = FITC.Inputs.calc inducing training_inputs in

  let model = FITC.Model.calc inputs ~sigma2:noise_sigma2 in
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

(* let () = main () *)

module Gauss_deriv = struct
  include Make_FITC_deriv (Gauss_deriv_all_vec)

  let k = Gauss.k
end

open Gauss_deriv

let with_sigma2 sigma2 =
  let inducing = Deriv.Inducing.calc k inducing_inputs in
  let inputs = Deriv.Inputs.calc inducing training_inputs in
  let model = Deriv.Model.calc ~sigma2 inputs in
  let trained = Deriv.Trained.calc model ~targets:training_targets in
  let eval_trained = Deriv.Trained.calc_eval trained in
  let eval_trained_evidence = Eval.Trained.calc_evidence eval_trained in
  printf "eval trained evidence: %.15f@." eval_trained_evidence;
  let _, deriv_evidence = Deriv.Model.calc_evidence_sigma2 model in
  let trained_deriv_evidence =
    Deriv.Trained.calc_evidence_sigma2 trained deriv_evidence
  in
  printf "deriv trained evidence: %.15f@." trained_deriv_evidence;
  eval_trained_evidence

let main () =
  let e1 = with_sigma2 noise_sigma2 in
  let epsilon = 10e-10 in
  let e2 = with_sigma2 (noise_sigma2 +. epsilon) in
  printf "finite: %.15f\n" ((e2 -. e1) /. epsilon)

let find_sigma2 () =
  let inducing = Deriv.Inducing.calc k inducing_inputs in
  let inputs = Deriv.Inputs.calc inducing training_inputs in
  let eval_inputs = Deriv.Inputs.calc_eval inputs in

  let model_ref = ref None in

  let multim_f ~x =
    let sigma2 = exp x.{0} in
    let model =
      match !model_ref with
      | None ->
          let model = Eval.Model.calc ~sigma2 eval_inputs in
          model_ref := Some model;
          model
      | Some model -> Eval.Model.update_sigma2 model sigma2
    in
    let trained = Eval.Trained.calc model ~targets:training_targets in
    let evidence = Eval.Trained.calc_evidence trained in
    -. evidence
  in

  let dmodel_ref = ref None in

  let multim_dcommon ~x ~g =
    let sigma2 = exp x.{0} in
    let dmodel =
      match !dmodel_ref with
      | None ->
          let dmodel = Deriv.Model.calc ~sigma2 inputs in
          dmodel_ref := Some dmodel;
          dmodel
      | Some dmodel -> Deriv.Model.update_sigma2 dmodel sigma2
    in
    let _, model_deriv_evidence = Deriv.Model.calc_evidence_sigma2 dmodel in
    let trained = Deriv.Trained.calc dmodel ~targets:training_targets in
    let devidence =
      Deriv.Trained.calc_evidence_sigma2 trained model_deriv_evidence
    in
    g.{0} <- -. devidence *. sigma2;
    trained
  in

  let multim_df ~x ~g =
    ignore (multim_dcommon ~x ~g)
  in

  let multim_fdf ~x ~g =
    let trained = multim_dcommon ~x ~g in
    let evidence =
      Eval.Trained.calc_evidence (Deriv.Trained.calc_eval trained)
    in
    -. evidence
  in

  let multim_fun_fdf =
    {
      Gsl_fun.
      multim_f = multim_f;
      multim_df = multim_df;
      multim_fdf = multim_fdf;
    }
  in

  let init = Gsl_vector.create ~init:(log 1.) 1 in

  let module Gd = Gsl_multimin.Deriv in
  let mumin =
    Gd.make Gd.VECTOR_BFGS2 1
      multim_fun_fdf ~x:init ~step:10e-2 ~tol:10e-4
  in
  let x = Gsl_vector.create 1 in
  let rec loop last_evidence =
    let nll = Gd.minimum ~x mumin in
    let evidence = -. nll in
    let diff = abs_float (1. -. (evidence /. last_evidence)) in
    printf "diff: %f\n%!" diff;
    if diff < 0.001 then nll
    else (
      printf "evidence: %f\n%!" evidence;
      Gd.iterate mumin;
      loop evidence)
  in
  let nll = loop neg_infinity in
  -. nll, exp x.{0}

let main () =
  let evidence, sigma2 = find_sigma2 () in
  printf "evidence: %.15f  sigma2: %.15f\n" evidence sigma2

let main () =
  let sigma2 = noise_sigma2 in

  let epsilon = 10e-6 in

  let inducing = Deriv.Inducing.calc k inducing_inputs in
  let inputs = Deriv.Inputs.calc inducing training_inputs in
  let model = Deriv.Model.calc ~sigma2 inputs in

  let new_k =
    {
      k with
      Gauss_all_vec_spec.Kernel.a =
        k.Gauss_all_vec_spec.Kernel.a +. epsilon
    }
  in
  let inducing = Deriv.Inducing.calc new_k inducing_inputs in
  let inputs = Deriv.Inputs.calc inducing training_inputs in
  let model2 = Deriv.Model.calc ~sigma2 inputs in

  let hyper_model = Deriv.Model.prepare_hyper model in
  let mev, model_evidence =
    Deriv.Model.calc_evidence hyper_model
      Gauss_deriv_all_vec.Deriv_spec.Hyper.A
  in

  let mf1 = Eval.Model.calc_evidence (Deriv.Model.calc_eval model) in
  let mf2 = Eval.Model.calc_evidence (Deriv.Model.calc_eval model2) in

  printf "mdevidence: %f\n%!" mev;
  printf "mdfinite:   %f\n%!" ((mf2 -. mf1) /. epsilon);

  let trained = Deriv.Trained.calc model ~targets:training_targets in
  let trained2 = Deriv.Trained.calc model2 ~targets:training_targets in

  let hyper_trained = Deriv.Trained.prepare_hyper trained hyper_model in
  let deriv = Deriv.Trained.calc_evidence hyper_trained model_evidence in

  let f1 = Eval.Trained.calc_evidence (Deriv.Trained.calc_eval trained) in
  let f2 = Eval.Trained.calc_evidence (Deriv.Trained.calc_eval trained2) in

  printf "evidence: %f\n%!" f1;
  printf "devidence: %f\n%!" deriv;
  printf "dfinite:   %f\n%!" ((f2 -. f1) /. epsilon)

let () = main ()
