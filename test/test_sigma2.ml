open Format

open Lacaml.Impl.D
open Lacaml.Io

open Gpr
open Utils

open Test_kernels.SE_iso
open Gen_data

let with_sigma2 ?(print = false) sigma2 =
  let module Eval = FITC.Eval in
  let module Deriv = FITC.Deriv in
  let eval_prep_inducing = Eval.Inducing.Prepared.calc inducing_inputs in
  let deriv_prep_inducing = Deriv.Inducing.Prepared.calc eval_prep_inducing in
  let inducing = Deriv.Inducing.calc kernel deriv_prep_inducing in
  let eval_prep_inputs =
    Eval.Inputs.Prepared.calc eval_prep_inducing training_inputs
  in
  let deriv_prep_inputs =
    Deriv.Inputs.Prepared.calc deriv_prep_inducing eval_prep_inputs
  in
  let inputs = Deriv.Inputs.calc inducing deriv_prep_inputs in

  let model = Deriv.Model.calc ~sigma2 inputs in
  let eval_model = Deriv.Model.calc_eval model in
  let eval_model_log_evidence = Eval.Model.calc_log_evidence eval_model in
  if print then
    printf "eval model log evidence: %.15f@." eval_model_log_evidence;
  let model_deriv_log_evidence = Deriv.Model.calc_log_evidence_sigma2 model in
  if print then
    printf "deriv model log evidence: %.15f@.@." model_deriv_log_evidence;

  let trained = Deriv.Trained.calc model ~targets:training_targets in
  let eval_trained = Deriv.Trained.calc_eval trained in
  let eval_trained_log_evidence = Eval.Trained.calc_log_evidence eval_trained in
  if print then
    printf "eval trained log evidence: %.15f@." eval_trained_log_evidence;
  let trained_deriv_log_evidence =
    Deriv.Trained.calc_log_evidence_sigma2 trained
  in
  if print then
    printf "deriv trained log evidence: %.15f@.@." trained_deriv_log_evidence;

  eval_model_log_evidence, eval_trained_log_evidence

let main () =
  let me1, te1 = with_sigma2 ~print:true noise_sigma2 in
  let epsilon = 1e-6 in
  let me2, te2 = with_sigma2 (noise_sigma2 +. epsilon) in
  printf "finite diff model: %.15f@\n" ((me2 -. me1) /. epsilon);
  printf "finite diff trained: %.15f@\n" ((te2 -. te1) /. epsilon)

let () = main ()
