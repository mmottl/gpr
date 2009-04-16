open Format

open Lacaml.Impl.D
open Lacaml.Io

open Gpr
open Utils

open Test_kernels.SE_iso
open Gen_data

let main () =
  let sigma2 = noise_sigma2 in

  let epsilon = 1e-6 in

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

  let new_kernel =
    let params = Eval.Spec.Kernel.get_params kernel in
    let new_log_ell = params.Cov_se_iso.Params.log_ell +. epsilon in
    let new_params =
      { params with Cov_se_iso.Params.log_ell = new_log_ell }
    in
    Eval.Spec.Kernel.create new_params
  in

  let inducing = Deriv.Inducing.calc new_kernel deriv_prep_inducing in
  let inputs = Deriv.Inputs.calc inducing deriv_prep_inputs in
  let model2 = Deriv.Model.calc ~sigma2 inputs in

  let hyper_model = Deriv.Model.prepare_hyper model in
  let mev, model_log_evidence =
    Deriv.Model.calc_log_evidence hyper_model `Log_ell
  in

  let mf1 = Eval.Model.calc_log_evidence (Deriv.Model.calc_eval model) in
  let mf2 = Eval.Model.calc_log_evidence (Deriv.Model.calc_eval model2) in

  print_float "mdlog_evidence" mev;
  print_float "mdfinite" ((mf2 -. mf1) /. epsilon);

  let trained = Deriv.Trained.calc model ~targets:training_targets in
  let trained2 = Deriv.Trained.calc model2 ~targets:training_targets in

  let hyper_trained = Deriv.Trained.prepare_hyper trained hyper_model in
  let deriv =
    Deriv.Trained.calc_log_evidence hyper_trained model_log_evidence
  in

  let f1 = Eval.Trained.calc_log_evidence (Deriv.Trained.calc_eval trained) in
  let f2 = Eval.Trained.calc_log_evidence (Deriv.Trained.calc_eval trained2) in

  print_float "log evidence" f1;
  print_float "dlog_evidence" deriv;
  print_float "dfinite" ((f2 -. f1) /. epsilon)

(*
let main () =
  Lacaml.Io.Context.set_dim_defaults (Some (Context.create 5));

  let sigma2 = noise_sigma2 in

  let epsilon = 1e-6 in

  let module Eval = FITC.Eval in
  let module Deriv = FITC.Deriv in

  Utils.print_mat "inducing_inputs" inducing_inputs;
  let run () =
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

    let hyper_model = Deriv.Model.prepare_hyper model in
    let mev, _model_log_evidence =
      Deriv.Model.calc_log_evidence hyper_model (`Inducing_hyper {
        Cov_se_iso.ind = 3; dim = 1 })
    in

    let mf = Eval.Model.calc_log_evidence (Deriv.Model.calc_eval model) in
    mev, mf
  in

  let mev, mf1 = run () in
  inducing_inputs.{1, 3} <- inducing_inputs.{1, 3} +. epsilon;
  let _, mf2 = run () in

  print_float "mdlog_evidence" mev;
  print_float "mdfinite" ((mf2 -. mf1) /. epsilon)
*)

let () = main ()
