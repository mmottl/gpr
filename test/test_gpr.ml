open Format

open Lacaml.Impl.D
open Lacaml.Io

open Gpr

module FIC = Fitc.Make_FIC (Kernel.Gauss)

let timing name f =
  let t1 = Unix.times () in
  let res = f () in
  let t2 = Unix.times () in
  printf "%s %.2f@." name (t2.Unix.tms_utime -. t1.Unix.tms_utime);
  res

let rng = Gsl_rng.make (Gsl_rng.default ())

let noise_sigma = 0.5
let noise_sigma2 = noise_sigma *. noise_sigma

let f ?(noise = false) x =
  sin x /. x +. if noise then Gsl_randist.gaussian rng ~sigma:noise_sigma else 0.

let get_data ?noise n =
  let inputs = Mat.create 1 n in
  let targets = Vec.create n in
  for i = 1 to n do
    let x = Random.float 10. -. 5. in
    inputs.{1, i} <- x;
    targets.{i} <- f ?noise x;
  done;
  inputs, targets

let dir = ref "/"

let gen_write pp file obj =
  let oc = open_out (Filename.concat !dir file) in
  fprintf (formatter_of_out_channel oc) "%a@." pp obj;
  close_out oc

let write_float = gen_write pp_print_float
let write_vec = gen_write pp_fvec
let write_mat = gen_write pp_fmat

let get_training () =
(*
  Random.self_init ();
  Gsl_rng.set_default_seed (Random.nativeint Nativeint.max_int);
*)
  Random.init 0;
  Gsl_rng.set rng 0n;
  let n_inputs = 1000 in
  let training_inputs, training_targets = get_data ~noise:true n_inputs in
  write_float "sigma2" noise_sigma2;
  write_mat "inputs" training_inputs;
  write_vec "targets" training_targets;
  let n_inducing_inputs = 10 in
  let inducing_inputs = Mat.create 1 n_inducing_inputs in
  for i = 1 to n_inducing_inputs do
    inducing_inputs.{1, i} <- Random.float 10. -. 5.;
  done;
  write_mat "inducing_inputs" inducing_inputs;
  training_inputs, training_targets, inducing_inputs

let main1 () =
  dir := "data1";
(*   timing "all" (fun () -> *)
  let training_inputs, training_targets, inducing_inputs = get_training () in
  let common =
    {
      FIC.Spec.
      kernel = ();
      sigma2 = noise_sigma2;
      inducing_inputs = inducing_inputs;
    }
  in
(*   timing "all" (fun () -> *)
  let trained = FIC.Trained.train common ~inputs:training_inputs ~targets:training_targets in
  let model = FIC.Full_predictor.of_trained trained in
  let means, variances =
    FIC.Full_predictor.means_variances model training_inputs
  in
  printf "neg_log_likelihood: %.3f@." (FIC.Trained.neg_log_likelihood trained);
  write_vec "means" means;
  write_vec "variances" variances

open Experimental

module FITC_spec = struct
  module Kernel = Kernel.Gauss

  let get_sigma2 _ = noise_sigma2
  let jitter = 10e-9
end


module FITC = Make_FITC (FITC_spec)

open FITC

let main2 () =
  dir := "data2";
  let training_inputs, training_targets, inducing_inputs = get_training () in
(*   timing "all" (fun () -> *)
  let inducing = Inducing.calc () inducing_inputs in
  let reduceds = Reduceds.calc inducing training_inputs in
  let model = Model.calc reduceds in
  let trained = Trained.calc model ~targets:training_targets in
  printf "neg_log_likelihood: %.3f@." (Trained.neg_log_likelihood trained);
  let weights = Weights.calc trained in
  let means = Means.calc_inputs weights model in
  let inducing_means = Means.copy (Means.calc_inducing weights model) in
  let means_vec = Means.copy means in
  let variances = Variances.copy ~predictive:false (Variances.calc_inputs model) in
  write_vec "means" means_vec;
  write_vec "inducing_means" inducing_means;
  write_vec "variances" variances;
  let covariances = Covariances.calc_inputs model in
  let samplers = Samplers.calc ~predictive:false means covariances in
  write_vec "sample1" (Samplers.sample samplers);
  write_vec "sample2" (Samplers.sample samplers);
  write_vec "sample3" (Samplers.sample samplers);
  write_vec "sample4" (Samplers.sample samplers);
  write_vec "sample5" (Samplers.sample samplers)

(*
let main () =
  let d = 50 in
  let n = 50000 in
  let mat = Mat.random d n in
  let vec = Vec.random n in
  let common =
    let n_inducing_inputs = 100 in
    let inducing_inputs = Mat.random d n_inducing_inputs in
    {
      FIC.Common.
      kernel = ();
      sigma2 = 0.1;
      inducing_inputs = inducing_inputs;
    }
  in
  timing "all" (fun () ->
    let trained = FIC.Trained.train common ~inputs:mat ~targets:vec in
    let model = FIC.Mean_predictor.of_trained trained in
    let input = Vec.random d in
    let mean = FIC.Mean_predictor.mean model input in
    printf "%f@." mean;
    ());
  ()
*)

let () =
  main1 ();
  main2 ()
