open Format

open Lacaml.Impl.D
open Lacaml.Io

module FIC = Make_FIC (Kernel)

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
    let x = Random.float 20. -. 10. in
    inputs.{1, i} <- x;
    targets.{i} <- f ?noise x;
  done;
  inputs, targets

let gen_write pp file obj =
  let oc = open_out (Filename.concat "data" file) in
  fprintf (formatter_of_out_channel oc) "%a@." pp obj;
  close_out oc

let write_mat = gen_write pp_fmat
let write_vec = gen_write pp_fvec

let main () =
(*
  Random.self_init ();
  Gsl_rng.set_default_seed (Random.nativeint Nativeint.max_int);
*)
  let n_inputs = 1000 in
  let training_inputs, training_targets = get_data ~noise:true n_inputs in
  write_mat "inputs" training_inputs;
  write_vec "targets" training_targets;
  let common =
    let n_inducing_points = 50 in
    let inducing_points = Mat.create 1 n_inducing_points in
    for i = 1 to n_inducing_points do
      inducing_points.{1, i} <- Random.float 20. -. 10.;
    done;
    write_mat "inducing_points" inducing_points;
    {
      FIC.Common.
      kernel = ();
      sigma2 = noise_sigma2;
      inducing_points = inducing_points;
    }
  in
(*   timing "all" (fun () -> *)
    let trained = FIC.Trained.train common ~inputs:training_inputs ~targets:training_targets in
    printf "neg_log_likelihood: %.3f@." (FIC.Trained.neg_log_likelihood trained);
    let model = FIC.Full_predictor.of_trained `Fitc ~fmt:`Chols trained in
(*     let means = Vec.create n_inputs in *)
    let means, variances =
      FIC.Full_predictor.means_variances model ~inputs:training_inputs
    in
    write_vec "means" means;
    write_vec "variances" variances

(*
let main () =
  let d = 50 in
  let n = 50000 in
  let mat = Mat.random d n in
  let vec = Vec.random n in
  let common =
    let n_inducing_points = 100 in
    let inducing_points = Mat.random d n_inducing_points in
    {
      FIC.Common.
      kernel = ();
      sigma2 = 0.1;
      inducing_points = inducing_points;
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

let () = main ()
