(* File: ocaml_gpr.ml

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

module Args = struct
  type cmd = [ `Train | `Test ]

  type t = {
    cmd : cmd;
    model_file : string;
    with_stddev : bool;
    predictive : bool;
    max_iter : int option;
    n_inducing : int;
    sigma2 : float;
    amplitude : float;
    dim_red : int option;
    log_het_sked : float option;
    multiscale : bool;
    tol : float;
    step : float;
    eps : float;
    verbose : bool;
  }

  let cmd : cmd ref = ref `Train
  let model_file = ref None
  let with_stddev = ref false
  let predictive = ref false
  let max_iter = ref None
  let n_inducing = ref 10
  let sigma2 = ref 1.
  let amplitude = ref 1.
  let dim_red = ref None
  let log_het_sked = ref None
  let multiscale = ref false
  let tol = ref 0.1
  let step = ref 0.1
  let eps = ref 0.1
  let verbose = ref false

  let set_some n_ref n = n_ref := Some n

  let args =
    Arg.align
      [
        (
          "-cmd",
          Arg.Symbol ([ "train"; "test" ], function
            | "train" -> cmd := `Train
            | "test" -> cmd := `Test
            | _ -> assert false  (* impossible *)
          ),
          " train (default) or test model"
        );(
          "-model",
          Arg.String (fun str -> model_file := Some str),
          " model file to use"
        );(
          "-with-stddev",
          Arg.Set with_stddev,
          " make predictions with both mean and variance"
        );(
          "-predictive",
          Arg.Set predictive,
          " standard deviation includes noise level (predictive distribution)"
        );(
          "-max-iter",
          Arg.Int (set_some max_iter),
          " maximum number of optimization steps (default: limitless)"
        );(
          "-n-inducing",
          Arg.Set_int n_inducing,
          sprintf
            " sets number of randomly initialized inducing inputs (default: %d)"
            !n_inducing
        );(
          "-sigma2",
          Arg.Set_float sigma2,
          sprintf " sets initial noise level (default: %f)" !sigma2
        );(
          "-amplitude",
          Arg.Set_float amplitude,
          sprintf " sets initial amplitude level (default: %f)" !amplitude
        );(
          "-dim-red",
          Arg.Int (set_some dim_red),
          " sets dimensionality reduction (default: none)"
        );(
          "-log-het-sked",
          Arg.Float (set_some log_het_sked),
          " turns on / sets log-heteroskedastic \
          noise (may require negative values)"
        );(
          "-multiscale",
          Arg.Set multiscale,
          " turns on multiscale approximation"
        );(
          "-tol",
          Arg.Set_float tol,
          sprintf " sets tolerance for gradient descent (default: %f)" !tol
        );(
          "-step",
          Arg.Set_float step,
          sprintf " sets step size for gradient descent (default: %f)" !step
        );(
          "-eps",
          Arg.Set_float eps,
          sprintf " sets epsilon for gradient descent (default: %f)" !eps
        );(
          "-verbose",
          Arg.Set verbose,
          " prints information while training"
        );
      ]

  let usage_msg = sprintf "%s: -cmd [ train | test ] -model file" Sys.argv.(0)

  let anon_fun _ = failwith "no anonymous arguments allowed"

  let some name opt_ref =
    match !opt_ref with
    | Some v -> v
    | None ->
        eprintf "command line option %s not provided\n\n%!" name;
        prerr_endline usage_msg;
        exit 1

  let get () =
    Arg.parse args anon_fun usage_msg;
    {
      cmd = !cmd;
      model_file = some "model" model_file;
      with_stddev = !with_stddev;
      predictive = !predictive;
      max_iter = !max_iter;
      n_inducing = !n_inducing;
      sigma2 = !sigma2;
      amplitude = !amplitude;
      dim_red = !dim_red;
      log_het_sked = !log_het_sked;
      multiscale = !multiscale;
      tol = !tol;
      step = !step;
      eps = !eps;
      verbose = !verbose;
    }
end

let read_samples () =
  let rex = Str.regexp "," in
  let split str = Array.of_list (Str.split rex str) in
  match try Some (read_line ()) with _ -> None with
  | None -> failwith "no data"
  | Some line ->
      let conv_line line =
        try Array.map ~f:Float.of_string (split line)
        with exc -> Exn.reraisef exc "failure '%s' converting sample" line ()
      in
      let sample = conv_line line in
      let d = Array.length sample in
      let rec loop samples =
        match try Some (read_line ()) with _ -> None with
        | Some line ->
            let floats = conv_line line in
            if Array.length floats <> d then
              failwithf
                "incompatible dimension of sample in line %d: %s"
                (List.length samples + 1) line ()
            else loop (floats :: samples)
        | None -> Array.of_list (List.rev samples)
      in
      loop [sample]

open Lacaml.D

open Gpr

module GP = Fitc_gp.Make_deriv (Cov_se_fat.Deriv)
module FIC = GP.Variational_FIC.Eval

module Model = struct
  type t = {
    sigma2 : float;
    target_mean : float;
    input_means : vec;
    input_stddevs : vec;
    kernel : Cov_se_fat.Eval.Kernel.t;
    inducing_points : FIC.Spec.Inducing.t;
    coeffs : vec;
    co_variance_coeffs : FIC.Model.co_variance_coeffs;
  }
end

let read_training_samples () =
  let samples = read_samples () in
  let n = Array.length samples in
  let d = Array.length samples.(0) - 1 in
  let inputs = Mat.create d n in
  let targets = Vec.create n in
  Array.iteri samples ~f:(fun c0 sample ->
    for r1 = 1 to d do inputs.{r1, c0 + 1} <- sample.(r1 - 1) done;
    targets.{c0 + 1} <- sample.(d));
  inputs, targets

let write_model model_file ~target_mean ~input_means ~input_stddevs trained =
  let oc = open_out model_file in
  let model =
    let model = FIC.Trained.get_model trained in
    let sigma2 = FIC.Model.get_sigma2 model in
    let kernel = FIC.Model.get_kernel model in
    let inducing = FIC.Model.get_inducing model in
    let inducing_points = FIC.Inducing.get_points inducing in
    let mean_predictor = FIC.Mean_predictor.calc_trained trained in
    let coeffs = FIC.Mean_predictor.get_coeffs mean_predictor in
    let co_variance_coeffs = FIC.Model.calc_co_variance_coeffs model in
    {
      Model.
      sigma2; target_mean; input_means; input_stddevs; kernel;
      inducing_points; coeffs; co_variance_coeffs;
    }
  in
  Marshal.to_channel oc model [];
  Out_channel.close oc

exception Bailout

let train args =
  let
    {
      Args.
      model_file; max_iter; n_inducing; sigma2; amplitude; dim_red;
      log_het_sked; multiscale; tol; step; eps = epsabs; verbose
    } = args
  in
  let inputs, targets = read_training_samples () in
  let big_dim = Mat.dim1 inputs in
  let n_inputs = Mat.dim2 inputs in
  let f_inputs = float n_inputs in
  let calc_mean vec = Vec.sum vec /. float (Vec.dim vec) in
  let target_mean = calc_mean targets in
  let targets = Vec.map (fun n -> n -. target_mean) targets in
  let target_variance = Vec.sqr_nrm2 targets /. f_inputs in
  if verbose then eprintf "target variance: %.5f\n%!" target_variance;
  let input_means = Vec.create big_dim in
  let input_stddevs = Vec.create big_dim in
  for i = 1 to big_dim do
    let input = Mat.copy_row inputs i in
    let mean = calc_mean input in
    input_means.{i} <- mean;
    let stddev = sqrt (Vec.ssqr ~c:mean input) in
    input_stddevs.{i} <- stddev;
    for j = 1 to n_inputs do
      inputs.{i, j} <- (inputs.{i, j} -. mean) /. stddev;
    done;
  done;
  let n_inducing = min n_inducing (Vec.dim targets) in
  Random.self_init ();
  let params =
    let log_sf2 = 2. *. log amplitude in
    let d, tproj =
      match dim_red with
      | None -> big_dim, None
      | Some small_dim ->
          let small_dim = min big_dim small_dim in
          let tproj = Mat.random big_dim small_dim in
          Mat.scal (1. /. float big_dim) tproj;
          small_dim, Some tproj
    in
    let log_hetero_skedasticity =
      match log_het_sked with
      | Some log_het_sked -> Some (Vec.make n_inducing log_het_sked)
      | None -> None
    in
    let log_multiscales_m05 =
      if multiscale then Some (Mat.make0 d n_inducing)
      else None
    in
    Cov_se_fat.Params.create
      {
        Cov_se_fat.Params.
        d; log_sf2; tproj; log_hetero_skedasticity; log_multiscales_m05
      }
  in
  let kernel = Cov_se_fat.Eval.Kernel.create params in
  let get_trained_stats trained =
    let { FIC.Stats.smse; msll; mad; maxad } = FIC.Stats.calc trained in
    sprintf
      "MSLL=%7.7f SMSE=%7.7f MAD=%7.7f MAXAD=%7.7f"
      msll smse mad maxad
  in
  let best_trained = ref None in
  let report_trained_model, report_gradient_norm =
    let got_signal = ref false in
    Signal.Expert.set Signal.int (`Handle (fun _ -> got_signal := true));
    let bailout ~iter _ =
      if !got_signal then raise Bailout;
      match max_iter with
      | Some max_iter when iter > max_iter -> raise Bailout
      | _ -> ()
    in
    if verbose then
      let last_eval_time = ref 0. in
      let last_deriv_time = ref 0. in
      let maybe_print last_time line =
        let now = Unix.gettimeofday () in
        if !last_time +. 1. < now then begin
          last_time := now;
          prerr_endline line;
        end
      in
      Some (fun ~iter trained ->
        best_trained := Some trained;
        bailout ~iter ();
        maybe_print last_eval_time
          (sprintf "iter %4d: %s" iter (get_trained_stats trained))),
      Some (fun ~iter norm ->
        bailout ~iter ();
        maybe_print last_deriv_time
          (sprintf "iter %4d: |gradient|=%.5f" iter norm))
    else Some bailout, None
  in
  match
    try
      Some (
        GP.Variational_FIC.Deriv.Optim.Gsl.train
          ?report_trained_model ?report_gradient_norm
          ~kernel ~sigma2 ~n_rand_inducing:n_inducing
          ~tol ~step ~epsabs ~inputs ~targets ())
    with GP.FIC.Deriv.Optim.Gsl.Optim_exception Bailout -> !best_trained
  with
  | None -> ()
  | Some trained ->
      if verbose then eprintf "result: %s\n%!" (get_trained_stats trained);
      write_model model_file ~target_mean ~input_means ~input_stddevs trained

let read_test_samples big_dim =
  let samples = read_samples () in
  let n = Array.length samples in
  if n = 0 then Mat.empty
  else begin
    let input_dim = Array.length samples.(0) in
    if input_dim <> big_dim then
      failwithf
        "incompatible dimension of inputs (%d), expected %d"
        input_dim big_dim ();
    let inputs = Mat.create big_dim n in
    Array.iteri samples ~f:(fun c0 sample ->
      for r1 = 1 to big_dim do inputs.{r1, c0 + 1} <- sample.(r1 - 1) done);
    inputs
  end

let read_model model_file : Model.t =
  let ic = open_in model_file in
  let model = Marshal.from_channel ic in
  In_channel.close ic;
  model

let test args =
  let { Args.model_file; with_stddev; predictive } = args in
  let
    {
      Model.
      sigma2; target_mean; input_means; input_stddevs; kernel;
      inducing_points; coeffs; co_variance_coeffs
    } = read_model model_file
  in
  let big_dim = Vec.dim input_stddevs in
  let inputs = read_test_samples big_dim in
  let n_inputs = Mat.dim2 inputs in
  for i = 1 to big_dim do
    let mean = input_means.{i} in
    let stddev = input_stddevs.{i} in
    for j = 1 to n_inputs do
      inputs.{i, j} <- (inputs.{i, j} -. mean) /. stddev;
    done;
  done;
  let mean_predictor = FIC.Mean_predictor.calc inducing_points ~coeffs in
  let inducing = FIC.Inducing.calc kernel inducing_points in
  let inputs = FIC.Inputs.calc inputs inducing in
  let means = FIC.Means.get (FIC.Means.calc mean_predictor inputs) in
  let renorm_mean mean = mean +. target_mean in
  if with_stddev then
    let co_variance_predictor =
      FIC.Co_variance_predictor.calc kernel inducing_points co_variance_coeffs
    in
    let vars = FIC.Variances.calc co_variance_predictor ~sigma2 inputs in
    let vars = FIC.Variances.get ~predictive vars in
    Vec.iteri (fun i pre_mean ->
      let mean = renorm_mean pre_mean in
      printf "%f,%f\n" mean (sqrt vars.{i})) means
  else Vec.iter (fun mean -> printf "%f\n" (renorm_mean mean)) means

let main () =
  let args = Args.get () in
  match args.Args.cmd with
  | `Train -> train args
  | `Test -> test args

let () = main ()
