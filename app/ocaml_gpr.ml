open Printf

module Args = struct
  type cmd = [ `train | `test ]

  type t = {
    cmd : cmd;
    model_file : string;
    with_stddev : bool;
    max_iter : int option;
    n_inducing : int;
    sigma2 : float;
    amplitude : float;
    dim_red : int option;
    het_sked : bool;
    multiscale : bool;
    tol : float;
    step : float;
    eps : float;
    verbose : bool;
  }

  let cmd : cmd ref = ref `train
  let model_file = ref None
  let with_stddev = ref false
  let max_iter = ref None
  let n_inducing = ref 10
  let sigma2 = ref 1.
  let amplitude = ref 1.
  let dim_red = ref None
  let het_sked = ref false
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
            | "train" -> cmd := `train
            | "test" -> cmd := `test
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
          "-max-iter",
          Arg.Int (set_some max_iter),
          " make predictions with 95%-confidence interval (default: limitless)"
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
          "-het-sked",
          Arg.Set het_sked,
          " turns on learning of heteroskedastic noise"
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
      max_iter = !max_iter;
      n_inducing = !n_inducing;
      sigma2 = !sigma2;
      amplitude = !amplitude;
      dim_red = !dim_red;
      het_sked = !het_sked;
      multiscale = !multiscale;
      tol = !tol;
      step = !step;
      eps = !eps;
      verbose = !verbose;
    }
end

let read_samples () =
  let rex = Pcre.regexp "," in
  match try Some (read_line ()) with _ -> None with
  | None -> failwith "no data"
  | Some line ->
      let sample = Array.map float_of_string (Pcre.asplit ~rex line) in
      let d = Array.length sample in
      let rec loop samples =
        match try Some (read_line ()) with _ -> None with
        | Some line ->
            let sample = Pcre.asplit ~rex line in
            if Array.length sample <> d then
              failwith (
                sprintf "incompatible dimension of sample in line %d: %s"
                  (List.length samples + 1) line)
            else loop (Array.map float_of_string sample :: samples)
        | None -> Array.of_list (List.rev samples)
      in
      loop [sample]

open Lacaml.Impl.D

open Gpr

module GP = Fitc_gp.Make_deriv (Cov_se_fat.Deriv)
module FIC = GP.FIC.Eval

module Model = struct
  type t = {
    mean : float;
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
  let fill c0 sample =
    for r1 = 1 to d do inputs.{r1, c0 + 1} <- sample.(r1 - 1) done;
    targets.{c0 + 1} <- sample.(d)
  in
  Array.iteri fill samples;
  inputs, targets

let write_model model_file mean trained =
  let oc = open_out model_file in
  let model =
    let model = FIC.Trained.get_model trained in
    let kernel = FIC.Model.get_kernel model in
    let inducing = FIC.Model.get_inducing model in
    let inducing_points = FIC.Inducing.get_points inducing in
    let mean_predictor = FIC.Mean_predictor.calc_trained trained in
    let coeffs = FIC.Mean_predictor.get_coeffs mean_predictor in
    let co_variance_coeffs = FIC.Model.calc_co_variance_coeffs model in
    {
      Model.
      mean = mean;
      kernel = kernel;
      inducing_points = inducing_points;
      coeffs = coeffs;
      co_variance_coeffs = co_variance_coeffs;
    }
  in
  Marshal.to_channel oc model [];
  close_out oc

let train args =
  let
    {
      Args.
      model_file = model_file;
      n_inducing = n_inducing;
      sigma2 = sigma2;
      amplitude = amplitude;
      dim_red = dim_red;
      het_sked = het_sked;
      multiscale = multiscale;
      tol = tol;
      step = step;
      eps = epsabs;
      verbose = verbose;
    } = args
  in
  let inputs, targets = read_training_samples () in
  let big_dim = Mat.dim1 inputs in
  let mean = Vec.sum targets /. float (Vec.dim targets) in
  let targets = Vec.map (fun n -> n -. mean) targets in
  let n_inducing = min n_inducing (Vec.dim targets) in
  Random.self_init ();
  let trained =
    let params =
      let log_sf2 = 2. *. log amplitude in
      let d, tproj =
        match dim_red with
        | None -> big_dim, None
        | Some small_dim ->
            let small_dim = min big_dim small_dim in
            small_dim, Some (Mat.make big_dim small_dim (1. /. float big_dim))
      in
      let log_hetero_skedasticity =
        if het_sked then Some (Vec.make0 n_inducing)
        else None
      in
      let log_multiscales_m05 =
        if multiscale then Some (Mat.make0 d n_inducing)
        else None
      in
      Cov_se_fat.Params.create
        {
          Cov_se_fat.Params.
          d = d;
          log_sf2 = log_sf2;
          tproj = tproj;
          log_hetero_skedasticity = log_hetero_skedasticity;
          log_multiscales_m05 = log_multiscales_m05;
        }
    in
    let kernel = Cov_se_fat.Eval.Kernel.create params in
    let report_trained_model, report_gradient_norm =
      if verbose then
        Some (fun ~iter trained ->
          let le = FIC.Trained.calc_log_evidence trained in
          let rmse = FIC.Trained.calc_rmse trained in
          printf "iter %4d:  log evidence: %.5f  rmse: %.5f\n%!" iter le rmse),
        Some (fun ~iter norm ->
          printf "iter %4d:  |gradient| = %.5f\n%!" iter norm)
      else None, None
    in
    GP.FIC.Deriv.Optim.Gsl.train
      ?report_trained_model ?report_gradient_norm
      ~kernel ~sigma2 ~n_rand_inducing:n_inducing
      ~tol ~step ~epsabs ~inputs ~targets ()
  in
  write_model model_file mean trained

let read_test_samples () =
  let samples = read_samples () in
  let n = Array.length samples in
  let d = Array.length samples.(0) - 1 in
  let inputs = Mat.create d n in
  let fill c0 sample =
    for r1 = 1 to d do inputs.{r1, c0 + 1} <- sample.(r1 - 1) done
  in
  Array.iteri fill samples;
  inputs

let read_model model_file : Model.t =
  let ic = open_in model_file in
  let model = Marshal.from_channel ic in
  close_in ic;
  model

let test args =
  let { Args.model_file = model_file; with_stddev = with_stddev } = args in
  let
    {
      Model.
      mean = mean;
      kernel = kernel;
      inducing_points = inducing_points;
      coeffs = coeffs;
      co_variance_coeffs = co_variance_coeffs;
    } = read_model model_file
  in
  let inputs = read_test_samples () in
  let mean_predictor = FIC.Mean_predictor.calc inducing_points ~coeffs in
  let inducing = FIC.Inducing.calc kernel inducing_points in
  let inputs = FIC.Inputs.calc inducing inputs in
  let means = FIC.Means.get (FIC.Means.calc mean_predictor inputs) in
  let means = Vec.map (fun n -> n +. mean) means in
  if with_stddev then
    let co_variance_predictor =
      FIC.Co_variance_predictor.calc kernel inducing_points co_variance_coeffs
    in
    let variances =
      FIC.Variances.calc co_variance_predictor ~sigma2:0. inputs
    in
    let variances = FIC.Variances.get ~predictive:false variances in
    Vec.iteri (fun i mean -> printf "%f,%f\n" mean variances.{i}) means
  else
    Vec.iter (fun n -> printf "%f\n" n) means

let main () =
  let args = Args.get () in
  match args.Args.cmd with
  | `train -> train args
  | `test -> test args

let () = main ()
