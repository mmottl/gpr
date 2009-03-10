open Lacaml.Impl.D

open Gpr
open Utils

let n_inputs = 1000
let n_inducing_inputs = 10
let noise_sigma = 1.5
let noise_sigma2 = noise_sigma *. noise_sigma

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
  let inducing_inputs = Mat.create 1 n_inducing_inputs in
  for i = 1 to n_inducing_inputs do
    inducing_inputs.{1, i} <-
      -5. +. float i *. 10. /. float (n_inducing_inputs + 1)
  done;
  training_inputs, training_targets, inducing_inputs

let training_inputs, training_targets, inducing_inputs = get_training ()
