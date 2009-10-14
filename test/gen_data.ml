open Lacaml.Impl.D

open Gpr
open Utils

let n_inputs = 1000
let n_inducing = 10
let noise_sigma = 0.5
let noise_sigma2 = noise_sigma *. noise_sigma

let f ?(with_noise = false) x =
  let v = (sin (3. *. x)) /. x +. abs_float (x -. 3.) /. (x *. x +. 1.) in
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

let training_inputs, training_targets = get_data ~with_noise:true n_inputs
