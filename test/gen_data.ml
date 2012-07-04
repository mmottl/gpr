(* File: gen_data.ml

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
open Lacaml.D

open Gpr

let n_inputs = 1000
let n_inducing = 10
let noise_sigma = 0.7
let noise_sigma2 = noise_sigma *. noise_sigma

let f ?(with_noise = false) x =
  let v = (sin (3. *. x)) /. x +. abs_float (x -. 3.) /. (x *. x +. 1.) in
  if with_noise then
    v +. Gsl.Randist.gaussian Utils.default_rng ~sigma:noise_sigma
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

let gen_write pp file obj =
  let oc = open_out (Filename.concat "test/data" file) in
  Format.fprintf (Format.formatter_of_out_channel oc) "%a@." pp obj;
  close_out oc

let write_float file = gen_write Format.pp_print_float file
let write_vec file = gen_write pp_vec file
let write_mat file = gen_write pp_mat file
