(* File: cov_se_iso.mli

   OCaml-GPR - Gaussian Processes for OCaml

     Copyright (C) 2009-  Markus Mottl
     email: markus.mottl@gmail.com
     WWW:   http://www.ocaml.info

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*)

(** {6 Isotropic squared exponential covariance} *)

(** The covariance is defined as:

    [k(x, y) = sf^2 * exp(-1/2 * |1/ell*(x-y)|^2)]

    where [sf^2] is the amplitude, and [ell] is the length scale.
*)

open Lacaml.D

open Interfaces.Specs

module Params : sig type t = { log_ell : float; log_sf2 : float } end

type inducing_hyper = { ind : int; dim : int }

module Eval :
  Eval
    with type Kernel.params = Params.t
    with type Inducing.t = mat
    with type Input.t = vec
    with type Inputs.t = mat

module Deriv :
  Deriv
    with module Eval = Eval
    with type Hyper.t =
      [ `Log_ell | `Log_sf2 | `Inducing_hyper of inducing_hyper ]
