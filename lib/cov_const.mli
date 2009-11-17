(* File: cov_const.mli

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

(** {6 Covariance of a constant function} *)

(** The covariance is defined as:

    [k(x, y) = 1/s^2]
    [logtheta = log(s)]
*)

open Interfaces.Specs

module Params : sig type t = { log_theta : float } end

module Eval :
  Eval
    with type Kernel.params = Params.t
    with type Inducing.t = int
    with type Input.t = unit
    with type Inputs.t = int

module Deriv :
  Deriv
    with module Eval = Eval
    with type Hyper.t = [ `Log_theta ]
