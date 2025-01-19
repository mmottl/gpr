(* OCaml-GPR - Gaussian Processes for OCaml

   Copyright © 2009- Markus Mottl <markus.mottl@gmail.com>

   This library is free software; you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the Free
   Software Foundation; either version 2.1 of the License, or (at your option)
   any later version.

   This library is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
   details.

   You should have received a copy of the GNU Lesser General Public License
   along with this library; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA *)

(** {2 Covariance of linear functions with Automatic Relevance Determination} *)

(** The covariance is defined as:

    [k(x, y) = x*inv(P)*y]

    where P is a diagonal matrix containing ARD parameters ell_1^2,...,ell_D^2,
    and D is the dimensionality of the input space. *)

open Lacaml.D
open Interfaces.Specs

module Params : sig
  type t = { log_ells : vec }
end

module Eval :
  Eval
    with type Kernel.params = Params.t
    with type Inducing.t = mat
    with type Input.t = vec
    with type Inputs.t = mat

module Deriv :
  Deriv with module Eval = Eval with type Hyper.t = [ `Log_ell of int ]
