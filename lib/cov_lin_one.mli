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

(** {2 Covariance of linear functions with one hyperparameter} *)

(** The covariance is defined as:

    [k(x, y) = x*inv(P)*y + 1/t^2] [logtheta = log(t)]

    where P is a diagonal matrix containing [t^2] along the diagonal. *)

open Lacaml.D
open Interfaces.Specs

module Params : sig
  type t = { log_theta : float }
end

module Eval :
  Eval
    with type Kernel.params = Params.t
    with type Inducing.t = mat
    with type Input.t = vec
    with type Inputs.t = mat

module Deriv : Deriv with module Eval = Eval with type Hyper.t = [ `Log_theta ]
