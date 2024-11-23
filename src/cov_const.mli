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

(** {2 Covariance of a constant function} *)

(** The covariance is defined as:

    [k(x, y) = 1/s^2] [logtheta = log(s)] *)

open Interfaces.Specs

module Params : sig
  type t = { log_theta : float }
end

module Eval :
  Eval
    with type Kernel.params = Params.t
    with type Inducing.t = int
    with type Input.t = unit
    with type Inputs.t = int

module Deriv : Deriv with module Eval = Eval with type Hyper.t = [ `Log_theta ]
