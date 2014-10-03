(* File: block_diag.mli

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

open Lacaml.D

(** Type of block diagonal matrices *)
type t = private { data : mat array; n : int }

(** [create mats] @return a block diagonal matrix whose block elements are made
    of the matrices in [mats]. *)
val create : mat array -> t

(** [copy bm] @return a copy of block diagonal matrix [bm]. *)
val copy : t -> t

(** [potrf ?jitter bm] perform Cholesky factorization on block diagonal matrix
    [bm] using Cholesky [jitter] if given.

    @param jitter default = no jitter
*)
val potrf : ?jitter : float -> t -> unit

(** [potri ?jitter ?factorize bm] invert block diagonal matrix [bm] using
    its Cholesky factor.  If [factorize] is [false], it is assumed that the
    factorization has already been performed, otherwise it will be calculated
    using Cholesky [jitter] if given.

    @param jitter default = no jitter
    @param factorize default = [true]
*)
val potri : ?jitter : float -> ?factorize : bool -> t -> unit
