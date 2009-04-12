open Lacaml.Impl.D

type t = private { data : mat array; n : int }

val create : mat array -> t

val copy : t -> t

val potrf : ?jitter : float -> t -> unit

val potri : ?jitter : float -> ?factorize : bool -> t -> unit
