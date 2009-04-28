open Printf

open Lacaml.Impl.D

type t = { data : mat array; n : int }

let check_square (i, size) mat =
  let m = Mat.dim1 mat in
  let n = Mat.dim2 mat in
  if m = n then i + 1, size + n
  else
    failwith (
      sprintf
        "Block_diag.check_square: matrix at index %d not square: m = %d, n = %d"
        i m n)

let create mats =
  { data = mats; n = snd (Array.fold_left check_square (0, 0) mats) }

let copy t = { t with data = Array.map (fun mat -> lacpy mat) t.data }

let reraise_exc loc i exc =
  failwith (
    sprintf "Block_diag.%s: failed at index %d: %s"
      loc i (Printexc.to_string exc))

let potrf ?jitter t =
  let acti i mat =
    try potrf ?jitter mat with exc -> reraise_exc "potrf" i exc
  in
  Array.iteri acti t.data

let potri ?jitter ?factorize t =
  let acti i mat =
    try potri ?jitter ?factorize mat with exc -> reraise_exc "potri" i exc
  in
  Array.iteri acti t.data
