(lang dune 1.10)
(name gpr)

(generate_opam_files true)

(source (github mmottl/gpr))
(license "LGPL-2.1+ with OCaml linking exception")
(homepage "https://mmottl.github.io/gpr")
(documentation "https://mmottl.github.io/gpr/api")

(maintainers "Markus Mottl <markus.mottl@gmail.com>")
(authors "Markus Mottl <markus.mottl@gmail.com>")

(package
  (name gpr)
  (synopsis "GPR - Library and Application for Gaussian Process Regression")
  (description "\
Gaussian process regression is a modern Bayesian approach to machine learning,
and GPR implements some of the latest advances in this field.")
  (depends
    (ocaml (>= 4.08))
    (dune (>= 1.10))
    base-threads
    (core (>= 0.9.1))
    (lacaml (>= 11.0.0))
    (gsl (>= 1.24.0))
  )
)
