.PHONY: all clean doc lint

all:
	dune build @install

clean:
	dune clean

doc:
	dune build @doc

lint:
	@opam lint gpr.opam
