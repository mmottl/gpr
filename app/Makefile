TARGETS = $(addsuffix .bc, ocaml_gpr)

.PHONY: all clean

all:
	@echo jbuilder: Entering directory \`$(abspath $(dir $(lastword $(PWD))))\'
	@jbuilder build --dev $(TARGETS)

clean:
	@jbuilder clean
