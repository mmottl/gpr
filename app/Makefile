TARGETS = $(addsuffix .bc, ocaml_gpr)

.PHONY: all clean

all:
	@jbuilder build --dev $(TARGETS)

clean:
	@jbuilder clean
