NAME := gpr_manual

PDFLATEX := pdflatex
BIBTEX := bibtex
MAKEINDEX := makeindex

all: $(NAME).pdf

%.pdf: %.tex
	$(PDFLATEX) < /dev/null $<
	-$(BIBTEX) $(basename $<)
	-$(MAKEINDEX) $(basename $<)
	$(PDFLATEX) < /dev/null $<
	$(PDFLATEX) < /dev/null $<

.PHONY: clean
clean:
	@rm -f $(NAME).dvi $(NAME).log $(NAME).aux $(NAME).idx $(NAME).toc \
	       $(NAME).bbl $(NAME).blg $(NAME).pdf missfont.log $(NAME).fls
