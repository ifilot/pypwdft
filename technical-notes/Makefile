LATEXMK ?= latexmk

SOURCE := pseudopotentials.tex
PDF := $(SOURCE:.tex=.pdf)

.PHONY: all clean

all: $(PDF)

$(PDF): $(SOURCE)
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(SOURCE)

clean:
	$(LATEXMK) -C $(SOURCE)
