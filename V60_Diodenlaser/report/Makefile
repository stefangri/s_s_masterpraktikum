all: build/main.pdf

build/main.pdf: content/* lit.bib | build
			lualatex --output-directory=build main.tex

full: content/* lit.bib | build
		lualatex --output-directory=build main.tex
		biber build/main.bcf
		lualatex --output-directory=build main.tex


build:
		mkdir -p build

clean:
		rm -rf build
.PHONY: clean all results
