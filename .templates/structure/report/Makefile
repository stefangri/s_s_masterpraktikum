all: build/main.pdf

debug: content/* lit.bib build/header.fmt| build
	 lualatex --output-directory=build --draft -shell-escape "&./build/header main.tex"

build/main.pdf: content/*  lit.bib | build
			lualatex --output-directory=build main.tex

full: content/* lit.bib | build
		lualatex --output-directory=build main.tex
		biber build/main.bcf
		lualatex --output-directory=build main.tex

results: ../analysis/results 		

../analysis/results: ../analysis/analysis.py 
	python ../analysis/analysis.py	

build:
		mkdir -p build

clean:
		rm -rf build
.PHONY: clean all results
