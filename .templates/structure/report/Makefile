all: build/main.pdf

build/main.pdf: content/* lit.bib ../analysis/results/* ../analysis/analysis.py| build
			lualatex --output-directory=build main.tex

full: content/* ../analysis/results/* lit.bib ../analysis/analysis.py| build
		lualatex --output-directory=build main.tex
		biber build/main.bcf
		lualatex --output-directory=build main.tex

results: ../analysis/results/* 

../analysis/results/*: ../analysis/analysis.py | ../analysis/results
	python ../analysis/analysis.py	

build:
		mkdir -p build

clean:
		rm -rf build
.PHONY: clean all results
