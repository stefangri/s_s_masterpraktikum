all: build/main.pdf

build/main.pdf: content/* results lit.bib | build
			lualatex --output-directory=build main.tex

full: content/* results lit.bib | build
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
