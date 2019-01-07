SHELL=/bin/bash
filename=jpaper
OPTS = --shell-escape

.PHONY: clean clean_inkscape clean figures clean_all pdf


pdf:${filename}.pdf

${filename}.pdf:${filename}.tex
	@echo "Building PDF from"
	pdflatex ${OPTS} ${^:.tex=}
	bibtex $^||true
	pdflatex ${OPTS} ${^:.tex=}
	pdflatex ${OPTS} ${^:.tex=}


clean:
	@rm -f ${filename}.{ps,pdf,log,aux,out,dvi,bbl,blg}
	@rm -f tmp1.qdf tmp2.qdf
	@rm -rf auto/
	@rm -rf ltximg/
	@rm -rf jpaper.prv
	@rm -f *.log

clean_inkscape:
	@rm -rf figures/svg-inkscape/*

clean_figures:
	@rm -f figures/*.svg

clean_all:clean clean_inkscape clean_figures
	@echo "Cleaning all figures"
