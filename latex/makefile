SHELL=/bin/bash
filename=jpaper
OPTS = --shell-escape

pdf:
	echo "in PDF-latex"
	pdflatex ${OPTS} ${filename}.tex
	bibtex ${filename}||true
	pdflatex ${OPTS} ${filename}.tex
	pdflatex ${OPTS} ${filename}.tex


clean:
	rm -f ${filename}.{ps,pdf,log,aux,out,dvi,bbl,blg}
	rm -f tmp1.qdf tmp2.qdf
	rm -rf auto/
	rm -rf ltximg/
	rm -rf jpaper.prv
	rm -f *.log
clean_inkscape:
	rm -f -r figures/svg-inkscape/*

clean_figures:
	rm -f figures/drift_fit.svg

clean_all:clean clean_inkscape clean_figures
	echo "stuff"
