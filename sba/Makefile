all: sba.pdf

R_OPTS=--no-save --no-restore --no-init-file --no-site-file

sba.pdf: sba.Rnw
	R ${R_OPTS} -e "library(knitr);knit('sba.Rnw')"
	latexmk sba.tex

clean:
	rm -rf auto cache figures/* sba.tex sba.aux sba.bbl sba.blg \
	sba.fdb_latexmk sba.fls sba.log sba.pdf sba.tex
