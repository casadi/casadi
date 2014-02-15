REVISION=$(shell git rev-parse HEAD)

default: html pdf

html: pdf
	latex2html -local_icons -info "" -init_file customlatex2html.perl -bottom_navigation -no_math -show_section_numbers -scalable_fonts -split 3 casadi-users_guide.tex
	
snippets:
	rm -f pytex_*
	pdflatex -interaction=nonstopmode casadi-users_guide.tex
	
pdf: casadi-users_guide.tex
	rm -f pytex_*
	pdflatex -interaction=nonstopmode $<
	python pytex.py
	pdflatex -interaction=nonstopmode $<
	rm -f *.aux *.out *.toc #pytex_*
