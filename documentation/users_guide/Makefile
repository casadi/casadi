SVNREVISION=$(shell svn info | grep Revision | cut -d " " -f 2)$(shell git svn info | grep Revision | cut -d " " -f 2)

default: html pdf

html: pdf
	latex2html -local_icons -info "" -bottom_navigation -no_math -show_section_numbers -scalable_fonts -split 3 casadi-users_guide.tex
	
	
pdf: casadi-users_guide.tex
	pdflatex -interaction=nonstopmode $<
	python pytex.py
	pdflatex -interaction=nonstopmode $<
	rm *.aux *.out *.toc #pytex_*
