PYTHON=python
ifdef WITH_PYTHON3
PYTHON=python3
endif
REVISION=$(shell $(PYTHON) -c "import casadi;print(casadi.__version__)")

default: html pdf

html: pdf
	latex2html -local_icons -info "" -init_file customlatex2html.perl -bottom_navigation -no_math -show_section_numbers -scalable_fonts -split 3 casadi-users_guide.tex
	find -name \*.html -exec sed -i "s/currentversionplaceholder/$(REVISION)/g" {} \;

snippets:
	rm -f pytex_*
	pdflatex -interaction=nonstopmode "\newcommand{\currentversion}{$(REVISION)\xspace} \input{casadi-users_guide.tex}" -j    obname=casadi-users_guide

pdf: casadi-users_guide.tex
	rm -f pytex_*
	pdflatex -interaction=nonstopmode "\newcommand{\currentversion}{$(REVISION)\xspace} \input{casadi-users_guide.tex}" -jobname=casadi-users_guide
	python pytex.py $(PYTHON)
	pdflatex -interaction=nonstopmode "\newcommand{\currentversion}{$(REVISION)\xspace} \input{casadi-users_guide.tex}" -jobname=casadi-users_guide
	rm -f *.out *.toc #pytex_*
