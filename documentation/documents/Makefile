TEX = $(wildcard *.tex)
PDF = $(TEX:.tex=.pdf)

all: ${PDF}

%.pdf : %.tex
	pdflatex -interaction=nonstopmode $<
	pdflatex -interaction=nonstopmode $<
