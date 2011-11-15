default: html pdf

html: casadi-doc.tex
	latex2html -local_icons -info "" -bottom_navigation -no_math -show_section_numbers -scalable_fonts -split 3 $<
	
pdf: casadi-doc.tex
	pdflatex -interaction=nonstopmode $<
	pdflatex -interaction=nonstopmode $<
	rm *.aux *.log *.out *.toc
