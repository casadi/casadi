default: html pdf

html: casadi-users_guide.tex
	latex2html -local_icons -info "" -bottom_navigation -no_math -show_section_numbers -scalable_fonts -split 3 $<
	
pdf: casadi-users_guide.tex
	pdflatex -interaction=nonstopmode $<
	pdflatex -interaction=nonstopmode $<
	rm *.aux *.log *.out *.toc
