all: dox

options:
	python extra/generate_options_doc.py
	
examples:
	python extra/fetch_examples_pdf.py

dox:
	doxygen Doxyfile.in
	cp style/casadi.png html/casadi.png

full: dox
	python extra/generate_options_doc.py
	python extra/fetch_examples_pdf.py
	doxygen Doxyfile.in

publish: 
	./publish.sh
	
clean: 
	rm -rf XML/ html/ 
