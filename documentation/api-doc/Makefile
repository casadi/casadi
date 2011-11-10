all: dox

options:
	python extra/generate_options_doc.py

tutorials:
	cd ../tutorials/python && make

examples:
	cd ../examples && make
	python extra/fetch_examples_pdf.py

dox:
	doxygen Doxyfile.in
	cp style/casadi.png html/casadi.png

swig:
	sh extra/doc2swig.sh

full: dox
	cd ../examples && make
	cd ../tutorials/python && python make_tutorials.py
	cd ../cheatsheet && make
	python extra/generate_options_doc.py
	python extra/fetch_examples_pdf.py
	sh extra/doc2swig.sh
	doxygen Doxyfile.in

publish: 
	./publish.sh
	
clean: 
	rm -rf XML/ html/ 
