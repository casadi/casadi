all: doxpublish

options:
	python extra/generate_options_doc.py

tutorials:
	cd ../tutorials/python && make

examples:
	cd ../examples && make
	python extra/fetch_examples_pdf.py

doxlocal:
	doxygen Doxygen.in
	cp style/casadi.png html/casadi.png
	
doxpublish:
	patch -u Doxyfile.in publish.patch -o - | doxygen -
	cp style/casadi.png html/casadi.png

swig:
	sh extra/doc2swig.sh

full: doxpublish
	cd ../examples && make
	cd ../tutorials/python && python make_tutorials.py
	cd ../tutorials/cpp && make
	cd ../cheatsheet && make
	python extra/generate_options_doc.py
	python extra/fetch_examples_pdf.py
	sh extra/doc2swig.sh
	patch -u Doxyfile.in publish.patch -o - | doxygen -

publish: 
	./publish.sh
	
clean: 
	rm -rf XML/ html/ 
