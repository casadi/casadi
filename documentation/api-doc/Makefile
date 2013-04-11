REVISION=$(shell git rev-parse HEAD)


all: doxpublish

options:
	python extra/generate_options_doc.py

tutorials:
	cd ../tutorials/python && make

documents:
	cd ../documents && make

examples:
	cd ../examples && make
	python extra/fetch_examples_pdf.py

doxlocal:
	doxygen Doxygen.in
	cp style/casadi.png html/casadi.png
	
doxpublish:
	patch -u Doxyfile.in publish.patch -o - | sed 's/PROJECT_NUMBER         = /PROJECT_NUMBER         = $(REVISION)/' | doxygen -
	cp style/casadi.png html/casadi.png
	find -name *-example.html | xargs perl -n -i~ -w -e 'print unless /<div class="line"><span class="comment">#[^\!].*?\/span>\n?/;'
	find -name *-example.html | xargs perl -p -i~ -w -e 's/<span class="comment">#[^\!].*?\/span>\n?//g;'

swig:
	sh extra/doc2swig.sh

full: doxpublish
	cd ../examples && make
	cd ../documents && make
	cd ../tutorials/python && python make_tutorials.py
	#cd ../tutorials/cpp && make
	cd ../users_guide && make
	cd ../cheatsheet && make
	python extra/generate_options_doc.py
	python extra/fetch_examples_pdf.py
	sh extra/doc2swig.sh
	patch -u Doxyfile.in publish.patch -o -  | sed 's/PROJECT_NUMBER         = /PROJECT_NUMBER         = $(REVISION)/' | doxygen -
	find -name *-example.html | xargs perl -n -i~ -w -e 'print unless /<div class="line"><span class="comment">#[^\!].*?\/span>\n?/;'
	find -name *-example.html | xargs perl -p -i~ -w -e 's/<span class="comment">#[^\!].*?\/span>\n?//g;'

publish: 
	./publish.sh
	
clean: 
	rm -rf XML/ html/ 
