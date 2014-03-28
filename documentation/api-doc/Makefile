REVISION=$(shell git rev-parse HEAD)
READABLEREVISION=$(shell python -c "import casadi;print casadi.__version__")

all: doxpublish

options:
	python extra/generate_options_doc.py

tutorials:
	cd ../tutorials/python && make

example_pack:
	cd ../example_pack && make

documents:
	cd ../documents && make

examples:
	cd ../examples && make
	python extra/fetch_examples_pdf.py

doxlocal:
	doxygen Doxygen.in
	cp style/casadi.png html/casadi.png
	
doxpublish:
	patch -u Doxyfile.in publish.patch -o - | sed 's/PROJECT_NUMBER         = /PROJECT_NUMBER         = $(REVISION)/' |    sed 's/PROJECT_BRIEF          = /PROJECT_BRIEF          = $(READABLEREVISION)/' | doxygen -
	sed -i.bak "s/q:query/q:query,db:dbLocation/g" html/search/search.js 
	doxyindexer searchdata.xml -o html
	cp style/casadi.png html/casadi.png
	patch -u Doxyfile.in publish.patch -o - | sed 's/PROJECT_NUMBER         = /PROJECT_NUMBER         = $(REVISION)/' |    sed 's/PROJECT_BRIEF          = /PROJECT_BRIEF          = $(READABLEREVISION)/' |    sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = INTERNAL/' |   sed 's/HTML_OUTPUT            = html/HTML_OUTPUT            = internal/' | sed 's/COLLABORATION_GRAPH    = NO/COLLABORATION_GRAPH    = YES/' | doxygen -
	sed -i.bak "s/q:query/q:query,db:dbLocation/g" internal/search/search.js 
	doxyindexer searchdata.xml -o internal
	cp style/casadi.png internal/casadi.png
	find -name *-example.html | xargs perl -n -i~ -w -e 'print unless /<div class="line"><span class="comment">#[^\!].*?\/span>\n?/;'
	find -name *-example.html | xargs perl -p -i~ -w -e 's/<span class="comment">#[^\!].*?\/span>\n?//g;'

swig:
	sh extra/doc2swig.sh

full: 
	cd ../examples && make
	cd ../documents && make
	cd ../tutorials/python && python make_tutorials.py
	#cd ../tutorials/cpp && make
	cd ../users_guide && make
	cd ../cheatsheet && make
	patch -u Doxyfile.in publish.patch -o - | sed 's/PROJECT_NUMBER         = /PROJECT_NUMBER         = $(REVISION)/' |    sed 's/PROJECT_BRIEF          = /PROJECT_BRIEF          = $(READABLEREVISION)/' |    sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = INTERNAL/' |   sed 's/HTML_OUTPUT            = html/HTML_OUTPUT            = internal/' | doxygen -
	python extra/generate_options_doc.py
	python extra/fetch_examples_pdf.py
	sh extra/doc2swig.sh
	make doxpublish
	cd ../example_pack && make

clean: 
	rm -rf XML/ html/ 
