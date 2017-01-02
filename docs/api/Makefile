REVISION=$(shell git rev-parse HEAD)
READABLEREVISION=$(shell python -c "import casadi;print casadi.__version__")

.PHONY: export

all: doxpublish

options:
	python extra/generate_options_doc.py

example_pack:
	cd ../example_pack && make

documents:
	cd ../documents && make

examples:
	cd ../examples && make
	python extra/fetch_examples_pdf.py

doxlocal:
	doxygen Doxyfile.in
	cp style/casadi.png html/casadi.png
	
doxpublish:
	rm -rf XML XML_internal XML_clutter
	patch -u Doxyfile.in publish.patch -o - | sed 's/PROJECT_NUMBER         = /PROJECT_NUMBER         = $(REVISION)/' |    sed 's/PROJECT_BRIEF          = /PROJECT_BRIEF          = $(READABLEREVISION)/' | doxygen -
	sed -i.bak "s/q:query/q:query,db:dbLocation/g" html/search/search.js
	find html/ -name "*html" -type f -exec sed -i -e 's/CASADI_[A-Z_]*_EXPORT//g' {} \;
	doxyindexer searchdata.xml -o html
	cp style/casadi.png html/casadi.png
	python extra/generate_overview.py
	perl -n -i~ -w -e 'if (/<div class="contents">/) { print qq|<div class="contents"> <img src="overview.png" border="0" alt="" usemap="#overview"/>|; open(IN, "<overview.map"); print while(<IN>); } else {print;}' html/index.html
	cp overview.png html/overview.png
	patch -u Doxyfile.in publish.patch -o - | sed 's/PROJECT_NUMBER         = /PROJECT_NUMBER         = $(REVISION)/' |    sed 's/PROJECT_BRIEF          = /PROJECT_BRIEF          = $(READABLEREVISION)/' |    sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = INTERNAL DEPRECATED CLUTTER UNSAFE SWIGINTERNAL/' |   sed 's/HTML_OUTPUT            = html/HTML_OUTPUT            = internal/' | sed 's/COLLABORATION_GRAPH    = NO/COLLABORATION_GRAPH    = YES/' | sed 's/INLINE_SOURCES         = NO/INLINE_SOURCES         = YES/' | sed 's/REFERENCED_BY_RELATION = NO/REFERENCED_BY_RELATION = YES/' | sed 's/REFERENCES_RELATION    = NO/REFERENCES_RELATION    = YES/'| sed 's/                         \*\.dox/                         *.dox *.cpp/' | sed 's/style\/header.html/style\/header_internal.html/'| sed 's/XML_OUTPUT             = XML/XML_OUTPUT             = XML_internal/'  | sed 's/PREDEFINED+=             SWIG//' | doxygen -
	sed -i.bak "s/q:query/q:query,db:dbLocation/g" internal/search/search.js 
	doxyindexer searchdata.xml -o internal
	cp style/casadi.png internal/casadi.png
	python extra/generate_overview.py
	perl -n -i~ -w -e 'if (/<div class="contents">/) { print qq|<div class="contents"> <img src="overview.png" border="0" alt="" usemap="#overview"/>|; open(IN, "<overview.map"); print while(<IN>); } else {print;}' internal/index.html
	cp overview.png internal/overview.png
	find -name *-example.html | xargs perl -n -i~ -w -e 'print unless /<div class="line"><span class="comment">#[^\!].*?\/span>\n?/;'
	find -name *-example.html | xargs perl -p -i~ -w -e 's/<span class="comment">#[^\!].*?\/span>\n?//g;'

doxpublishlocal:
	rm -rf XML XML_internal XML_clutter
	cat Doxyfile.in | sed 's/PROJECT_NUMBER         = /PROJECT_NUMBER         = $(REVISION)/' |    sed 's/PROJECT_BRIEF          = /PROJECT_BRIEF          = $(READABLEREVISION)/' | doxygen -
	find html/ -name "*html" -type f -exec sed -i -e 's/CASADI_[A-Z_]*_EXPORT//g' {} \;
	cp style/casadi.png html/casadi.png
	python extra/generate_overview.py
	perl -n -i~ -w -e 'if (/<div class="contents">/) { print qq|<div class="contents"> <img src="overview.png" border="0" alt="" usemap="#overview"/>|; open(IN, "<overview.map"); print while(<IN>); } else {print;}' html/index.html
	cp overview.png html/overview.png
	cat Doxyfile.in | sed 's/PROJECT_NUMBER         = /PROJECT_NUMBER         = $(REVISION)/' |    sed 's/PROJECT_BRIEF          = /PROJECT_BRIEF          = $(READABLEREVISION)/' |    sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = INTERNAL DEPRECATED CLUTTER UNSAFE SWIGINTERNAL/' |   sed 's/HTML_OUTPUT            = html/HTML_OUTPUT            = internal/' | sed 's/COLLABORATION_GRAPH    = NO/COLLABORATION_GRAPH    = YES/' | sed 's/INLINE_SOURCES         = NO/INLINE_SOURCES         = YES/' | sed 's/REFERENCED_BY_RELATION = NO/REFERENCED_BY_RELATION = YES/' | sed 's/REFERENCES_RELATION    = NO/REFERENCES_RELATION    = YES/'| sed 's/                         \*\.dox/                         *.dox *.cpp/' | sed 's/style\/header.html/style\/header_internal.html/'| sed 's/XML_OUTPUT             = XML/XML_OUTPUT             = XML_internal/'  | sed 's/PREDEFINED+=             SWIG//' | doxygen -
	cp style/casadi.png internal/casadi.png
	python extra/generate_overview.py
	perl -n -i~ -w -e 'if (/<div class="contents">/) { print qq|<div class="contents"> <img src="overview.png" border="0" alt="" usemap="#overview"/>|; open(IN, "<overview.map"); print while(<IN>); } else {print;}' internal/index.html
	cp overview.png internal/overview.png
	find -name *-example.html | xargs perl -n -i~ -w -e 'print unless /<div class="line"><span class="comment">#[^\!].*?\/span>\n?/;'
	find -name *-example.html | xargs perl -p -i~ -w -e 's/<span class="comment">#[^\!].*?\/span>\n?//g;'

swig:
	sh extra/doc2swig.sh

doxtest:
	patch -u Doxyfile.in publish.patch -o - |  sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = INTERNAL DEPRECATED CLUTTER UNSAFE SWIGINTERNAL/' | doxygen - 2>&1 | if grep -e "File ended in the middle of a comment block" -e "Found unknown command" -e "warning: expected whitespace" -e "warning: Conditional section"; then exit 1; else exit 0;fi


doxexport:
	rm -rf XML XML_internal XML_clutter
	patch -u Doxyfile.in publish.patch -o - |  sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = INTERNAL DEPRECATED CLUTTER UNSAFE SWIGINTERNAL/' |   sed 's/GENERATE_HTML          = YES/GENERATE_HTML          = NO/' | sed 's/                         \*\.dox/                         *.dox *.cpp/' | sed 's/XML_OUTPUT             = XML/XML_OUTPUT             = XML_internal/'  | sed 's/PREDEFINED+=             SWIG//' | doxygen -
	python extra/generate_options_doc.py
	rm -rf XML XML_internal XML_clutter
	patch -u Doxyfile.in publish.patch -o - |  sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = DEPRECATED CLUTTER UNSAFE SWIGINTERNAL/' |   sed 's/GENERATE_HTML          = YES/GENERATE_HTML          = NO/' | sed 's/                         \*\.dox/                         *.dox *.cpp/' | sed 's/XML_OUTPUT             = XML/XML_OUTPUT             = XML_cluttered /' | doxygen -
	patch -u Doxyfile.in publish.patch -o - |  sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = INTERNAL DEPRECATED CLUTTER UNSAFE SWIGINTERNAL/' |   sed 's/GENERATE_HTML          = YES/GENERATE_HTML          = NO/' | sed 's/                         \*\.dox/                         *.dox *.cpp/' | sed 's/XML_OUTPUT             = XML/XML_OUTPUT             = XML_internal/'  | sed 's/PREDEFINED+=             SWIG//' | doxygen -
	sh extra/doc2swig.sh

doxexport_old:
	rm -rf XML XML_internal XML_clutter
	patch -u Doxyfile.in publish.patch -o - |  sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = INTERNAL DEPRECATED CLUTTER UNSAFE SWIGINTERNAL/' |   sed 's/GENERATE_HTML          = YES/GENERATE_HTML          = NO/' | sed 's/                         \*\.dox/                         *.dox *.cpp/' | sed 's/XML_OUTPUT             = XML/XML_OUTPUT             = XML_internal/' | sed 's/PREDEFINED+=             SWIG//' | doxygen -
	python extra/generate_options_doc.py
	rm -rf XML XML_internal XML_clutter
	patch -u Doxyfile.in publish.patch -o - |  sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = DEPRECATED CLUTTER UNSAFE SWIGINTERNAL EXPANDED/' |   sed 's/GENERATE_HTML          = YES/GENERATE_HTML          = NO/' | sed 's/XML_OUTPUT             = XML/XML_OUTPUT             = XML_cluttered /' | sed 's/PREDEFINED+=             WITH_DEPRECATED_FEATURES//' | sed 's/EXPAND_ONLY_PREDEF     = YES/EXPAND_ONLY_PREDEF     = NO/' | doxygen -
	patch -u Doxyfile.in publish.patch -o - |  sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = INTERNAL DEPRECATED CLUTTER UNSAFE SWIGINTERNAL EXPANDED/' |   sed 's/GENERATE_HTML          = YES/GENERATE_HTML          = NO/' | sed 's/XML_OUTPUT             = XML/XML_OUTPUT             = XML_internal/' |  sed 's/EXPAND_ONLY_PREDEF     = YES/EXPAND_ONLY_PREDEF     = NO/' | doxygen -
	sh extra/doc2swig.sh

full: 
	cd examples && make
	cd ../documents && make
	cd ../users_guide && make
	cd ../cheatsheet && make
	patch -u Doxyfile.in publish.patch -o - |  sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = INTERNAL DEPRECATED CLUTTER UNSAFE SWIGINTERNAL/' |   sed 's/GENERATE_HTML          = YES/GENERATE_HTML          = NO/' | sed 's/                         \*\.dox/                         *.dox *.cpp/' | sed 's/XML_OUTPUT             = XML/XML_OUTPUT             = XML_internal/'  | sed 's/PREDEFINED+=             SWIG//' | doxygen -
	python extra/fetch_examples_pdf.py
	make doxpublish
	cd ../example_pack && make

fulllocal: 
	cd examples && make
	cd ../documents && make
	cd ../users_guide && make
	cd ../cheatsheet && make
	patch -u Doxyfile.in publish.patch -o - |  sed 's/ENABLED_SECTIONS       = /ENABLED_SECTIONS       = INTERNAL DEPRECATED CLUTTER UNSAFE SWIGINTERNAL/' |   sed 's/GENERATE_HTML          = YES/GENERATE_HTML          = NO/' | sed 's/                         \*\.dox/                         *.dox *.cpp/' | sed 's/XML_OUTPUT             = XML/XML_OUTPUT             = XML_internal/'  | sed 's/PREDEFINED+=             SWIG//' | doxygen -
	python extra/fetch_examples_pdf.py
	make doxpublishlocal
	cd ../example_pack && make

export:
	mkdir -p export/api export/users_guide export/cheatsheets/ export/documents/
	cp -r html/ export/api/
	cp -r internal/ export/api/
	cp  ../users_guide/*.pdf export/users_guide/
	cp -r ../users_guide/casadi-users_guide export/users_guide/html/
	cp  ../cheatsheet/*.pdf export/cheatsheets/
	cp  ../documents/*.pdf export/documents/

clean: 
	rm -rf XML/ html/ 

