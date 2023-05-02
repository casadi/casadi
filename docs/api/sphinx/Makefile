# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = CasADi
SOURCEDIR     = source
BUILDDIR      = build

EXPRESSIONS != grep DECL ../../../swig/casadi.i | grep -o -P "(?<=casadi_)\w*?(?=\()" | sort -u | paste -sd "," -

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile bare_body bare_toc

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


bare_body:
	sed 's/@@EXPRESSIONS@@/  :members: $(EXPRESSIONS)/g' source/api.rst.in > source/api.rst
	rm -rf build
	make singlehtml SPHINXOPTS="-D html_theme=bare_body "
	sed -i 's/{{</{\&lbrace;</g' build/singlehtml/index.html
	sed -i 's/>}}/>}\&rbrace;/g' build/singlehtml/index.html
	mkdir -p ../../../web/content/python-api
	echo "---\n---\n" > ../../../web/content/python-api/_index.html
	cat build/singlehtml/index.html >> ../../../web/content/python-api/_index.html

bare_toc:
	sed 's/@@EXPRESSIONS@@//g' source/api.rst.in > source/api.rst
	make singlehtml SPHINXOPTS="-D html_theme=bare_body -D 'html_theme_options.mode=toc'"
	sed -i 's/index\.html//g' build/singlehtml/index.html
	mkdir -p ../../../web/themes/casadi-theme/layouts/partials/
	cat build/singlehtml/index.html > ../../../web/themes/casadi-theme/layouts/partials/python-api-sidebar.html
