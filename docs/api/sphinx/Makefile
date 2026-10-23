# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = CasADi
SOURCEDIR     = source
BUILDDIR      = build

# Module-level functions with a docstring, except the solver factories documented in their own sections
EXPRESSIONS != python -c "import casadi, inspect; print(','.join(sorted(k for k, v in vars(casadi).items() if inspect.isfunction(v) and not k.startswith('_') and v.__doc__ and k not in ('nlpsol', 'qpsol', 'conic', 'integrator', 'rootfinder'))))"

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile bare_body bare_toc

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


bare_body:
	python ../extra/doctest_sidecar.py ../../../casadi source/_doctest_sidecar.json
	sed 's/@@EXPRESSIONS@@/  :members: $(EXPRESSIONS)/g' source/api.rst.in > source/api.rst
	rm -rf build
	make singlehtml SPHINXOPTS="-D html_theme=bare_body "
	sed -i 's/{{</{\&lbrace;</g' build/singlehtml/index.html
	sed -i 's/>}}/>}\&rbrace;/g' build/singlehtml/index.html
	mkdir -p ../../../web/content/python-api
	echo "---\n---\n" > ../../../web/content/python-api/_index.html
	cat build/singlehtml/index.html >> ../../../web/content/python-api/_index.html

bare_toc:
	python ../extra/doctest_sidecar.py ../../../casadi source/_doctest_sidecar.json
	sed 's/@@EXPRESSIONS@@//g' source/api.rst.in > source/api.rst
	make singlehtml SPHINXOPTS="-D html_theme=bare_body -D 'html_theme_options.mode=toc'"
	sed -i 's/index\.html//g' build/singlehtml/index.html
	mkdir -p ../../../web/themes/casadi-theme/layouts/partials/
	cat build/singlehtml/index.html > ../../../web/themes/casadi-theme/layouts/partials/python-api-sidebar.html
