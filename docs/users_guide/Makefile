# Minimal makefile for Sphinx documentation
#


# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = CasADi
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile bare_body bare_toc users_guide.pdf

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.DEFAULT_GOAL := users_guide.pdf

all: users_guide.pdf bare_toc bare_body

bare_body:
	make singlehtml SPHINXOPTS="-D html_theme=bare_body "
	sed -i 's/{{</{\&lbrace;</g' build/singlehtml/index.html
	sed -i 's/>}}/>}\&rbrace;/g' build/singlehtml/index.html
	echo "---\n---\n" > build/_index.html
	cat build/singlehtml/index.html >> build/_index.html

bare_toc:
	make singlehtml SPHINXOPTS="-D html_theme=bare_body -D 'html_theme_options.mode=toc'"
	sed -i 's/index\.html//g' build/singlehtml/index.html
	cat addendum.html build/singlehtml/index.html > build/docs-sidebar.html

users_guide.pdf:
	make latexpdf
	cp $(BUILDDIR)/latex/CasADi.pdf $@

