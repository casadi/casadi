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

bare_body:
	make singlehtml SPHINXOPTS="-D html_theme=bare_body "

bare_toc:
	make singlehtml SPHINXOPTS="-D html_theme=bare_body -D 'html_theme_options.mode=toc'"

users_guide.pdf:
	make latexpdf
	cp $(BUILDDIR)/latex/CasADi.pdf $@
