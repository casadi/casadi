PY=$(shell find -name "*.py"  -not -name "_temp.py" -not -name "*ctemplate*")
CPP=$(shell find -name "*.cc" -not -name "_temp.cc" -not -name "*ctemplate*")
EXE=$(CPP:.cc=.run)
LOG=$(CPP:.cc=.log)
PDF=$(PY:.py=.pdf)
BUILDDIR=../../../build

OCTAVE_PATH=$(BUILDDIR)/lib
all: cpp cpprun python

.PHONY = python cpp cpprun

python: $(PDF)
cpp: $(EXE)
cpprun: $(LOG)

clean: 
	rm -f _temp.* $(PDF) $(EXE)

$(PDF): %.pdf : %.py
	rm -rf pylab.py pylab.pyc matplotlibrc
	jupytext --to notebook $<
	jupyter nbconvert --to pdf --execute $(<:.py=.ipynb)
	#rm $(<:.py=.ipynb)
	
$(EXE): %.run : %.cc ctemplate/compiler.sh
	ctemplate/compiler.sh $< $@
	
$(LOG): %.log : %.run
	$< > $@ 2>&1
	
ctemplate/compiler.sh: ctemplate/ctemplate.py
	cd ctemplate && python ctemplate.py && chmod +x compiler.sh && cd ..
