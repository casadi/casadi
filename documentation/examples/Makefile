PY=$(shell find -name "*.py"  -not -name "_temp.py" -not -name "*ctemplate*")
CPP=$(shell find -name "*.cc" -not -name "_temp.cc" -not -name "*ctemplate*")
OCT=$(shell find -name "*.m" -not -name "_temp.m" -not -name "*ctemplate*")
EXE=$(CPP:.cc=.run)
LOG=$(CPP:.cc=.log)
OCTLOG=$(OCT:.m=.out)
PDF=$(PY:.py=.pdf)
BUILDDIR=../../build

CASADI_OBJECT=$(BUILDDIR)/python/casadi/_casadi.so
OCTAVE_PATH=$(BUILDDIR)/lib
all: python cpp cpprun octave

.PHONY = python cpp cpprun

python: $(PDF)
cpp: $(EXE)
cpprun: $(LOG)
octave: $(OCTLOG)

clean: 
	rm -f _temp.* $(PDF) $(EXE)

$(PDF): %.pdf : %.py $(CASADI_OBJECT)
	rm -rf pylab.py pylab.pyc matplotlibrc
	grep -v "^#[^\!]" $< > _temp.py
	pyreport -d -l _temp.py -o _temp.pdf
	rm _temp.py
	mv _temp.pdf $@
	
$(EXE): %.run : %.cc $(BUILDDIR)/bin/ctemplate ctemplate/compiler.sh ctemplate/linker.sh
	cp $< _temp.cc
	ctemplate/compiler.sh _temp
	ctemplate/linker.sh _temp
	rm _temp.cc _temp.o
	mv _temp.run $@
	
$(LOG): %.log : %.run
	$< > $@ 2>&1
	
$(OCTLOG): %.out : %.m
	octave --no-init-file -p $(OCTAVE_PATH) -q $< 2>&1 > $@ 
	
ctemplate/compiler.sh ctemplate/linker.sh: $(BUILDDIR)/documentation/examples/ctemplate/CMakeFiles/ctemplate.dir/link.txt  $(BUILDDIR)/documentation/examples/ctemplate/CMakeFiles/ctemplate.dir/flags.make ctemplate/ctemplate.py
	cd ctemplate && python ctemplate.py && chmod +x compiler.sh linker.sh && cd ..
