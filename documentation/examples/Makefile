PY=$(shell find -name "*.py"  -not -name "_temp.py" -not -name "*ctemplate*")
CPP=$(shell find -name "*.cc" -not -name "_temp.cc" -not -name "*ctemplate*")
EXE=$(CPP:.cc=.run)
LOG=$(CPP:.cc=.log)
PDF=$(PY:.py=.pdf)
BUILDDIR=../../build

CASADI_OBJECT=$(BUILDDIR)/swig/python/casadi/_casadi.so
all: python cpp cpprun

.PHONY = python cpp cpprun

python: $(PDF)
cpp: $(EXE)
cpprun: $(LOG)

clean: 
	rm -f _temp.* $(PDF) $(EXE)

$(PDF): %.pdf : %.py $(CASADI_OBJECT)
	cp $< _temp.py
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
	
	
ctemplate/compiler.sh ctemplate/linker.sh: $(BUILDDIR)/documentation/examples/ctemplate/CMakeFiles/ctemplate.dir/link.txt  $(BUILDDIR)/documentation/examples/ctemplate/CMakeFiles/ctemplate.dir/flags.make
	cd ctemplate && python ctemplate.py && chmod +x compiler.sh linker.sh && cd ..
