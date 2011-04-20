PY=$(shell find -name *.py)
PDF=$(PY:.py=.pdf)
CASADI_OBJECT=../../../build/swig/python/casadi/_casadi.so
all: $(PDF)

clean: 
	rm -f py.temp $(PDF)

$(PDF): %.pdf : %.py $(CASADI_OBJECT)
	cp $< py.temp
	pyreport -l py.temp -o temp.pdf
	rm py.temp
	mv temp.pdf $@

