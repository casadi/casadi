all: python

python:
	cd python && python alltests.py && cd..

octave:
	cd octave && octave test.m && cd ..

full: python octave examples tutorials
	
examples:
	python examples_py.py '../documentation/examples'
	python examples_cpp.py

tutorials:
	python examples_py.py '../documentation/tutorials/python/src'
