all: unittests

full: unittests examples tutorials benchmarks

benchmarks:
	cd python && python complexity.py && cd ..

unittests: unittests_py unittests_oct
	
unittests_py:
	cd python && python alltests.py && cd ..

unittests_oct:
	python internal/test_oct.py octave

examples: examples_indoc examples_code

examples_indoc: examples_indoc_py examples_indoc_cpp examples_indoc_oct

examples_code: examples_code_py examples_code_cpp

examples_indoc_py:
	python internal/test_py.py '../documentation/examples'
	
examples_code_py:
	python internal/test_py.py '../examples/python' -skipfiles="c_code_generation.py"
	
examples_indoc_cpp:
	python internal/test_cpp.py '../documentation/examples'
	
examples_code_cpp:
	python internal/test_cppcmake.py '../examples' -skipfiles="test_liftopt.cpp casadi_error_handling.cpp"

examples_indoc_oct:
	python internal/test_oct.py '../documentation/examples'

tutorials:
	python internal/test_py.py '../documentation/tutorials/python/src'
