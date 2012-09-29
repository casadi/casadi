all: unittests





# Each of these targets will be tested by CasADi's trunktesterbot as individual tests.
# These targets must all go on one line
trunktesterbot: unittests_py unittests_oct examples_indoc_py examples_indoc_oct examples_indoc_cpp tutorials examples_code_py examples_code_cpp users_guide user_guide_snippets_py

# Set to -memcheck if you want to include a check for memory leaks
MEMCHECK = 

full: unittests examples tutorials benchmarks

benchmarks:
	cd python && python complexity.py && cd ..

unittests: unittests_py unittests_oct
	
unittests_py:
	python internal/test_py.py python -skipfiles="alltests.py helpers.py complexity.py speed.py" $(MEMCHECK)

unittests_oct:
	python internal/test_oct.py octave

examples: examples_indoc examples_code

examples_indoc: examples_indoc_py examples_indoc_cpp examples_indoc_oct

examples_code: examples_code_py examples_code_cpp

# No memchecks because there is a memory leak in SuperLU
examples_indoc_py:
	python internal/test_py.py '../documentation/examples'
	
# No memchecks because it leaks like hell, and is slow
examples_code_py:
	python internal/test_py.py '../examples/python' -skipfiles="c_code_generation.py parametric_sensitivities.py reduced_hessian.py"
	
examples_indoc_cpp:
	python internal/test_cpp.py '../documentation/examples' $(MEMCHECK)
	
# No memchecks because it's slow
examples_code_cpp:
	python internal/test_cppcmake.py '../examples/cplusplus' -skipfiles="test_liftopt.cpp casadi_error_handling.cpp parametric_sensitivities.cpp reduced_hessian.cpp"

examples_indoc_oct:
	python internal/test_oct.py '../documentation/examples'

tutorials:
	python internal/test_py.py '../documentation/tutorials/python/src' $(MEMCHECK)

users_guide:
	cd '../documentation/users_guide/' && make snippets && cd ..

user_guide_snippets_py: users_guide
	python internal/test_py.py '../documentation/users_guide' -skipfiles="pytex.py temp_pytex.py"
