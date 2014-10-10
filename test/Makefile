all: unittests





# Each of these targets will be tested by CasADi's trunktesterbot as individual tests.
# These targets must all go on one line
trunktesterbot: unittests_py examples_indoc_py examples_indoc_cpp tutorials examples_code_py examples_code_cpp users_guide user_guide_snippets_py lint

trunktesterbot_knownbugs: unittests_py_knownbugs examples_indoc_py examples_indoc_cpp tutorials examples_code_py examples_code_cpp users_guide user_guide_snippets_py

# Set to -memcheck if you want to include a check for memory leaks
MEMCHECK = 

lint :
	cd ../build && make lint; cd ../test

full: unittests examples tutorials benchmarks

full_knownbugs: unittests_knownbugs examples tutorials benchmarks

benchmarks:
	cd python && python complexity.py; cd ..

python: unittests_py examples_indoc_py examples_code_py user_guide_snippets_py tutorials

python_knownbugs: unittests_py_knownbugs examples_indoc_py examples_code_py user_guide_snippets_py tutorials

cpp: examples_indoc_cpp  examples_code_cpp

unittests: unittests_py

unittests_knownbugs: unittests_py_knownbugs

unittests_py:
ifdef MEMCHECK
	python internal/test_py.py python -skipfiles="alltests.py helpers.py complexity.py speed.py" -passoptions="--ignore_memory_light --run_slow"
	python internal/test_py.py python -skipfiles="alltests.py helpers.py complexity.py speed.py" $(MEMCHECK) -passoptions="--ignore_memory_heavy --run_slow"
else
	python internal/test_py.py python -skipfiles="alltests.py helpers.py complexity.py speed.py" # No memory checks
endif

unittests_py_knownbugs:
	python internal/test_py.py python -skipfiles="alltests.py helpers.py complexity.py speed.py" -passoptions="--known_bugs"

examples: examples_indoc examples_code

examples_indoc: examples_indoc_py examples_indoc_cpp

examples_code: examples_code_py examples_code_cpp

unittests_matlab:
	python internal/test_matlab.py matlab

# No memchecks because there is a memory leak in SuperLU
examples_indoc_py:
	python internal/test_py.py '../docs/api/examples'

# No memchecks because it leaks like hell, and is slow
examples_code_py:
	python internal/test_py.py '../docs/examples/python' -skipfiles="c_code_generation.py parametric_sensitivities.py reduced_hessian.py"

examples_indoc_cpp:
	python internal/test_cpp.py '../docs/api/examples' $(MEMCHECK)

# No memchecks because it's slow
examples_code_cpp:
	python internal/test_cppcmake.py '../docs/examples/cplusplus' -skipfiles="test_liftopt.cpp casadi_error_handling.cpp parametric_sensitivities.cpp reduced_hessian.cpp test_opencl.cpp nlp_codegen.cpp codegen_usage.cpp stabilized_sqpmethod_nl.cpp"

tutorials:
	python internal/test_py.py '../docs/tutorials/python/src' $(MEMCHECK)

users_guide:
	cd '../docs/users_guide/' && make snippets && cd ..

user_guide_snippets_py: users_guide
	python internal/test_py.py '../docs/users_guide' -skipfiles="pytex.py temp_pytex.py"
