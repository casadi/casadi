all: unittests





# Each of these targets will be tested by CasADi's trunktesterbot as individual tests.
# These targets must all go on one line
trunktesterbot: unittests_py examples_indoc_py examples_indoc_cpp examples_code_py examples_code_cpp users_guide user_guide_snippets_py

trunktesterbot_no_unittests_py: examples_indoc_py examples_indoc_cpp examples_code_py examples_code_cpp users_guide user_guide_snippets_py

trunktesterbot_knownbugs: unittests_py_knownbugs examples_indoc_py examples_indoc_cpp examples_code_py examples_code_cpp users_guide user_guide_snippets_py

# Set to -memcheck if you want to include a check for memory leaks
MEMCHECK =

lint :
	cd ../build && make lint; cd ../test

full: unittests examples benchmarks

full_knownbugs: unittests_knownbugs examples benchmarks

benchmarks:
	cd python && python complexity.py; cd ..

python: unittests_py examples_indoc_py examples_code_py user_guide_snippets_py

matlab: unittests_matlab examples_matlab
octave: unittests_octave examples_octave

python_knownbugs: unittests_py_knownbugs examples_indoc_py examples_code_py user_guide_snippets_py

cpp: examples_indoc_cpp  examples_code_cpp

unittests: unittests_py

unittests_knownbugs: unittests_py_knownbugs

unittests_py_slow:
	python internal/test_py.py python -skipfiles="alltests.py helpers.py complexity.py speed.py" -passoptions="--run_slow"

unittests_py_valgrind:
	python internal/test_py.py python -skipfiles="alltests.py helpers.py complexity.py speed.py" -memcheck -passoptions="--ignore_memory_heavy"

unittests_py:
ifdef MEMCHECK
	python internal/test_py.py python -skipfiles="alltests.py helpers.py complexity.py speed.py" -passoptions="--ignore_memory_light --run_slow"
	python internal/test_py.py python -skipfiles="alltests.py helpers.py complexity.py speed.py" $(MEMCHECK) -passoptions="--ignore_memory_heavy  --run_slow"
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

unittests_octave:
	python internal/test_octave.py octave
	python internal/test_octave.py matlab -skipfiles="callback.m nlp_callback.m"

examples_matlab:
	python internal/test_matlab.py '../docs/examples/matlab' -skipfiles="MyCallback.m NlpCallback.m"

examples_octave:
	python internal/test_octave.py '../docs/examples/matlab' -skipfiles="MyCallback.m NlpCallback.m callback.m nlp_callback.m"

# No memchecks because there is a memory leak in SuperLU
examples_indoc_py:
	python internal/test_py.py '../docs/api/examples'

# No memchecks because it leaks like hell, and is slow
examples_code_py:
	python internal/test_py.py '../docs/examples/python' -skipfiles="c_code_generation.py"

examples_indoc_cpp:
	python internal/test_cpp.py '../docs/api/examples' $(MEMCHECK)

doxygen:
	cd '../docs/api/' ; make doxtest

# No memchecks because it's slow
examples_code_cpp:
	python internal/test_cppcmake.py '../docs/examples/cplusplus' -skipfiles="test_liftopt.cpp casadi_error_handling.cpp test_opencl.cpp nlp_codegen.cpp codegen_usage.cpp stabilized_sqpmethod_nl.cpp knitro_nl.cpp worhp_nl.cpp snopt_nl.cpp rocket_snopt.cpp test_blocksqp.cpp" -skipdirs="cmake_find_package"
	cd "../docs/examples/cplusplus/cmake_find_package" && ./run.sh && cd ../../../../test

users_guide:
	cd '../docs/users_guide/' && make snippets && cd ..

user_guide_snippets_py: users_guide
	python internal/test_py.py '../docs/users_guide' -skipfiles="pytex.py temp_pytex.py"
