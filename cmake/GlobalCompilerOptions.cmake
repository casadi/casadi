if (WIN32)

    # Enable compiler warnings globally, interpret source code as UTF-8, and 
    # define the "not" alternative operator
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /utf-8 /Dnot=!")

    # Build-type specific flags
    set(CMAKE_CXX_FLAGS_DEBUG "/DDEBUG \
        /DEIGEN_INITIALIZE_MATRICES_BY_NAN")
    set(CMAKE_C_FLAGS_DEBUG   "/DDEBUG")
    set(CMAKE_DEBUG_POSTFIX "-debug")

else()

    # Enable compiler warnings globally
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
        -Wall -Wextra -Werror -pedantic -pedantic-errors -fdiagnostics-show-option \
        -Wmissing-include-dirs \
        -Wno-unknown-pragmas \
        -Wno-error=unused-parameter -Wno-error=unused-variable")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
        -Wall -Wextra -Werror -pedantic -pedantic-errors -fdiagnostics-show-option")

    # Build-type specific flags
    set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -O0 -DDEBUG \
        -DEIGEN_INITIALIZE_MATRICES_BY_NAN")
    set(CMAKE_C_FLAGS_DEBUG   "-ggdb3 -O0 -DDEBUG")
    set(CMAKE_DEBUG_POSTFIX "-debug")

    set(CMAKE_CXX_FLAGS_ASAN "-g3 -O0 \
        -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract \
        -fsanitize-address-use-after-scope \
        -fstack-protector-all \
        -DEIGEN_INITIALIZE_MATRICES_BY_NAN")
    set(CMAKE_C_FLAGS_ASAN   "-g3 -O0 \
        -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract")
    set(CMAKE_ASAN_POSTFIX "-asan")

    set(CMAKE_CXX_FLAGS_TSAN "-g3 -O0 \
        -fsanitize=thread,undefined")
    set(CMAKE_C_FLAGS_TSAN   "-g3 -O0 \
        -fsanitize=thread,undefined")
    set(CMAKE_TSAN_POSTFIX "-tsan")

    # Coverage
    set(CMAKE_C_FLAGS_COVERAGE "${CMAKE_C_FLAGS_ASAN} \
        --coverage -fno-inline")
    set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_ASAN} \
        --coverage -fno-inline")
    set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "--coverage")
    set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)

endif()