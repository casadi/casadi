function(add_warnings_target tgt_name warnings_as_errors)

    # GCC, Clang, AppleClang
    set(COMMON_WARNINGS_CXX
        -fdiagnostics-show-option
        -Wall
        -Wextra
        -pedantic
        -Wpedantic
        -pedantic-errors
        -Wdouble-promotion
        -Wswitch-default
        -Wswitch-enum
        -Wno-missing-braces
        -Wno-psabi
        -Wimplicit-fallthrough
        -Wuninitialized
        -Wconversion
    )
    set(COMMON_WARNINGS_Fortran
        -fdiagnostics-show-option
        -Wall
        -Wextra
        -Wpedantic
        -fimplicit-none
    )
    set(COMMON_WARNINGS_C
        -fdiagnostics-show-option
        -Wall
        -Wextra
        -pedantic
        -Wpedantic
        -pedantic-errors
        -Wimplicit-fallthrough
        -Wuninitialized
        -Wconversion
    )
    # GCC
    set(GCC_WARNINGS_CXX
        -Wno-error=unused-but-set-variable
        -Wsuggest-override
        -Wno-error=attributes
        -ftemplate-backtrace-limit=99
        -fconcepts-diagnostics-depth=99
    )
    # Clang, AppleClang
    set(CLANG_WARNINGS_CXX
        -Wno-error=unknown-warning-option
        -Wno-newline-eof
        -Wno-error=unused-but-set-variable
        -Winconsistent-missing-override
        -Wno-gnu-zero-variadic-macro-arguments
        -Wunreachable-code
        -Wunreachable-code-break
        -Wunreachable-code-fallthrough
        -Wunreachable-code-return
        -Wunreachable-code-aggressive
        -Wno-error=self-assign-overloaded
        -Wno-non-c-typedef-for-linkage
    )
    # Flang
    set(FLANG_WARNINGS_Fortran
        -fimplicit-none
    )
    # MSVC (Microsoft)
    set(MSVC_WARNINGS_CXX
        /W3
        /wd4127 # conditional expression is constant
        /wd4458 # declaration of 'x' hides class member
        /wd4251 # 'x' needs to have dll-interface to be used by clients of 'y'
        /wd4305 # 'initializing': truncation from 'double' to 'float'
        /wd4661 # no suitable definition provided for explicit template instantiation request
        /wd5030 # attribute 'attribute-name' is not recognized
    )
    set(MSVC_WARNINGS_C
        /W4
    )
    # Intel ICC
    set(INTEL_WARNINGS_CXX
        -Wall
        -Wextra
    )
    set(INTEL_WARNINGS_Fortran
        -Wall
        -Wextra
    )
    set(INTEL_WARNINGS_C
        -Wall
        -Wextra
    )

    # Add target that defines all the warning options in its interface.
    add_library(${tgt_name} INTERFACE)
    target_compile_options(${tgt_name} INTERFACE ${${PROJECT}_WARNNGS})
    
    foreach (LANG C CXX Fortran)
        # Enable warnings as errors
        if (warnings_as_errors)
            list(APPEND COMMON_WARNINGS_${LANG} -Werror)
            list(APPEND MSVC_WARNINGS_${LANG} /WX)
            list(APPEND INTEL_WARNINGS_${LANG} -Werror)
            list(APPEND FLANG_WARNINGS_${LANG} -Werror)
        endif()
        if (CMAKE_${LANG}_COMPILER_ID MATCHES "GNU")
            target_compile_options(${tgt_name} INTERFACE
            $<$<COMPILE_LANGUAGE:${LANG}>:$<BUILD_INTERFACE:${COMMON_WARNINGS_${LANG}} ${GCC_WARNINGS_${LANG}}>>)
        elseif (CMAKE_${LANG}_COMPILER_ID MATCHES ".*Clang")
            target_compile_options(${tgt_name} INTERFACE
            $<$<COMPILE_LANGUAGE:${LANG}>:$<BUILD_INTERFACE:${COMMON_WARNINGS_${LANG}} ${CLANG_WARNINGS_${LANG}}>>)
        elseif (CMAKE_${LANG}_COMPILER_ID MATCHES "Flang")
            target_compile_options(${tgt_name} INTERFACE
            $<$<COMPILE_LANGUAGE:${LANG}>:$<BUILD_INTERFACE:${FLANG_WARNINGS_${LANG}}>>)
        elseif (CMAKE_${LANG}_COMPILER_ID MATCHES "MSVC")
            target_compile_options(${tgt_name} INTERFACE
                $<$<COMPILE_LANGUAGE:${LANG}>:$<BUILD_INTERFACE:${MSVC_WARNINGS_${LANG}}>>)
        elseif (CMAKE_${LANG}_COMPILER_ID MATCHES "Intel")
            target_compile_options(${tgt_name} INTERFACE
            $<$<COMPILE_LANGUAGE:${LANG}>:$<BUILD_INTERFACE:${INTEL_WARNINGS_${LANG}}>>)
        else()
            message(WARNING "No known warnings for this ${LANG} compiler")
        endif()
    endforeach()

endfunction()