function(add_warnings_target tgt_name warnings_as_errors)

    # GCC, Clang, AppleClang
    set(COMMON_WARNINGS
        -fdiagnostics-show-option
        -Wall
        -Wextra
        -pedantic
        -Wpedantic
        -pedantic-errors
        -Wdouble-promotion
        -Wswitch-default
        -Wswitch-enum
        -Wimplicit-fallthrough
        -Wuninitialized
        -Wno-missing-braces
    )
    # GCC
    set(GCC_WARNINGS
        -Wno-error=unused-but-set-variable
        -Wsuggest-override
        -Wno-error=attributes
    )
    # Clang, AppleClang
    set(CLANG_WARNINGS
        -Wno-error=unknown-warning-option
        -Wno-newline-eof
        -Wno-error=unused-but-set-variable
        -Winconsistent-missing-override
        -Wno-gnu-zero-variadic-macro-arguments
    )
    # MSVC (Microsoft)
    set(MSVC_WARNINGS
        /W3
        /wd4127 # conditional expression is constant
        /wd4458 # declaration of 'x' hides class member
        /wd4251 # 'x' needs to have dll-interface to be used by clients of 'y'
        /wd4305 # 'initializing': truncation from 'double' to 'float'
        /wd4661 # no suitable definition provided for explicit template instantiation request
        /permissive-
    )
    # Intel ICC
    set(INTEL_WARNINGS 
        -Wall
        -Wextra
    )

    # Enable warnings as errors
    if (warnings_as_errors)
        if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
            list(APPEND MSVC_WARNINGS /WX)
        elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
            list(APPEND INTEL_WARNINGS -Werror)
        else()
            list(APPEND COMMON_WARNINGS -Werror)
        endif()
    endif()

    # Add target that defines all the warning options in its interface.
    add_library(${tgt_name} INTERFACE)
    target_compile_options(${tgt_name} INTERFACE ${${PROJECT}_WARNNGS})
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(${tgt_name} INTERFACE
            $<BUILD_INTERFACE:${COMMON_WARNINGS} ${GCC_WARNINGS}>)
    elseif (CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
        target_compile_options(${tgt_name} INTERFACE
            $<BUILD_INTERFACE:${COMMON_WARNINGS} ${CLANG_WARNINGS}>)
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        target_compile_options(${tgt_name} INTERFACE
            $<BUILD_INTERFACE:${MSVC_WARNINGS}>)
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        target_compile_options(${tgt_name} INTERFACE
            $<BUILD_INTERFACE:${INTEL_WARNINGS}>)
    else()
        message(WARNING "No known warnings for this compiler")
    endif()

endfunction()