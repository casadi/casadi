function(add_clang_tidy_target TARGET_NAME SOURCES_LIST0)

  list(REMOVE_DUPLICATES SOURCES_LIST0)
  list(SORT SOURCES_LIST0)

  # filtering out unwanted files
  set(SOURCES_LIST)
  foreach(item ${SOURCES_LIST0})
    string(REGEX MATCH ".*_runtime_str\\.h" dummy ${item})
    if(NOT dummy)
      string(REGEX MATCH ".*meta\\.cpp" dummy ${item})
    endif()
    if(NOT dummy)
      string(REGEX MATCH ".*feasiblesqpmethod.*" dummy ${item})
    endif()
    if(NOT dummy)
      list(APPEND SOURCES_LIST ${item})
    endif()
  endforeach()

  list(LENGTH SOURCES_LIST LISTCOUNT) 

  if(LISTCOUNT)
    add_custom_target(clang_tidy_${TARGET_NAME}
      COMMAND clang-tidy
              -p "${PROJECT_BINARY_DIR}/compile_commands.json"
              -warnings-as-errors=*
              -checks=readability-*,-readability-implicit-bool-cast,-readability-implicit-bool-conversion,-readability-braces-around-statements,-readability-else-after-return,readability-inconsistent-declaration-parameter-name,-readability-inconsistent-declaration-parameter-name,-readability-identifier-length,-readability-function-cognitive-complexity,-readability-isolate-declaration,-readability-use-anyofallof,-readability-magic-numbers,mbugprone-*,-bugprone-macro-parentheses,cppcoreguidelines-*,-cppcoreguidelines-owning-memory,-cppcoreguidelines-special-member-functions,-cppcoreguidelines-pro-bounds-pointer-arithmetic,-cppcoreguidelines-pro-type-member-init,-cppcoreguidelines-interfaces-global-init,-cppcoreguidelines-pro-type-cstyle-cast,-cppcoreguidelines-pro-type-vararg,-cppcoreguidelines-pro-type-static-cast-downcast,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-cppcoreguidelines-c-copy-assignment-signature,-cppcoreguidelines-pro-type-union-access,-cppcoreguidelines-pro-type-reinterpret-cast,-cppcoreguidelines-pro-type-const-cast,-cppcoreguidelines-pro-bounds-constant-array-index,-cppcoreguidelines-macro-usage,-cppcoreguidelines-avoid-magic-numbers,-cppcoreguidelines-init-variables,-clang-diagnostic-error
              "-header-filter=${PROJECT_SOURCE_DIR}/casadi/.*"
              ${SOURCES_LIST}
      DEPENDS ${SOURCES_LIST}
      COMMENT "Clang_tidy ${TARGET_NAME}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      VERBATIM)
  else()
    add_custom_target(clang_tidy_${TARGET_NAME}
    COMMAND echo "Nothing to do"
    COMMENT "Clang_tidy ${TARGET_NAME}"
    VERBATIM)
  endif()

  if(LISTCOUNT)
  add_custom_target(clang_tidy_fix_${TARGET_NAME}
    COMMAND clang-tidy
            -fix
            -checks=llvm-*,-llvm-header-guard
            -p "${PROJECT_BINARY_DIR}/compile_commands.json"
            "-header-filter=${PROJECT_SOURCE_DIR}/casadi/.*"
            ${SOURCES_LIST}
    DEPENDS ${SOURCES_LIST}
    COMMENT "Clang_tidy ${TARGET_NAME}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    VERBATIM)
  else()
    add_custom_target(clang_tidy_fix_${TARGET_NAME}
    COMMAND echo "Nothing to do"
    COMMENT "Clang_tidy ${TARGET_NAME}"
    VERBATIM)
  endif()
endfunction()# portability-*,
