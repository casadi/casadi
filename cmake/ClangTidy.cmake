function(add_clang_tidy_target TARGET_NAME SOURCES_LIST0)

  list(REMOVE_DUPLICATES SOURCES_LIST0)
  list(SORT SOURCES_LIST0)

  # filtering out unwanted files
  set(SOURCES_LIST)
  foreach(item ${SOURCES_LIST0})
    string(REGEX MATCH ".*\\.hpp" dummy ${item})
    if(NOT dummy)
      string(REGEX MATCH ".*meta\\.cpp" dummy ${item})
      if(NOT dummy)
        list(APPEND SOURCES_LIST ${item})
      endif()
    endif()
  endforeach()

  add_custom_target(clang_tidy_${TARGET_NAME}
    COMMAND clang-tidy
            -p "${PROJECT_BINARY_DIR}/compile_commands.json"
            -warnings-as-errors=*
            -checks=readability-*,-readability-implicit-bool-cast,-readability-implicit-bool-conversion,-readability-braces-around-statements,-readability-else-after-return,readability-inconsistent-declaration-parameter-name,-readability-inconsistent-declaration-parameter-name,bugprone-*,-bugprone-macro-parentheses,cppcoreguidelines-*,-cppcoreguidelines-owning-memory,-cppcoreguidelines-special-member-functions,-cppcoreguidelines-pro-bounds-pointer-arithmetic,-cppcoreguidelines-pro-type-member-init,-cppcoreguidelines-interfaces-global-init,-cppcoreguidelines-pro-type-cstyle-cast,-cppcoreguidelines-pro-type-vararg,-cppcoreguidelines-pro-type-static-cast-downcast,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-cppcoreguidelines-c-copy-assignment-signature,-cppcoreguidelines-pro-type-union-access,-cppcoreguidelines-pro-type-reinterpret-cast,-cppcoreguidelines-pro-type-const-cast,-cppcoreguidelines-pro-bounds-constant-array-index
            "-header-filter=${PROJECT_SOURCE_DIR}/casadi/.*"
            ${SOURCES_LIST}
    DEPENDS ${SOURCES_LIST}
    COMMENT "Clang_tidy ${TARGET_NAME}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    VERBATIM)

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
endfunction()# portability-*,
