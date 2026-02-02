# Creates dSYM for versioned library if LIB_PATH is a symlink
# Usage: cmake -DLIB_PATH=/path/to/lib.dylib -P create_versioned_dsym.cmake

if(IS_SYMLINK "${LIB_PATH}")
  file(READ_SYMLINK "${LIB_PATH}" _link_target)
  get_filename_component(_lib_dir "${LIB_PATH}" DIRECTORY)
  if(NOT IS_ABSOLUTE "${_link_target}")
    set(_real_path "${_lib_dir}/${_link_target}")
  else()
    set(_real_path "${_link_target}")
  endif()
  execute_process(COMMAND dsymutil "${_real_path}" -o "${_real_path}.dSYM")
endif()
