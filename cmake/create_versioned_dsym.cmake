# Creates dSYM symlink for versioned library if LIB_PATH is a symlink
# Usage: cmake -DLIB_PATH=/path/to/lib.dylib -P create_versioned_dsym.cmake

if(IS_SYMLINK "${LIB_PATH}")
  file(READ_SYMLINK "${LIB_PATH}" _link_target)
  get_filename_component(_lib_dir "${LIB_PATH}" DIRECTORY)
  get_filename_component(_dsym_name "${LIB_PATH}" NAME)
  # Create symlink: libfoo.2.dylib.dSYM -> libfoo.dylib.dSYM
  file(CREATE_LINK "${_dsym_name}.dSYM" "${_lib_dir}/${_link_target}.dSYM" SYMBOLIC)
endif()
