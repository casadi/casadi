# Turn a runtime file into a file with strings
macro(CASADI_STRINGIFY STRFILE)
  # Start with an empty file
  file(WRITE ${STRFILE} "")
  foreach(FILE ${ARGN})
    # Add declaration of string
    get_filename_component(FILENAME ${FILE} NAME_WE)
    file(APPEND ${STRFILE} "const char* ${FILENAME}_str =")
    # Append file as strings.  file(READ) instead of file(STRINGS): the latter
    # merges a line ending in backslash with the next, injecting a stray ';'.
    file(READ ${FILE} FILE_RAW)
    string(REPLACE "\r" "" FILE_RAW "${FILE_RAW}")
    # Note: square brackets, backslashes and semicolons corrupt the line splitting
    # Trick by Petr Kmoch
    string(REPLACE "^" "^!" FILE_RAW "${FILE_RAW}")
    string(REPLACE "[" "^a" FILE_RAW "${FILE_RAW}")
    string(REPLACE "]" "^b" FILE_RAW "${FILE_RAW}")
    string(REPLACE "\\" "^c" FILE_RAW "${FILE_RAW}")
    string(REPLACE ";" "^d" FILE_RAW "${FILE_RAW}")
    string(REPLACE "\n" ";" FILE_CONTENTS "${FILE_RAW}")
    foreach(LINE ${FILE_CONTENTS})
      string(REPLACE "\"" "\\\"" LINE "${LINE}") # Replace " with \"
      string(REPLACE "^d" ";" LINE "${LINE}")
      string(REPLACE "^c" "\\\\" LINE "${LINE}") # Replace \ with \\
      string(REPLACE "^b" "]" LINE "${LINE}")
      string(REPLACE "^a" "[" LINE "${LINE}")
      string(REPLACE "^!" "^" LINE "${LINE}")
      file(APPEND ${STRFILE} "\n  \"${LINE}\\n\"")
    endforeach()
    # End declaration
    file(APPEND ${STRFILE} ";\n\n")
  endforeach()
endmacro()

# Stringify C runtime
CASADI_STRINGIFY(${OUTPUT}
  ${SOURCES}
)
