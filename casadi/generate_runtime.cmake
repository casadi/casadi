# Turn a runtime file into a file with strings
macro(CASADI_STRINGIFY STRFILE)
  # Start with an empty file
  file(WRITE ${STRFILE} "")
  foreach(FILE ${ARGN})
    # Add declaration of string
    get_filename_component(FILENAME ${FILE} NAME_WE)
    file(APPEND ${STRFILE} "const char* ${FILENAME}_str =")
    # Append file as strings
    file(STRINGS ${FILE} FILE_CONTENTS)
    foreach(LINE ${FILE_CONTENTS})
      string(REPLACE "\\" "\\\\" LINE "${LINE}") # Replace \ with \\
      string(REPLACE "\"" "\\\"" LINE "${LINE}") # Replace " with \"
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
