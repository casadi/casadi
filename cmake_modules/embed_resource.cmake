file(STRINGS "${INPUT}" FILE_CONTENTS)
file(WRITE ${OUTPUT} "namespace CasADi {\nconst char * resource_${SYMBOL} =" )

foreach(line ${FILE_CONTENTS})
  string(REGEX REPLACE "\"" "\\\\\"" line "${line}")
  file(APPEND ${OUTPUT} "\n  \"${line}\"" )
endforeach()

file(APPEND ${OUTPUT} ";\n}" )
