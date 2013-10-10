file(STRINGS "${INPUT}" FILE_CONTENTS)
file(WRITE ${OUTPUT}.hpp "#ifndef RESOURCE_${SYMBOL}_HPP\n#define RESOURCE_${SYMBOL}_HPP\nnamespace CasADi {\nextern const char * resource_${SYMBOL};" )

file(APPEND ${OUTPUT}.hpp "\n}\n#endif // RESOURCE_${SYMBOL}_HPP" )

file(WRITE ${OUTPUT}.cpp "#include \"resource_${SYMBOL}.hpp\"\n namespace CasADi {\nconst char * resource_${SYMBOL} =" )

foreach(line ${FILE_CONTENTS})
  string(REGEX REPLACE "\"" "\\\\\"" line "${line}")
  file(APPEND ${OUTPUT}.cpp "\n  \"${line}\\n\"" )
endforeach()

file(APPEND ${OUTPUT}.cpp ";\n}\n" )
