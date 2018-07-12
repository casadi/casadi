file(GLOB panoc_SRC
    "*.c"
)
file(GLOB constraint_SRC
    "constraints/*.c"
)

add_library(PANOC_lib STATIC ${panoc_SRC} ${constraint_SRC})