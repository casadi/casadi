include_directories(../../)

embed_resource("boilerplate" "profilereport.boilerplate" profilereport_boilerplate)

if(WITH_PROFILING)
add_executable(profilereport profilereport.cpp ${profilereport_boilerplate})
endif()
