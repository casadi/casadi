# Abseil Dependencies DSL
#
# This file contains the dependency graph for abseil libraries extracted from abslTargets.cmake.
# Format: library_name [dep1, dep2, ...]
#
# To regenerate when upgrading abseil:
#   cd build
#   python3 ../misc/generate_absl_deps_dsl.py > /tmp/new_dsl.txt 2>&1
#   # Then update the ABS_LIBS_WITH_DEPS variable below with the new DSL
#
# This DSL captures the INTERFACE_LINK_LIBRARIES relationships between abseil libraries,
# allowing CMake to properly resolve dependencies without needing --whole-archive linker flags.

set(ABS_LIBS_WITH_DEPS "
int128 []
log_severity []
random_internal_platform []
exponential_biased []
spinlock_wait []
strerror []
flags_commandlineflag_internal []
leak_check []
time_zone []
city []
civil_time []
utf8_for_code_point []
log_internal_nullguard []
low_level_hash [int128]
raw_logging_internal [log_severity]
random_internal_randen_hwaes_impl [random_internal_platform]
random_internal_randen_slow [random_internal_platform]
periodic_sampler [exponential_biased]
decode_rust_punycode [utf8_for_code_point]
debugging_internal [raw_logging_internal]
bad_optional_access [raw_logging_internal]
cordz_functions [exponential_biased, raw_logging_internal]
bad_any_cast_impl [raw_logging_internal]
base [log_severity, raw_logging_internal, spinlock_wait]
random_seed_gen_exception [raw_logging_internal]
strings_internal [raw_logging_internal]
scoped_set_env [raw_logging_internal]
bad_variant_access [raw_logging_internal]
throw_delegate [raw_logging_internal]
random_internal_randen_hwaes [random_internal_platform, random_internal_randen_hwaes_impl]
demangle_rust [decode_rust_punycode]
stacktrace [debugging_internal, raw_logging_internal]
log_internal_conditions [base]
crc_cpu_detect [base]
tracing_internal [base]
malloc_internal [base, raw_logging_internal]
string_view [base, throw_delegate]
random_internal_randen [random_internal_platform, random_internal_randen_hwaes, random_internal_randen_slow]
demangle_internal [demangle_rust]
crc_internal [crc_cpu_detect, raw_logging_internal]
graphcycles_internal [base, malloc_internal, raw_logging_internal]
poison [malloc_internal]
strings [string_view, strings_internal, base, int128, raw_logging_internal, throw_delegate]
flags_commandlineflag [flags_commandlineflag_internal, strings]
time [base, civil_time, int128, raw_logging_internal, strings, time_zone]
log_internal_fnmatch [strings]
log_internal_check_op [base, leak_check, log_internal_nullguard, strings]
random_distributions [strings]
flags_marshalling [log_severity, int128, strings]
hash [city, int128, strings, low_level_hash]
crc32c [crc_cpu_detect, crc_internal, strings]
random_internal_seed_material [raw_logging_internal, strings]
random_internal_distribution_test_util [raw_logging_internal, strings]
symbolize [debugging_internal, demangle_internal, base, malloc_internal, raw_logging_internal, strings]
die_if_null [strings]
log_internal_proto [base, strings]
str_format_internal [strings, int128]
flags_private_handle_accessor [flags_commandlineflag, flags_commandlineflag_internal, strings]
log_entry [log_severity, strings, time]
kernel_timeout_internal [base, raw_logging_internal, time]
log_internal_globals [log_severity, raw_logging_internal, strings, time]
crc_cord_state [crc32c, strings]
random_internal_pool_urbg [base, random_internal_randen, random_internal_seed_material, random_seed_gen_exception, raw_logging_internal]
examine_stack [stacktrace, symbolize, raw_logging_internal]
log_internal_structured_proto [log_internal_proto, strings]
log_sink [log_entry]
log_internal_format [log_internal_globals, log_severity, strings, time]
cord_internal [crc_cord_state, raw_logging_internal, strings, throw_delegate]
random_seed_sequences [random_internal_pool_urbg, random_internal_seed_material, random_seed_gen_exception, string_view]
failure_signal_handler [examine_stack, stacktrace, base, raw_logging_internal]
cord [base, cord_internal, cordz_functions, cordz_info, crc32c, crc_cord_state, raw_logging_internal, strings]
cordz_handle [base, raw_logging_internal, synchronization]
cordz_info [base, cord_internal, cordz_functions, cordz_handle, raw_logging_internal, stacktrace, synchronization, time]
cordz_sample_token [cordz_handle, cordz_info]
flags_config [flags_program_name, strings, synchronization]
flags_internal [base, flags_commandlineflag, flags_commandlineflag_internal, flags_config, flags_marshalling, synchronization]
flags_parse [flags_config, flags_commandlineflag, flags_commandlineflag_internal, flags_internal, flags_private_handle_accessor, flags_program_name, flags_reflection, flags_usage, strings, synchronization]
flags_program_name [strings, synchronization]
flags_reflection [flags_commandlineflag, flags_private_handle_accessor, flags_config, strings, synchronization]
flags_usage [flags_usage_internal, raw_logging_internal, strings, synchronization]
flags_usage_internal [flags_config, flags_commandlineflag, flags_internal, flags_private_handle_accessor, flags_program_name, flags_reflection, strings, synchronization]
hashtablez_sampler [base, exponential_biased, raw_logging_internal, synchronization, time]
log_flags [log_globals, log_severity, flags_marshalling, strings, vlog_config_internal]
log_globals [hash, log_severity, raw_logging_internal, strings, vlog_config_internal]
log_initialize [log_globals, log_internal_globals, time]
log_internal_log_sink_set [base, log_internal_globals, log_globals, log_entry, log_severity, log_sink, raw_logging_internal, synchronization, strings]
log_internal_message [base, examine_stack, log_internal_format, log_internal_globals, log_internal_proto, log_internal_log_sink_set, log_internal_nullguard, log_internal_structured_proto, log_globals, log_entry, log_severity, log_sink, raw_logging_internal, strerror, strings, time]
raw_hash_set [hash, hashtablez_sampler, raw_logging_internal]
status [cord, leak_check, raw_logging_internal, stacktrace, strerror, strings, symbolize]
statusor [base, raw_logging_internal, status, strings]
synchronization [graphcycles_internal, kernel_timeout_internal, base, malloc_internal, raw_logging_internal, stacktrace, symbolize, tracing_internal, time, tracing_internal]
vlog_config_internal [base, log_internal_fnmatch, strings, synchronization]
")

# Parse the DSL and create abseil targets with full dependency graph
#
# This function reads the ABS_LIBS_WITH_DEPS DSL and creates IMPORTED STATIC library
# targets for each abseil library, with proper INTERFACE_LINK_LIBRARIES relationships.
#
# The targets are created as absl-external::library_name and linked to the provided
# parent target (typically protobuf::libprotobuf).
#
# Arguments:
#   PARENT_TARGET - The target to link all abseil libraries to (e.g., protobuf::libprotobuf)
#
function(setup_abseil_dependencies PARENT_TARGET)
  # Parse the DSL and create targets with full dependency graph
  string(REPLACE "\n" ";" ABS_LIBS_LINES "${ABS_LIBS_WITH_DEPS}")

  foreach(line IN LISTS ABS_LIBS_LINES)
    string(STRIP "${line}" line)
    if(line STREQUAL "" OR line MATCHES "^#")
      continue()
    endif()

    # Parse: "library_name [dep1, dep2]"
    if(line MATCHES "^([a-z_0-9]+) \\[(.*)\\]$")
      set(lib_name "${CMAKE_MATCH_1}")
      set(deps_str "${CMAKE_MATCH_2}")

      # Create the imported library target
      add_library(absl-external::${lib_name} STATIC IMPORTED)
      set_target_properties(absl-external::${lib_name} PROPERTIES
        IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/external_projects/lib/${CMAKE_STATIC_LIBRARY_PREFIX}absl_${lib_name}${CMAKE_STATIC_LIBRARY_SUFFIX}"
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/external_projects/include")

      # Link to parent target (flat list for backward compatibility)
      target_link_libraries(${PARENT_TARGET} INTERFACE absl-external::${lib_name})

      # Add dependencies between abseil libraries
      # This recreates the INTERFACE_LINK_LIBRARIES graph from abseil's CMake config,
      # ensuring that when a library is linked, all its dependencies are automatically
      # included by CMake, which allows the linker to properly resolve all symbols
      # including RTTI typeinfo symbols.
      if(NOT deps_str STREQUAL "")
        string(REPLACE "," ";" dep_list "${deps_str}")
        foreach(dep IN LISTS dep_list)
          string(STRIP "${dep}" dep)
          if(NOT dep STREQUAL "")
            # Link this library to its dependencies
            target_link_libraries(absl-external::${lib_name} INTERFACE absl-external::${dep})
          endif()
        endforeach()
      endif()
    else()
      message(WARNING "Could not parse abseil dependency line: ${line}")
    endif()
  endforeach()
endfunction()
