#pragma once

#include <alpaqa/export.hpp>

#include <string>
#include <typeinfo>

/// Get the pretty name of the given type as a string.
ALPAQA_EXPORT std::string demangled_typename(const std::type_info &t);
