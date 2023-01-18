#pragma once

#include <alpaqa/casadi-loader-export.h>

#ifdef _WIN32
#define CASADI_LOADER_EXPORT_EXTERN_TEMPLATE(strcls, name, ...)                \
    extern template strcls name<__VA_ARGS__>
#define CASADI_LOADER_EXPORT_TEMPLATE(strcls, name, ...)                       \
    template strcls CASADI_LOADER_EXPORT name<__VA_ARGS__>
#else
#define CASADI_LOADER_EXPORT_EXTERN_TEMPLATE(strcls, name, ...)                \
    extern template strcls CASADI_LOADER_EXPORT name<__VA_ARGS__>
#define CASADI_LOADER_EXPORT_TEMPLATE(strcls, name, ...)                       \
    template strcls name<__VA_ARGS__>
#endif