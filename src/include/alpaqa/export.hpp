#pragma once

#include <alpaqa/export.h>

#ifdef _WIN32
#define ALPAQA_EXPORT_EXTERN_TEMPLATE(strcls, name, ...)                       \
    extern template strcls name<__VA_ARGS__>
#define ALPAQA_EXPORT_TEMPLATE(strcls, name, ...)                              \
    template strcls ALPAQA_EXPORT name<__VA_ARGS__>
#else
#define ALPAQA_EXPORT_EXTERN_TEMPLATE(strcls, name, ...)                       \
    extern template strcls ALPAQA_EXPORT name<__VA_ARGS__>
#define ALPAQA_EXPORT_TEMPLATE(strcls, name, ...)                              \
    template strcls name<__VA_ARGS__>
#endif