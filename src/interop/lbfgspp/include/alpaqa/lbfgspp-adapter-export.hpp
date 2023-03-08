#pragma once

#include <alpaqa/lbfgspp-adapter-export.h>

#ifndef DOXYGEN

#ifdef _WIN32
#define ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(strcls, name, ...)               \
    extern template strcls name<__VA_ARGS__>
#define ALPAQA_LBFGSPP_EXPORT_TEMPLATE(strcls, name, ...)                      \
    template strcls LBFGSPP_ADAPTER_EXPORT name<__VA_ARGS__>
#else
#define ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(strcls, name, ...)               \
    extern template strcls LBFGSPP_ADAPTER_EXPORT name<__VA_ARGS__>
#define ALPAQA_LBFGSPP_EXPORT_TEMPLATE(strcls, name, ...)                      \
    template strcls name<__VA_ARGS__>
#endif

#else // DOXYGEN

#define ALPAQA_LBFGSPP_EXPORT_EXTERN_TEMPLATE(...)
#define ALPAQA_LBFGSPP_EXPORT_TEMPLATE(...)

#endif // DOXYGEN