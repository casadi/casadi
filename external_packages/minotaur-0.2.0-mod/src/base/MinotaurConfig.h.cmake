#ifdef MINOTAURCONFIG_H
#error "Do not include MinotaurConfig.h in a header file."
#else
#define MINOTAURCONFIG_H
#endif

#include "MinotaurCFortran.h"

#ifdef F77_GLOBAL
#define F77_FUNC F77_GLOBAL
#define F77_FUNC_ F77_GLOBAL_
#endif

/* Define to 1 if you have the getrusage() function. */
#cmakedefine MINOTAUR_RUSAGE

