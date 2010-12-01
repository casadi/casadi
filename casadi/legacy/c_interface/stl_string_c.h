#ifndef STL_STRING_C_H
#define STL_STRING_C_H

#ifdef __cplusplus
  extern "C" {
#endif 

#include <stdio.h>

/* STL string class */
typedef void* string_ptr;

/* Allocate STL string */
string_ptr casadi_string_new(void);

/* Delete STL string */
int casadi_string_delete(string_ptr str);

/* Assign STL string */
int casadi_string_assign(string_ptr str, const char* s);

/* Get char array */
const char* casadi_string_get(string_ptr str);

#ifdef __cplusplus
  }
#endif 

#endif /* STL_STRING_C_H */
