#ifndef fmi3PlatformTypes_h
#define fmi3PlatformTypes_h

/*
This header file defines the data types of FMI 3.0.
It must be used by both FMU and FMI master.

Copyright (C) 2008-2011 MODELISAR consortium,
              2012-2022 Modelica Association Project "FMI"
              All rights reserved.

This file is licensed by the copyright holders under the 2-Clause BSD License
(https://opensource.org/licenses/BSD-2-Clause):

----------------------------------------------------------------------------
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
----------------------------------------------------------------------------
*/

/* Include the integer and boolean type definitions */
#include <stdint.h>
#include <stdbool.h>


/* tag::Component[] */
typedef           void* fmi3Instance;             /* Pointer to the FMU instance */
/* end::Component[] */

/* tag::ComponentEnvironment[] */
typedef           void* fmi3InstanceEnvironment;  /* Pointer to the FMU environment */
/* end::ComponentEnvironment[] */

/* tag::FMUState[] */
typedef           void* fmi3FMUState;             /* Pointer to the internal FMU state */
/* end::FMUState[] */

/* tag::ValueReference[] */
typedef        uint32_t fmi3ValueReference;       /* Handle to the value of a variable */
/* end::ValueReference[] */

/* tag::VariableTypes[] */
typedef           float fmi3Float32;  /* Single precision floating point (32-bit) */
/* tag::fmi3Float64[] */
typedef          double fmi3Float64;  /* Double precision floating point (64-bit) */
/* end::fmi3Float64[] */
typedef          int8_t fmi3Int8;     /* 8-bit signed integer */
typedef         uint8_t fmi3UInt8;    /* 8-bit unsigned integer */
typedef         int16_t fmi3Int16;    /* 16-bit signed integer */
typedef        uint16_t fmi3UInt16;   /* 16-bit unsigned integer */
typedef         int32_t fmi3Int32;    /* 32-bit signed integer */
typedef        uint32_t fmi3UInt32;   /* 32-bit unsigned integer */
typedef         int64_t fmi3Int64;    /* 64-bit signed integer */
typedef        uint64_t fmi3UInt64;   /* 64-bit unsigned integer */
typedef            bool fmi3Boolean;  /* Data type to be used with fmi3True and fmi3False */
typedef            char fmi3Char;     /* Data type for one character */
typedef const fmi3Char* fmi3String;   /* Data type for character strings
                                         ('\0' terminated, UTF-8 encoded) */
typedef         uint8_t fmi3Byte;     /* Smallest addressable unit of the machine
                                         (typically one byte) */
typedef const fmi3Byte* fmi3Binary;   /* Data type for binary data
                                         (out-of-band length terminated) */
typedef            bool fmi3Clock;    /* Data type to be used with fmi3ClockActive and
                                         fmi3ClockInactive */

/* Values for fmi3Boolean */
#define fmi3True  true
#define fmi3False false

/* Values for fmi3Clock */
#define fmi3ClockActive   true
#define fmi3ClockInactive false
/* end::VariableTypes[] */

#endif /* fmi3PlatformTypes_h */
