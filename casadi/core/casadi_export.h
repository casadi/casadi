/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef CASADI_EXPORT_H
#define CASADI_EXPORT_H

#ifdef CASADI_STATIC_DEFINE
#define CASADI_EXPORT
#define CASADI_NO_EXPORT
#else // CASADI_STATIC_DEFINE
#ifndef CASADI_EXPORT
#ifdef casadi_EXPORTS
// We are building this library
#define CASADI_EXPORT __attribute__((visibility("default")))
#else // casadi_EXPORTS
// We are using this library
#define CASADI_EXPORT __attribute__((visibility("default")))
#endif // casadi_EXPORTS
#endif // CASADI_EXPORT

#ifndef CASADI_NO_EXPORT
#define CASADI_NO_EXPORT __attribute__((visibility("hidden")))
#endif // CASADI_NO_EXPORT
#endif // CASADI_STATIC_DEFINE

#ifndef CASADI_DEPRECATED
#define CASADI_DEPRECATED __attribute__ ((__deprecated__))
#endif // CASADI_DEPRECATED

#ifndef CASADI_DEPRECATED_EXPORT
#define CASADI_DEPRECATED_EXPORT CASADI_EXPORT CASADI_DEPRECATED
#endif // CASADI_DEPRECATED_EXPORT

#ifndef CASADI_DEPRECATED_NO_EXPORT
#define CASADI_DEPRECATED_NO_EXPORT CASADI_NO_EXPORT CASADI_DEPRECATED
#endif // CASADI_DEPRECATED_NO_EXPORT

#endif
