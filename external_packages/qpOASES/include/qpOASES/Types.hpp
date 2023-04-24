/*
 *  This file is part of qpOASES.
 *
 *  qpOASES -- An Implementation of the Online Active Set Strategy.
 *  Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
 *  Christian Kirches et al. All rights reserved.
 *
 *  qpOASES is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  qpOASES is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with qpOASES; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *  \file include/qpOASES/Types.hpp
 *  \author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *  \version 3.2
 *  \date 2007-2015
 *
 *  Declaration of all non-built-in types (except for classes).
 */


#ifndef QPOASES_TYPES_HPP
#define QPOASES_TYPES_HPP


/* If your compiler does not support the snprintf() function,
 * uncomment the following line and try to compile again. */
/* #define __NO_SNPRINTF__ */


/* Uncomment the following line for setting the __DSPACE__ flag. */
/* #define __DSPACE__ */

/* Uncomment the following line for setting the __XPCTARGET__ flag. */
/* #define __XPCTARGET__ */


/* Uncomment the following line for setting the __NO_FMATH__ flag. */
/* #define __NO_FMATH__ */

/* Uncomment the following line to enable debug information. */
/* #define __DEBUG__ */

/* Uncomment the following line to enable suppress any kind of console output. */
/* #define __SUPPRESSANYOUTPUT__ */


/** Forces to always include all implicitly fixed bounds and all equality constraints
 *  into the initial working set when setting up an auxiliary QP. */
#define __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__


/* Uncomment the following line to activate the use of an alternative Givens
 * plane rotation requiring only three multiplications. */
/* #define __USE_THREE_MULTS_GIVENS__ */

/* Uncomment the following line to activate the use of single precision arithmetic. */
/* #define __USE_SINGLE_PRECISION__ */



/* Work-around for Borland BCC 5.5 compiler. */
#ifdef __BORLANDC__
#if __BORLANDC__ < 0x0561
  #define __STDC__ 1
#endif
#endif


/* Work-around for Microsoft compilers. */
#ifdef _MSC_VER
  #define __NO_SNPRINTF__
  #pragma warning( disable : 4061 4100 4250 4514 4996 )
#endif


#ifdef __DSPACE__

    /** Macro for switching on/off the beginning of the qpOASES namespace definition. */
    #define BEGIN_NAMESPACE_QPOASES

    /** Macro for switching on/off the end of the qpOASES namespace definition. */
    #define END_NAMESPACE_QPOASES

    /** Macro for switching on/off the use of the qpOASES namespace. */
    #define USING_NAMESPACE_QPOASES

    /** Macro for switching on/off references to the qpOASES namespace. */
    #define REFER_NAMESPACE_QPOASES ::

#else

    /** Macro for switching on/off the beginning of the qpOASES namespace definition. */
    #define BEGIN_NAMESPACE_QPOASES  namespace casadi_qpOASES {

    /** Macro for switching on/off the end of the qpOASES namespace definition. */
    #define END_NAMESPACE_QPOASES    }

    /** Macro for switching on/off the use of the qpOASES namespace. */
    #define USING_NAMESPACE_QPOASES  using namespace casadi_qpOASES;

    /** Macro for switching on/off references to the qpOASES namespace. */
    #define REFER_NAMESPACE_QPOASES  casadi_qpOASES::

#endif


/* Avoid any printing on embedded platforms. */
#if defined(__DSPACE__) || defined(__XPCTARGET__)
  #define __SUPPRESSANYOUTPUT__
  #define __NO_SNPRINTF__
#endif


#ifdef __NO_SNPRINTF__
  #if (!defined(_MSC_VER)) || defined(__DSPACE__) || defined(__XPCTARGET__)
    /* If snprintf is not available, provide an empty implementation... */
    int snprintf( char* s, size_t n, const char* format, ... );
  #else
    /* ... or substitute snprintf by _snprintf for Microsoft compilers. */
    #define snprintf _snprintf
  #endif
#endif /* __NO_SNPRINTF__ */



/** Macro for accessing the Cholesky factor R. */
#define RR( I,J )  R[(I)+nV*(J)]

/** Macro for accessing the orthonormal matrix Q of the QT factorisation. */
#define QQ( I,J )  Q[(I)+nV*(J)]

/** Macro for accessing the triangular matrix T of the QT factorisation. */
#define TT( I,J )  T[(I)*sizeT+(J)]


/* If neither MA57 nor MA27 are selected, activate the dummy solver */
#if !defined(SOLVER_MA27) && !defined(SOLVER_MA57) && !defined(SOLVER_NONE)
#define SOLVER_NONE
#endif

BEGIN_NAMESPACE_QPOASES


/** Defines real_t for facilitating switching between double and float. */
#ifdef __USE_SINGLE_PRECISION__
typedef float real_t;
#else
typedef double real_t;
#endif /* __USE_SINGLE_PRECISION__ */


/** Defines int_t for facilitating switching between int and long int. */
#ifdef __USE_LONG_INTEGERS__
typedef long int_t;
typedef unsigned long uint_t;
#else
typedef int int_t;
typedef unsigned int uint_t;
#endif /* __USE_LONG_INTEGERS__ */


/** typedef for Fortran INTEGER type. Might be platform dependent! */
typedef int fint;


/**
 * Integer type for sparse matrix row/column entries. Make this "int"
 * for 32 bit entries, and "long" for 64-bit entries on x86_64 platform.
 *
 * Most sparse codes still assume 32-bit entries here (HSL, BQPD, ...)
 */
typedef int_t sparse_int_t;


/** Summarises all possible logical values. */
enum BooleanType
{
    BT_FALSE,                   /**< Logical value for "false". */
    BT_TRUE                     /**< Logical value for "true". */
};


/** Summarises all possible print levels. Print levels are used to describe
 *  the desired amount of output during runtime of qpOASES. */
enum PrintLevel
{
    PL_DEBUG_ITER = -2,         /**< Full tabular debugging output. */
    PL_TABULAR,                 /**< Normal tabular output. */
    PL_NONE,                    /**< No output. */
    PL_LOW,                     /**< Print error messages only. */
    PL_MEDIUM,                  /**< Print error and warning messages as well as concise info messages. */
    PL_HIGH                     /**< Print all messages with full details. */
};


/** Defines visibility status of a message. */
enum VisibilityStatus
{
    VS_HIDDEN,                  /**< Message not visible. */
    VS_VISIBLE                  /**< Message visible. */
};


/** Summarises all possible states of the (S)QProblem(B) object during the
solution process of a QP sequence. */
enum QProblemStatus
{
    QPS_NOTINITIALISED,         /**< QProblem object is freshly instantiated or reset. */
    QPS_PREPARINGAUXILIARYQP,   /**< An auxiliary problem is currently setup, either at the very beginning
                                 *   via an initial homotopy or after changing the QP matrices. */
    QPS_AUXILIARYQPSOLVED,      /**< An auxilary problem was solved, either at the very beginning
                                 *   via an initial homotopy or after changing the QP matrices. */
    QPS_PERFORMINGHOMOTOPY,     /**< A homotopy according to the main idea of the online active
                                 *   set strategy is performed. */
    QPS_HOMOTOPYQPSOLVED,       /**< An intermediate QP along the homotopy path was solved. */
    QPS_SOLVED                  /**< The solution of the actual QP was found. */
};


/** Summarises all possible types of the QP's Hessian matrix. */
enum HessianType
{
    HST_ZERO,               /**< Hessian is zero matrix (i.e. LP formulation). */
    HST_IDENTITY,           /**< Hessian is identity matrix. */
    HST_POSDEF,             /**< Hessian is (strictly) positive definite. */
    HST_POSDEF_NULLSPACE,   /**< Hessian is positive definite on null space of active bounds/constraints. */
    HST_SEMIDEF,            /**< Hessian is positive semi-definite. */
    HST_INDEF,              /**< Hessian is indefinite. */
    HST_UNKNOWN             /**< Hessian type is unknown. */
};


/** Summarises all possible types of bounds and constraints. */
enum SubjectToType
{
    ST_UNBOUNDED,       /**< Bound/constraint is unbounded. */
    ST_BOUNDED,         /**< Bound/constraint is bounded but not fixed. */
    ST_EQUALITY,        /**< Bound/constraint is fixed (implicit equality bound/constraint). */
    ST_DISABLED,        /**< Bound/constraint is disabled (i.e. ignored when solving QP). */
    ST_UNKNOWN          /**< Type of bound/constraint unknown. */
};


/** Summarises all possible states of bounds and constraints. */
enum SubjectToStatus
{
    ST_LOWER = -1,          /**< Bound/constraint is at its lower bound. */
    ST_INACTIVE,            /**< Bound/constraint is inactive. */
    ST_UPPER,               /**< Bound/constraint is at its upper bound. */
    ST_INFEASIBLE_LOWER,    /**< (to be documented) */
    ST_INFEASIBLE_UPPER,    /**< (to be documented) */
    ST_UNDEFINED            /**< Status of bound/constraint undefined. */
};

/** Flag indicating which type of update generated column in Schur complement. */
enum SchurUpdateType
{
    SUT_VarFixed,           /**< Free variable gets fixed. */
    SUT_VarFreed,           /**< Fixed variable gets freed. */
    SUT_ConAdded,           /**< Constraint becomes active. */
    SUT_ConRemoved,         /**< Constraint becomes inactive. */
    SUT_UNDEFINED           /**< Type of Schur update is undefined. */
};

/**
 *  \brief Stores internal information for tabular (debugging) output.
 *
 *  Struct storing internal information for tabular (debugging) output
 *  when using the (S)QProblem(B) objects.
 *
 *  \author Hans Joachim Ferreau
 *  \version 3.2
 *  \date 2013-2015
 */
struct TabularOutput {
    int_t idxAddB;      /**< Index of bound that has been added to working set. */
    int_t idxRemB;      /**< Index of bound that has been removed from working set. */
    int_t idxAddC;      /**< Index of constraint that has been added to working set. */
    int_t idxRemC;      /**< Index of constraint that has been removed from working set. */
    int_t excAddB;      /**< Flag indicating whether a bound has been added to working set to keep a regular projected Hessian. */
    int_t excRemB;      /**< Flag indicating whether a bound has been removed from working set to keep a regular projected Hessian. */
    int_t excAddC;      /**< Flag indicating whether a constraint has been added to working set to keep a regular projected Hessian. */
    int_t excRemC;      /**< Flag indicating whether a constraint has been removed from working set to keep a regular projected Hessian. */
};



/**
 *  \brief Struct containing the variable header for mat file.
 *
 *  Struct storing the header of a variable to be stored in
 *  Matlab's binary format (using the outdated Level 4 variant
 *  for simplictiy).
 *
 *  Note, this code snippet has been inspired from the document
 *  "Matlab(R) MAT-file Format, R2013b" by MathWorks
 *
 *  \author Hans Joachim Ferreau
 *  \version 3.2
 *  \date 2013-2015
 */
typedef struct {
    long numericFormat;     /**< Flag indicating numerical format. */
    long nRows;             /**< Number of rows. */
    long nCols;             /**< Number of rows. */
    long imaginaryPart;     /**< (to be documented) */
    long nCharName;         /**< Number of character in name. */
} MatMatrixHeader;

/** User-defined linear solver memory type */
typedef void* linsol_memory_t;

/** Initialization function for a user-defined linear solver function
 * Sparsity pattern in sparse triplet format (0-based)
 */
typedef int_t (*linsol_init_t)(linsol_memory_t mem,
                               int_t dim,
                               int_t nnz,
                               const int_t* row,
                               const int_t* col);

/** Symbolical factorization function for a user-defined linear solver function
  * Requires knowledge if the numerical values in order to perform e.g.
  * partial pivoting or to eliminate numerical zeros from the sparsity pattern
  */
typedef int_t (*linsol_sfact_t)(linsol_memory_t mem,
                                const real_t* vals);

/** Numerical factorization function for a user-defined linear solver function
  * Assumes a (not necessarily up-to-date) symbolic factorization is available.
  * The routine must calculate the number of negative eigenvalues as well as rank.
  */
typedef int_t (*linsol_nfact_t)(linsol_memory_t mem,
                                const real_t* vals, int_t* neig, int_t* rank);

/** Solve a factorized linear system for a user-defined linear solver function
  * Multiple right-hand-sides. The solution overwrites the right-hand-side.
  */
typedef int_t (*linsol_solve_t)(linsol_memory_t mem,
                                int_t nrhs, real_t* rhs);

/** Custom printing function.
  */
typedef void (*printf_t)(const char* s);

END_NAMESPACE_QPOASES


#endif  /* QPOASES_TYPES_HPP */


/*
 *  end of file
 */
