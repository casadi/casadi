/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2012 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file include/qpOASES/Types.hpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0beta
 *	\date 2007-2012
 *
 *	Declaration of all non-built-in types (except for classes).
 */


#ifndef QPOASES_TYPES_HPP
#define QPOASES_TYPES_HPP


/* If your compiler does not support the snprintf() function,
 * uncomment the following line and try to compile again. */
/* #define snprintf _snprintf */


/* Uncomment the following line for setting the __DSPACE__ flag. */
/* #define __DSPACE__ */

/* Uncomment the following line for setting the __XPCTARGET__ flag. */
/* #define __XPCTARGET__ */


/* Uncomment the following line to enable debug information. */
/* #define __DEBUG__ */

/* Uncomment the following line to enable debug iteration output. */
/* #define __DEBUG_ITER__ */

/* Uncomment the following line to enable suppress any kind of console output. */
/* #define __SUPPRESSANYOUTPUT__ */


/** Forces to always include all implicitly fixed bounds and all equality constraints
 *  into the initial working set when setting up an auxiliary QP. */
#define __ALWAYS_INITIALISE_WITH_ALL_EQUALITIES__

/* Uncomment the following line to activate the use of a special treatment of 
 * inactive constraints that is more efficient in case of QP formulations 
 * comprising many constraints. */
/* #define __MANY_CONSTRAINTS__ */

/* Uncomment the following line to activate the use of an alternative Givens 
 * plane rotation requiring only three multiplications. */
/* #define __USE_THREE_MULTS_GIVENS__ */

/* Uncomment the following line to activate the use of single precision arithmetic. */
/* #define __USE_SINGLE_PRECISION__ */



/* Work-around for Borland BCC 5.5 compiler. */
#ifdef __BORLANDC__
  #define __STDC__ 1
#endif


/* Work-around for Microsoft compilers. */
#ifdef _MSC_VER
  #define snprintf _snprintf
#endif


#ifdef __DSPACE__
  /* This fix ensures a compilable code only,
   * but all snprintf commands won't work. */
  #define snprintf printf
#endif

#ifdef __DSPACE__

	/** Macro for switching on/off the beginning of the qpOASES namespace definition. */
	#define BEGIN_NAMESPACE_QPOASES
    
	/** Macro for switching on/off the end of the qpOASES namespace definition. */
	#define END_NAMESPACE_QPOASES

	/** Macro for switching on/off the use of the qpOASES namespace. */
	#define USING_NAMESPACE_QPOASES

	/** Macro for switching on/off references to the qpOASES namespace. */
	#define REFER_NAMESPACE_QPOASES

#else

	/** Macro for switching on/off the beginning of the qpOASES namespace definition. */
	#define BEGIN_NAMESPACE_QPOASES  namespace qpOASES {

	/** Macro for switching on/off the end of the qpOASES namespace definition. */
	#define END_NAMESPACE_QPOASES    }
	
	/** Macro for switching on/off the use of the qpOASES namespace. */
	#define USING_NAMESPACE_QPOASES  using namespace qpOASES;
	
	/** Macro for switching on/off references to the qpOASES namespace. */
	#define REFER_NAMESPACE_QPOASES  qpOASES::

#endif


/** Macro for accessing the Cholesky factor R. */
#define RR( I,J )  R[(I)+nV*(J)]

/** Macro for accessing the orthonormal matrix Q of the QT factorisation. */
#define QQ( I,J )  Q[(I)+nV*(J)]

/** Macro for accessing the triangular matrix T of the QT factorisation. */
#define TT( I,J )  T[(I)*sizeT+(J)]



BEGIN_NAMESPACE_QPOASES


/** Defines real_t for facilitating switching between double and float. */
#ifdef __USE_SINGLE_PRECISION__
typedef float real_t;
#else
typedef double real_t;
#endif /* __USE_SINGLE_PRECISION__ */


/** Summarises all possible logical values. */
enum BooleanType
{
	BT_FALSE,					/**< Logical value for "false". */
	BT_TRUE						/**< Logical value for "true". */
};


/** Summarises all possible print levels. Print levels are used to describe
 *	the desired amount of output during runtime of qpOASES. */
enum PrintLevel
{
	PL_TABULAR = -1,			/**< Tabular output. */
	PL_NONE,					/**< No output. */
	PL_LOW,						/**< Print error messages only. */
	PL_MEDIUM,					/**< Print error and warning messages as well as concise info messages. */
	PL_HIGH						/**< Print all messages with full details. */
};


/** Defines visibility status of a message. */
enum VisibilityStatus
{
	VS_HIDDEN,					/**< Message not visible. */
	VS_VISIBLE					/**< Message visible. */
};


/** Summarises all possible states of the (S)QProblem(B) object during the
solution process of a QP sequence. */
enum QProblemStatus
{
	QPS_NOTINITIALISED,			/**< QProblem object is freshly instantiated or reset. */
	QPS_PREPARINGAUXILIARYQP,	/**< An auxiliary problem is currently setup, either at the very beginning
								 *   via an initial homotopy or after changing the QP matrices. */
	QPS_AUXILIARYQPSOLVED,		/**< An auxilary problem was solved, either at the very beginning
								 *   via an initial homotopy or after changing the QP matrices. */
	QPS_PERFORMINGHOMOTOPY,		/**< A homotopy according to the main idea of the online active
								 *   set strategy is performed. */
	QPS_HOMOTOPYQPSOLVED,		/**< An intermediate QP along the homotopy path was solved. */
	QPS_SOLVED					/**< The solution of the actual QP was found. */
};


/** Summarises all possible types of the QP's Hessian matrix. */
enum HessianType
{
	HST_ZERO,					/**< Hessian is zero matrix (i.e. LP formulation). */
	HST_IDENTITY,				/**< Hessian is identity matrix. */
	HST_POSDEF,					/**< Hessian is (strictly) positive definite. */
	HST_POSDEF_NULLSPACE,		/**< Hessian is positive definite on null space of active bounds/constraints. */
	HST_SEMIDEF,				/**< Hessian is positive semi-definite. */
	HST_UNKNOWN					/**< Hessian type is unknown. */
};


/** Summarises all possible types of bounds and constraints. */
enum SubjectToType
{
	ST_UNBOUNDED,		/**< Bound/constraint is unbounded. */
	ST_BOUNDED,			/**< Bound/constraint is bounded but not fixed. */
	ST_EQUALITY,		/**< Bound/constraint is fixed (implicit equality bound/constraint). */
	ST_DISABLED,
	ST_UNKNOWN			/**< Type of bound/constraint unknown. */
};


/** Summarises all possible states of bounds and constraints. */
enum SubjectToStatus
{
	ST_INACTIVE,			/**< Bound/constraint is inactive. */
	ST_LOWER,				/**< Bound/constraint is at its lower bound. */
	ST_UPPER,				/**< Bound/constraint is at its upper bound. */
	ST_INFEASIBLE_LOWER,
	ST_INFEASIBLE_UPPER,
	ST_UNDEFINED			/**< Status of bound/constraint undefined. */
};



END_NAMESPACE_QPOASES


#endif	/* QPOASES_TYPES_HPP */


/*
 *	end of file
 */
