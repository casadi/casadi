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
 *	\file include/qpOASES/Utils.hpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0beta
 *	\date 2007-2012
 *
 *	Declaration of some utilities for working with the different QProblem classes.
 */


#ifndef QPOASES_UTILS_HPP
#define QPOASES_UTILS_HPP


#include <qpOASES/MessageHandling.hpp>



BEGIN_NAMESPACE_QPOASES


/** Prints a vector.
 * \return SUCCESSFUL_RETURN */
returnValue print(	const real_t* const v,	/**< Vector to be printed. */
					int n					/**< Length of vector. */
					);

/** Prints a permuted vector.
 * \return SUCCESSFUL_RETURN */
returnValue print(	const real_t* const v,		/**< Vector to be printed. */
					int n,						/**< Length of vector. */
					const int* const V_idx		/**< Pemutation vector. */
					);

/** Prints a named vector.
 * \return SUCCESSFUL_RETURN */
returnValue print(	const real_t* const v,	/**< Vector to be printed. */
					int n,					/**< Length of vector. */
					const char* name		/** Name of vector. */
					);

/** Prints a matrix.
 * \return SUCCESSFUL_RETURN */
returnValue print(	const real_t* const M,	/**< Matrix to be printed. */
					int nrow,				/**< Row number of matrix. */
					int ncol				/**< Column number of matrix. */
					);

/** Prints a permuted matrix.
 * \return SUCCESSFUL_RETURN */
returnValue print(	const real_t* const M,		/**< Matrix to be printed. */
					int nrow,					/**< Row number of matrix. */
					int ncol	,				/**< Column number of matrix. */
					const int* const ROW_idx,	/**< Row pemutation vector. */
					const int* const COL_idx	/**< Column pemutation vector. */
					);

/** Prints a named matrix.
 * \return SUCCESSFUL_RETURN */
returnValue print(	const real_t* const M,	/**< Matrix to be printed. */
					int nrow,				/**< Row number of matrix. */
					int ncol,				/**< Column number of matrix. */
					const char* name		/** Name of matrix. */
					);

/** Prints an index array.
 * \return SUCCESSFUL_RETURN */
returnValue print(	const int* const index,	/**< Index array to be printed. */
					int n					/**< Length of index array. */
					);

/** Prints a named index array.
 * \return SUCCESSFUL_RETURN */
returnValue print(	const int* const index,	/**< Index array to be printed. */
					int n,					/**< Length of index array. */
					const char* name		/**< Name of index array. */
					);


/** Prints a string to desired output target (useful also for MATLAB output!).
 * \return SUCCESSFUL_RETURN */
returnValue myPrintf(	const char* s	/**< String to be written. */
						);


/** Prints qpOASES copyright notice.
 * \return SUCCESSFUL_RETURN */
returnValue printCopyrightNotice( );


/** Reads a real_t matrix from file.
 * \return SUCCESSFUL_RETURN \n
 		   RET_UNABLE_TO_OPEN_FILE \n
		   RET_UNABLE_TO_READ_FILE */
returnValue readFromFile(	real_t* data,				/**< Matrix to be read from file. */
							int nrow,					/**< Row number of matrix. */
							int ncol,					/**< Column number of matrix. */
							const char* datafilename	/**< Data file name. */
							);

/** Reads a real_t vector from file.
 * \return SUCCESSFUL_RETURN \n
 		   RET_UNABLE_TO_OPEN_FILE \n
		   RET_UNABLE_TO_READ_FILE */
returnValue readFromFile(	real_t* data,				/**< Vector to be read from file. */
							int n,						/**< Length of vector. */
							const char* datafilename	/**< Data file name. */
							);

/** Reads an integer (column) vector from file.
 * \return SUCCESSFUL_RETURN \n
 		   RET_UNABLE_TO_OPEN_FILE \n
		   RET_UNABLE_TO_READ_FILE */
returnValue readFromFile(	int* data,					/**< Vector to be read from file. */
							int n,						/**< Length of vector. */
							const char* datafilename	/**< Data file name. */
							);


/** Writes a real_t matrix into a file.
 * \return SUCCESSFUL_RETURN \n
 		   RET_UNABLE_TO_OPEN_FILE  */
returnValue writeIntoFile(	const real_t* const data,	/**< Matrix to be written into file. */
							int nrow,					/**< Row number of matrix. */
							int ncol,					/**< Column number of matrix. */
							const char* datafilename,	/**< Data file name. */
							BooleanType append			/**< Indicates if data shall be appended if the file already exists (otherwise it is overwritten). */
							);

/** Writes a real_t vector into a file.
 * \return SUCCESSFUL_RETURN \n
 		   RET_UNABLE_TO_OPEN_FILE  */
returnValue writeIntoFile(	const real_t* const data,	/**< Vector to be written into file. */
							int n,						/**< Length of vector. */
							const char* datafilename,	/**< Data file name. */
							BooleanType append			/**< Indicates if data shall be appended if the file already exists (otherwise it is overwritten). */
							);

/** Writes an integer (column) vector into a file.
 * \return SUCCESSFUL_RETURN \n
 		   RET_UNABLE_TO_OPEN_FILE */
returnValue writeIntoFile(	const int* const integer,	/**< Integer vector to be written into file. */
							int n,						/**< Length of vector. */
							const char* datafilename,	/**< Data file name. */
							BooleanType append			/**< Indicates if integer shall be appended if the file already exists (otherwise it is overwritten). */
							);


/** Returns the current system time.
 * \return current system time */
real_t getCPUtime( );


/** Returns the Euclidean norm of a vector.
 * \return 0: successful */
real_t getNorm(	const real_t* const v,	/**< Vector. */
				int n					/**< Vector's dimension. */
				);


/** Returns sign of a real_t-valued argument.
 * \return	 1.0: argument is non-negative \n
		 	-1.0: argument is negative */
inline real_t getSign(	real_t arg	/**< Double-valued argument whose sign is to be determined. */
						);


/** Returns maximum of two integers.
 * \return	Maximum of two integers */
inline int getMax(	int x,	/**< First integer. */
					int y	/**< Second integer. */
					);
					
/** Returns minimum of two integers.
 * \return	Minimum of two integers */
inline int getMin(	int x,	/**< First integer. */
					int y	/**< Second integer. */
					);

					
/** Returns maximum of two reals.
 * \return	Maximum of two reals */
inline real_t getMax(	real_t x,	/**< First real_t. */
						real_t y	/**< Second real_t. */
						);

/** Returns minimum of two reals.
 * \return	Minimum of two reals */
inline real_t getMin(	real_t x,	/**< First real_t. */
						real_t y	/**< Second real_t. */
						);

/** Returns the absolute value of a real number.
 * \return	Absolute value of a real number */
inline real_t getAbs(	real_t x	/**< Real number. */
						);


/** Computes "residual" of KKT system.  */
void getKKTResidual(	int nV,						/**< Number of variables. */
						int nC,						/**< Number of constraints. */
						const real_t* const H,		/**< Hessian matrix. */
						const real_t* const g,		/**< Sequence of gradient vectors. */
						const real_t* const A,		/**< Constraint matrix. */
						const real_t* const lb,		/**< Sequence of lower bound vectors (on variables). */
						const real_t* const ub,		/**< Sequence of upper bound vectors (on variables). */
						const real_t* const lbA,	/**< Sequence of lower constraints' bound vectors. */
						const real_t* const ubA,	/**< Sequence of upper constraints' bound vectors. */
						const real_t* const x,		/**< Sequence of primal trial vectors. */
						const real_t* const y,		/**< Sequence of dual trial vectors. */
						real_t& stat,				/**< Maximum value of stationarity condition residual. */
						real_t& feas,				/**< Maximum value of primal feasibility violation. */
						real_t& cmpl				/**< Maximum value of complementarity residual. */
						);


/** Writes a value of BooleanType into a string.
 * \return SUCCESSFUL_RETURN */
returnValue convertBooleanTypeToString(	BooleanType value, 		/**< Value to be written. */
										char* const string		/**< Input: String of sufficient size, \n
																	 Output: String containing value. */
										);

/** Writes a value of SubjectToStatus into a string.
 * \return SUCCESSFUL_RETURN */
returnValue convertSubjectToStatusToString(	SubjectToStatus value,	/**< Value to be written. */
											char* const string		/**< Input: String of sufficient size, \n
																		 Output: String containing value. */
											);

/** Writes a value of PrintLevel into a string.
 * \return SUCCESSFUL_RETURN */
returnValue convertPrintLevelToString(	PrintLevel value, 		/**< Value to be written. */
										char* const string		/**< Input: String of sufficient size, \n
																	 Output: String containing value. */
										);


#ifdef __DEBUG__
/** Writes matrix with given dimension into specified file. */
extern "C" void gdb_printmat(	const char *fname,			/**< File name. */
								real_t *M,					/**< Matrix to be written. */
								int n,						/**< Number of rows. */
								int m,						/**< Number of columns. */
								int ldim					/**< Leading dimension. */
								);
#endif /* __DEBUG__ */


END_NAMESPACE_QPOASES


#include <qpOASES/Utils.ipp>

#endif	/* QPOASES_UTILS_HPP */


/*
 *	end of file
 */
