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
 *	\file include/qpOASES/Matrices.hpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.0beta
 *	\date 2009-2012
 *
 *  Various matrix classes: Abstract base matrix class, dense and sparse matrices,
 *  including symmetry exploiting specializations.
 */



#ifndef QPOASES_MATRICES_HPP
#define QPOASES_MATRICES_HPP


#ifdef __USE_SINGLE_PRECISION__

	/** Macro for calling level 3 BLAS operation in single precision. */
	#define GEMM sgemm_
	/** Macro for calling level 3 BLAS operation in single precision. */
	#define SYR ssyr_
	/** Macro for calling level 3 BLAS operation in single precision. */
	#define SYR2 ssyr2_
	/** Macro for calling level 3 BLAS operation in single precision. */
	#define POTRF spotrf_

#else

	/** Macro for calling level 3 BLAS operation in double precision. */
	#define GEMM dgemm_
	/** Macro for calling level 3 BLAS operation in double precision. */
	#define SYR  dsyr_
	/** Macro for calling level 3 BLAS operation in double precision. */
	#define SYR2 dsyr2_
	/** Macro for calling level 3 BLAS operation in double precision. */
	#define POTRF dpotrf_

#endif /* __USE_SINGLE_PRECISION__ */


extern "C"
{
	/** Performs one of the matrix-matrix operation in double precision. */
	void dgemm_ ( const char*, const char*, const unsigned long*, const unsigned long*, const unsigned long*,
			const double*, const double*, const unsigned long*, const double*, const unsigned long*,
			const double*, double*, const unsigned long* );
	/** Performs one of the matrix-matrix operation in single precision. */
	void sgemm_ ( const char*, const char*, const unsigned long*, const unsigned long*, const unsigned long*,
			const float*, const float*, const unsigned long*, const float*, const unsigned long*,
			const float*, float*, const unsigned long* );

	/** Performs a symmetric rank 1 operation in double precision. */
	void dsyr_ ( const char *, const unsigned long *, const double *, const double *,
				 const unsigned long *, double *, const unsigned long *);
	/** Performs a symmetric rank 1 operation in single precision. */
	void ssyr_ ( const char *, const unsigned long *, const float *, const float *,
				 const unsigned long *, float *, const unsigned long *);

	/** Performs a symmetric rank 2 operation in double precision. */
	void dsyr2_ ( const char *, const unsigned long *, const double *, const double *,
				  const unsigned long *, const double *, const unsigned long *, double *, const unsigned long *);
	/** Performs a symmetric rank 2 operation in single precision. */
	void ssyr2_ ( const char *, const unsigned long *, const float *, const float *,
				  const unsigned long *, const float *, const unsigned long *, float *, const unsigned long *);

	/** Calculates the Cholesky factorization of a real symmetric positive definite matrix in double precision. */
	void dpotrf_ ( const char *, const unsigned long *, double *, const unsigned long *, long * );
	/** Calculates the Cholesky factorization of a real symmetric positive definite matrix in single precision. */
	void spotrf_ ( const char *, const unsigned long *, float *, const unsigned long *, long * );
}


BEGIN_NAMESPACE_QPOASES

/*
 * integer type for sparse matrix row/column entries. Make this "int"
 * for 32 bit entries, and "long" for 64-bit entries on x86_64 platform.
 *
 * Most sparse codes still assume 32-bit entries here (HSL, BQPD, ...)
 */
typedef int sparse_int_t;


/**
 *	\brief Abstract base class for interfacing tailored matrix-vector operations.
 *
 *	Abstract base matrix class. Supplies interface to matrix vector
 *  products, including products with submatrices given by (ordered) working set
 *  index lists (see \a SubjectTo).
 *
 *	\author Andreas Potschka, Christian Kirches, Hans Joachim Ferreau
 *	\version 3.0beta
 *	\date 2011
 */
class Matrix
{
	public:
		/** Default constructor. */
		Matrix( ) { doNotFreeMemory(); };

		/** Destructor. */
		virtual ~Matrix( ) { };

		/** Frees all internal memory. */
		virtual void free( ) = 0;

		/** Returns a deep-copy of the Matrix object.
		 *	\return Deep-copy of Matrix object */
		virtual Matrix* duplicate( ) const = 0;

		/** Returns i-th diagonal entry.
		 *	\return i-th diagonal entry */
		virtual real_t diag(	int i			/**< Index. */
								) const = 0;

		/** Checks whether matrix is square and diagonal.
		 *	\return BT_TRUE  iff matrix is square and diagonal; \n
		 *	        BT_FALSE otherwise. */
		virtual BooleanType isDiag( ) const = 0;

        /** Get the two-norm of the matrix
         *  \return Two-norm of the matrix
         */
        virtual real_t getNorm () const = 0;

        /** Get the two-norm of a row
         *  \return Two-norm of row \a rNum
         */
        virtual real_t getRowNorm (int rNum) const = 0;

		/** Retrieve indexed entries of matrix row multiplied by alpha.
		 *	\return SUCCESSFUL_RETURN */
		virtual returnValue getRow(
				int rNum,						/**< Row number. */
				const Indexlist* const icols,	/**< Index list specifying columns. */
				real_t alpha,					/**< Scalar factor. */
				real_t *row						/**< Output row vector. */
				) const = 0;

		/** Retrieve indexed entries of matrix column multiplied by alpha.
		 *	\return SUCCESSFUL_RETURN */
		virtual returnValue getCol(
				int cNum,						/**< Column number. */
				const Indexlist* const irows,	/**< Index list specifying rows. */
				real_t alpha,					/**< Scalar factor. */
				real_t *col						/**< Output column vector. */
				) const = 0;

		/** Evaluate Y=alpha*A*X + beta*Y.
		 *	\return SUCCESSFUL_RETURN */
		virtual returnValue times (	int xN,					/**< Number of vectors to multiply. */
									real_t alpha,			/**< Scalar factor for matrix vector product. */
									const real_t *x,		/**< Input vector to be multiplied. */
									int xLD,				/**< Leading dimension of input x. */
									real_t beta,			/**< Scalar factor for y. */
									real_t *y,				/**< Output vector of results. */
									int yLD					/**< Leading dimension of output y. */
									) const = 0;

		/** Evaluate Y=alpha*A'*X + beta*Y.
		 *	\return SUCCESSFUL_RETURN */
		virtual returnValue transTimes (	int xN,				/**< Number of vectors to multiply. */
											real_t alpha,		/**< Scalar factor for matrix vector product. */
											const real_t *x,	/**< Input vector to be multiplied. */
											int xLD,			/**< Leading dimension of input x. */
											real_t beta,		/**< Scalar factor for y. */
											real_t *y,			/**< Output vector of results. */
											int yLD				/**< Leading dimension of output y. */
											) const = 0;

		/** Evaluate matrix vector product with submatrix given by Indexlist.
		 *	\return SUCCESSFUL_RETURN */
		virtual returnValue times (	const Indexlist* const irows,	/**< Index list specifying rows. */
									const Indexlist* const icols,	/**< Index list specifying columns. */
									int xN,							/**< Number of vectors to multiply. */
									real_t alpha,					/**< Scalar factor for matrix vector product. */
									const real_t *x,				/**< Input vector to be multiplied. */
									int xLD,						/**< Leading dimension of input x. */
									real_t beta,					/**< Scalar factor for y. */
									real_t *y,						/**< Output vector of results. */
									int yLD,						/**< Leading dimension of output y. */
									BooleanType yCompr = BT_TRUE	/**< Compressed storage for y. */
									) const = 0;

		/** Evaluate matrix transpose vector product.
		 *	\return SUCCESSFUL_RETURN */
		virtual returnValue transTimes (	const Indexlist* const irows,	/**< Index list specifying rows. */
											const Indexlist* const icols,	/**< Index list specifying columns. */
											int xN,							/**< Number of vectors to multiply. */
											real_t alpha,					/**< Scalar factor for matrix vector product. */
											const real_t *x,				/**< Input vector to be multiplied. */
											int xLD,						/**< Leading dimension of input x. */
											real_t beta,					/**< Scalar factor for y. */
											real_t *y,						/**< Output vector of results. */
											int yLD							/**< Leading dimension of output y. */
											) const = 0;

		/** Adds given offset to diagonal of matrix.
		 *	\return SUCCESSFUL_RETURN \n
		 			RET_NO_DIAGONAL_AVAILABLE */
		virtual returnValue addToDiag(	real_t alpha		/**< Diagonal offset. */
										) = 0;


		/** Prints matrix to screen.
		 *	\return SUCCESSFUL_RETURN \n
		 *	        RET_NOT_YET_IMPLEMENTED */
		virtual returnValue print( ) const = 0;


		/** Returns whether internal memory needs to be de-allocated.
		 *	\return BT_TRUE  iff internal memory needs to be de-allocated, \n
		 			BT_FALSE otherwise */
		BooleanType needToFreeMemory( ) const { return freeMemory; };

		/** Enables de-allocation of internal memory. */
		void doFreeMemory( ) { freeMemory = BT_TRUE; };

		/** Disables de-allocation of internal memory. */
		void doNotFreeMemory( ) { freeMemory = BT_FALSE; };


	protected:
			BooleanType freeMemory;				/**< Indicating whether internal memory needs to be de-allocated. */

};


/**
 *	\brief Abstract base class for interfacing matrix-vector operations tailored to symmetric matrices.
 *
 *	Abstract base class for symmetric matrices. Extends Matrix interface with
 *  bilinear form evaluation.
 *
 *	\author Andreas Potschka, Christian Kirches, Hans Joachim Ferreau
 *	\version 3.0beta
 *	\date 2011
 */
class SymmetricMatrix : public virtual Matrix
{
	public:

		/** Compute bilinear form y = x'*H*x using submatrix given by index list.
		 *	\return SUCCESSFUL_RETURN */
		virtual returnValue bilinear(
				const Indexlist* const icols,	/**< Index list specifying columns of x. */
				int xN,							/**< Number of vectors to multiply. */
				const real_t *x,				/**< Input vector to be multiplied (uncompressed). */
				int xLD,						/**< Leading dimension of input x. */
				real_t *y,						/**< Output vector of results (compressed). */
				int yLD							/**< Leading dimension of output y. */
				) const = 0;

};


/**
 *	\brief Interfaces matrix-vector operations tailored to general dense matrices.
 *
 *	Dense matrix class (row major format).
 *
 *	\author Andreas Potschka, Christian Kirches, Hans Joachim Ferreau
 *	\version 3.0beta
 *	\date 2011
 */
class DenseMatrix : public virtual Matrix
{
	public:
		/** Default constructor. */
		DenseMatrix( ) : nRows(0), nCols(0), leaDim(0), val(0) {};

		/** Constructor from vector of values.
		 *  Caution: Data pointer must be valid throughout lifetime
		 */
		DenseMatrix(	int m,			/**< Number of rows. */
						int n,			/**< Number of columns. */
						int lD,			/**< Leading dimension. */
						real_t *v		/**< Values. */
						) : nRows(m), nCols(n), leaDim(lD), val(v) {}


		/** Destructor. */
		virtual ~DenseMatrix();

		/** Frees all internal memory. */
		virtual void free( );

		/** Returns a deep-copy of the Matrix object.
		 *	\return Deep-copy of Matrix object */
		virtual Matrix *duplicate( ) const;

		/** Returns i-th diagonal entry.
		 *	\return i-th diagonal entry */
		virtual real_t diag(	int i			/**< Index. */
								) const;

		/** Checks whether matrix is square and diagonal.
		 *	\return BT_TRUE  iff matrix is square and diagonal; \n
		 *	        BT_FALSE otherwise. */
		virtual BooleanType isDiag( ) const;

        /** Get the two-norm of the matrix
         *  \return Two-norm of the matrix
         */
        virtual real_t getNorm () const;
		
        /** Get the two-norm of a row
         *  \return Two-norm of row \a rNum
         */
        virtual real_t getRowNorm (int rNum) const;

        /** Retrieve indexed entries of matrix row multiplied by alpha.
		 *  \return SUCCESSFUL_RETURN */
		virtual returnValue getRow(
				int rNum,						/**< Row number. */
				const Indexlist* const icols,	/**< Index list specifying columns. */
				real_t alpha,					/**< Scalar factor. */
				real_t *row						/**< Output row vector. */
				) const;

		/** Retrieve indexed entries of matrix column multiplied by alpha.
		 *  \return SUCCESSFUL_RETURN */
		virtual returnValue getCol(
				int cNum,						/**< Column number. */
				const Indexlist* const irows,	/**< Index list specifying rows. */
				real_t alpha,					/**< Scalar factor. */
				real_t *col						/**< Output column vector. */
				) const;

		/** Evaluate Y=alpha*A*X + beta*Y.
		 *  \return SUCCESSFUL_RETURN. */
		returnValue times (	int xN,					/**< Number of vectors to multiply. */
							real_t alpha,			/**< Scalar factor for matrix vector product. */
							const real_t *x,		/**< Input vector to be multiplied. */
							int xLD,				/**< Leading dimension of input x. */
							real_t beta,			/**< Scalar factor for y. */
							real_t *y,				/**< Output vector of results. */
							int yLD					/**< Leading dimension of output y. */
							) const;

		/** Evaluate Y=alpha*A'*X + beta*Y.
		 *  \return SUCCESSFUL_RETURN. */
		returnValue transTimes (	int xN,				/**< Number of vectors to multiply. */
									real_t alpha,		/**< Scalar factor for matrix vector product. */
									const real_t *x,	/**< Input vector to be multiplied. */
									int xLD,			/**< Leading dimension of input x. */
									real_t beta,		/**< Scalar factor for y. */
									real_t *y,			/**< Output vector of results. */
									int yLD				/**< Leading dimension of output y. */
									) const;

		/** Evaluate matrix vector product with submatrix given by Indexlist.
		 *	\return SUCCESSFUL_RETURN */
		virtual returnValue times (	const Indexlist* const irows,	/**< Index list specifying rows. */
									const Indexlist* const icols,	/**< Index list specifying columns. */
									int xN,							/**< Number of vectors to multiply. */
									real_t alpha,					/**< Scalar factor for matrix vector product. */
									const real_t *x,				/**< Input vector to be multiplied. */
									int xLD,						/**< Leading dimension of input x. */
									real_t beta,					/**< Scalar factor for y. */
									real_t *y,						/**< Output vector of results. */
									int yLD,						/**< Leading dimension of output y. */
									BooleanType yCompr = BT_TRUE	/**< Compressed storage for y. */
									) const;

		/** Evaluate matrix transpose vector product.
		 *	\return SUCCESSFUL_RETURN */
		virtual returnValue transTimes (	const Indexlist* const irows,	/**< Index list specifying rows. */
											const Indexlist* const icols,	/**< Index list specifying columns. */
											int xN,							/**< Number of vectors to multiply. */
											real_t alpha,					/**< Scalar factor for matrix vector product. */
											const real_t *x,				/**< Input vector to be multiplied. */
											int xLD,						/**< Leading dimension of input x. */
											real_t beta,					/**< Scalar factor for y. */
											real_t *y,						/**< Output vector of results. */
											int yLD							/**< Leading dimension of output y. */
											) const;

		/** Adds given offset to diagonal of matrix.
		 *	\return SUCCESSFUL_RETURN \n
		 			RET_NO_DIAGONAL_AVAILABLE */
		virtual returnValue addToDiag(	real_t alpha		/**< Diagonal offset. */
										);

		/** Prints matrix to screen.
		 *	\return SUCCESSFUL_RETURN */
		virtual returnValue print( ) const;


	protected:
		int nRows;			/**< Number of rows. */
		int nCols;			/**< Number of columns. */
		int leaDim;			/**< Leading dimension. */
		real_t *val;		/**< Vector of entries. */
};


/**
 *	\brief Interfaces matrix-vector operations tailored to symmetric dense matrices.
 *
 *	Symmetric dense matrix class.
 *
 *	\author Andreas Potschka, Christian Kirches, Hans Joachim Ferreau
 *	\version 3.0beta
 *	\date 2011
 */
class SymDenseMat : public DenseMatrix, public SymmetricMatrix
{
	public:
		/** Default constructor. */
		SymDenseMat() : DenseMatrix() {}

		/** Constructor from vector of values. */
		SymDenseMat(	int m,			/**< Number of rows. */
						int n,			/**< Number of columns. */
						int lD,			/**< Leading dimension. */
						real_t *v		/**< Values. */
						) : DenseMatrix(m, n, lD, v) {}

		/** Returns a deep-copy of the Matrix object.
		 *	\return Deep-copy of Matrix object */
		virtual Matrix *duplicate( ) const;

		/** Compute bilinear form y = x'*H*x using submatrix given by index list.
		 *	\return SUCCESSFUL_RETURN */
		virtual returnValue bilinear(
				const Indexlist* const icols,	/**< Index list specifying columns of x. */
				int xN,							/**< Number of vectors to multiply. */
				const real_t *x,				/**< Input vector to be multiplied (uncompressed). */
				int xLD,						/**< Leading dimension of input x. */
				real_t *y,						/**< Output vector of results (compressed). */
				int yLD							/**< Leading dimension of output y. */
				) const;
};


/**
 *	\brief Interfaces matrix-vector operations tailored to general sparse matrices.
 *
 *	Sparse matrix class (col compressed format).
 *
 *	\author Andreas Potschka, Christian Kirches, Hans Joachim Ferreau
 *	\version 3.0beta
 *	\date 2011
 */
class SparseMatrix : public virtual Matrix
{
	public:
		/** Default constructor. */
		SparseMatrix();

		/** Constructor with arguments. */
		SparseMatrix(
				int nr, 		/**< Number of rows. */
				int nc, 		/**< Number of columns. */
				sparse_int_t *r, 		/**< Row indices (length). */
				sparse_int_t *c, 		/**< Indices to first entry of columns (nCols+1). */
				real_t *v,		/**< Vector of entries (length). */
				sparse_int_t *d = 0		/**< Indices to first entry of lower triangle (including diagonal) (nCols). */
				);

		/** Constructor from dense matrix. */
		SparseMatrix(
				int nr, 				/**< Number of rows. */
				int nc,			 		/**< Number of columns. */
				int ld,					/**< Leading dimension. */
				const real_t * const v	/**< Row major stored matrix elements. */
				);

		/** Destructor. */
		virtual ~SparseMatrix();

		/** Frees all internal memory. */
		virtual void free( );

		/** Returns a deep-copy of the Matrix object.
		 *	\return Deep-copy of Matrix object */
		virtual Matrix *duplicate( ) const;

		/** Returns i-th diagonal entry.
		 *	\return i-th diagonal entry */
		virtual real_t diag(	int i			/**< Index. */
								) const;

		/** Checks whether matrix is square and diagonal.
		 *	\return BT_TRUE  iff matrix is square and diagonal; \n
		 *	        BT_FALSE otherwise. */
		virtual BooleanType isDiag( ) const;


        /** Get the two-norm of the matrix
         *  \return Two-norm of the matrix
         */
        virtual real_t getNorm () const;

        /** Get the two-norm of a row
         *  \return Two-norm of row \a rNum
         */
        virtual real_t getRowNorm (int rNum) const;

		/** Retrieve indexed entries of matrix row multiplied by alpha. */
		virtual returnValue getRow (
				int rNum,						/**< Row number. */
				const Indexlist* const icols,	/**< Index list specifying columns. */
				real_t alpha,					/**< Scalar factor. */
				real_t *row						/**< Output row vector. */
				) const;

		/** Retrieve indexed entries of matrix column multiplied by alpha. */
		virtual returnValue getCol (
				int cNum,						/**< Column number. */
				const Indexlist* const irows,	/**< Index list specifying rows. */
				real_t alpha,					/**< Scalar factor. */
				real_t *col						/**< Output column vector. */
				) const;

		/** Evaluate Y=alpha*A*X + beta*Y. */
		virtual returnValue times (	int xN,					/**< Number of vectors to multiply. */
									real_t alpha,			/**< Scalar factor for matrix vector product. */
									const real_t *x,		/**< Input vector to be multiplied. */
									int xLD,				/**< Leading dimension of input x. */
									real_t beta,			/**< Scalar factor for y. */
									real_t *y,				/**< Output vector of results. */
									int yLD					/**< Leading dimension of output y. */
									) const;

		/** Evaluate Y=alpha*A'*X + beta*Y. */
		virtual returnValue transTimes (	int xN,				/**< Number of vectors to multiply. */
											real_t alpha,		/**< Scalar factor for matrix vector product. */
											const real_t *x,	/**< Input vector to be multiplied. */
											int xLD,			/**< Leading dimension of input x. */
											real_t beta,		/**< Scalar factor for y. */
											real_t *y,			/**< Output vector of results. */
											int yLD				/**< Leading dimension of output y. */
											) const;

		/** Evaluate matrix vector product with submatrix given by Indexlist. */
		virtual returnValue times (	const Indexlist* const irows,	/**< Index list specifying rows. */
									const Indexlist* const icols,	/**< Index list specifying columns. */
									int xN,							/**< Number of vectors to multiply. */
									real_t alpha,					/**< Scalar factor for matrix vector product. */
									const real_t *x,				/**< Input vector to be multiplied. */
									int xLD,						/**< Leading dimension of input x. */
									real_t beta,					/**< Scalar factor for y. */
									real_t *y,						/**< Output vector of results. */
									int yLD,						/**< Leading dimension of output y. */
									BooleanType yCompr = BT_TRUE	/**< Compressed storage for y. */
									) const;

		/** Evaluate matrix transpose vector product. */
		virtual returnValue transTimes (	const Indexlist* const irows,	/**< Index list specifying rows. */
											const Indexlist* const icols,	/**< Index list specifying columns. */
											int xN,							/**< Number of vectors to multiply. */
											real_t alpha,					/**< Scalar factor for matrix vector product. */
											const real_t *x,				/**< Input vector to be multiplied. */
											int xLD,						/**< Leading dimension of input x. */
											real_t beta,					/**< Scalar factor for y. */
											real_t *y,						/**< Output vector of results. */
											int yLD							/**< Leading dimension of output y. */
											) const;

		/** Adds given offset to diagonal of matrix.
		 *	\return SUCCESSFUL_RETURN \n
		 			RET_NO_DIAGONAL_AVAILABLE */
		virtual returnValue addToDiag(	real_t alpha		/**< Diagonal offset. */
										);

		/** Create jd field from ir and jc.
		 *  \return Pointer to jd. */
		sparse_int_t *createDiagInfo();

		/** Allocate and create dense matrix in row major format.
		 *  \return Pointer to matrix array. */
		real_t *full() const;

		/** Prints matrix to screen.
		 *	\return RET_NOT_YET_IMPLEMENTED */
		virtual returnValue print( ) const;


	protected:
		int nRows;			/**< Number of rows. */
		int nCols;			/**< Number of columns. */
		sparse_int_t *ir;			/**< Row indices (length). */
		sparse_int_t *jc;			/**< Indices to first entry of columns (nCols+1). */
		sparse_int_t *jd;			/**< Indices to first entry of lower triangle (including diagonal) (nCols). */
		real_t *val;		/**< Vector of entries (length). */
};


/**
 *	\brief Interfaces matrix-vector operations tailored to general sparse matrices.
 *
 *	Sparse matrix class (row compressed format).
 *
 *	\author Andreas Potschka, Christian Kirches, Hans Joachim Ferreau
 *	\version 3.0beta
 *	\date 2011
 */
class SparseMatrixRow : public virtual Matrix
{
	public:
		/** Default constructor. */
		SparseMatrixRow();

		/** Constructor with arguments. */
		SparseMatrixRow(
				int nr, 		/**< Number of rows. */
				int nc, 		/**< Number of columns. */
				sparse_int_t *r, 		/**< Indices to first entry of rows (nRows+1). */
				sparse_int_t *c, 		/**< Column indices (length). */
				real_t *v,		/**< Vector of entries (length). */
				sparse_int_t *d = 0		/**< Indices to first entry of upper triangle (including diagonal) (nRows). */
				);

		/** Constructor from dense matrix. */
		SparseMatrixRow(
				int nr, 				/**< Number of rows. */
				int nc,			 		/**< Number of columns. */
				int ld,					/**< Leading dimension. */
				const real_t * const v	/**< Row major stored matrix elements. */
				);

		/** Destructor. */
		virtual ~SparseMatrixRow();

		/** Frees all internal memory. */
		virtual void free( );

		/** Returns a deep-copy of the Matrix object.
		 *	\return Deep-copy of Matrix object */
		virtual Matrix *duplicate( ) const;

		/** Returns i-th diagonal entry.
		 *	\return i-th diagonal entry */
		virtual real_t diag(	int i			/**< Index. */
								) const;

		/** Checks whether matrix is square and diagonal.
		 *	\return BT_TRUE  iff matrix is square and diagonal; \n
		 *	        BT_FALSE otherwise. */
		virtual BooleanType isDiag( ) const;

		
        /** Get the two-norm of the matrix
         *  \return Two-norm of the matrix
         */
        virtual real_t getNorm () const;

        /** Get the two-norm of a row
         *  \return Two-norm of row \a rNum
         */
        virtual real_t getRowNorm (int rNum) const;

		/** Retrieve indexed entries of matrix row multiplied by alpha. */
		virtual returnValue getRow (
				int rNum,						/**< Row number. */
				const Indexlist* const icols,	/**< Index list specifying columns. */
				real_t alpha,					/**< Scalar factor. */
				real_t *row						/**< Output row vector. */
				) const;

		/** Retrieve indexed entries of matrix column multiplied by alpha. */
		virtual returnValue getCol (
				int cNum,						/**< Column number. */
				const Indexlist* const irows,	/**< Index list specifying rows. */
				real_t alpha,					/**< Scalar factor. */
				real_t *col						/**< Output column vector. */
				) const;

		/** Evaluate Y=alpha*A*X + beta*Y. */
		virtual returnValue times (	int xN,					/**< Number of vectors to multiply. */
									real_t alpha,			/**< Scalar factor for matrix vector product. */
									const real_t *x,		/**< Input vector to be multiplied. */
									int xLD,				/**< Leading dimension of input x. */
									real_t beta,			/**< Scalar factor for y. */
									real_t *y,				/**< Output vector of results. */
									int yLD					/**< Leading dimension of output y. */
									) const;

		/** Evaluate Y=alpha*A'*X + beta*Y. */
		virtual returnValue transTimes (	int xN,				/**< Number of vectors to multiply. */
											real_t alpha,		/**< Scalar factor for matrix vector product. */
											const real_t *x,	/**< Input vector to be multiplied. */
											int xLD,			/**< Leading dimension of input x. */
											real_t beta,		/**< Scalar factor for y. */
											real_t *y,			/**< Output vector of results. */
											int yLD				/**< Leading dimension of output y. */
											) const;

		/** Evaluate matrix vector product with submatrix given by Indexlist. */
		virtual returnValue times (	const Indexlist* const irows,	/**< Index list specifying rows. */
									const Indexlist* const icols,	/**< Index list specifying columns. */
									int xN,							/**< Number of vectors to multiply. */
									real_t alpha,					/**< Scalar factor for matrix vector product. */
									const real_t *x,				/**< Input vector to be multiplied. */
									int xLD,						/**< Leading dimension of input x. */
									real_t beta,					/**< Scalar factor for y. */
									real_t *y,						/**< Output vector of results. */
									int yLD,						/**< Leading dimension of output y. */
									BooleanType yCompr = BT_TRUE	/**< Compressed storage for y. */
									) const;

		/** Evaluate matrix transpose vector product. */
		virtual returnValue transTimes (	const Indexlist* const irows,	/**< Index list specifying rows. */
											const Indexlist* const icols,	/**< Index list specifying columns. */
											int xN,							/**< Number of vectors to multiply. */
											real_t alpha,					/**< Scalar factor for matrix vector product. */
											const real_t *x,				/**< Input vector to be multiplied. */
											int xLD,						/**< Leading dimension of input x. */
											real_t beta,					/**< Scalar factor for y. */
											real_t *y,						/**< Output vector of results. */
											int yLD							/**< Leading dimension of output y. */
											) const;

		/** Adds given offset to diagonal of matrix.
		 *	\return SUCCESSFUL_RETURN \n
		 			RET_NO_DIAGONAL_AVAILABLE */
		virtual returnValue addToDiag(	real_t alpha		/**< Diagonal offset. */
										);

		/** Create jd field from ir and jc.
		 *  \return Pointer to jd. */
		sparse_int_t *createDiagInfo();

		/** Allocate and create dense matrix in row major format.
		 *  \return Pointer to matrix array. */
		real_t *full() const;

		/** Prints matrix to screen.
		 *	\return RET_NOT_YET_IMPLEMENTED */
		virtual returnValue print( ) const;


	protected:
		int nRows;			/**< Number of rows. */
		int nCols;			/**< Number of columns. */
		sparse_int_t *jr;			/**< Indices to first entry of row (nRows+1). */
		sparse_int_t *ic;			/**< Column indices (length). */
		sparse_int_t *jd;			/**< Indices to first entry of upper triangle (including diagonal) (nRows). */
		real_t *val;		/**< Vector of entries (length). */
};


/**
 *	\brief Interfaces matrix-vector operations tailored to symmetric sparse matrices.
 *
 *	Symmetric sparse matrix class (column compressed format).
 *
 *	\author Andreas Potschka, Christian Kirches, Hans Joachim Ferreau
 *	\version 3.0beta
 *	\date 2011
 */
class SymSparseMat : public SymmetricMatrix, public SparseMatrix
{
	public:
		/** Default constructor. */
		SymSparseMat() : SparseMatrix() {}

		/** Constructor with arguments. */
		SymSparseMat(
				int nr, 		/**< Number of rows. */
				int nc, 		/**< Number of columns. */
				sparse_int_t *r, 		/**< Row indices (length). */
				sparse_int_t *c, 		/**< Indices to first entry of columns (nCols+1). */
				real_t *v,		/**< Vector of entries (length). */
				sparse_int_t *d = 0		/**< Indices to first entry of lower triangle (including diagonal) (nCols). */
				) : SparseMatrix(nr, nc, r, c, v, d) {}

		/** Constructor from dense matrix. */
		SymSparseMat(
				int nr, 				/**< Number of rows. */
				int nc,			 		/**< Number of columns. */
				int ld,					/**< Leading dimension. */
				const real_t * const v	/**< Row major stored matrix elements. */
				) : SparseMatrix(nr, nc, ld, v) {}

		/** Returns a deep-copy of the Matrix object.
		 *	\return Deep-copy of Matrix object */
		virtual Matrix *duplicate( ) const;


		/** Compute bilinear form y = x'*H*x using submatrix given by index list.
		*	\return SUCCESSFUL_RETURN */
		virtual returnValue bilinear(
				const Indexlist* const icols,	/**< Index list specifying columns of x. */
				int xN,							/**< Number of vectors to multiply. */
				const real_t *x,				/**< Input vector to be multiplied (uncompressed). */
				int xLD,						/**< Leading dimension of input x. */
				real_t *y,						/**< Output vector of results (compressed). */
				int yLD							/**< Leading dimension of output y. */
				) const;
};


END_NAMESPACE_QPOASES


#endif	/* QPOASES_MATRICES_HPP */
