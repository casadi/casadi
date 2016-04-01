// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 


#ifndef MINOTAUREIGEN_H
#define MINOTAUREIGEN_H

#include "LinearFunction.h"
#include "QuadraticFunction.h"

namespace Minotaur {

  class QuadraticFunction;
  typedef boost::shared_ptr<QuadraticFunction> QuadraticFunctionPtr;
  typedef boost::shared_ptr<const QuadraticFunction> ConstQuadraticFunctionPtr;
  //class EigenVector;
  //typedef boost::shared_ptr<const EigenVector> EigenVectorPtr;
  //typedef LinearFunction EigenVector;
  //typedef LinearFunctionPtr EigenVectorPtr;
  typedef std::pair<double, LinearFunctionPtr> EigenPair;
  typedef std::vector<EigenPair>::const_iterator EigenPairConstIterator;

  class Eigen;
  typedef boost::shared_ptr<Eigen> EigenPtr;
  
  // /**
  // In Minotaur, EigenVectors are used only in the context of quadratic 
  // functions. Thus an Eigen Vector is just a linear function.
  // */
  // class EigenVector {
  //   public:
  //     // /**
  //     // Default constructor
  //     // */
  //     EigenVector();

  //     // /**
  //     // Destroy
  //     // */
  //     ~EigenVector() {}

  //     // /**
  //     // Display the nonzero coefficients of the eigen vector.
  //     // */
  //     void write(std::ostream &s) const;

  //   private:
  //     // /**
  //     // The vector
  //     // */
  //     LinearFunctionPtr lf_;
  // };


  // /**
  // The EigenCalculator class is used to calculate eigen vectors and values
  // for matrices and other objects that have an associated matrix (like a
  // quadratic).
  // */
  class EigenCalculator {
  public:

    // /**
    // Construct using a quadratic function.
    // */
    EigenCalculator();

    // /**
    // Destroy
    // */
    ~EigenCalculator() {};

    // /**
    // Calculate EigenValues only
    // */
    EigenPtr findValues(ConstQuadraticFunctionPtr qf);

    // /**
    // Calculate EigenVectors as well
    // */
    EigenPtr findVectors(ConstQuadraticFunctionPtr qf); 

    // /**
    // Let qf = x'Ax, lf = cx. First find eigen vectors of the hessian of
    // qf. Then, x'Ax = x'QRERQ'x, where Q is orthogonal (QQ' = I). R is a
    // diagonal matrix with i-th diagonal element being the square root of
    // i-th eigen value. E has entries 1,-1 along the diagonal. Let y =
    // RQ'x. Then c'x = c'QR^(-1)y = b'y.  b = R^(-1)Q'c. This function
    // calculates Q, RER, b. 
    //
    // The "x" vectors in qf and lf are not the same. e.g.  x0*x0 + x1*x1 +
    // 2x0*x1 + 2x0 + x2 + x3.  here lf does not contain x1 and has extra
    // variables x2 and x3. The vector 'c' thus corresponds to only 2x0 and
    // we get: 1*(x0 + x1)^2 + 0*(x0 - x1)^2 + 1*(x0+x1) + 1*(x0-x1) + x2 +
    // x3.  p_terms will thus have (x0 + x1 + 0.5)^2, n_terms will have
    // terms that have negative eigen values.  lin_terms will have
    // x0-x1+x2+x3 and cb = -0.25.
    // */
    void getSumOfSquares (
                          std::vector<LinearFunctionPtr> & p_terms, 
                          std::vector<LinearFunctionPtr> & n_terms,
                          std::vector<double> & p_const,
                          std::vector<double> & n_const,
                          LinearFunctionPtr & lin_terms, double & c,
                          ConstQuadraticFunctionPtr qf, ConstLinearFunctionPtr lf);

  private:

    // /**
    // The quadratic function for whose Hessian we wish to find the eigen values.
    // */
    ConstQuadraticFunctionPtr qf_;

    // /**
    // Dimension of the square matrix
    // */
    UInt n_;

    // /**
    // The square matrix is stored as a single array. The element A[i,j] can
    // be accessed at A_[i+j*n_]. And element A_[i] = A[i mod n_, i/n_].
    // */
    double *A_;

    // /**
    // Number of eigen values found.
    // */
    int m_;

    // /**
    // The absolute error tolerance for the eigenvalues.  An approximate
    // eigenvalue is accepted as converged when it is determined to lie in an
    // interval [a,b] of width less than or equal to 
    // ABSTOL + EPS *   max(|a|,|b|)
    // */
    double abstol_;

    // /**
    // N for eigen values only, V for values and vectors.
    // */
    char findVectors_; 

    // /**
    // Array where the eigen vectors are stored by LAPACK.
    // */
    double *w_;

    // /**
    // the first M columns of z contain the orthonormal eigenvectors of the
    // matrix A corresponding to the selected eigenvalues, with the i-th
    // column of Z holding the eigenvector associated with W(i). If JOBZ =
    // 'N', then Z is not referenced.  Note: the user must ensure that at
    // least max(1,M) columns are supplied in the array Z; if RANGE = 'V',
    // the exact value of M is not known in advance and an upper bound must
    // be used.  Supplying N columns is always safe.
    // */
    double *z_;        

    // /**
    // The i-th eigenvector is nonzero only in elements 
    // ISUPPZ( 2*i-1 ) through ISUPPZ( 2*i ).
    // */
    int *isuppz_;    

    // /**
    // A map of what is the  index in matrix A_ of a variable
    // */
    std::map<ConstVariablePtr, UInt, CompareVariablePtr> indices_;

    // /**
    // What variable does column 'i' of matrix A_ represent. 
    // */
    std::vector<ConstVariablePtr> vars_;

    // /**
    // Dot product of coefficients of two linear functions
    // */
    double getDotProduct_(ConstLinearFunctionPtr lf1, 
                          ConstLinearFunctionPtr lf2);

    // /**
    // Call lapack routines to calculate the values.
    // */
    void calculate_();

    // /**
    // Allocate the A_ matrix and fill in the values from the quadratic
    // function qf.
    // */
    void fillA_();

    // /**
    // Get eigen values and (if calculated) eigen vectors from the
    // calculator.
    // */
    EigenPtr getEigen_();

    // /**
    // Construct a linear function based on the eigen vectors of A_.
    // */
    LinearFunctionPtr getLinearFunction_(const int i);
  };

  class Eigen {
  public:
    // /**
    // Default constructor
    // */
    Eigen();

    // /**
    // Add an eigen value and an eigen vector to the current list.
    // */
    void add(double value, LinearFunctionPtr e_vector);

    // /**
    // Get the number of negative eigen values
    // */
    UInt numNegative() const;

    // /**
    // Get the number of zero eigen values
    // */
    UInt numZero() const;

    // /**
    // Get the number of positive eigen values
    // */
    UInt numPositive() const;

    // /**
    // Get the first evPair
    // */
    EigenPairConstIterator begin() const;

    // /**
    // Get the last evPair
    // */
    EigenPairConstIterator end() const;

    // /**
    // Display the values and the vectors.
    // */
    void write(std::ostream &s) const;

  private:
    // /**
    // Each item in this vector is a pair of an eigen value and the
    // corresponding vector.
    // */
    std::vector<EigenPair> evPairs_;

    // /**
    // Number of negative eigen values.
    // */
    UInt neg_;

    // /**
    // Number of eigen values that are zero.
    // */
    UInt zero_;

    // /**
    // Number of positive eigen values.
    // */
    UInt pos_;
  };

}
#endif

// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
