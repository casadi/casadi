// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"
#include "Eigen.h"
#include "LinearFunction.h"
#include "Variable.h"

using namespace Minotaur;

#ifdef F77_FUNC

extern "C"
{
  void F77_FUNC(dsyevr, DSYEVR)(char *jobz, char *range, char *uplo, int *n,
                                double *a, int *lda, int *vl, int *vu,
                                int *il, int *iu, double *abstol, int *m,
                                double *w, double *z, int *ldz, int *isuppz,
                                double *work, int *lwork, int *iwork,
                                int *liwork, int *info );
}


EigenCalculator::EigenCalculator()
  :n_(0),
   A_(0),
   abstol_(1e-6)
{
}


EigenPtr EigenCalculator::findValues(ConstQuadraticFunctionPtr qf)
{
  EigenPtr ePtr = EigenPtr(); // NULL
  if (qf) {
    qf_ = qf;
    fillA_();
    findVectors_ = 'N'; // N for eigen values only, V for values and vectors.
    w_ = new double[n_];
    z_ = 0; 
    isuppz_ = 0;

    calculate_();
    ePtr = getEigen_();
    delete [] w_;
    delete [] A_;
    vars_.clear();
    indices_.clear();
  }
  return ePtr;
}


EigenPtr EigenCalculator::findVectors(ConstQuadraticFunctionPtr qf)
{
  EigenPtr ePtr = EigenPtr(); // NULL
  if (qf) {
    qf_ = qf;
    fillA_();
    findVectors_ = 'V'; // N for eigen values only, V for values and vectors.
    w_ = new double[n_];
    z_ = new double [n_*n_];
    isuppz_ = new int [2*n_];

    // fill z_ with 0. necessary because of a bug in lapack for n_ = 2
    // also bug in lapack for n_ = 1
    // We are ignore isuppz_ altogether now. lapack doesnt fill it correctly.
    // if (n_<=2) {
    //   std::fill(z_, z_+(n_*n_), 0.0);
    //   std::fill(isuppz_, isuppz_+2*n_, 1);
    // }

    calculate_();
    ePtr = getEigen_();
    delete [] w_;
    delete [] A_;
    delete [] z_;
    delete [] isuppz_;
    w_ = 0;
    A_ = 0;
    z_ = 0;
    isuppz_ = 0;
  }
  return ePtr;
}


void EigenCalculator::fillA_()
{
  UInt i,j;
  VarCountConstMap * qf_map = qf_->getVarMap();
  n_ = qf_map->size();

  // allocate space and create a map of what row/column corresponds to each
  // variable in qf.
  A_ = new double [n_*n_];
  std::fill(A_, A_+(n_*n_), 0);
  i=0;
  indices_.clear();
  vars_.clear();
  for (VarCountConstMap::const_iterator it = qf_map->begin();
      it != qf_map->end(); ++it, ++i) {
    indices_[it->first] = i;
    vars_.push_back(it->first);
  }

  // visit each term in qf and populate A_
  for (VariablePairGroupConstIterator it = qf_->begin(); it != qf_->end(); 
      ++it) {
    i = indices_[it->first.first];
    j = indices_[it->first.second];
    A_[i+j*n_] += 0.5*it->second;
    A_[j+i*n_] += 0.5*it->second;
  }
}


void EigenCalculator::getSumOfSquares (
    std::vector<LinearFunctionPtr> & p_terms, 
    std::vector<LinearFunctionPtr> & n_terms,
    std::vector<double> & p_const,
    std::vector<double> & n_const,
    LinearFunctionPtr & lin_terms, double & c,
    ConstQuadraticFunctionPtr qf, ConstLinearFunctionPtr lf)
{
  EigenPtr ePtr;
  LinearFunctionPtr lf2;
  ConstVariablePtr v_ptr;
  VarCountConstMap * qf_map;
  VarCountConstMap::const_iterator qf_it;
  double evalue, coeff;
  LinearFunctionPtr evector;

  assert(qf);

  c = 0;
  p_terms.clear();
  n_terms.clear();
  p_const.clear();
  n_const.clear();
  lin_terms.reset();
  ePtr = findVectors(qf);
  //ePtr->write(std::cout);

  // visit each term of lf and check if it exists in qf, if it does, add it to
  // a new linear function lf2.
  qf_map = qf_->getVarMap();
  lf2 = (LinearFunctionPtr) new LinearFunction();
  if (lf) {
    for (VariableGroupConstIterator it=lf->termsBegin(); it!=lf->termsEnd(); 
        ++it) {
      v_ptr = it->first;
      qf_it = qf_map->find(v_ptr);
      if (qf_it != qf_map->end()) {
        lf2->addTerm(v_ptr, it->second);
      }
    }

    lin_terms = lf->clone();
  } else {
    lin_terms = (LinearFunctionPtr) new LinearFunction();
  }

  // lf2 only has variables that occur in qf. it may not have all of qf's
  // variables.
  // now calculate b = R^(-1)Q'c. b_i = 
  for (EigenPairConstIterator it=ePtr->begin(); it!=ePtr->end(); ++it) {
    evector = it->second;
    evalue = it->first;
    if (fabs(evalue) > abstol_) {
      coeff = getDotProduct_(lf2, evector);
      //std::cout << "value = " << evalue << std::endl;
      //std::cout << "eigen vector = ";
      //evector->write(std::cout);
      //std::cout << "\ncoeff = " << coeff << std::endl;
      //std::cout << "final expression = " << evalue << "(";
      //evector->write(std::cout);
      //std::cout << " + " << 0.5/evalue*coeff << ")^2" << std::endl;
      // subtract from lin_terms
      (*lin_terms) += -1*coeff*evector;
      //std::cout << "remaining linear terms = ";
      //lin_terms->write(std::cout);
      //std::cout << "\n";

      // constant term
      c -= 0.25*coeff*coeff/evalue;

      if (evalue > 0) {
        (*evector) *= sqrt(evalue);
        p_terms.push_back(evector);
        p_const.push_back(0.5/sqrt(evalue)*coeff);
      } else {
        (*evector) *= sqrt(-evalue);
        n_terms.push_back(evector);
        n_const.push_back(-0.5/sqrt(-evalue)*coeff);
      }
    } else {
    }
  }
  //std::cout << "c = " << c << std::endl;
}


double EigenCalculator::getDotProduct_(ConstLinearFunctionPtr lf1, 
    ConstLinearFunctionPtr lf2)
{
  ConstLinearFunctionPtr tmp; // for exchange
  double d_product = 0;   // the dot product of coefficients of lf1 and lf2.
  ConstVariablePtr v_ptr;

  if (!lf1 || !lf2) {
    return 0;
  }

  if (lf1->getNumTerms() > lf2->getNumTerms()) {
    tmp = lf1;
    lf1 = lf2;
    lf2 = tmp;
  }

  for (VariableGroupConstIterator it1 = lf1->termsBegin(); 
      it1 != lf1->termsEnd(); ++it1) {
    v_ptr = it1->first;
    d_product += it1->second * lf2->getWeight(v_ptr);
  }

  return d_product;
}


void EigenCalculator::calculate_()
{
  char range = 'A'; // A for all eigen values/vectors, V for values in range (vl, vu]
                    // I for il-th through iu-th eigen value.
  char uplo = 'L';  // L for storing only lower triangular part of the matrix,
                    // U for upper.
  int lda = n_;     // The leading dimension of A, lda >= max(1,n)
  int vl=0, vu=0;   // Not used when range='A'
  int il=0, iu=0;   // Not used when ='A'
  int n=n_;
  int ldz = 1;      // The leading dimension of the array z_.
  double *work;     // Work space.
  int lwork;        // length of the array 'work'
  int *iwork;       // Work space.
  int liwork;       // length of the array 'iwork'
  int info = 0;     //  0 => successful exit
                    // -i => i-th argument has some problems
                    //  i => i off-diagonal elements of an intermediate
                    //       tridiagonal form did not converge to zero.
  m_=0;             // number of eigen values found.
  assert(n_);
  assert(A_);

  if (findVectors_ == 'V') {
    ldz = n_;
  }

  // get the required size.
  liwork = lwork = -1;
  work = new double[1];
  iwork = new int[1];
  work[0] = 0.0;
  iwork[0] = 0;

  F77_FUNC(dsyevr,DSYEVR)(&findVectors_, &range, &uplo, &n, A_, &lda, &vl, &vu,
      &il, &iu, &abstol_, &m_, w_, z_, &ldz, isuppz_, work, &lwork, iwork,
      &liwork, &info);
  assert(info==0);

  // allocate memory
  lwork = (int) work[0];
  liwork = iwork[0];
  delete [] work; 
  delete [] iwork;
  work = new double[lwork];
  iwork = new int[liwork];
  std::fill(work, work+lwork, 0.0);
  std::fill(iwork, iwork+liwork, 0);

  // do actual evaluation.
  F77_FUNC(dsyevr,DSYEVR)(&findVectors_, &range, &uplo, &n, A_, &lda, &vl,
      &vu, &il, &iu, &abstol_, &m_, w_, z_, &ldz, isuppz_, work, &lwork, iwork,
      &liwork, &info);
  assert(info==0);
  // free
  delete [] work; 
  delete [] iwork;
}


// If we calculated eigen vectors, then copy values and vectors, otherwise
// just copy values.
EigenPtr EigenCalculator::getEigen_()
{
  LinearFunctionPtr null_ptr = LinearFunctionPtr();
  LinearFunctionPtr lf;
  EigenPtr eigen = (EigenPtr) new Eigen();
  if (findVectors_=='V') {
    for (int i=0; i<m_; ++i) {
      lf = getLinearFunction_(i);
      if (fabs(w_[i]) < abstol_) {
        eigen->add(0, lf);
      } else {
        eigen->add(w_[i], lf);
      }
    }
  } else if (findVectors_=='N') {
    for (int i=0; i<m_; ++i) {
      if (fabs(w_[i]) < abstol_) {
        eigen->add(0, null_ptr);
      } else {
        eigen->add(w_[i], null_ptr);
      }
    }
  } else {
    assert (!"findVectors_ value is not supported!");
  }
  return eigen;
}


LinearFunctionPtr EigenCalculator::getLinearFunction_(const int i)
{
  LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
  // lapack has a bug in dstemr for when n_ = 2.
  // It also does not give correct values in isuppz_. We will ignore isuppz_
  // altogether.
  for (UInt j=0; j<n_; ++j) {
    lf->addTerm(vars_[j], z_[i*n_+j]);
  }
  return lf;
}


// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //


Eigen::Eigen()
  : neg_(0),
    zero_(0),
    pos_(0)
{
}


void Eigen::add(double value, LinearFunctionPtr e_vector)
{
  evPairs_.push_back(std::make_pair(value, e_vector));
  if (fabs(value) < 1e-7) {
    ++zero_;
  } else if (value < 0) {
    ++neg_;
  } else {
    ++pos_;
  }
}


UInt Eigen::numNegative() const
{
  return neg_;
}


UInt Eigen::numZero() const
{
  return zero_;
}


UInt Eigen::numPositive() const
{
  return pos_;
}


EigenPairConstIterator Eigen::begin() const
{
  return evPairs_.begin();
}


EigenPairConstIterator Eigen::end() const
{
  return evPairs_.end();
}


void Eigen::write(std::ostream &out) const
{
  for (EigenPairConstIterator it=evPairs_.begin(); it!=evPairs_.end(); ++it) { 
    out << "eigen value = " << it->first << std::endl;
    out << "eigen vector = ";
    if (it->second) {
      it->second->write(out);
    }
    out << std::endl;
  }
}


// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //

// EigenVector::EigenVector()
// {
//   lf_ = LinearFunctionPtr(); // NULL
// }

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
