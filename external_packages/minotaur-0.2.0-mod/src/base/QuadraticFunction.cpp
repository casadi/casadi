// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file QuadraticFunction.cpp
 * \brief Define class QuadraticFunction for storing quadratic functions.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <queue>

#include "MinotaurConfig.h"
#include "HessianOfLag.h"
#include "LinearFunction.h"
#include "QuadraticFunction.h"
#include "Variable.h"

using namespace Minotaur;


QuadraticFunction::QuadraticFunction() 
  : etol_(1e-8),
    hCoeffs_(0),
    hFirst_(0),
    hOff_(0),
    hSecond_(0),
    terms_(), 
    varFreq_()
{
}


QuadraticFunction::QuadraticFunction(UInt nz, double *vals, UInt *irow,
                                     UInt *jcol, VariableConstIterator vbeg)
: etol_(1e-8),
  hCoeffs_(0),
  hFirst_(0),
  hOff_(0),
  hSecond_(0),
  terms_(), 
  varFreq_()
{
  VariablePtr v1, v2;
  for (UInt i=0; i<nz; ++i) {
    v1 = *(vbeg+irow[i]); 
    v2 = *(vbeg+jcol[i]); 
    if (v1==v2) {
      addTerm(v1, v2, 0.5*vals[i]);
    } else {
      addTerm(v1, v2, vals[i]);
    }
  }
}


QuadraticFunction::~QuadraticFunction() 
{
  if (hCoeffs_) {
    delete [] hCoeffs_;
    delete [] hOff_;
    delete [] hFirst_;
    delete [] hSecond_;
  }
  terms_.clear();
  varFreq_.clear();
}


QuadraticFunctionPtr QuadraticFunction::clone() const
{
   QuadraticFunctionPtr qf = (QuadraticFunctionPtr) new QuadraticFunction();
   qf->terms_.insert(terms_.begin(), terms_.end());
   qf->varFreq_.insert(varFreq_.begin(), varFreq_.end());
   return qf;
}


QuadraticFunctionPtr QuadraticFunction::cloneWithVars(VariableConstIterator vbeg) 
  const
{
  QuadraticFunctionPtr qf = (QuadraticFunctionPtr) new QuadraticFunction();
  for(VariablePairGroupConstIterator it = begin(); it != end(); ++it) {
    qf->addTerm(*(vbeg+it->first.first->getIndex()), 
        *(vbeg+it->first.second->getIndex()), it->second);
  }
  return qf;
}


double QuadraticFunction::eval(const std::vector<double> &x) const
{
   double sum = 0.0;
   for(VariablePairGroupConstIterator it = begin(); it != end(); ++it) {
      sum += it->second * x[it->first.first->getIndex()] * 
        x[it->first.second->getIndex()];
   }
   return sum;
}


double QuadraticFunction::eval(const double *x) const
{
   double sum = 0.0;
   for(VariablePairGroupConstIterator it = begin(); it != end(); ++it) {
      sum += it->second * x[it->first.first->getIndex()] * 
        x[it->first.second->getIndex()];
   }
   return sum;
}


void QuadraticFunction::evalGradient(const double *x, double *grad_f)
{
  assert (grad_f);
  if (x) {
    for(VariablePairGroupConstIterator it = terms_.begin(); it != terms_.end(); ++it) {
      grad_f[it->first.first->getIndex()] +=  it->second * 
        x[it->first.second->getIndex()];
      grad_f[it->first.second->getIndex()] +=  it->second * 
        x[it->first.first->getIndex()];
    }
  }

}


void QuadraticFunction::evalGradient(const std::vector<double> & x, 
    std::vector<double> & grad_f)
{
  for(VariablePairGroupConstIterator it = terms_.begin(); it != terms_.end(); ++it) {
    grad_f[it->first.first->getIndex()] +=  it->second * 
      x[it->first.second->getIndex()];
    grad_f[it->first.second->getIndex()] +=  it->second * 
      x[it->first.first->getIndex()];
  }
}


void QuadraticFunction::evalHessian(const double mul, const double *, 
                                    const LTHessStor *, double *values, int *)
{
  for (UInt i=0; i<terms_.size(); ++i) {
    values[hOff_[i]] += mul*hCoeffs_[i];
  }
}


void  QuadraticFunction::fillHessStor(LTHessStor *stor)
{
  VariablePtr v;
  UInt vindex;
  UInt j;
  std::deque<UInt> *inds;
  std::deque<UInt>::iterator it;

  prepHess();

  j = 0;
  for (UInt i=0; i<hStarts_.size()-1; ++i) {
    vindex = hSecond_[hStarts_[i]];
    while (vindex!=stor->rows[j]->getIndex()) {
      ++j;
    }
    inds = stor->colQs+j;

    it = inds->begin();
    for (UInt i2=hStarts_[i]; i2<hStarts_[i+1]; ++i2) {
      while (it!=inds->end() && (*it)<hFirst_[i2]) {
        ++it;
      }
      if (it==inds->end()) {
        inds->push_back(hFirst_[i2]);
      } else if ((*it)!=hFirst_[i2]) {
        it = inds->insert(it,hFirst_[i2]);
      } else {
      }
    }
  }
}


void QuadraticFunction::finalHessStor(const LTHessStor *stor)
{
  UInt vindex;
  UInt *inds;
  UInt nz = 0;
  UInt j;
  UInt off;

  j = 0;
  for (UInt i=0; i<hStarts_.size()-1; ++i) {
    vindex = hSecond_[hStarts_[i]];
    while (vindex!=stor->rows[j]->getIndex()) {
      ++j;
    }
    inds = stor->cols+stor->starts[j];
    off  = stor->starts[j];

    for (UInt i2=hStarts_[i]; i2<hStarts_[i+1]; ++i2) {
      while (*inds != hFirst_[i2]) {
        ++inds;
        ++off;
        assert(off < stor->starts[j+1]);
      }
      hOff_[nz] = off;
      ++nz;
    }
  }
  assert(nz==terms_.size());
}


void QuadraticFunction::fillJac(const double *x, double *values, int *) 
{
  UInt i=0;
  for (VariablePairGroupConstIterator it = terms_.begin(); it != terms_.end(); 
      ++it) {
    values[jacOff_[i]] += it->second * x[jacInd_[i]];
    ++i;
    values[jacOff_[i]] += it->second * x[jacInd_[i]];
    ++i;
  }
}


void QuadraticFunction::getVars(VariableSet *vars)
{
  for (VarIntMap::const_iterator it=varFreq_.begin(); it!= varFreq_.end();
       ++it) {
    vars->insert(it->first);
  }
}


void QuadraticFunction::mult(double c) {
  if (fabs(c) < 1e-7) {
    terms_.clear();
    varFreq_.clear();
  } else {
    for (VariablePairGroupIterator it = terms_.begin(); it != terms_.end(); 
        ++it) {
      it->second *= c;
    }
  }
}


void QuadraticFunction::addTerm(VariablePair vp, const double weight)
{
  assert (vp.first->getId() <= vp.second->getId());
  if (fabs(weight) >= etol_) {
    terms_.insert(std::make_pair(vp, weight));
    varFreq_[vp.first] += 1;
    varFreq_[vp.second] += 1;
  }
}


void QuadraticFunction::addTerm(ConstVariablePtr v1, ConstVariablePtr v2, 
    const double weight)
{
  VariablePair vp = (v1->getId() < v2->getId()) ? std::make_pair(v1, v2) :
    std::make_pair(v2, v1);
  addTerm(vp, weight);
}


void QuadraticFunction::incTerm(ConstVariablePair vp, const double a)
{
  if (fabs(a) > etol_) {
    VariablePairGroupIterator it = terms_.find(vp);
    if (it == terms_.end()) {
      varFreq_[vp.first] += 1;
      varFreq_[vp.second] += 1;
      terms_.insert(std::make_pair(vp, a));
    } else {
      double nv = (it->second + a);
      if (fabs(nv) < etol_) {
        terms_.erase(vp);
        std::map<ConstVariablePtr, UInt>::iterator vit;
        vit = varFreq_.find(vp.first);
        vit->second -= 1;
        if (vit->second < 1) {
          varFreq_.erase(vit);
        }
        vit = varFreq_.find(vp.second);
        vit->second -= 1;
        if (vit->second < 1) {
          varFreq_.erase(vit);
        }
      } else {
        it->second = nv;
      }
    }
  }
}


void QuadraticFunction::incTerm(ConstVariablePtr v1, ConstVariablePtr v2,
    const double a)
{
  VariablePair vp = (v1->getId() < v2->getId()) ? std::make_pair(v1, v2)
    : std::make_pair(v2, v1);
  incTerm(vp, a);
}


UInt QuadraticFunction::getNumTerms() const
{ 
  return terms_.size(); 
}


UInt QuadraticFunction::getNumVars() const
{ 
  return varFreq_.size(); 
}


VarCountConstMap * QuadraticFunction::getVarMap() const
{
  return &varFreq_;
}


double QuadraticFunction::getWeight(ConstVariablePair & vp) 
{ 
  VariablePairGroupConstIterator it = terms_.find(vp);
  if (it == terms_.end()) {
    return 0.0;
  }
  return it->second;
}


double QuadraticFunction::getWeight(ConstVariablePtr v1, ConstVariablePtr v2) 
{ 
  ConstVariablePair vp = (v1->getId() < v2->getId()) ? std::make_pair(v1, v2)
    : std::make_pair(v2, v1);
  return getWeight(vp);
}


int QuadraticFunction::getFreq(ConstVariablePtr v1)
{
  return varFreq_[v1];
}


bool QuadraticFunction::hasVar(ConstVariablePtr v) const
{
  return (varFreq_.find(v) != varFreq_.end());
}


void QuadraticFunction::removeVar(VariablePtr v, double val, 
    LinearFunctionPtr lf) 
{
  for (VariablePairGroupIterator it = terms_.begin(); it != terms_.end();) {
    if (it->first.first == v && it->first.first == it->first.second) {
      terms_.erase(it++);
    } else if (it->first.first == v) {
      lf->incTerm(it->first.second, it->second * val);
      terms_.erase(it++);
    } else if (it->first.second == v) {
      lf->incTerm(it->first.first, it->second * val);
      terms_.erase(it++);
    } else {
      ++it;
    }
  }
}


void QuadraticFunction::prepJac(VarSetConstIter vbeg, VarSetConstIter vend)
{
  UInt i;
  VarIntMap omap;
  std::map<ConstVariablePtr, UInt>::iterator vit;

  i=0;
  for (VarSetConstIter it=vbeg; it!=vend; ++it, ++i) {
    vit = varFreq_.find(*it);
    if (vit!=varFreq_.end()) {
      omap[*it] = i;
    }
  }
  jacOff_.resize(2*terms_.size());
  jacInd_.resize(2*terms_.size());

  i=0;
  for (VariablePairGroupIterator it = terms_.begin(); it != terms_.end(); 
      ++it) {
    jacOff_[i] = omap[it->first.first];
    jacInd_[i] = it->first.second->getIndex();
    ++i;
    jacOff_[i] = omap[it->first.second];
    jacInd_[i] = it->first.first->getIndex();
    ++i;
  }
}


void QuadraticFunction::prepHess()
{
  UInt nterms    = terms_.size();
  UInt *first    = new UInt[nterms];
  UInt *f        = first;
  UInt *second   = new UInt[nterms];
  UInt *s        = second;
  double *coeffs = new double[nterms];
  double *c      = coeffs;
  UInt prev;

  if (hFirst_) {
    delete [] hCoeffs_;
    delete [] hFirst_;
    delete [] hOff_;
    delete [] hSecond_;
    hFirst_ = 0;
    hSecond_ = 0;
    hCoeffs_ = 0;
  }

  hStarts_.clear();
  hStarts_.push_back(0);

  if (0==nterms) {
    // delete not necessary, but valgrind complains.
    delete [] first;
    delete [] second;
    delete [] coeffs;
    return;
  }

  for (VariablePairGroupIterator it = terms_.begin(); it != terms_.end();
       ++it, ++f, ++s, ++c) {
    *f = it->first.first->getIndex();
    *s = it->first.second->getIndex();
    if (*f == *s) {
      *c = 2.0*it->second;
    } else {
      *c = it->second;
    }
  }

  // remember, we need lower triangle.
  sortLT_(nterms, first, second, coeffs);

  prev = second[0];
#if DEBUG
  for (UInt i=0; i+1<nterms; ++i) {
    assert(second[i]<second[i+1] || (second[i] == second[i+1] && first[i]<first[i+1]));
  }
#endif 
  for (UInt i=0; i<nterms; ++i) {
    if (second[i]>prev) {
      prev = second[i];
      hStarts_.push_back(i);
    }
  }
  hStarts_.push_back(nterms);

  hFirst_  = first;
  hSecond_ = second;
  hCoeffs_ = coeffs;
  hOff_    = new UInt[nterms];
}


void QuadraticFunction::sortLT_(UInt n, UInt *f, UInt *s, double *c)
{
  UInt l = 0;
  UInt r = n-1;
  UInt pivot = r/2;
  UInt sp, fp;
  double dtmp;
  UInt itmp;

  if (n<2) {
    return;
  }

  sp = s[pivot];
  fp = f[pivot];
  while (l<r) {
    while (l<=pivot) {
      if (s[l] > sp) {
        break;
      } else if (s[l] == sp && f[l] >= fp) {
        break;
      }
      ++l;
    }
    while (r>=pivot) {
      if (s[r] < sp) {
        break;
      } else if (s[r] == sp && f[r] <= fp) {
        break;
      }
      --r;
    }
    if (l<r) {
      itmp = s[l];
      s[l] = s[r];
      s[r] = itmp;

      itmp = f[l];
      f[l] = f[r];
      f[r] = itmp;

      dtmp = c[l];
      c[l] = c[r];
      c[r] = dtmp;
      if (l==pivot) {
        pivot=r;
        ++l;
      } else if (r==pivot) {
        pivot=l;
        --r;
      }
    }
  }
  sortLT_(pivot, f, s, c);
  sortLT_(n-pivot-1, f+pivot+1, s+pivot+1, c+pivot+1);
}


void QuadraticFunction::subst(VariablePtr out, VariablePtr in, double rat)
{
  std::queue <std::pair<VariablePair, double> > newterms;
  std::pair <VariablePair, double> vpg;
  std::map<ConstVariablePtr, UInt>::iterator vit;

  vit = varFreq_.find(out);
  if (vit==varFreq_.end()) {
    return;
  }

  for (VariablePairGroupIterator it = terms_.begin(); it != terms_.end();){
    if (it->first.first == out || it->first.second==out) {
      vpg = *it;
      if (vpg.first.first==out) {
        vpg.first.first=in;
      }
      if (vpg.first.second==out) {
        vpg.first.second=in;
      }
      vpg.second *= rat;
      newterms.push(vpg);
      varFreq_[it->first.first]  -= 1;
      varFreq_[it->first.second] -= 1;
      terms_.erase(it++);
    } else {
      ++it;
    }
  }
  varFreq_.erase(vit);

  while (!newterms.empty()) {
    vpg = newterms.front();
    incTerm(vpg.first, vpg.second);
    newterms.pop();
  }
}


double QuadraticFunction::getFixVarOffset(VariablePtr v, double val)
{
  if (fabs(val)>1e-7) {
    ConstVariablePair vp = std::make_pair(v, v);
    return val * val * getWeight(vp);
  } else {
    return 0.0;
  }
}


VarIntMapConstIterator QuadraticFunction::varsBegin() const
{
  return varFreq_.begin();
}


VarIntMapConstIterator QuadraticFunction::varsEnd() const
{
  return varFreq_.end();
}


void QuadraticFunction::operator+=(ConstQuadraticFunctionPtr q2)
{
  if (q2) {
    for (VariablePairGroupConstIterator it = q2->terms_.begin(); 
        it != q2->terms_.end(); ++it) {
      incTerm(it->first, it->second);
    }
  }
}


void QuadraticFunction::operator*=(const double c)
{
  mult(c);
}


void QuadraticFunction::write(std::ostream &s) const
{
   for (VariablePairGroupConstIterator it = terms_.begin(); it != terms_.end(); 
       ++it) {
     if (it->second > 0.0) {
       s << "+ ";
     } 
     else {
       s << "- ";
     }
     s << std::abs(it->second) << "*" << it->first.first->getName();
     s << "*" << it->first.second->getName() << " ";
   }
}


namespace Minotaur {
QuadraticFunctionPtr operator-(ConstQuadraticFunctionPtr q1, 
     ConstQuadraticFunctionPtr q2)
{
  QuadraticFunctionPtr qf = QuadraticFunctionPtr();  //NULL
  if (!q1 && !q2) {
    // do nothing.
  } else if (!q1) {
    qf = -1. * q2;
  } else if (!q2) {
    qf = q1->clone();
  } else {
    qf = q1->clone();
    for (VariablePairGroupConstIterator it = q2->terms_.begin(); 
        it != q2->terms_.end(); it++) {
      qf->incTerm(it->first, -1*it->second);
    }
  }
  return qf;
}


QuadraticFunctionPtr operator+(ConstQuadraticFunctionPtr q1, 
     ConstQuadraticFunctionPtr q2)
{
  QuadraticFunctionPtr qf = QuadraticFunctionPtr();  //NULL
  if (!q1 && !q2) {
    // do nothing.
  } else if (!q1) {
    qf = q2->clone();
  } else if (!q2) {
    qf = q1->clone();
  } else {
    qf = q1->clone();
    for (VariablePairGroupConstIterator it = q2->terms_.begin(); 
        it != q2->terms_.end(); it++) {
      qf->incTerm(it->first, it->second);
    }
  }
  return qf;
}


QuadraticFunctionPtr operator*(const double c, ConstQuadraticFunctionPtr q2)
{
  // creates a linear function even when c = 0.
  // Returns NULL if q2 is NULL.
  QuadraticFunctionPtr qf = QuadraticFunctionPtr();  //NULL
  if (q2) {
    qf = (QuadraticFunctionPtr) new QuadraticFunction();
    for (VariablePairGroupConstIterator it = q2->terms_.begin(); 
        it != q2->terms_.end(); it++) {
      qf->addTerm(it->first, c*it->second);
    }
  } else {
    // do nothing
  }
  return qf;
}

}


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
