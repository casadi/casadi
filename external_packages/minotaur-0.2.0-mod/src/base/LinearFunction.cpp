// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 


#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"
#include "LinearFunction.h"
#include "Operations.h"
#include "QuadraticFunction.h"
#include "Variable.h"

using namespace Minotaur;


LinearFunction::LinearFunction()
  : hasChanged_(true),
    tol_(1e-9)
{
  terms_.clear();
}


LinearFunction::LinearFunction(const double tol)
  : hasChanged_(true),
    tol_(tol)
{
  terms_.clear();
}


LinearFunction::LinearFunction(double *a, VariableConstIterator vbeg, 
    VariableConstIterator vend, double tol)
  : hasChanged_(true),
    tol_(tol)
{
  VariablePtr v;
  for (VariableConstIterator it=vbeg; it!=vend; ++it) {
    v = *it;
    addTerm(v, *a);
    ++a;
  }
}


LinearFunction::~LinearFunction()
{
  terms_.clear();
}


void LinearFunction::add(LinearFunctionPtr lf)
{
  if (lf) {
    for (VariableGroupConstIterator it = lf->terms_.begin(); it != lf->terms_.end(); 
         ++it) {
      incTerm(it->first, it->second);
    }
    hasChanged_ = true;
  }
}


void LinearFunction::addTerm(ConstVariablePtr var, const double a) 
{
  if (fabs(a) > tol_) {
    terms_.insert(std::make_pair(var, a));
    hasChanged_ = true;
  }
}


LinearFunctionPtr LinearFunction::clone() const
{
   LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
   for (VariableGroupConstIterator it = terms_.begin(); it != terms_.end(); 
       ++it) {
      lf->addTerm(it->first, it->second);
   }
   return lf;
}


LinearFunctionPtr LinearFunction::cloneWithVars(VariableConstIterator vbeg) 
  const
{
   LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
   for (VariableGroupConstIterator it = terms_.begin(); it != terms_.end(); 
       ++it) {
      lf->addTerm(*(vbeg+it->first->getIndex()), it->second);
   }
   return lf;
}


void LinearFunction::incTerm(ConstVariablePtr var, const double a)
{
  if (fabs(a) > tol_) {
    double nv = (terms_[var] += a);
    if (fabs(nv) < tol_) {
      terms_.erase(var);
    } 
    hasChanged_ = true;
  }
}


double LinearFunction::eval(const std::vector<double> &x) const
{
   return(InnerProduct(x, terms_));
}


double LinearFunction::eval(const double *x) const
{
  double value = 0;
  for (VariableGroupConstIterator it=terms_.begin(); it!=terms_.end(); ++it) {
    value += x[it->first->getIndex()] * it->second;
  }
  return value;
}


void LinearFunction::evalGradient(double *grad_f) const
{
  for (VariableGroupConstIterator it=terms_.begin(); it!=terms_.end(); ++it) {
    grad_f[it->first->getIndex()] += it->second;
  }
}


void LinearFunction::fillJac(double *values, int *) 
{
  double *v = values;
  for (UInt i=0; i<off_.size(); ++i, ++v) {
    *v += off_[i];
  }
}


void LinearFunction::computeBounds(double *l, double *u)
{
  double lb = 0.0;
  double ub = 0.0;
  double a;

  for (VariableGroupConstIterator it=terms_.begin(); it!=terms_.end(); ++it) {
    a = it->second;
    if (a>0) {
      lb += a*it->first->getLb();
      ub += a*it->first->getUb();
    } else {
      lb += a*it->first->getUb();
      ub += a*it->first->getLb();
    }
  }
  *l = lb;
  *u = ub;
}


void LinearFunction::getVars(VariableSet *vars)
{
  for (VariableGroupConstIterator it=terms_.begin(); it!=terms_.end(); ++it) {
    vars->insert(it->first);
  }
}


double LinearFunction::getWeight(ConstVariablePtr var) const
{
   VariableGroupConstIterator it = terms_.find(var);
   if (it == terms_.end()) {
     return 0.0;
   }
   return it->second;
}


bool LinearFunction::hasVar(ConstVariablePtr v) const
{
  return (terms_.find(v) != terms_.end());
}


void LinearFunction::multiply(double d)
{ 
  if (fabs(d) < 1e-7) {
    terms_.clear();
    hasChanged_ = true;
  } else {
    for (VariableGroupIterator it = terms_.begin(); it != terms_.end(); ++it) {
      it->second *= d;
    }
  }
  hasChanged_ = true;
}


VariableGroupConstIterator LinearFunction::termsBegin() const
{
  return terms_.begin();
}


VariableGroupConstIterator LinearFunction::termsEnd() const
{
  return terms_.end();
}


QuadraticFunctionPtr LinearFunction::square()
{
  QuadraticFunctionPtr qf = (QuadraticFunctionPtr) new QuadraticFunction();
  double d, dd;
  ConstVariablePtr v1;
  VariableGroupConstIterator it1, it2;

  // assume that no variable repeats in the linear terms.
  for (it1 = terms_.begin(); it1 != terms_.end(); ++it1) {
    d = it1->second;
    dd = 2.0*d;
    v1 = it1->first;
    qf->addTerm(v1, v1, d*d);
    it2 = it1;
    ++it2;
    for (; it2 != terms_.end(); ++it2) {
      qf->addTerm(v1, it2->first, dd*it2->second);
    }
  }
  return qf;
}


void LinearFunction::removeVar(VariablePtr v, double )
{
  terms_.erase(v);
  hasChanged_ = true;
}

void LinearFunction::clearAll()
{
  terms_.clear();
  off_.clear();
  hasChanged_ = true;
}


double LinearFunction::getFixVarOffset(VariablePtr v, double val)
{
  return val*getWeight(v);
}


void LinearFunction::prepJac(VarSetConstIter vbeg, VarSetConstIter vend)
{
  if (hasChanged_) {
    UInt i;
    VariableGroupConstIterator vit;

    i = 0;
    for (VarSetConstIter it=vbeg; it!=vend; ++it, ++i) {}
    off_.resize(i);

    i = 0;
    for (VarSetConstIter it=vbeg; it!=vend; ++it, ++i) {
      off_[i] = getWeight(*it);
    }
    hasChanged_ = false;
  }
}


void LinearFunction::operator+=(ConstLinearFunctionPtr l2)
{
  if (l2) {
    for (VariableGroupConstIterator it = l2->terms_.begin(); 
        it != l2->terms_.end(); ++it) {
      incTerm(it->first, it->second);
    }
    hasChanged_ = true;
  }
}


void LinearFunction::operator-=(ConstLinearFunctionPtr l2)
{
  if (l2) {
    for (VariableGroupConstIterator it = l2->terms_.begin(); 
        it != l2->terms_.end(); ++it) {
      incTerm(it->first, -1.*it->second);
    }
    hasChanged_ = true;
  }
}


void LinearFunction::operator*=(const double c)
{
  if (fabs(c) < 1e-7) {
    terms_.clear();
  } else {
    for (VariableGroupIterator it = terms_.begin(); it != terms_.end(); ++it) {
      it->second *= c;
    }
  }
  hasChanged_ = true;
}


void LinearFunction::write(std::ostream &s) const
{
   for (VariableGroupConstIterator it = terms_.begin(); it != terms_.end(); 
       ++it) {
     if (it->second >= 0.0) {
       s << "+ ";
     } else {
       s << "- ";
     }
     s << std::abs(it->second) << "*" << it->first->getName();
     s << " ";
   }
}


namespace Minotaur {
LinearFunctionPtr operator-(ConstLinearFunctionPtr l1, 
     ConstLinearFunctionPtr l2)
{
  LinearFunctionPtr lf = LinearFunctionPtr();  //NULL
  if (!l1 && !l2) {
    // do nothing.
  } else if (!l1) {
    lf = l2->clone();
    (*lf) *= -1.;
  } else if (!l2) {
    lf = l1->clone();
  } else {
    lf = l1->clone();
    for (VariableGroupConstIterator it = l2->terms_.begin(); 
        it != l2->terms_.end(); ++it) {
      lf->incTerm(it->first, -1*it->second);
    }
  }
  return lf;
}


LinearFunctionPtr operator+(ConstLinearFunctionPtr l1, 
     ConstLinearFunctionPtr l2)
{
  LinearFunctionPtr lf = LinearFunctionPtr();  //NULL
  if (!l1 && !l2) {
    // do nothing.
  } else if (!l1) {
    lf = l2->clone();
  } else if (!l2) {
    lf = l1->clone();
  } else {
    lf = l1->clone();
    for (VariableGroupConstIterator it = l2->terms_.begin(); 
        it != l2->terms_.end(); ++it) {
      lf->incTerm(it->first, it->second);
    }
  }
  return lf;
}


LinearFunctionPtr operator*(const double c, ConstLinearFunctionPtr l2)
{
  // creates a linear function even when c = 0.
  // Returns NULL if l2 is NULL.
  LinearFunctionPtr lf = LinearFunctionPtr();  //NULL
  if (l2 && fabs(c)>1e-7) {
    lf = (LinearFunctionPtr) new LinearFunction();
    for (VariableGroupConstIterator it = l2->terms_.begin(); 
        it != l2->terms_.end(); ++it) {
      lf->addTerm(it->first, c*it->second);
    }
  } else {
    // do nothing
  }
  return lf;
}


QuadraticFunctionPtr operator*(ConstLinearFunctionPtr l1, 
    ConstLinearFunctionPtr l2)
{
  QuadraticFunctionPtr qf = QuadraticFunctionPtr(); // NULL
  if (l1 && l2) {
    qf = (QuadraticFunctionPtr) new QuadraticFunction();
    for (VariableGroupConstIterator it1 = l1->terms_.begin(); 
        it1 != l1->terms_.end(); ++it1) {
      for (VariableGroupConstIterator it2 = l2->terms_.begin(); 
          it2 != l2->terms_.end(); ++it2) {
        qf->incTerm(it1->first, it2->first, it1->second * it2->second);
      }
    }
  } else {
    // do nothing
  }
  return qf;
} 

} // namespace Minotaur

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
