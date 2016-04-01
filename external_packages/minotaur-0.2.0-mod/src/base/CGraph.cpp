//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file CGraph.cpp
 * \brief Define class CGraph for storing computational graph of a nonlinear
 * function.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>
#include <stack>

#include "MinotaurConfig.h"
#include "CGraph.h"
#include "CNode.h"
#include "HessianOfLag.h"
#include "VarBoundMod.h"
#include "Variable.h"

using namespace Minotaur;

CGraph::CGraph()
  : aNodes_(0),
    changed_(false),
    hInds_(0),
    hNnz_(0),
    hOffs_(0),
    hStarts_(0),
    gOffs_(0),
    oNode_(0)
{
  dq_.clear();
  varNode_.clear();
  vq_.clear();
}


CGraph::~CGraph()
{
  varNode_.clear();
  vq_.clear();
  dq_.clear();
  for (UInt i=0; i<aNodes_.size(); ++i) {
    delete aNodes_[i];
  }
  aNodes_.clear();
}


void CGraph::addConst(const double eps, int &)
{
  CNode *n = newNode(eps);

  if (oNode_) {
    oNode_ = newNode(OpPlus, oNode_, n);
  } else {
    oNode_ = n;
  }
  changed_ = true;
  finalize();
}


NonlinearFunctionPtr CGraph::clone(int *err) const
{
  return clone_(err);
}


CGraphPtr CGraph::clone_(int *err) const
{
  CGraphPtr cg = (CGraphPtr) new CGraph();
  std::map<const CNode*, CNode*>nnmap;
  std::map<const CNode*, CNode*>::iterator mit;
  const CNode *const_node;
  CNode *node = 0;
  VariablePtr v;

  for (UInt i=0; i<aNodes_.size(); ++i) {
    const_node = aNodes_[i];
#if DEBUG
    mit = nnmap.find(const_node);
    assert (mit==nnmap.end());
#endif
    node = const_node->clone();
    cg->aNodes_.push_back(node);
    nnmap.insert(std::pair<const CNode*, CNode*>(const_node, node));
    if (OpVar==node->getOp()) {
      v= getVar(const_node);
      node->setV(v);
      cg->varNode_.insert(std::pair<ConstVariablePtr, CNode*> (v, node));
      cg->vars_.insert(v);
    }
  }
  for (UInt i=0; i<aNodes_.size(); ++i) {
    aNodes_[i]->copyParChild(cg->aNodes_[i], &nnmap);
  }

  mit = nnmap.find(oNode_);
  assert(mit!=nnmap.end());
  cg->oNode_ = mit->second;

  cg->finalize();

  cg->hInds_ = hInds_;
  cg->hNnz_ = hNnz_;
  cg->hOffs_ = hOffs_;
  cg->hStarts_ = hStarts_;
  cg->gOffs_ = gOffs_;

  cg->changed_ = true;
  *err = 0;
  return cg;
}


NonlinearFunctionPtr CGraph::cloneWithVars(VariableConstIterator vbeg,
                                           int *) const
{
  CGraphPtr cg = (CGraphPtr) new CGraph();
  std::map<const CNode*, CNode*>nnmap;
  std::map<const CNode*, CNode*>::iterator mit;
  const CNode *const_node;
  CNode *node = 0;
  VariablePtr v;

  for (UInt i=0; i<aNodes_.size(); ++i) {
    const_node = aNodes_[i];
#if DEBUG
    mit = nnmap.find(const_node);
    assert (mit==nnmap.end());
#endif
    node = const_node->clone();
    cg->aNodes_.push_back(node);
    nnmap.insert(std::pair<const CNode*, CNode*>(const_node, node));
    if (OpVar==node->getOp()) {
      v = *(vbeg+const_node->getV()->getIndex());
      node->setV(v);
      cg->varNode_.insert(std::pair<ConstVariablePtr, CNode*> (v, node));
      cg->vars_.insert(v);
    }
  }
  for (UInt i=0; i<aNodes_.size(); ++i) {
    aNodes_[i]->copyParChild(cg->aNodes_[i], &nnmap);
  }

  mit = nnmap.find(oNode_);
  assert(mit!=nnmap.end());
  cg->oNode_ = mit->second;

  cg->finalize();

  cg->hInds_ = hInds_;
  cg->hNnz_ = hNnz_;
  cg->hOffs_ = hOffs_;
  cg->hStarts_ = hStarts_;
  cg->gOffs_ = gOffs_;

  return cg;
}


void CGraph::computeBounds(double *lb, double *ub, int *error)
{
  *error = 0;
  for (CNodeQ::iterator it=vq_.begin(); it!=vq_.end(); ++it) {
    (*it)->updateBnd(error);
  }
  for (CNodeQ::iterator it=dq_.begin(); it!=dq_.end(); ++it) {
    (*it)->updateBnd(error);
  }
  *lb = oNode_->getLb();
  *ub = oNode_->getUb();
}


double CGraph::eval(const double *x, int *error)
{
  for (CNodeQ::iterator it=vq_.begin(); it!=vq_.end(); ++it) {
    (*it)->eval(x, error);
  }

  for (CNodeQ::iterator it=dq_.begin(); it!=dq_.end(); ++it) {
    (*it)->eval(x, error);
    if (0!=*error) {
      break;
    }
  }
  return oNode_->getVal();
}


void CGraph::evalGradient(const double *x, double *grad_f, int *error)
{
  eval(x, error);
  if (*error>0) {
    return;
  }
  grad_(error);
  if (*error>0) {
    return;
  }
  for (CNodeQ::iterator it=vq_.begin(); it!=vq_.end(); ++it) {
    grad_f[(*it)->getV()->getIndex()] += (*it)->getG();
  }
}


void CGraph::evalHessian(double mult, const double *x, 
                         const LTHessStor *, double *values, int *error)
{
  UInt i = 0;
  UInt vind;
  VariablePtr v;
  UInt nz = 0;
  bool use2 = true;
  std::stack<CNode *> st2;
  double thresh = vq_.size();

  thresh = thresh*(thresh-1.0)/2.0 * 0.75;

  if (hNnz_>thresh) {
    use2 = false;
  }

  // always eval. We do not assume that evaluations of x are already
  // available. It creates a big mess and doesn't save much.
  eval(x, error);
  if (true==use2) {
    for (CNodeQ::iterator it=dq_.begin(); it!=dq_.end(); ++it) {
      (*it)->setB(false);
      (*it)->setGi(0.0);
      (*it)->setH(0.0);
      (*it)->setTempI(0);
    }
    for (CNodeQ::iterator it=vq_.begin(); it!=vq_.end(); ++it) {
      (*it)->setB(false);
      (*it)->setGi(0.0);
      (*it)->setH(0.0);
      (*it)->setTempI(0);
    }
    grad_(error);
  }

  //std::cout << std::endl << " evaling hessian ";
  //write(std::cout);
  //std::cout << std::endl;
  for (VarNodeMap::iterator it=varNode_.begin(); it!=varNode_.end(); ++it, 
       ++i) {
    if (hStarts_[i]<hStarts_[i+1]) {
      v = it->first;
      vind = v->getIndex();
      if (true==use2) {
        fwdGrad2_(&st2, it->second); 
        revHess2_(&st2, mult, vind, values, &nz, error);
      } else {
        fwdGrad_(it->second); 
        revHess_(error);
        for (CNodeQ::iterator it2=vq_.begin(); 
             it2!=vq_.end() && (*it2)->getV()->getIndex() <= vind; ++it2) {
          //std::cout << std::endl << v->getName() << " " 
          //<< (*it2)->getV()->getName() << " " << (*it2)->getH() 
          //<< " hinds = " << hInds_[nz] << " my index = " 
          //<< (*it2)->getV()->getIndex() << std::endl;
          if (hInds_[nz] == (*it2)->getV()->getIndex()) {
            values[hOffs_[nz]] += mult * (*it2)->getH();
            ++nz;
            if (nz == hStarts_[i+1]) {
              break;
            }
          }
        }
      }
    }
  }
  //for (UInt i=0; i<hNnz_; ++i) {
  //  std::cout << "h[" << i << "] = " << values[i] << std::endl;
  //}
  //exit(0);
}


void CGraph::fillHessInds_(CNode *node, UIntQ *inds)
{
  CNode **c1=0, **c2=0;

  for (CNodeQ::iterator it=dq_.begin(); it!=dq_.end(); ++it) {
    (*it)->setTempI(0);
    (*it)->setB(false);
  }
  for (CNodeQ::iterator it=vq_.begin(); it!=vq_.end(); ++it) {
    (*it)->setTempI(0);
    (*it)->setB(false);
  }
  node->setTempI(1);
  //std::cout << node->getV()->getName() << ": ";

  for (CNodeQ::iterator it=dq_.begin(); it!=dq_.end(); ++it) {
    switch((*it)->numChild()) {
    case (0):
      break;
    case (1):
      (*it)->setTempI((*it)->getL()->getTempI());
      break;
    case (2):
      if ((*it)->getL()->getTempI()>0 || (*it)->getR()->getTempI()>0) {
        (*it)->setTempI(1);
      }
      break;
    default:
      c1 = (*it)->getListL();
      c2 = (*it)->getListR();
      while (c1<c2) {
        if ((*c1)->getTempI()>0) {
          (*it)->setTempI(1);
          break;
        }
        ++c1;
      }
    }
  }


  for (CNodeQ::reverse_iterator it=dq_.rbegin(); it!=dq_.rend(); ++it) {
    (*it)->propHessSpa();
  }

  for (VarNodeMap::iterator it=varNode_.begin(); it!=varNode_.end(); ++it) {
    if (it->second->getB()==true) {
      //std::cout << (it)->second->getV()->getName() << " ";
      inds->push_back(it->first->getIndex());
    }
  }
  //std::cout << "\n\n"; 
}


void CGraph::fillHessInds2_(CNode *node, UIntQ *inds)
{
  std::stack<CNode *>st, st2;
  CNode *n, *n2;

  node->setTempI(1);
  st.push(node);
  st2.push(node);
  while (!st.empty()) {
    n = st.top();
    st.pop();
    switch(n->numPar()) {
    case 0:
      assert(n==oNode_);
      break;
    case 1:
      n2 = n->getUPar();
      if (0==n2->getTempI()) {
        n2->setTempI(1);
        st.push(n2);
        st2.push(n2);
      }
    default:
      {
      CQIter2 *it = n->getParB();
      while (it) {
        n2 = it->node;
        if (0==n2->getTempI()) {
          n2->setTempI(1);
          st.push(n2);
          st2.push(n2);
        }
        it = it->next;
      }
      }
    }
  }

  CNodeRSet nset;
  CNodeRSet::iterator sit; 
  nset.insert(oNode_);
  while (!nset.empty()) {
    sit = nset.begin();
    n = *sit;
    nset.erase(sit);
    n->propHessSpa2(&nset);
    if (OpVar==n->getOp() && n->getB()==true) {
      inds->push_back(n->getV()->getIndex());
    }
    n->setB(false);
    n->setTempI(0);
    //std::cout << nset.size() << std::endl;
  }
  while (!st2.empty()) {
    n = st2.top();
    st2.pop();
    n->setB(false);
    n->setTempI(0);
  }
}


void CGraph::fillHessStor(LTHessStor *stor)
{
  // for each variable in this nonlinear function, find sparsity pattern of 
  // \grad (e_i . \grad(f)).

  VariablePtr v;
  UIntQ *inds = new UIntQ();
  UIntQ::iterator it2, it_st;
  VariablePtr *stor_rows = stor->rows;
  UIntQ *st_inds = stor->colQs;
  UInt vind;
  bool use2 = true;

  hInds_.clear();
  hOffs_.clear();
  hStarts_.clear();
  hStarts_.reserve(varNode_.size()+1);
  hStarts_.push_back(0);
  hNnz_ = 0;

  if (true == changed_) {
    simplifyDq_();
    changed_ = false;
  }

  if (use2) {
    for (CNodeQ::iterator it=dq_.begin(); it!=dq_.end(); ++it) {
      (*it)->setTempI(0);
      (*it)->setB(false);
      if ((*it)->numPar()==1) {
        assert((*it)->getUPar());
      }
    }
    for (CNodeQ::iterator it=vq_.begin(); it!=vq_.end(); ++it) {
      (*it)->setTempI(0);
      (*it)->setB(false);
    }
  }


  for (VarNodeMap::iterator it=varNode_.begin(); it!=varNode_.end(); ++it) {
    v = it->first;
    vind = v->getIndex();
    inds->clear();
    //std::cout << "variable " << v->getName() << std::endl;
    if (use2) {
      fillHessInds2_(it->second, inds);
    } else {
      fillHessInds_(it->second, inds);
    }
    //std::cout << "size = " << inds->size() << std::endl;

    while (*stor_rows != v) {
      ++stor_rows;
      ++st_inds;
    }
    it_st = st_inds->begin();

    // copy the indices from inds into st_inds
    for (it2 = inds->begin(); it2 != inds->end() && vind >= *it2; ++it2) {
      while (true) {
        if (it_st == st_inds->end()) {
          st_inds->push_back(*it2);
          it_st = st_inds->end();
          break;
        } else if (*it_st > *it2) {
          it_st = st_inds->insert(it_st, *it2);
          break;
        } else if (*it_st == *it2) {
          break;
        } else {
          ++it_st;
        }
      }
      hInds_.push_back(*it2);
      ++hNnz_;
    }
    hStarts_.push_back(hNnz_);
  }
  inds->clear();
  delete inds;
}


void CGraph::fillJac(const double *x, double *values, int *error)
{
  UInt *goff = &gOffs_[0];

  *error = 0;
  eval(x, error);
  if (*error>0) {
    return;
  }
  grad_(error);
  if (*error>0) {
    return;
  }

  for (VarNodeMap::iterator it=varNode_.begin(); it!=varNode_.end(); 
       ++it,++goff) {
    values[*goff] += it->second->getG();
  }
}


void CGraph::finalHessStor(const LTHessStor *stor)
{
  UInt *st_cols;
  UInt *st_starts = stor->starts;
  VariablePtr *stor_rows = stor->rows;
  VariablePtr v;
  UInt i, j, ind2, off;

  // we need to fill offsets.
  hOffs_.clear();
  hOffs_.reserve(hNnz_);

  // visit all indices in hInds_ and see what position (i) do they appear in
  // stor. Then put i in hOffs_.
  i = 0;
  off = 0;
  for (VarNodeMap::iterator it=varNode_.begin(); it!=varNode_.end();
       ++it, ++i) {
    // find v in stor_rows.
    v = it->first;
    while (*stor_rows != v) {
      ++st_starts;
      ++stor_rows;
    }

    off = *st_starts;
    st_cols = stor->cols+off;
    for (j=hStarts_[i]; j!=hStarts_[i+1]; ++j) {
      ind2 = hInds_[j];
      // now find ind1, ind2 in stor.
      while (*st_cols != ind2) {
        ++st_cols;
        ++off;
      }
      hOffs_.push_back(off);
    }
  }
  assert (hNnz_ == hOffs_.size());
}


void CGraph::finalize()
{
  std::stack<CNode *>st;
  CNode *n1, *lchild, *rchild;
  CNode **p1, **p2;
  UInt id = 0;
  assert(oNode_);
  st.push(oNode_);


  vq_.clear();
  dq_.clear();

  for (VarNodeMap::iterator it=varNode_.begin(); it!=varNode_.end(); ++it) {
    vq_.push_back(it->second);
    it->second->setId(0);
  }

  for (UInt i=0; i<aNodes_.size(); ++i) {
    aNodes_[i]->setTempI(aNodes_[i]->numChild()+1);
  }

  id = 1;
  while(!st.empty()) {
    n1 = st.top();
    if (0==n1->numChild()) {
      n1->setTempI(0); 
      st.pop();
    } else if (1==n1->numChild()) {
      lchild = n1->getL();
      if (0==lchild->getTempI()) {
        n1->setTempI(0); 
        st.pop();
        dq_.push_back(n1);
        n1->setId(id); 
        ++id;
      } else {
        st.push(lchild);
      }
    } else if (2==n1->numChild()) {
      lchild = n1->getL();
      if (0==lchild->getTempI()) {
        rchild = n1->getR();
        if (0==rchild->getTempI()) {
          n1->setTempI(0); 
          st.pop();
          dq_.push_back(n1);
          n1->setId(id);
          ++id;
        } else {
          st.push(rchild);
        }
      } else {
        st.push(lchild);
      }
    } else {
      UInt i = n1->numChild()+1-n1->getTempI();
      p1 = n1->getListL()+i;
      p2 = n1->getListR();
      while (p1<p2) {
        if (0!=(*p1)->getTempI()) {
          st.push(*p1);
          n1->setTempI(n1->numChild()-i);
          break;
        }
        ++i;
        ++p1;
      }
      if (p1==p2) {
        n1->setTempI(0); 
        st.pop();
        dq_.push_back(n1);
        n1->setId(id);
        ++id;
      }
    }
  }
  //std::cout << std::endl
  //          << "dq_ has " << dq_.size() << " elements." << std::endl
  //          << "vq_ has " << vq_.size() << " elements." << std::endl
  //          << "aNodes_ has " << aNodes_.size() << " elements." << std::endl;
  //std::cout << std::endl << "dq_:" << std::endl;
  //for (CNodeQ::iterator it=dq_.begin(); it!=dq_.end(); ++it) {
  //  (*it)->writeSubExp(std::cout);
  //  std::cout << std::endl;
  //}

}


void CGraph::fwdGrad_(CNode *node) 
{
  for (CNodeQ::iterator it=dq_.begin(); it!=dq_.end(); ++it) {
    (*it)->setGi(0.0);
    (*it)->setG(0.0);
    (*it)->setH(0.0);
  }
  for (CNodeQ::iterator it=vq_.begin(); it!=vq_.end(); ++it) {
    (*it)->setGi(0.0);
    (*it)->setG(0.0);
    (*it)->setH(0.0);
  }
  node->setGi(1.0);

  for (CNodeQ::iterator it=dq_.begin(); it!=dq_.end(); ++it) {
    (*it)->fwdGrad();
  }
}


void CGraph::fwdGrad2_(std::stack<CNode *> *st2, CNode *node) 
{
  
  CNode *n, *n2;
  CNodeSet nset; // different compare method.
  CNodeSet::iterator sit; 
  CNodeQ qq;

  while(!(st2->empty())) {
    st2->pop();
  }

  node->setGi(1.0);

  node->setTempI(1);
  nset.insert(node);
  while (!nset.empty()) {
    sit = nset.begin();
    n = *sit;
    nset.erase(sit);
    qq.push_back(n);
    switch(n->numPar()) {
    case 0:
      assert(n==oNode_);
      break;
    case 1:
      n2 = n->getUPar();
      if (0==n2->getTempI()) {
        n2->setTempI(1);
        nset.insert(n2);
      }
    default:
      {
      CQIter2 *it = n->getParB();
      while (it) {
        n2 = it->node;
        if (0==n2->getTempI()) {
          n2->setTempI(1);
          nset.insert(n2);
        }
        it = it->next;
      }
      }
    }
  }

  for (CNodeQ::iterator it=qq.begin(); it!=qq.end(); ++it) {
    n = *it;
    n->fwdGrad();
    st2->push(n);
  }
  qq.clear();
}


double CGraph::getFixVarOffset(VariablePtr, double)
{
  return 0.0;
}


UInt CGraph::getNumNodes()
{
  return aNodes_.size();
}


NonlinearFunctionPtr CGraph::getPersp(VariablePtr z, double eps, int *err) const
{
  CNode *znode = 0;
  CNode *dnode = 0;
  CNode *vnode = 0;
  CNode *anode = 0;
  VarNodeMap::iterator mit;
  VariablePtr v;
  CGraphPtr nlf = clone_(err);
  CQIter2 *cqit2;

  if (*err) {
    return NonlinearFunctionPtr();
  }

  if (nlf->hasVar(z)) {
    *err = 1;
    return NonlinearFunctionPtr();
  }

  // remove all nodes that have a variable from aNodes_ 
  for (CNodeVector::iterator it2=nlf->aNodes_.begin();
       it2!=nlf->aNodes_.end();) {
    if (OpVar == (*it2)->getOp()) {
      it2 = nlf->aNodes_.erase(it2);
    } else {
      ++it2;
    }
  }

  znode = nlf->newNode(z);
  if (eps>0.0) {
    anode = nlf->newNode(eps);
    znode = nlf->newNode(OpPlus, anode, znode);
  } 

  // visit all nodes that have variables in them
  for (VarSetConstIter it = nlf->vars_.begin(); it!=nlf->vars_.end(); ++it) {
    v = *it;
    if (v != z) {
      mit = nlf->varNode_.find(v);
      dnode = mit->second;

      nlf->varNode_.erase(mit);

      vnode = nlf->newNode(v);
      anode = nlf->newNode(OpDiv, vnode, znode);

      // set parents of anode.
      switch(dnode->numPar()) {
      case 0:
        break;
      case 1: 
        anode->addPar(dnode->getUPar());
        dnode->getUPar()->changeChild(dnode, anode);
        break;
      default:
        cqit2 = dnode->getParB();
        while (cqit2) {
          anode->addPar(cqit2->node);
          cqit2->node->changeChild(dnode, anode);
        }
        break;
      }
      delete dnode;
    } 
  }

  nlf->oNode_ = nlf->newNode(OpMult, nlf->oNode_, znode);
  nlf->changed_ = true;
  nlf->finalize();
  //nlf->write(std::cout);

  return nlf;
}


const CNode* CGraph::getOut() const
{
  return oNode_;
}


UInt CGraph::getHessNz()
{
  return hNnz_;
}


FunctionType CGraph::getType() const
{ 
  if (vars_.empty()) {
    return Constant;
  } else if (oNode_) {
    return oNode_->getType();
  } 
  return Constant;
}


void CGraph::grad_(int *error)
{
  for (CNodeQ::iterator it=dq_.begin(); it!=dq_.end(); ++it) {
    (*it)->setG(0.0);
  }
  for (CNodeQ::iterator it=vq_.begin(); it!=vq_.end(); ++it) {
    (*it)->setG(0.0);
  }

  oNode_->setG(1.0);
  for (CNodeQ::reverse_iterator it=dq_.rbegin(); it!=dq_.rend(); ++it) {
    (*it)->grad(error);
  }
}


VariablePtr CGraph::getVar(const CNode *cnode) const
{
  //ugly and stupid.
  for (VarNodeMap::const_iterator it=varNode_.begin(); it!=varNode_.end();
       ++it) {
    if (it->second==cnode) {
      return it->first;
    }
  }
  assert(!"should never happen!!");
  return VariablePtr();
}


void CGraph::getVars(VariableSet *vars)
{
  for (VarNodeMap::iterator it=varNode_.begin(); it!=varNode_.end(); ++it) {
    vars->insert(it->first);
  }
}


bool CGraph::isIdenticalTo(CGraphPtr cg)
{
  CNodeVector::iterator it1, it2;
  CNode *n1, *n2;
  OpCode o1, o2;
  if (!cg) {
    return false;
  }

  if (cg->aNodes_.size()!=aNodes_.size() ||
      cg->varNode_.size()!=varNode_.size() ||
      cg->dq_.size()!=dq_.size()) {
    return false;
  }

  it1 = cg->aNodes_.begin();
  it2 = aNodes_.begin();
  for (; it2 != aNodes_.end(); ++it1, ++it2) {
    n1 = *it1;
    n2 = *it2;
    o1 = (*it1)->getOp();
    o2 = (*it2)->getOp();
    if (o1!=o2) {
      return false;
    } else if ((OpNum==o1 || OpInt==o1) && 
               fabs(n1->getVal()-n2->getVal())>1e-12) {
      return false;
    } else if (OpVar==o1 && n1->getV() != n2->getV()) {
      return false;
    }
  }
  return true;
}


bool CGraph::isSumOfSquares() const
{
  return isSOSRec_(oNode_);
}


bool CGraph::isSOSRec_(CNode *node) const
{
  if (OpSqr==node->getOp()) {
    return true;
  } else if (OpPlus==node->getOp()) {
    return (isSOSRec_(node->getL()) && isSOSRec_(node->getR()));
  } else if (OpSumList==node->getOp()) {
    CNode** c1 = (node)->getListL();
    CNode** c2 = (node)->getListR();
    while (c1<c2) {
      if (!isSOSRec_(*c1)) {
        return false;
      }
      ++c1;
    }
    return true;
  } else if (OpNum==node->getOp() && node->getVal()>=0.0) {
    return true;
  }
  return false;
}

void CGraph::multiply(double c)
{
  if (fabs(c+1.0)<1e-12) {
    CNode *node = new CNode(OpUMinus, oNode_, 0);
    oNode_ = node;
    aNodes_.push_back(node);
    if (dq_.empty()==false) {
      dq_.push_back(node);
    }
  } else {
    assert(!"cannot multiply in cgraph!");
  }
}


CNode* CGraph::newNode(OpCode op, CNode *lnode, CNode *rnode)
{
  CNode *node = new CNode(op, lnode, rnode);
  aNodes_.push_back(node);
  return node;
}


CNode* CGraph::newNode(OpCode op, CNode **child, UInt n)
{
  CNode *node = new CNode(op, child, n);
  aNodes_.push_back(node);
  return node;
}


CNode* CGraph::newNode(double d)
{
  CNode *z = 0;
  CNode *node = new CNode(OpNum, z, z);
  node->setDouble(d);
  node->setVal(d);
  aNodes_.push_back(node);
  return node;
}


CNode* CGraph::newNode(int i)
{
  CNode *z = 0;
  CNode *node = new CNode(OpInt, z, z);
  node->setInt(i);
  node->setVal(i);
  aNodes_.push_back(node);
  return node;
}


CNode* CGraph::newNode(VariablePtr v)
{
  CNode *z = 0;
  VarNodeMap::iterator it = varNode_.find(v);
  if (it==varNode_.end()) {
    CNode *node = new CNode(OpVar, z, z);
    node->setV(v);
    varNode_.insert(std::pair<ConstVariablePtr, CNode*>(v, node));
    aNodes_.push_back(node);
    vars_.insert(v);
    return node;
  } else {
    return it->second;
  }
  return 0;
}


#ifdef NDEBUG
void CGraph::prepJac(VarSetConstIter vb, VarSetConstIter )
#else
void CGraph::prepJac(VarSetConstIter vb, VarSetConstIter ve)
#endif
{
  VarSetConstIter it = vb;
  UInt i=0;

  gOffs_.reserve(varNode_.size());
  for (VarSetConstIter it2 = vars_.begin(); it2!=vars_.end(); ++it2) {
    while (*it != *it2) {
      assert (it!=ve);
      ++it;
      ++i;
    }
    gOffs_.push_back(i);
  }
}


void CGraph::removeVar(VariablePtr v, double val)
{
  VarNodeMap::iterator it = varNode_.find(v);

  if (it!=varNode_.end()) {
    CNode *cnode = it->second;
    for (CNodeQ::iterator it2=vq_.begin(); it2!=vq_.end(); ++it2) {
      if (*it2 == cnode) {
        vq_.erase(it2);
        break;
      }
    }
    cnode->setVal(val);
    cnode->setDouble(val);
    cnode->setV(VariablePtr());
    cnode->setType(Constant);
    cnode->setOp(OpNum);
    vars_.erase(v);
    varNode_.erase(it);
    changed_ = true;
  }
}


void CGraph::revHess_(int *error)
{
  oNode_->setG(1.0);
  for (CNodeQ::reverse_iterator it=dq_.rbegin(); it!=dq_.rend(); ++it) {
    (*it)->grad(error); // reverse mode gradient.
    (*it)->hess(error);
  }
}


void CGraph::revHess2_(std::stack<CNode *> *st2, double mult, UInt vind,
                       double *values, UInt *nz, int *error)
{
  CNodeRSet nset;
  CNodeRSet::iterator sit; 
  CNode *n;

  nset.insert(oNode_);
  oNode_->setG(1.0);
  while (!nset.empty()) {
    sit = nset.begin();
    n = *sit;
    nset.erase(sit);
    n->hess2(&nset, error);
    if (OpVar==n->getOp() && n->getB()==true && vind >= n->getV()->getIndex()) {
      values[hOffs_[*nz]] += mult * n->getH();
      (*nz) += 1;
    }
    n->setB(false);
    n->setGi(0.0);
    n->setH(0.0);
    n->setTempI(0);
    //std::cout << nset.size() << std::endl;
  }
  while (!st2->empty()) {
    n = st2->top();
    st2->pop();
    n->setB(false);
    n->setGi(0.0);
    n->setH(0.0);
    n->setTempI(0);
  }
}


void CGraph::simplifyDq_()
{
  UInt id = 1;
  for (CNodeQ::iterator it=dq_.begin(); it!=dq_.end();) {
    if (Constant==(*it)->findFType()) {
      it = dq_.erase(it);
    } else {
      (*it)->setId(id); 
      ++id;
      ++it;
    }
  }
}


void CGraph::sqrRoot(int &err)
{
  CNode *n = 0;
  if (!oNode_) {
    err = 1;
    return;
  }

  oNode_ = newNode(OpSqrt, oNode_, n);
  changed_ = true;
  finalize();
}


void CGraph::subst(VariablePtr out, VariablePtr in, double rat)
{
  CNode *nin = 0, *nout = 0, *nmult = 0;
  VarNodeMap::iterator it;
  UInt minid;
  bool btmp;

  if (vars_.find(out)==vars_.end()) {
    return;
  }
  //std::cout << "substituting variable " << out->getName() << " by "
  //  << rat << " " << in->getName() << "\n";
  vars_.erase(out);

  it = varNode_.find(out);
  nout = it->second;
  varNode_.erase(it);


  it = varNode_.find(in);
  if (it==varNode_.end()) {
    nin = new CNode(OpVar, nin, nin);
    nin->setV(in);
    varNode_.insert(std::pair<ConstVariablePtr, CNode*>(in, nin));
    vars_.insert(in);
    for (CNodeVector::iterator it2=aNodes_.begin(); it2!=aNodes_.end(); ++it2) {
      if ((*it2)==nout) {
        *it2 = nin;
        break;
      }
    }
    btmp = false;
    for (CNodeQ::iterator it2=vq_.begin(); it2!=vq_.end(); ++it2) {
      if ((*it2)->getV()->getIndex()>in->getIndex()) {
        vq_.insert(it2, nin);
        btmp = true;
        break;
      }
    }
    if (!btmp) {
      vq_.push_back(nin);
    }
  } else {
    CNodeVector::iterator itin, itout;
    nin = it->second;
    for (CNodeVector::iterator it2=aNodes_.begin(); it2!=aNodes_.end(); ++it2) {
      if ((*it2)==nout) {
        itout = it2;
      } else if ((*it2)==nin) {
        itin = it2;
      }
    }
    if (itin<itout) {
      aNodes_.erase(itout);
    } else {
      *itout = *itin;
      aNodes_.erase(itin);
    }
  }
  if (rat!=1.0) {
    nmult = newNode(rat);
    nmult = newNode(OpMult, nin, nmult);
  }

  minid = aNodes_.size();
  for (CNodeQ::iterator it2=dq_.begin(); it2!=dq_.end(); ++it2) {
    switch ((*it2)->numChild()) {
    case (0):
      break;
    case (1):
      if ((*it2)->getL()==nout) {
        if (rat==1.0) {
          (*it2)->setL(nin);
          nin->addPar(*it2);
        } else {
          (*it2)->setL(nmult);
          nmult->addPar(*it2);
          if ((*it2)->getId()<minid) {
            minid = (*it2)->getId();
          }
        }
      }
      break;
    case (2):
      if ((*it2)->getL()==nout) {
        if (rat==1.0) {
          (*it2)->setL(nin);
          nin->addPar(*it2);
        } else {
          (*it2)->setL(nmult);
          nmult->addPar(*it2);
          if ((*it2)->getId()<minid) {
            minid = (*it2)->getId();
          }
        }
      }
      if ((*it2)->getR()==nout) {
        if (rat==1.0) {
          (*it2)->setR(nin);
          nin->addPar(*it2);
        } else {
          (*it2)->setR(nmult);
          nmult->addPar(*it2);
          if ((*it2)->getId()<minid) {
            minid = (*it2)->getId();
          }
        }
      }
      break;
    default:
      {
        CNode** c1 = (*it2)->getListL();
        CNode** c2 = (*it2)->getListR();
        while (c1<c2) {
          if ((*c1)==nout) {
            if (rat==1.0) {
              (*c1) = nin;
              nin->addPar(*c1);
            } else {
              (*c1) = nmult;
              nmult->addPar(*it2);
              if ((*it2)->getId()<minid) {
                minid = (*it2)->getId();
              }
            }
          }
          ++c1;
        }
      }
    }
  }

  for (CNodeQ::iterator it2=vq_.begin(); it2!=vq_.end(); ++it2) {
    if ((*it2)==nout) {
      vq_.erase(it2);
      break;
    }
  }
  for (CNodeQ::iterator it2=dq_.begin(); it2!=dq_.end(); ++it2) {
    if ((*it2)->getId()==minid) {
      (*it2)->setId(minid+1);
      nmult->setId(minid);
      it2 = dq_.insert(it2, nmult);
    } else if ((*it2)->getId()>minid+1) {
      (*it2)->setId((*it2)->getId()+1);
    }
  }

  if (oNode_==nout) {
    oNode_ = nin;
  }
  delete nout;
  changed_ = true;
}


void CGraph::setOut(CNode *node)
{
  oNode_ = node;
}


void CGraph::varBoundMods(double lb, double ub, VarBoundModVector &mods,
                          SolveStatus *status)
{
  double lb2 = -INFINITY;
  double ub2 = INFINITY;
  int error = 0;
  bool is_inf = false;

  computeBounds(&lb2, &ub2, &error);
  if (error>0) {
    *status = SolveError;
    return;
  }

  oNode_->setBounds(fmax(lb,lb2),fmin(ub,ub2));
  for (CNodeQ::reverse_iterator it=dq_.rbegin(); it!=dq_.rend(); ++it) {
    (*it)->propBounds(&is_inf, &error);
    if (true == is_inf) {
      *status = SolvedInfeasible;
      return;
    }
    if (error>0) {
      *status = SolveError;
      return;
    }
  }
  for (CNodeQ::iterator it=vq_.begin(); it!=vq_.end(); ++it) {
    if ((*it)->getLb()>(*it)->getV()->getLb()+1e-5) {
      mods.push_back((VarBoundModPtr) new VarBoundMod(getVar(*it), Lower,
                                                      (*it)->getLb()));
    }
    if ((*it)->getUb()<(*it)->getV()->getUb()-1e-5) {
      mods.push_back((VarBoundModPtr) new VarBoundMod(getVar(*it), Upper,
                                                      (*it)->getUb()));
    }
  }
}


void CGraph::write(std::ostream &out) const
{
  if (oNode_) {
    oNode_->writeSubExp(out);
  } else {
    out << 0 << std::endl;
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
