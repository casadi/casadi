/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include "kinvec.hpp"
#include <math.h>

using namespace std;
using namespace CasADi;

namespace KINEMATICS{

KinVec::KinVec()  {

}

KinVec::KinVec(const SXMatrix &x,const SXMatrix &y, const SXMatrix &z, bool type_,const Frame& ref_): type(type_), ref(ref_)   {
  v = SXMatrix(3,1);
  v[0]=x;
  v[1]=y;
  v[2]=z;
  order=0;
}

KinVec::KinVec(const KinVec &v): type(v.type), ref(v.ref), v(v.v), J(v.J) , c(v.c), q(v.q), dq(v.dq),ddq(v.ddq), order(v.order)    {
}
  
  
KinVec::KinVec(const SXMatrix& xyz,bool type_,const Frame& ref_): type(type_), ref(ref_), v(xyz)  {
  order=0;
}

KinVec::KinVec(const SXMatrix &xyz, bool type_,const Frame& ref_,const SXMatrix &J_,const SXMatrix & q_,const SXMatrix & dq_,const SXMatrix & ddq_,int order_): type(type_), ref(ref_), v(xyz), J(J_) , q(q_), dq(dq_) , ddq(ddq_) , order(order_) {

}

KinVec::KinVec(const SXMatrix &xyz, bool type_,const Frame& ref_,const SXMatrix &J_,const SXMatrix &c_,const SXMatrix & q_,const SXMatrix & dq_,const SXMatrix & ddq_,int order_): type(type_), ref(ref_), v(xyz), J(J_) , c(c_), q(q_), dq(dq_),  ddq(ddq_), order(order_)  {

}

KinVec::KinVec(const SXMatrix &xyz, bool type_,const Frame& ref_,const SXMatrix & q_,const SXMatrix & dq_,const SXMatrix & ddq_,int order_): type(type_), ref(ref_), v(xyz), q(q_), dq(dq_),ddq(ddq_), order(order_)  {

}


std::ostream& operator<<(std::ostream &stream, const KinVec &vec){

  SXMatrix v = vec.v;
  simplify(v);
  
 stream << "KinVec" << endl;
 stream << "  Components: " << v << endl;
 stream << "  type: "<< vec.type << endl;
 stream << "  Expressed in: " << vec.ref << endl;
 stream << "  J: " << vec.J  << endl;
 stream << "  c: " << vec.c  << endl;
 stream << "  order: " << vec.order << endl;
 stream << "  q: " << vec.q << endl;
 stream << "  dq: " << vec.dq << endl;
 stream << "  ddq: " << vec.ddq << endl;
  return stream;
}

const SXMatrix& KinVec::getQ() const{
  return q;
}

const SXMatrix& KinVec::getDQ() const{
  return dq;
}

const SXMatrix& KinVec::getDDQ() const{
  return ddq;
}

void KinVec::setDDQ(const SXMatrix& ddq_){
  ddq=ddq_;
  v=J*ddq+c;
}


SXMatrix KinVec::explicitize(const SXMatrix & ddq) const {
  if (q.size() != ddq.size()) throw("ddq argument must be of same size as interal q and dq");
  if (abs(order)!=2) throw("can only solve 2nd order KinVecs");
  SXMatrix Jc;
  SXMatrix cc=c;
  for (int i=0;i<q.size();i++) {
    if (ddq[i]->isNan()) {
     Jc << getColumn(J,i);
    } else {
      cc+=getColumn(J,i)*ddq[i];
    }
  }
  return -inv(Jc)*cc;
}

SXMatrix KinVec::explicitize(std::vector<int> &di_) const {  
    std::map<int,int> di;
    for (int i = 0; i < di_.size(); i++ ){
	  di[di_[i]]=i;
    }
    return explicitize(di);
}

SXMatrix KinVec::explicitize(std::vector<int> &di_,std::vector<int> &ri_) const {  
    std::map<int,int> di;
    std::map<int,int> ri;
    for (int i = 0; i < di_.size(); i++ ){
	  di[di_[i]]=i;
    }
    for (int i = 0; i < ri_.size(); i++ ){
	  ri[ri_[i]]=i;
    }
    return explicitize(di,ri);
}

SXMatrix KinVec::explicitize(std::map<int,int> &di) const {
    std::map<int,int> ri;
    ri[0]=0;ri[1]=1;ri[2]=2;
    return explicitize(di,ri);
}

SXMatrix KinVec::explicitize(std::map<int,int> &di,std::map<int,int> &ri) const {

    //if ( abs(order) !=2  ){ THROWERROR( RET_CAN_ONLY_SOLVE_2ND_ORDER_KINVECS ); ASSERT( 1 == 0 ); }
    if ( J.size()==0) throw("KinVec jacobian information was lost. Some operators destroy this information. (e.g. cross)");

    //Expression Jc;               // SHOULD BE INTERMEDIATE STATE.
    SXMatrix Jc=zeros(ri.size(),di.size());       
    SXMatrix cc=zeros(ri.size(),1);
    
    for (int i = 0; i < 3; i++ ){
	  if (ri.count(i)) {
	    cc[ri[i]]=c[i];
	  }
    }
    
    for (int j = 0; j < ddq.size(); j++ ){
      for (int i = 0; i < 3; i++ ){
	if (di.count(j)) {
	  if (ri.count(i)) {
	    Jc(ri[i],di[j])=J(i,j);
	  }
	} else {
	  if (ri.count(i)) {
	    cc[ri[i]]+=J(i,j)*ddq[j];
	  }
	}
      }
    }

    return -inv(Jc)*cc;
}


KinVec KinVec::expressedIn(const Frame& f) {
  SXMatrix p(3,1);
  return KinVec(ref.chain(this->v,f,this->type),this->type,f,ref.chain(this->J,f,0),ref.chain(this->c,f,0),this->q,this->dq,this->ddq,this->order);
}




SXMatrix KinVec::getCoords() const{
  return v;
}

KinVec pos(const Frame& f,const Frame& ei) {
  SXMatrix p=zeros(3,1);
  return KinVec(f.chain(p,ei,1),1,ei);
}

KinVec vel(const Frame& f, const Frame& wt, const Frame& ei, const SXMatrix & q, const SXMatrix & dq,const SXMatrix & ddq) {
   SXMatrix p=zeros(3,1);
   p=f.chain(p,wt,1);
   SXMatrix J=jacobian(p,q);
   J=wt.chain(J,ei,0);
 
   return KinVec(J*dq,0,ei,J,q,dq,ddq,1);
}


KinVec acc(const Frame& f, const Frame& wt, const Frame& ei, const SXMatrix & q, const SXMatrix & dq,const SXMatrix & ddq) {
   SXMatrix p=zeros(3,1);
   p=f.chain(p,wt,1);
   SXMatrix J=jacobian(p,q);
   SXMatrix acnst=jacobian(J*dq,q)*dq;
   
   
   J=wt.chain(J,ei,0);
   acnst=wt.chain(acnst,ei,0);
   return KinVec(J*ddq+acnst,0,ei,J,acnst,q,dq,ddq,2);
}


KinVec rotVel(const Frame& f, const Frame& wt, const Frame& ei, const SXMatrix & q, const SXMatrix & dq,const SXMatrix & ddq) {
   SXMatrix R=eye(3);
   R=f.chain(R,wt,0);

   SXMatrix J;
   
   J << getRow(R,1)*jacobian(trans(getRow(R,2)),q);
   J << getRow(R,2)*jacobian(trans(getRow(R,0)),q);
   J << getRow(R,0)*jacobian(trans(getRow(R,1)),q);
   
   J=wt.chain(J,ei,0);

   return KinVec(J*dq,0,ei,J,q,dq,ddq,-1);
}


KinVec rotAcc(const Frame& f, const Frame& wt, const Frame& ei, const SXMatrix & q, const SXMatrix & dq,const SXMatrix & ddq) {
   SXMatrix R=eye(3);
   R=f.chain(R,wt,0);

   SXMatrix J;
   
   J << getRow(R,1)*jacobian(trans(getRow(R,2)),q);
   J << getRow(R,2)*jacobian(trans(getRow(R,0)),q);
   J << getRow(R,0)*jacobian(trans(getRow(R,1)),q);
   
   J=wt.chain(J,ei,0);
   
   SXMatrix acnst=jacobian(J*dq,q)*dq;
   
   return KinVec(J*ddq+acnst,0,ei,J,acnst,q,dq,ddq,-2);
}

KinVec vel(const Frame& f, const Frame& wt, const Frame& ei) {return vel(f,wt,ei,f.getQ(),f.getDQ(),f.getDDQ());}
KinVec acc(const Frame& f, const Frame& wt, const Frame& ei) {return acc(f,wt,ei,f.getQ(),f.getDQ(),f.getDDQ());}
KinVec rotVel(const Frame& f, const Frame& wt, const Frame& ei) {return rotVel(f,wt,ei,f.getQ(),f.getDQ(),f.getDDQ());}
KinVec rotAcc(const Frame& f, const Frame& wt, const Frame& ei) {return rotAcc(f,wt,ei,f.getQ(),f.getDQ(),f.getDDQ());}

KinVec ex(const Frame& f) {return KinVec(1,0,0,0,f);}
KinVec ey(const Frame& f) {return KinVec(0,1,0,0,f);}
KinVec ez(const Frame& f) {return KinVec(0,0,1,0,f);}
KinVec ex(const Frame& f,const Frame& ei) {return ex(f).expressedIn(ei);}
KinVec ey(const Frame& f,const Frame& ei) {return ey(f).expressedIn(ei);}
KinVec ez(const Frame& f,const Frame& ei) {return ez(f).expressedIn(ei);}

// Operator overloads

Frame expressCommon(KinVec &a,KinVec &b) {
   Frame c=a.ref.getCommonFrame(b.ref); 
   a=a.expressedIn(c);
   b=b.expressedIn(c);
   return c;
}

KinVec operator+(const KinVec &a,const KinVec &b) {
 KinVec A=a; 
 KinVec B=b; 
 
 SXMatrix q,dq,ddq;
 if (A.q.empty() && !B.q.empty()) {q=B.q;dq=B.dq;ddq=B.ddq;};
 if (!A.q.empty() && B.q.empty()) {q=A.q;dq=A.dq;ddq=A.ddq;};
 if (!A.q.empty() && !B.q.empty()) {
  if (!isEqual(A.dq,B.dq)) throw("Vectors were not constructed with same states (q,dq,ddq)");
  q=A.q;dq=A.dq;ddq=A.ddq;
 }

 //cout << A << endl;
 //cout << B << endl;
 Frame c=expressCommon(A,B);
 //cout << A << endl;
 //cout << B << endl;
 if (A.type && B.type) throw("Cannot add 2 position vectors/ 1-vectors");
 if (A.order && B.order && A.order != B.order) throw("Cannot add vectors of different order");
 SXMatrix J,cn;
 if (!A.J.empty() && !B.J.empty()) {
   J=A.J+B.J;
    if (!A.c.empty() && !B.c.empty()) {
    cn=A.c+B.c;
    } else if (!A.c.empty()) {
	cn=A.c;
    } else {
	cn=B.c;
    }
 } else if (!A.J.empty()) {
   J=A.J;
   if (!A.c.empty()) cn=A.c + B.v;
   if (A.c.empty()) cn=B.v;
 } else {
   J=B.J;
   if (!B.c.empty()) cn=B.c + A.v;
   if (B.c.empty()) cn=A.v;
 }

 return KinVec(A.v+B.v,A.type || B.type,c,J,cn,q,dq,ddq,A.order);
}

KinVec& KinVec::operator+=(const KinVec &b) {
 return *this = *this + b;
}

KinVec KinVec::operator-(const KinVec &b) {
 return *this+(-b);
}

KinVec& KinVec::operator-=(const KinVec &b) {
 return *this = *this - b;
}

SXMatrix norm(const KinVec &v) {return sqrt(v.v[0]*v.v[0]+v.v[1]*v.v[1]+v.v[2]*v.v[2]);}
    
KinVec operator*(const KinVec &a,const SXMatrix &b) {return KinVec(a.v*b,a.type,a.ref,a.q,a.dq,a.ddq,a.order);}
KinVec operator*(const SXMatrix &a,const KinVec &b) {
  SXMatrix J,c;
  if (!b.J.empty()) J = a*b.J;
  if (!b.c.empty()) c = a*b.c;
  return KinVec(a*b.v,b.type,b.ref,J,c,b.q,b.dq,b.ddq,b.order);
}
KinVec operator/(const KinVec &a,const SXMatrix &b) {return KinVec(a.v/b,a.type,a.ref,a.q,a.dq,a.ddq,a.order);}

KinVec KinVec::operator-() const{
  return KinVec(-v,type,ref,-J,-c,q,dq,ddq,order);
}


KinVec cross(const KinVec &a,const KinVec &b) {
 KinVec A=a; 
 KinVec B=b; 
 
 Frame c=expressCommon(A,B);
 if (A.type && B.type) throw("Cannot multiply 2 position vectors/ 1-vectors");
 //if (!A.q.isEqual(B.q)) throw("Vectors were not constructed with same states (q)");
 //if (!A.dq.isEqual(B.dq)) throw("Vectors were not constructed with same states (dq)");

 return KinVec(A.v[1]*B.v[2]-B.v[1]*A.v[2],B.v[0]*A.v[2]-A.v[0]*B.v[2],A.v[0]*B.v[1]-B.v[0]*A.v[1],0,c);
}

SXMatrix operator*(const KinVec &a,const KinVec &b) {
 KinVec A=a; 
 KinVec B=b; 
 
 Frame c=expressCommon(A,B);
 if (A.type && B.type) throw("Cannot multiply 2 position vectors/ 1-vectors");
 if (A.order != B.order ) throw("Cannot multiply vectors of different order");
 if (!isEqual(A.q,B.q)) throw("Vectors were not constructed with same states (q)");
 if (!isEqual(A.dq,B.dq)) throw("Vectors were not constructed with same states (dq)");

 return A.v[0]*B.v[0]+A.v[1]*B.v[1]+A.v[2]*B.v[2];
}

SXMatrix::Element KinVec::operator[](int i) {
  J=SXMatrix();
  c=SXMatrix();
  return v[i];
}

const SX KinVec::operator[](int i) const{
  return v[i];
}



KinVec KinVec::der() {
  SXMatrix q=ref.getQ();
  SXMatrix dq=ref.getDQ();
  SXMatrix ddq=ref.getDDQ();
  SXMatrix J=jacobian(v,dq);
  SXMatrix c=jacobian(v,q)*dq;
  
  return KinVec(J*ddq+c,0,ref,J,c,q,dq,ddq,copysign(abs(order)+1,order));
}

} // namespace KINEMATICS

