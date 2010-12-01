#include "frame.hpp"
#include "frame_node.hpp"

//#include <casadi/expression_tools.hpp>
using namespace std;
using namespace CasADi;
namespace KINEMATICS{

Frame::Frame(){
    node = 0;
    // Why not:
    //node = new FrameNode(name);
    //node->count++;
}

Frame::Frame(const std::string& name,const SXMatrix &q,const SXMatrix &dq,const SXMatrix &ddq) {
    node = new FrameNode(name,q,dq,ddq);
    node->count++;
}

SXMatrix norm(const SXMatrix &v) {return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);}

Frame::Frame(const std::string& name, const Frame& ref,const SXMatrix & T){
    //cout << "Constructor of frame" << endl;
    //cout << T << endl;
    //cout << norm(getColumn(T,0)) << endl;
    //cout << norm(getColumn(T,1)) << endl;
    //cout << norm(getColumn(T,2)) << endl;
    node = new FrameNode(name,ref,T);
    node->count++;
}

Frame::Frame(Frame const& frame){
  node = frame.node;
  node->count++;
}

Frame::Frame(FrameNode* ptr){
   node = ptr;
   node->count++;
}

Frame::~Frame(){
    if(node != 0){
      node->count--;
      if(node->count==0) delete node;
    }
}

Frame& Frame::operator=(const Frame &frame){
  // quick return if the old and new pointers point to the same object
  if(node == frame.node) return *this;

if(node){
  // decrease the object count for the old pointer
  node->count--;

  // delete if this was the last pointer	
  if(node->count == 0) delete node;
}

  // save the new pointer
  node = frame.node;
  if(node) node->count++;
  return *this;
}


std::ostream& operator<<(std::ostream &stream, const Frame &frame){
  if(frame.node) frame.node->print(stream);
  return stream;
}

const std::string& Frame::getName() const{
  return node->name;
}


const SXMatrix& Frame::getQ() const{
  return node->q;
}

const SXMatrix& Frame::getDQ() const{
  return node->dq;
}

const SXMatrix& Frame::getDDQ() const{
  return node->ddq;
}

// ------ Functions that deal with chaining

std::set<FrameNode*> Frame::getReferences() const {
  std::set<FrameNode*> refs;
  for(FrameNode* ptr = node; ptr != 0; ptr = ptr->ref.node){
      refs.insert(ptr);
  }
  return refs;
}

FrameNode* Frame::getCommonFramePtr(Frame other) const {
  std::set<FrameNode*> refs=getReferences();
  for(FrameNode* ptr = other.node; ptr != 0; ptr = ptr->ref.node){
    if(refs.count(ptr) > 0)
      return ptr;
  }
  throw "Could not find common frame";
}

Frame Frame::getCommonFrame(Frame other) const {
  return Frame(getCommonFramePtr(other));
}


std::vector<FrameNode*> Frame::getFrameChain(FrameNode* endptr) const {
  std::vector<FrameNode*> myChain;
  for(FrameNode* ptr = node; ptr != endptr; ptr = ptr->ref.node)
    myChain.push_back(ptr);
  return myChain;
}

// ---------

// The actual kinematics algorithms go below:

// Consider the tree  
//
//        ___0___
//       1       2
//     3   4       5

// (Frame 3).chain(e,(Frame 5))
//
// will return  inv(T_52)*inv(T_20)*T_10*T_31*e
// the third argument discriminates between position and velocity vector
// e can be a matrix as well

SXMatrix Frame::chain(const SXMatrix &e_,const Frame &ei,bool type) const {
  if (e_.size()==1 or e_.size()==0) return e_; // empty expression

  // why not const SXMatrix &e ?
  SXMatrix e=e_;
  
  FrameNode* common=getCommonFramePtr(ei);
  
  std::vector<FrameNode*> thisChain=     getFrameChain(common);
  std::vector<FrameNode*> eiChain  =  ei.getFrameChain(common);
  
  for (std::vector<FrameNode*>::iterator it=thisChain.begin() ; it != thisChain.end(); ++it ) {
    e=(*it)->R * e; // Going up
    if (type) e+=(*it)->p;
  }
  for (std::vector<FrameNode*>::reverse_iterator it=eiChain.rbegin() ; it != eiChain.rend(); ++it ) {
    e=trans((*it)->R)*e; // Going down
    if (type) e-=(trans((*it)->R) * (*it)->p);
  }

  return e;
}

// -----

} // namespace KINEMATICS

