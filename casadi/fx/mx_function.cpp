#include "mx_function.hpp"

#include "../mx/mx_node.hpp"
#include "../stl_vector_tools.hpp"

#include <stack>
#include <typeinfo>
#include <cassert>

using namespace std;

namespace CasADi{

MXFunctionNode::MXFunctionNode(const std::vector<MX>& inputv_, const std::vector<MX>& outputv_) : inputv(inputv_), outputv(outputv_){
  setOption("ad_order",1); // one by default

  // Allocate space for inputs
  input_.resize(inputv.size());
  for(int i=0; i<input_.size(); ++i)
    input_[i].setSize(inputv[i].size1(),inputv[i].size2());

  // Allocate space for outputs
  output_.resize(outputv.size());
  for(int i=0; i<output_.size(); ++i)
    output_[i].setSize(outputv[i].size1(),outputv[i].size2());

}

MXFunctionNode::~MXFunctionNode(){
}

void MXFunctionNode::makeAlgorithm(MXNode* root, vector<MXNode*> &nodes, map<const MXNode*,int>  &nodemap){

    // Stack
    stack<MXNode*> s;
    s.push(root);

    while(!s.empty()){

      // If the last element on the stack has not yet been added
      if (nodemap.find(s.top()) == nodemap.end()){


        // Loop over the children of the topmost element
        bool all_dependencies_added = true;
        for(vector<MX>::iterator it = s.top()->dep_.begin(); it!=s.top()->dep_.end(); ++it)
          if(nodemap.find(it->get()) == nodemap.end()){ // a dependency has not been added
            // Push the dependency to the top of the stack
            s.push(it->get());
            all_dependencies_added = false;
            break;
          }
  
        // Continue to the next loop if a dependency has been added
        if(!all_dependencies_added) continue;

        // Add to the vector of nodes
        nodes.push_back(s.top());

        // Save the index
        nodemap[s.top()] = nodes.size()-1;

        // Remove from stack
        s.pop();
    }
  }
}

void MXFunctionNode::init(){
  // Call the init function of the base class
  FXNode::init();

  // Clear the algorithm
  alg.clear();
  nodemap.clear();
  std::vector<MXNode*> nodes;

  // Order all nodes of a matrix syntax tree in the order of calculation
  for(vector<MX>::iterator it = outputv.begin(); it!=outputv.end(); ++it){
    makeAlgorithm(it->get(), nodes, nodemap);
  }

  // Create runtime elements for each node
  alg.resize(nodes.size());
  for(int i=0; i<alg.size(); ++i){
    // Save the node
    alg[i].assignNode(nodes[i]);

    // Set derivative order
    nodes[i]->maxord_ = 2;
    
    // Number of derivative directions
    nodes[i]->nfdir_ = nfdir_;
    nodes[i]->nadir_ = nadir_;

    nodes[i]->init();
  }
}

bool MXFunctionNode::hasEl(const MX& mx) const{
  // Locate the node:
  map<const MXNode*,int>::const_iterator it = nodemap.find(mx.get());

  // Make sure that the node was indeed found
  return it != nodemap.end();
}


int MXFunctionNode::findEl(const MX& mx) const{
  // Locate the node:
  map<const MXNode*,int>::const_iterator it = nodemap.find(mx.get());

  // Make sure that the node was indeed found
  if(it == nodemap.end()){
    stringstream ss;
    ss << "MXFunctionNode::findEl: Node \"" << mx << "\" not in the algorithm.";
    throw CasadiException(ss.str());
  }
  // Return the index
  return it->second;
}

void MXFunctionNode::setElement(int ind, const double* x, int ord){
  vector<double>& v = alg[ind]->val(ord);
  copy(x,x+v.size(),v.begin());
}

void MXFunctionNode::getElement(int ind, double *x, int ord) const{
  const vector<double>& v = alg[ind]->val(ord);
  copy(v.begin(),v.end(),x);
}

void MXFunctionNode::evaluate(int fsens_order, int asens_order){
  assert(fsens_order==0 || nfdir_==1);

  // Pass the inputs
  for(int ind=0; ind<input_.size(); ++ind)
    inputv[ind]->setOutput(input(ind).data());

  // Pass the inputs seeds
  if(fsens_order>0)
    for(int dir=0; dir<nfdir_; ++dir)
      for(int ind=0; ind<input_.size(); ++ind)
        inputv[ind]->setOutput(input(ind).dataF(dir),1); // TODO: add dir
  
  // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
  for(vector<MX>::iterator it=alg.begin(); it!=alg.end(); it++)
    (*it)->evaluate(fsens_order,0);

  // Get the outputs
  for(int ind=0; ind<outputv.size(); ++ind)
    outputv[ind]->getOutput(output(ind).data());

  // Get the output seeds
  if(fsens_order>0)
    for(int dir=0; dir<nfdir_; ++dir)
      for(int ind=0; ind<outputv.size(); ++ind)
        outputv[ind]->getOutput(output(ind).dataF(dir),1); // TODO: add dir
          
  if(asens_order>0){

    // Clear the seeds
    clear(1);

    // Pass the adjoint seeds
    for(int ind=0; ind<outputv.size(); ++ind)
      outputv[ind]->setOutput(output(ind).dataA(),1);

    // Evaluate all of the nodes of the algorithm: should only evaluate nodes that have not yet been calculated!
    for(vector<MX>::reverse_iterator it=alg.rbegin(); it!=alg.rend(); it++){
      (*it)->evaluate(0,1);
    }

    // Get the adjoint sensitivities
    for(int ind=0; ind<input_.size(); ++ind)
      inputv[ind]->getOutput(input(ind).dataA(),1);
    
    }
}

void MXFunctionNode::clear(int ord){
  // Loop over all nodes
  for(vector<MX>::iterator it=alg.begin(); it!=alg.end(); it++){
     // Loop over results
     for(vector<double>::iterator rit=(*it)->val(ord).begin(); rit !=(*it)->val(ord).end(); ++rit) 
       *rit = 0;  
  }
}

void MXFunctionNode::printAlgorithm(ostream &stream){
  int algcount = 0;
  for(vector<MX>::const_iterator it=alg.begin(); it!=alg.end(); ++it)
    cout << "[" << algcount++ << "]: " << *it << endl;
}

void MXFunctionNode::printValues(int ord, ostream &stream){
  int algcount = 0;
  for(vector<MX>::const_iterator it=alg.begin(); it!=alg.end(); ++it)
    cout << "[" << algcount++ << "]: " << (*it)->val_[ord] << endl;
}

MXFunction::MXFunction(){
}

MXFunction::MXFunction(const MX& inputm, const MX& outputm){
  vector<MX> inputv(1);
  inputv[0] = inputm;
  vector<MX> outputv(1);
  outputv[0] = outputm;
  assignNode(new MXFunctionNode(inputv,outputv));
}

MXFunction::MXFunction(const MX& inputm, const std::vector<MX>& outputv){
  vector<MX> inputv(1);
  inputv[0] = inputm;
  assignNode(new MXFunctionNode(inputv,outputv));
}

MXFunction::MXFunction(const std::vector<MX>& inputv, const MX& outputm){
  vector<MX> outputv(1);
  outputv[0] = outputm;
  assignNode(new MXFunctionNode(inputv,outputv));
}

MXFunction::MXFunction(const std::vector<MX>& inputv, const std::vector<MX>& outputv){
  assignNode(new MXFunctionNode(inputv,outputv));
}

const MXFunctionNode* MXFunction::operator->() const{
  return (const MXFunctionNode*)FX::operator->();
}

MXFunctionNode* MXFunction::operator->(){
  return (MXFunctionNode*)FX::operator->();
}

void MXFunctionNode::print(std::ostream &stream) const{
  FXNode::print(stream);
  stream << "input = " << inputv << endl;
  stream << "output = " << outputv << endl;  
}

MXFunctionNode* MXFunctionNode::clone() const{
  return new MXFunctionNode(*this);
}


} // namespace CasADi

