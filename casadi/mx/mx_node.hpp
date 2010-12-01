#ifndef MX_NODE_HPP
#define MX_NODE_HPP

#include "mx.hpp"
#include "../fx/mx_function.hpp"
#include <vector>

namespace CasADi{

/** \brief Node class for MX objects
  \author Joel Andersson 
  \date 2010
  A regular user is not supposed to work with this Node class. Use Option directly.
*/
class MXNode : public SharedObjectNode{
friend class MX;
friend class MXFunctionNode;

public:
  
/** \brief  Default constructor, multiple dependencies */
 explicit MXNode(const std::vector<MX>& dep=std::vector<MX>()); // multiple dependencies, default constructor
/** \brief Unary node constructor */
 explicit MXNode(const MX& dep); // unary node
/** \brief Binary node constructor */
 explicit MXNode(const MX& dep1, const MX& dep2); // binary node
/** \brief Ternary node constructor */
 explicit MXNode(const MX& dep1, const MX& dep2, const MX& dep3); // ternary node

/** \brief  Destructor */
  virtual ~MXNode();

/** \brief  Clone function */
//  virtual MXNode* clone() const;

/** \brief  Print description */
  virtual void print(std::ostream &stream) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order) = 0;

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
//  virtual void evaluateAdj();

/** \brief  Initialize */
  virtual void init();
  
/** \brief  Get the name */
  virtual const std::string& getName() const;
  
/** \brief  Check if symbolic */
  virtual bool isSymbolic() const;

/** \brief  Check if constant */
  virtual bool isConstant() const;
    
/** \brief  Set/get input/output */
  void setOutput(const vector<double>& val, int ord=0);
  void getOutput(vector<double>&, int ord=0) const;

  /** \brief  dependencies - functions that have to be evaluated before this one */
  MX& dep(int ind=0);
  const MX& dep(int ind=0) const;
  
  /** \brief  Number of dependencies */
  int ndep() const;
  
  /** \brief  Numerical value */
  const std::vector<double>& val(int order=0, int dir=0) const;
  std::vector<double>& val(int order=0, int dir=0);

  protected:
    //! Number of derivatives
    int maxord_;
    
    //! Number of derivative directions
    int nfdir_, nadir_;
    
    /** \brief  expression size */
    MatrixSize sz;
  
    /** \brief  dependencies - functions that have to be evaluated before this one */
    std::vector<MX> dep_;
    
  private:
  
    /** \brief  Numerical value */
    std::vector< std::vector<std::vector<double> > > val_;



};

} // namespace CasADi


#endif // MX_NODE_HPP
