#ifndef JACOBIAN_HPP
#define JACOBIAN_HPP

#include <set>
#include <map>
#include <vector>
#include <iostream>

#include "fx.hpp"

namespace CasADi{


  
  
/** \brief  Forward declaration of internal class */
class JacobianNode;

/** \brief Jacobian class	
  Universal Jacobian class, calculates the Jacobian of a function based on compression techniques 
  \author Joel Andersson 
  \date 2010
    joel.andersson@esat.kuleuven.be
*/ 
class Jacobian : public FX{
public:

/** \brief  Default constructor */
  Jacobian();

/** \brief  Copy constructor */
  Jacobian(const Jacobian& ref);

/** \brief  Create a Jacobian */
  explicit Jacobian(const FX& fcn, int iind=0, int oind=0);

/** \brief  Access functions of the node */
  JacobianNode* operator->();
  const JacobianNode* operator->() const;
};


/** \brief  Internal node class for Jacobian
  \author Joel Andersson 
  \date 2010
*/
class JacobianNode : public FXNode{
  friend class Jacobian;
  
  protected:
/** \brief  Multiple input, multiple output constructor, only to be accessed from Jacobian, therefore protected */
  JacobianNode(const FX& fcn, int iind, int oind);

  public:
  
/** \brief  Destructor */
  virtual ~JacobianNode();

/** \brief  Evaluate the algorithm */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Initialize */
  virtual void init();
  
  bool sparse_jac_, use_fd_, use_ad_fwd_;
  
  std::vector<int> elind_;
  std::vector<int> rr_, cc_, rind_;
  std::vector<double> epsilon_; // perturbations
  
  protected:
    
  // Dimensions
  int n_,m_;
  
/** \brief  Function to be differentiated */
  FX fcn_;
  
/** \brief  Output and input indices */
  int oind_, iind_;

};



} // namespace CasADi


#endif // JACOBIAN_HPP

