#ifndef MX_CONSTANT_HPP
#define MX_CONSTANT_HPP

#include "mx_node.hpp"

namespace CasADi{

/** \brief Represents an MX that is only composed of a constant.
	\author Joel Andersson 
	\date 2010

	A regular user is not supposed to work with this Node class.
	This user can call MX(double) directly, or even rely on implicit typecasting.
	\sa zeros , ones
*/
class MXConstant : public MXNode{
public:

/** \brief  Constructor */
  MXConstant(const double *x, int n, int m, char order ='R');
/** \brief Constructor that initialize to all zeros
     \sa zeros
*/
  MXConstant(int n, int m); // all zeros

/** \brief  Clone function */
  virtual MXConstant* clone() const;

/** \brief  Print */
  virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Evaluate the adjoint gradient and add the result in the dependency nodes */
/* virtual void evaluateAdj();*/
  
/** \brief  Evaluate the second order derivative (adjoint of the forward gradient) and store the result in the dependency nodes */
/** \brief    virtual void evaluateFoA(); */

/** \brief  Access an element matrix style */
   double operator()(int i, int j) const;
   double& operator()(int i, int j);

/** \brief  Access an element vector style */
   double operator[](int k) const;
   double& operator[](int k);

/** \brief  Print */
  friend std::ostream& operator<<(std::ostream &stream, const MXConstant &x);

  virtual bool isConstant() const;
  
/** \brief  data member */
  protected:
  std::vector<double> data;

};

} // namespace CasADi


#endif // MX_CONSTANT_HPP
