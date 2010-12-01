#ifndef EVALUATION_HPP
#define EVALUATION_HPP

#include "mx_node.hpp"
#include "../fx/fx.hpp"

namespace CasADi{

/** 
  \author Joel Andersson 
  \date 2010
*/
class Evaluation : public MXNode{
public:

/** \brief  Constructor */
explicit Evaluation(const FX& fcn, const std::vector<MX> &dep, int oind);

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
virtual void evaluate(int fsens_order, int asens_order);

protected:
  FX fcn_;
  int oind;

};

} // namespace CasADi

#endif // EVALUATION_HPP
