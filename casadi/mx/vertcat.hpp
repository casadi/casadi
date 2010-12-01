#ifndef VERTCAT_HPP
#define VERTCAT_HPP

#include "mx_node.hpp"

namespace CasADi{
/** \brief Represents a vertical concatenation of MXes
  \author Joel Andersson 
  \date 2010
*/
class Vertcat : public MXNode{
public:

/** \brief  Constructor */
explicit Vertcat(const std::vector<MX> &comp);

/** \brief  Clone function */
virtual Vertcat* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
virtual void evaluate(int fsens_order, int asens_order);

// virtual void setOutput(const std::vector<double>& val, int ord=0);
// virtual void getOutput(std::vector<double>& val, int ord=0) const;

protected:

};

} // namespace CasADi

#endif // VERTCAT_HPP
