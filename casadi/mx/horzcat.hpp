#ifndef HORZCAT_HPP
#define HORZCAT_HPP

#include "mx_node.hpp"

namespace CasADi{

/** \brief Represents a horizontal concatenation of MXes
  \author Joel Andersson 
  \date 2010
*/
class Horzcat : public MXNode{
public:

/** \brief  Constructor */
explicit Horzcat(const std::vector<MX> &comp);

/** \brief  Clone function */
virtual Horzcat* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;

/** \brief  Evaluate the function and store the result in the node */
virtual void evaluate(int fsens_order, int asens_order);

protected:

};

} // namespace CasADi

#endif // HORZCAT_HPP
