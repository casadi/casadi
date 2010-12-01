#ifndef NORM_HPP
#define NORM_HPP

#include "mx_node.hpp"

namespace CasADi{

/** \brief Represents any type of general norm
  \author Joel Andersson 
  \date 2010
*/
class Norm : public MXNode{
public:

/** \brief  Constructor */
Norm(const MX& x);

/** \brief  Constructor */
virtual void evaluate(int fsens_order, int asens_order);

};

/** \brief Represents a 2-norm operation on a MX
  \author Joel Andersson 
  \date 2010
*/
class Norm2 : public Norm{
public:

/** \brief  Constructor */
Norm2(const MX& x);

/** \brief  Clone function */
virtual Norm2* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;
};

/** \brief Represents a 1-norm operation on a MX
  \author Joel Andersson 
  \date 2010
*/
class Norm1 : public Norm{
public:

/** \brief  Constructor */
Norm1(const MX& x);

/** \brief  Clone function */
virtual Norm1* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;
};

/** \brief Represents an infinity-norm operation on a MX
  \author Joel Andersson 
  \date 2010
*/
class NormInf : public Norm{
public:

/** \brief  Constructor */
NormInf(const MX& x);

/** \brief  Clone function */
virtual NormInf* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream=std::cout) const;
};

} // namespace CasADi

#endif // NORM_HPP
