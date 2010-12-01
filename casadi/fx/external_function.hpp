#ifndef EXTERNAL_FUNCTION_HPP
#define EXTERNAL_FUNCTION_HPP

#include "c_function.hpp"
#include <string>

namespace CasADi{

  
/** \brief  Forward declaration of internal class */
class ExternalFunctionNode;

/** \brief  Interface for a function that is not implemented by CasADi symbolics
  \author Joel Andersson 
  \date 2010
	*/
class ExternalFunction : public CFunction{

public:

/** \brief  CONSTRUCTORS: */
/** \brief  default constructor */
  ExternalFunction();

/** \brief  Create an empty function */
  explicit ExternalFunction(const std::string& bin_name);

/** \brief  Access functions of the node */
  ExternalFunctionNode* operator->();
  const ExternalFunctionNode* operator->() const;
    
};
  
  
#if 0  
class ExternalFunctionNode : public FXNode{
  friend class ExternalFunction;
  public:
/** \brief  no public constructors */
  ~ExternalFunctionNode();

/** \brief  Set the input ind of derivative order ord */
  virtual void setArgument(const double *x, int ind=0, int ord=0);

/** \brief  Get the input */
  virtual void getArgument(double *x, int ind=0, int ord=0) const;

/** \brief  Set the output */
  virtual void setResult(const double *x, int ind=0, int ord=0);

/** \brief  Get the output */
  virtual void getResult(double *x, int ind=0, int ord=0) const;

/** \brief  initialize */
  virtual void init();

/** \brief  Clear the memory */
  virtual void clear(int ord=0);

/** \brief  Evaluate */
  virtual void evaluate(int tape_order=0);

/** \brief  Evaluate forward derivatives */
  virtual void evaluateFwd(bool use_tape=false);

/** \brief  Evaluate adjoint derivatives */
  virtual void evaluateAdj();

  protected:
/** \brief  constructor */
  explicit ExternalFunctionNode(const std::string& bin_name);

//@{
/** \brief  Function pointer types */
  typedef int (*evaluaterPtr)();
  typedef int (*clearerPtr)(int);
  typedef int (*setterPtr)(const double*, int, int);
  typedef int (*getterPtr)(double*, int, int);
//@}
  
//@{
/** \brief  Function pointers */
  setterPtr setArgument_ptr;
  getterPtr getArgument_ptr;
  setterPtr setResult_ptr;
  getterPtr getResult_ptr;
  clearerPtr clear_ptr;
  evaluaterPtr evaluate_ptr;
  evaluaterPtr evaluateFwd_ptr;
  evaluaterPtr evaluateAdj_ptr;
//@}
  
/** \brief  handle to the dll */
  void* handle;
};


#endif

} // namespace CasADi


#endif // EXTERNAL_FUNCTION_HPP
