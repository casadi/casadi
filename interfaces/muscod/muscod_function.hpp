#ifndef MUSCOD_FUNCTION_HPP
#define MUSCOD_FUNCTION_HPP

#include <symbolic/fx/fx.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include "muscod_interface.hpp"
#include <stack>

/** Comment:
The following function is depreciated, but it shows a very useful way of passing references to member functions of objects to c-function
where there is no way of passing and getting user_defined pointers.
Joel
*/


namespace CasADi{

   /** \brief  Input arguments of a MUSCOD function */
  enum MUSCOD_FCN_Input{
    MUSCOD_FCN_T, 
    MUSCOD_FCN_XD, 
    MUSCOD_FCN_XA,
    MUSCOD_FCN_U, 
    MUSCOD_FCN_P, 
    MUSCOD_FCN_NUM_IN
  };

  /** \brief  Output arguments of a MUSCOD function */
  enum MUSCOD_FCN_Output{
    MUSCOD_FCN_RHS,
    MUSCOD_FCN_RES,
    MUSCOD_FCN_NUM_OUT
  };
  
typedef void (*muscodFunctionPtr)(double  *t, double *xd, double *xa, double *u, double *p, double *rhs,double *rwh, long *iwh,  long *info);

/** This class interfaces between Muscods c-function passing and CasADi. 
    It gets extra complicated, since there is no way of passing information on which
    object a function corresponds to. All functions must be unique. To work around
    this, a solution with template meta programming has been implemented 
    */
  
/** Maximum number of instances allowed of the function */
const int MAX_NUM_INSTANCES = 100;
  
/** \brief  CasADi to MUSCOD function interface */
/// Note: only one instance is allowed for each template instaniation of this class
class MuscodFunction{
  public:
      /// Constructor (implicit typeconversion is allowed)
      MuscodFunction(const FX& f=FX());
      
      /// Destructor
      ~MuscodFunction();

      /// Copy constructor
      MuscodFunction(const MuscodFunction& fcn);

      /// Assignment
      MuscodFunction& operator=(const MuscodFunction& fcn);
      
      /// CasADi function
      FX f_;

      /// Get a function pointer to the current instance
      muscodFunctionPtr getPtr();
      
      /// Internal functions
      void fcn(double  *t, double *xd, double *xa, double *u, double *p, double *rhs,double *rwh, long *iwh,  long *info);

      /// Map the function pointers to functions
      static std::vector<muscodFunctionPtr> functions_;
      static std::vector<MuscodFunction*> instances_;
      static std::stack<int> free_instances_;
      
      /// Create a map from the function pointers to the place in the vector
      static std::vector<muscodFunctionPtr> generate_functions();

      /// The place place in the vector for the current instance
      int instance_no_;
};

/// Each class instance must have its own function pointer, it is not possible to pass pointers to member functions
template<int index>
void fcn_template(double  *t, double *xd, double *xa, double *u, double *p, double *rhs,double *rwh, long *iwh,  long *info){
  MuscodFunction::instances_[index]->fcn(t,xd,xa,u,p,rhs,rwh,iwh,info);
}


} // namespace CasADi

#endif //MUSCOD_FUNCTION_HPP
