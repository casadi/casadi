#ifndef ACADO_FUNCTION_HPP
#define ACADO_FUNCTION_HPP

#include <casadi/fx/fx.hpp>
#include <acado/function/c_function.hpp>

/** \brief  forward declarations */
namespace ACADO{
  class CFunction;
}

namespace CasADi{

/** \brief  CasADi to ACADO function interface */
class AcadoFunction{
  public:
      // Constructor
      AcadoFunction(const FX& f=FX());
      
      // Destructor
      ~AcadoFunction();
      
      // Initialize
      void init();

      // CasADi function
      FX f_;
      
      // ACADO c function pointer
      ACADO::CFunction *fcn_;

      // Input dimensions
      std::vector<int> dim_;

      // Number of equations
      int neq_;
      
      // C-callback functions called by ACADO
      static void fcn_wrapper( double *x, double *res, void *user_data );
      static void fcn_fwd_wrapper(int number, double *x, double *seed, double *f, double *df, void *user_data);
      static void fcn_adj_wrapper(int number, double *x, double *seed, double *f, double *df, void *user_data);

      // Internal functions
      void fcn( double *x, double *res);
      void fcn_fwd(int number, double *x, double *seed, double *f, double *df);
      void fcn_adj(int number, double *x, double *seed, double *f, double *df);
};

} // namespace CasADi

#endif //ACADO_FUNCTION_HPP
