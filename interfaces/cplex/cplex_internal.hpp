#ifndef CPLEX_INTERNAL_HPP
#define CPLEX_INTERNAL_HPP

#include "ilcplex/cplex.h"
#include "casadi/fx/nlp_solver_internal.hpp"

namespace CasADi{

/** \brief CplexMatrix is a class used to convert CasADi matrices to CPLEX format (similar to CSC).
  The class definition can be found in cplex_internal.cpp.
  
  @copydoc NLPSolver_doc
  \author Carlo Savorgnan
  \date 2011
*/
// The following class is just used to make the interfaced code cleaner. 
class CplexMatrix{
    bool symm_; // true if the matrix is symmetric
    CRSSparsity sparsity_;
    FX function_;
    int n_out_;
    std::vector<int> matcnt_;
    std::vector<double> data_; // used to store data for non symmetric matrices
    std::vector<int> mapping_; // used for non symmetric matrices
  public:
    /// reads matrix in casadi format
    void set(const FX& funct, int n_out, bool symm);
    /// returns non-zero values
    double* matval();
    /// returns indices of the beginning of columns
    int* matbeg();
    /// returns number of entries per column
    int* matcnt();
    /// returns row numbers
    int* matind();
};

class CplexInternal : public NLPSolverInternal{
  // TODO comment me!!!!
  public:
    explicit CplexInternal(const FX& F, const FX& G, const FX& H, const FX& J, const FX& GF);
    virtual ~CplexInternal();
    virtual CplexInternal* clone() const{ return new CplexInternal(*this);}
    void setX(const std::vector<double>& x);
    std::vector<double> getSol();
    virtual void init();
    virtual void evaluate(int nfdir, int nadir);
    
    /// point used for the linearization
    std::vector<double> x_;
    /// used to store the solution
    std::vector<double> sol_;
    /// Gradient of the objective function
    FX GF_;
    /// Hessian of the Lagrangian function (used for format conversion)
    CplexMatrix H_mat_;
    /// Jacobian of the constraint function (used for format conversion)
    CplexMatrix J_mat_; 
    
    // CPLEX environment pointer
    CPXENVptr env_;
    // CPLEX lp pointer
    CPXLPptr lp_;
    
    // CPLEX double parameter
    std::map<std::string, double> double_param_;
    
    // CPLEX int parameter
    std::map<std::string, int> int_param_;
    
    // sense of the optimization (min or max)
    int sense_;
};

} // namespace CasADi

#endif //CPLEX_INTERNAL_HPP
