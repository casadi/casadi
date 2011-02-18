#include "cplex_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <vector>

using namespace std;
namespace CasADi{

// The following class is just used to make the interfaced code cleaner. 
// ATTENTION! No checks are implemented!
enum CplexMatrixType{SYMMETRIC, NON_SYMMETRIC};
class CplexMatrix{
    MatType type_;
    CRSSparsity sparsity_;
    FX function_;
    int n_out_;
    vector<int> matcnt_;
    vector<double> data_; // used to store data for non symmetric matrices
    vector<int> mapping_; // used for non symmetric matrices
  public:
    /// reads matrix in casadi format
    void set(const& FX function, int n_out, CplexMatrixType type);
    /// returns non-zero values
    double* matval() const;
    /// returns indices of the beginning of columns
    int* matbeg() const;
    /// returns number of entries per column
    int* matcnt() const;
    /// returns row numbers
    int* matind() const;
};

void CplexMatrix::set(const& FX function, int n_out, CplexMatrixType type){
  function_ = function;
  n_out_ = n_out;
  type_ = type;
  
  if(type_ == SYMMETRIC){
    sparsity_ = function_.output(n_out).sparsity();
  } else { //type_ == NON_SYMMETRIC
    sparsity_ = function_.output(n_out).sparsity().transpose(mapping);
    data_.resize(sparsity_.size());
  }
  sparsity_ = mat.sparsity().transpose(mapping);
  
  matcnt_.resize(sparsity.size2());
  for(int ind=0; ind<matcnt_.size(); ++ind){
    matcnt_[ind] = sparsity.rowind()[ind+1]-sparsity.rowind()[ind];
  }
}

double* CplexMatrix::matval(){
  if(type_ == SYMMETRIC){
    return &function_.output(n_out_)[0];
  } else { //type_ == NON_SYMMETRIC
    for(int ind=0; ind<data_.size(); ++ind){
      data[ind] = function_.output(n_out_)[mapping_[ind]];
    }
    return &data[0];
  }
}

int* CplexMatrix::matbeg(){
  return &sparsity_.rowind()[0];
}

int* CplexMatrix::matcnt(){
  return &matcnt_[0];
}

int* CplexMatrix::matind(){
  return &sparsity_.col()[0];
}


// CPLEX interface

CplexInternal::CplexInternal(const FX& F, const FX& G, const FX& H, const FX& J, const FX& GF) : F_(F), G_(G), H_(H), J_(J), GF_(GF){
  casadi_warning("cplexInternal: the CPLEX interface is still experimental, more tests are needed");

  env_ = NULL;
  lp_ = NULL;
  
  addOption("objsense",               OT_INTEGER, CPX_MIN); // optimization sense (CPX_MIN or CPX_MAX)
//   addOption("reltol",                      OT_REAL,    1e-6); // relative tolerence for the IVP solution
//   addOption("linear_solver",               OT_STRING, "dense"); // "dense", "banded" or "iterative"
//   addOption("exact_jacobian",              OT_BOOLEAN,  false);
//   addOption("is_differential",             OT_INTEGERVECTOR,  Option());
}


CplexInternal::~CplexInternal(){
  int status;
  stringstream out;
  // Free problem memeory, if necessary
  if ( lp != NULL ) {
    status = CPXfreeprob (env, &lp);
      if ( status ) {
        out << status;
        throw CasadiException("CplexInternal::~CplexInternal: CPXfreeprob failed, error code" +  out.str() );
      }
   }
  
   // Free CPLEX environment, if necessary
   if ( env != NULL ) {
    status = CPXcloseCPLEX (&env);
    if ( status ) {
      throw CasadiException("CplexInternal::~CplexInternal: Could not close CPLEX environment");
      /* Note that CPXcloseCPLEX produces no output,
         so the only way to see the cause of the error is to use
         CPXgeterrorstring.  For other CPLEX routines, the errors will
         be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON. */
//       CPXgeterrorstring (env, status, errmsg);
//       fprintf (stderr, "%s", errmsg);
    }
  }
}

void CplexInternal::init(){
  // Initialize the functions
  F_.init();
  if(!G_.isNull()) G_.init();
  // n_ is the number of variables
  n_ = F_.input(0).numel();
  // m_ is the number of constraints
  m_ = G_.isNull() ? 0 : G_.output(0).numel();

  // Call the init method of the base class
  NLPSolverInternal::init();

  // Create a Jacobian if it does not already exists
  if(!G_.isNull() && J_.isNull()){
    J_ = G_.jacobian();
  }
  if(!J_.isNull()) J_.init();
  if(!H_.isNull()) H_.init();
  if(!GF_.isNull()) GF_.init();
  
  H_mat_.set(H_, 0, SYMMETRIC);
  J_mat_.set(J_, 0, NON_SYMMETRIC);
  
  // creating CPLEX environment and problem
  CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
  //   CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
  env = CPXopenCPLEX (&status);
  if (!env){
    throw CasadiException("CplexInternal::init: Cannot initialize CPLEX environment");
  }
//      if (debug)
//       CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
//    else
//       CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);

  lp = CPXcreateprob (env, &status, probname);
  if (!lp){
    throw CasadiException("CplexInternal::init: Cannot create CPLEX problem");
  }
}

void CplexInternal::evaluate(int fsens_order, int asens_order){
  vector<double> objGrad, rhs, rngval;
  vector<char> sense;
  objGrad.resize(n_);
  rhs.resize(m_);
  rngval.resize(m_);
  sense.resize(m_);
  
  // "Correct" upper and lower bounds
  for(vector<double>::iterator it=input(NLP_LBX).begin(); it!=input(NLP_LBX).end(); ++it)
    if(isinf(*it)) *it = -CPX_INFBOUND;
  for(vector<double>::iterator it=input(NLP_UBX).begin(); it!=input(NLP_UBX).end(); ++it)
    if(isinf(*it)) *it =  CPX_INFBOUND;
  // Construct upper end lower bounds in cplex format
  for(int ind=0; ind<m_; ++ind){
    double lo = input(NLP_LBG)[ind]; // lower bound
    double up = input(NLP_UBG)[ind]; // upper bound
    if ( isinf(lo) && isinf(up)){
      throw CasadiException("CplexInternal::evaluate: a constraint has no lower or upper bounds");
    }
    if ( up > lo ){ // NOTE what should be the behavior here?
      throw CasadiException("CplexInternal::evaluate: a constraint has a lower bound greater than the upper bound");
    }
    if ( lo == up ){
      // equality
      rhs[ind] = input(NLP_LBG)[ind];
      rngval[ind] = 0;
      sense[ind] = 'E';
    } else if isinf(lo){
      // less than
      rhs[ind] = up;
      rngval[ind] = 0;
      sense[ind] = 'L';
    } else if isinf(up) {
      // grater than
      rhs[ind] = lo;
      rngval[ind] = 0;
      sense[ind] = 'G';
    } else {
      // range
      rhs[ind] = lo;
      rngval[ind] = up-lo;
      sense[ind] = 'R';
    }
  }
  
  // Evaluation of all functions
  // Pass the argument to the function
  F_.setInput(x);
  F_.setAdjSeed(1.0);
  
  // Evaluate the function
  F_.evaluate(0,1);

  // Get the result
  F_.getOutput(obj);
  
  // Pass the argument to the function
  G_.setInput(x);
  
  // Pass the argument to the Jacobian function
  J_.setInput(x);
  
  // Evaluate the Jacobian function
  J_.evaluate();
  
  // The quadratic term is not used yet in the interface
//   H_.setInput(x);
//   H_.setInput(lambda,1);
//   H_.setInput(1.0,2);
//   // Evaluate
//   H_.evaluate();
  
  // Get the result
  F_.getAdjSens(objGrad);
  
  // setting problem data
  status = CPXcopylp (env, lp, n_, m_,
                      getOption("objsense").toInt(), // min or max
                      &objGrad[0],
                      &rhs[0],
                      &sense[0],
                      J_mat_.matbeg(),
                      J_mat_.matcnt(),
                      J_mat_.matind(),
                      J_mat_.matval(),
                      &input(NLP_LBX)[0],
                      &input(NLP_UBX)[0],
                      &rngval.begin());
  
  // This function will be needed to set the quadratic term in the cost
//   status = CPXcopyquad (env, lp, qmatbeg, qmatcnt, qmatind, qmatval);
  
//   status = CPXaddqconstr(env, lp, linnzcnt, quadnzcnt, rhs, sense, linind, linval, quadrow, quadcol, quadval, lname_str)

 
  // Lagrange multipliers
//   vector<double> lambda(n_+m_);
  


  status = CPXqpopt (env, lp);
  
  
  //////////////
//   std::cout << "Optimization(0): " << status << std::endl;
   if(status){
      std::cerr << "CPLEX: Failed to solve QP.\n";
   }
   status = CPXsolution (env, lp, &solstat, &objval, x.data().begin(), pi.data().begin(), slack, dj);  
   if(status){
      std::cerr << "CPLEX: Failed to get solution.\n";
   }
   

   if(debug){
      int solnstat = CPXgetstat (env, lp);
      if      ( solnstat == CPX_STAT_UNBOUNDED ) {
         printf ("Model is unbounded\n");
      }
      else if ( solnstat == CPX_STAT_INFEASIBLE ) {
         printf ("Model is infeasible\n");
      }
      else if ( solnstat == CPX_STAT_INForUNBD ) {
         printf ("Model is infeasible or unbounded\n");
      }
  }
}

} // namespace CasADi