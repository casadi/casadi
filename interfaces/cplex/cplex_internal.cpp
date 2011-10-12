#include "cplex_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <vector>

using namespace std;
namespace CasADi{

  CplexInternal::CplexInternal(const FX& F, const FX& G, const FX& H, const FX& J, const FX& GF) : NLPSolverInternal(F,G,H,J), GF_(GF){
  casadi_warning("cplexInternal: the CPLEX interface is still experimental, more tests are needed");

  env_ = NULL;
  lp_ = NULL;
  
  addOption("objsense",   OT_INTEGER, CPX_MIN); // optimization sense (CPX_MIN or CPX_MAX)
  addOption("name",       OT_STRING, "unnamed_cplex_problem");
//   addOption("reltol",                      OT_REAL,    1e-6); // relative tolerence for the IVP solution
//   addOption("linear_solver",               OT_STRING, "dense"); // "dense", "banded" or "iterative"
//   addOption("exact_jacobian",              OT_BOOLEAN,  false);
//   addOption("is_differential",             OT_INTEGERVECTOR,  Option());
}


CplexInternal::~CplexInternal(){
  int status;
  stringstream out;
  // Free problem memeory, if necessary
  if ( lp_ != NULL ) {
    status = CPXfreeprob (env_, &lp_);
      if ( status ) {
        out << status;
        throw CasadiException("CplexInternal::~CplexInternal: CPXfreeprob failed, error code" +  out.str() );
      }
   }
  
   // Free CPLEX environment, if necessary
   if ( env_ != NULL ) {
    status = CPXcloseCPLEX (&env_);
    if ( status ) {
      throw CasadiException("CplexInternal::~CplexInternal: Could not close CPLEX environment");
      /* Note that CPXcloseCPLEX produces no output,
         so the only way to see the cause of the error is to use
         CPXgeterrorstring.  For other CPLEX routines, the errors will
         be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON. */
//       CPXgeterrorstring (env_, status, errmsg);
//       fprintf (stderr, "%s", errmsg);
    }
  }
}

void CplexInternal::setX(const vector<double>& x){
  casadi_assert( x.size() == n_ );
  x_ = x;
}

vector<double> CplexInternal::getSol(){
  return sol_;
}

void CplexInternal::init(){
  int status;
  
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
  
  H_mat_.set(H_, 0, true);
  J_mat_.set(J_, 0, false);
  
  // creating CPLEX environment and problem
  CPXsetintparam (env_, CPX_PARAM_SCRIND, CPX_ON);
  //   CPXsetintparam (env_, CPX_PARAM_SCRIND, CPX_OFF);
  env_ = CPXopenCPLEX (&status);
  if (!env_){
    throw CasadiException("CplexInternal::init: Cannot initialize CPLEX environment");
  }
//      if (debug)
//       CPXsetintparam (env_, CPX_PARAM_SCRIND, CPX_ON);
//    else
//       CPXsetintparam (env_, CPX_PARAM_SCRIND, CPX_OFF);

  lp_ = CPXcreateprob (env_, &status, getOption("name").toString().c_str());
  if (!lp_){
    throw CasadiException("CplexInternal::init: Cannot create CPLEX problem");
  }
}

void CplexInternal::evaluate(int nfdir, int nadir){
  int status;
  vector<double> objGrad, rhs, rngval,lambda,slack,dj;
  vector<char> sense;
  objGrad.resize(n_);
  dj.resize(n_);
  rhs.resize(m_);
  lambda.resize(m_);
  slack.resize(m_);
  rngval.resize(m_);
  sense.resize(m_);
  sol_.resize(n_);
  double objval;
  int solstat;
  
  // "Correct" upper and lower bounds
  for(vector<double>::iterator it=input(NLP_LBX).begin(); it!=input(NLP_LBX).end(); ++it)
    if(isinf(*it)) *it = -CPX_INFBOUND;
  for(vector<double>::iterator it=input(NLP_UBX).begin(); it!=input(NLP_UBX).end(); ++it)
    if(isinf(*it)) *it =  CPX_INFBOUND;
  // Construct upper end lower bounds in cplex format
  for(int ind=0; ind<m_; ++ind){
       double lo = input(NLP_LBG).at(ind); // lower bound
       double up = input(NLP_UBG).at(ind); // upper bound
       if ( isinf(lo) && isinf(up)){
         throw CasadiException("CplexInternal::evaluate: a constraint has no lower or upper bounds");
       }
       if ( up > lo ){ // NOTE what should be the behavior here?
         throw CasadiException("CplexInternal::evaluate: a constraint has a lower bound greater than the upper bound");
       }
       if ( lo == up ){
         // equality
         rhs[ind] = input(NLP_LBG).at(ind);
         rngval[ind] = 0;
         sense[ind] = 'E';
       } else if ( isinf(lo) ){
         // less than
         rhs[ind] = up;
         rngval[ind] = 0;
         sense[ind] = 'L';
       } else if ( isinf(up) ){
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
     F_.setInput(x_);
     F_.setAdjSeed(1.0);
     
     // Evaluate the function
     F_.evaluate(0,1);

     // Get the result
   //   F_.getOutput(obj);
     
     // Pass the argument to the function
     G_.setInput(x_);
     
     // Pass the argument to the Jacobian function
     J_.setInput(x_);
     
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
     status = CPXcopylp (env_, lp_, n_, m_,
                         getOption("objsense").toInt(), // min or max
                         getPtr(objGrad),
                         getPtr(rhs),
                         getPtr(sense),
                         J_mat_.matbeg(),
                         J_mat_.matcnt(),
                         J_mat_.matind(),
                         J_mat_.matval(),
                         &input(NLP_LBX).front(),
                         &input(NLP_UBX).front(),
                         &rngval.front());
     
     // This function will be needed to set the quadratic term in the cost
   //   status = CPXcopyquad (env_, lp, qmatbeg, qmatcnt, qmatind, qmatval);
     
   //   status = CPXaddqconstr(env_, lp, linnzcnt, quadnzcnt, rhs, sense, linind, linval, quadrow, quadcol, quadval, lname_str)

    
     // Lagrange multipliers
   //   vector<double> lambda(n_+m_);
     


     status = CPXqpopt (env_, lp_);
     
     
     //////////////
   //   std::cout << "Optimization(0): " << status << std::endl;
      if(status){
         std::cerr << "CPLEX: Failed to solve QP.\n";
      }
      status = CPXsolution (env_, lp_,
                            &solstat,
                            &objval,
                            &sol_.front(),
                            &lambda.front(),
                            &slack.front(),
                            &dj.front());  
      if(status){
         std::cerr << "CPLEX: Failed to get solution.\n";
      }
   //    if(debug){
   //       int solnstat = CPXgetstat (env_, lp_);
   //       if      ( solnstat == CPX_STAT_UNBOUNDED ) {
   //          printf ("Model is unbounded\n");
   //       }
   //       else if ( solnstat == CPX_STAT_INFEASIBLE ) {
   //          printf ("Model is infeasible\n");
   //       }
   //       else if ( solnstat == CPX_STAT_INForUNBD ) {
   //          printf ("Model is infeasible or unbounded\n");
   //       }
   //   }
   }


   // implementation of CplexMatrix

   void CplexMatrix::set(const FX& funct, int n_out, bool symm){
     function_ = funct;
     n_out_ = n_out;
     symm_ = symm;
     
     if(symm_){
       sparsity_ = function_.output(n_out).sparsity();
     } else {
       sparsity_ = function_.output(n_out).sparsity().transpose(mapping_);
       data_.resize(sparsity_.size());
     }
     
     matcnt_.resize(sparsity_.size2());
     for(int ind=0; ind<matcnt_.size(); ++ind){
       matcnt_[ind] = sparsity_.rowind()[ind+1]-sparsity_.rowind()[ind];
     }
   }

   double* CplexMatrix::matval(){
     if(symm_){
       return &function_.output(n_out_).front();
     } else {
       for(int ind=0; ind<data_.size(); ++ind){
         data_.at(ind) = function_.output(n_out_).at(mapping_[ind]);
       }
       return &data_.front();
     }
   }

   int* CplexMatrix::matbeg(){
     return const_cast<int*>(&sparsity_.rowind().front());
   }

   int* CplexMatrix::matcnt(){
     return &matcnt_.front();
   }

   int* CplexMatrix::matind(){
     return const_cast<int*>(&sparsity_.col().front());
   }

} // namespace CasADi
