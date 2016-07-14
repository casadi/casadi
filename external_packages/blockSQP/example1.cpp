#include "blocksqp.hpp"
#include <limits>

static double const myInf = std::numeric_limits<double>::infinity();    ///< Used to mark sparse zeros in Jacobian

class MyProblem : public blockSQP::Problemspec {
public:
  blockSQP::Matrix xi0;                         ///< starting values for the optimization (dim nVar)

public:
  MyProblem( int nVar_,               ///< number of variables
             int nCon_,               ///< number of constraints
             int nBlocks_,            ///< number of diagonal blocks in the Hessian
             int *BlockIdx_,          ///< partition of the variable vector according to diagonal blocks (dim nBlocks+1)
             const blockSQP::Matrix &bl_,       ///< lower bounds for variables and constraints (dim nVar+nCon)
             const blockSQP::Matrix &bu_,       ///< upper bounds for variables and constraints (dim nVar+nCon)
             const blockSQP::Matrix &xi0_       ///< starting values for the optimization (dim nVar)
             );

  /// Set initial values for xi (and possibly lambda) and parts of the Jacobian that correspond to linear constraints (dense version).
  virtual void initialize( blockSQP::Matrix &xi,            ///< optimization variables
                           blockSQP::Matrix &lambda,        ///< Lagrange multipliers
                           blockSQP::Matrix &constrJac      ///< constraint Jacobian (dense)
                           );

  /// Set initial values for xi (and possibly lambda) and parts of the Jacobian that correspond to linear constraints (sparse version).
  virtual void initialize( blockSQP::Matrix &xi,            ///< optimization variables
                           blockSQP::Matrix &lambda,        ///< Lagrange multipliers
                           double *&jacNz,        ///< nonzero elements of constraint Jacobian
                           int *&jacIndRow,       ///< row indices of nonzero elements
                           int *&jacIndCol        ///< starting indices of columns
                           );

  /// Evaluate objective, constraints, and derivatives (dense version).
  virtual void evaluate( const blockSQP::Matrix &xi,        ///< optimization variables
                         const blockSQP::Matrix &lambda,    ///< Lagrange multipliers
                         double *objval,          ///< objective function value
                         blockSQP::Matrix &constr,          ///< constraint function values
                         blockSQP::Matrix &gradObj,         ///< gradient of objective
                         blockSQP::Matrix &constrJac,       ///< constraint Jacobian (dense)
                         blockSQP::SymMatrix *&hess,        ///< Hessian of the Lagrangian (blockwise)
                         int dmode,               ///< derivative mode
                         int *info                ///< error flag
                         );

  /// Evaluate objective, constraints, and derivatives (sparse version).
  virtual void evaluate( const blockSQP::Matrix &xi,        ///< optimization variables
                         const blockSQP::Matrix &lambda,    ///< Lagrange multipliers
                         double *objval,          ///< objective function value
                         blockSQP::Matrix &constr,          ///< constraint function values
                         blockSQP::Matrix &gradObj,         ///< gradient of objective
                         double *&jacNz,          ///< nonzero elements of constraint Jacobian
                         int *&jacIndRow,         ///< row indices of nonzero elements
                         int *&jacIndCol,         ///< starting indices of columns
                         blockSQP::SymMatrix *&hess,        ///< Hessian of the Lagrangian (blockwise)
                         int dmode,               ///< derivative mode
                         int *info                ///< error flag
                         );

  /// Generic method to convert dense constraint Jacobian to a sparse matrix in Harwell--Boeing (column compressed) format.
  virtual void convertJacobian( const blockSQP::Matrix &constrJac,  ///< constraint Jacobian (dense)
                                double *&jacNz,           ///< nonzero elements of constraint Jacobian
                                int *&jacIndRow,          ///< row indices of nonzero elements
                                int *&jacIndCol,          ///< starting indices of columns
                                bool firstCall = 0        ///< indicates if this method is called for the first time
                                );
};


MyProblem::MyProblem( int nVar_, int nCon_, int nBlocks_, int *blockIdx_, const blockSQP::Matrix &bl_, const blockSQP::Matrix &bu_, const blockSQP::Matrix &xi0_ ) {
  nVar = nVar_;
  nCon = nCon_;

  nBlocks = nBlocks_;
  blockIdx = new int[nBlocks+1];
  if (nBlocks == 1) {
    blockIdx[0] = 0;
    blockIdx[1] = nVar;
  } else {
    for (int i=0; i<nBlocks+1; i++ )
      blockIdx[i] = blockIdx_[i];
  }

  bl.Dimension( nVar + nCon ).Initialize( -myInf );
  bu.Dimension( nVar + nCon ).Initialize( myInf );

  for (int i=0; i<nVar+nCon; i++ ) {
    bl( i ) = bl_( i );
    bu( i ) = bu_( i );
  }

  objLo = -myInf;
  objUp = myInf;

  xi0.Dimension( nVar );
  for (int i=0; i<nVar; i++ )
    xi0( i ) = xi0_( i );
}


void MyProblem::convertJacobian( const blockSQP::Matrix &constrJac, double *&jacNz, int *&jacIndRow, int *&jacIndCol, bool firstCall ) {
  int nnz, count, i, j;

  if (firstCall ) {
    // 1st run: count nonzeros
    nnz = 0;
    for (j=0; j<nVar; j++ )
      for (i=0; i<nCon; i++ )
        if (fabs( constrJac( i, j ) < myInf )) nnz++;

    if (jacNz != 0 ) delete[] jacNz;
    if (jacIndRow != 0 ) delete[] jacIndRow;

    jacNz = new double[nnz];
    jacIndRow = new int[nnz + (nVar+1)];
    jacIndCol = jacIndRow + nnz;
  } else {
    nnz = jacIndCol[nVar];
    /* arrays jacInd* are already allocated! */
  }

  // 2nd run: store matrix entries columnwise in jacNz
  count = 0;
  for (j=0; j<nVar; j++ ) {
    jacIndCol[j] = count;
    for (i=0; i<nCon; i++ )
      if (fabs( constrJac( i, j ) < myInf ) ) {
        jacNz[count] = constrJac( i, j );
        jacIndRow[count] = i;
        count++;
      }
  }
  jacIndCol[nVar] = count;
  if (count != nnz )
    printf( "Error in convertJacobian: %i elements processed, should be %i elements!\n", count, nnz );
}


void MyProblem::initialize( blockSQP::Matrix &xi, blockSQP::Matrix &lambda, blockSQP::Matrix &constrJac ) {
  // set initial values for xi and lambda
  lambda.Initialize( 0.0 );
  for (int i=0; i<nVar; i++ )
    xi( i ) = xi0( i );
}


void MyProblem::initialize( blockSQP::Matrix &xi, blockSQP::Matrix &lambda, double *&jacNz, int *&jacIndRow, int *&jacIndCol ) {
  blockSQP::Matrix constrDummy, gradObjDummy, constrJac;
  blockSQP::SymMatrix *hessDummy;
  double objvalDummy;
  int info;

  // set initial values for xi and lambda
  lambda.Initialize( 0.0 );
  for (int i=0; i<nVar; i++ )
    xi( i ) = xi0( i );

  // find out Jacobian sparsity pattern by evaluating derivatives once
  constrDummy.Dimension( nCon ).Initialize( 0.0 );
  gradObjDummy.Dimension( nVar ).Initialize( 0.0 );
  constrJac.Dimension( nCon, nVar ).Initialize( myInf );
  evaluate( xi, lambda, &objvalDummy, constrDummy, gradObjDummy, constrJac, hessDummy, 1, &info );

  // allocate sparse Jacobian structures
  convertJacobian( constrJac, jacNz, jacIndRow, jacIndCol, 1 );
}

/*
 * PROBLEM-SPECIFIC PART STARTS HERE
 */

void MyProblem::evaluate( const blockSQP::Matrix &xi, const blockSQP::Matrix &lambda, double *objval, blockSQP::Matrix &constr,
                          blockSQP::Matrix &gradObj, double *&jacNz, int *&jacIndRow, int *&jacIndCol,
                          blockSQP::SymMatrix *&hess, int dmode, int *info ) {
  blockSQP::Matrix constrJac;

  constrJac.Dimension( nCon, nVar ).Initialize( myInf );
  evaluate( xi, lambda, objval, constr, gradObj, constrJac, hess, dmode, info );

  // Convert to sparse format
  if (dmode != 0 )
    convertJacobian( constrJac, jacNz, jacIndRow, jacIndCol, 0 );
}


void MyProblem::evaluate( const blockSQP::Matrix &xi, const blockSQP::Matrix &lambda, double *objval, blockSQP::Matrix &constr,
                          blockSQP::Matrix &gradObj, blockSQP::Matrix &constrJac, blockSQP::SymMatrix *&hess,
                          int dmode, int *info ) {
  *info = 0;

  /*
   * min   x1**2 - 0.5*x2**2
   * s.t.  x1 - x2 = 0
   */
  if (dmode >= 0 ) {
    *objval = xi(0)*xi(0) - 0.5*xi(1)*xi(1);
    constr( 0 ) = xi(0) - xi(1);
  }

  if (dmode > 0 ) {
    gradObj( 0 ) = 2.0 * xi(0);
    gradObj( 1 ) = -xi(1);

    constrJac( 0, 0 ) = 1.0;
    constrJac( 0, 1 ) = -1.0;
  }
}

int main() {
  int ret = 0;
  MyProblem *prob;
  blockSQP::SQPmethod *meth;
  blockSQP::SQPoptions *opts;
  blockSQP::SQPstats *stats;
  char outpath[255];
  strcpy( outpath, "./" );

  /*--------------------*/
  /* Setup problem data */
  /*--------------------*/
  int nVar = 2;
  int nCon = 1;

  // Initial values
  blockSQP::Matrix x0;
  x0.Dimension( nVar );
  x0( 0 ) = 10.0;
  x0( 1 ) = 10.0;

  // Variable bounds
  blockSQP::Matrix bl, bu;
  bl.Dimension( nVar+nCon ).Initialize( -myInf );
  bu.Dimension( nVar+nCon ).Initialize( myInf );

  // Constraint bounds
  bl( nVar ) = bu( nVar ) = 0.0;

  // Variable partition for block Hessian
  int nBlocks = nVar;
  int blockIdx[nBlocks+1];
  for (int i=0; i<nBlocks+1; i++ )
    blockIdx[i] = i;

  // Create problem evaluation object
  prob = new MyProblem( nVar, nCon, nBlocks, blockIdx, bl, bu, x0 );

  /*------------------------*/
  /* Options for SQP solver */
  /*------------------------*/
  opts = new blockSQP::SQPoptions();
  opts->opttol = 1.0e-12;
  opts->nlinfeastol = 1.0e-12;

  // 0: no globalization, 1: filter line search
  opts->globalization = 0;
  // 0: (scaled) identity, 1: SR1, 2: BFGS
  opts->hessUpdate = 0;
  // 0: initial Hessian is diagonal matrix, 1: scale initial Hessian according to Nocedal p.143,
  // 2: scale initial Hessian with Oren-Luenberger factor 3: scale initial Hessian with geometric mean of 1 and 2
  // 4: scale Hessian in every step with centered Oren-Luenberger sizing according to Tapia paper
  opts->hessScaling = 0;
  // scaling strategy for fallback BFGS update if SR1 and globalization is used
  opts->fallbackScaling = 0;
  // Size of limited memory
  opts->hessLimMem = 0;
  // If too many updates are skipped, reset Hessian
  opts->maxConsecSkippedUpdates = 200;
  // 0: full space Hessian approximation (ignore block structure), 1: blockwise updates
  opts->blockHess = 0;
  opts->whichSecondDerv = 0;
  opts->sparseQP = 1;
  opts->printLevel = 2;


  /*-------------------------------------------------*/
  /* Create blockSQP method object and run algorithm */
  /*-------------------------------------------------*/
  stats = new blockSQP::SQPstats( outpath );
  meth = new blockSQP::SQPmethod( prob, opts, stats );

  meth->init();
  ret = meth->run( 100 );
  meth->finish();
  if (ret == 1 )
    printf("\033[0;36m***Maximum number of iterations reached.***\n\033[0m");

  printf("\nPrimal solution:\n");
  meth->vars->xi.Print();
  printf("\nDual solution:\n");
  meth->vars->lambda.Print();
  printf("\nHessian approximation at the solution:\n");
  for (int i=0; i<meth->vars->nBlocks; i++ )
    meth->vars->hess[i].Print();
  //printf("\nFallback Hessian at the solution:\n");
  //for (int i=0; i<meth->vars->nBlocks; i++ )
  //meth->vars->hess2[i].Print();

  // Clean up
  delete prob;
  delete stats;
  delete opts;
  delete meth;
}
