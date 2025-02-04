//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

template<typename T1>
struct casadi_clarabel_prob {
  // The base QP problem (provided by CasADi)
  const casadi_qp_prob<T1>* qp;
  // CSC representation for the quadratic cost (P)
  const int *colindp, *rowp;
  // CSC representation for the constraint matrix (A)
  const int *colinda, *rowa;
};

template<typename T1>
void casadi_clarabel_setup(casadi_clarabel_prob<T1>* p) {
  // In this simple example the setup is trivial.
  // (If any preprocessing is needed on the CSC arrays, do it here.)
}

// The data structure that will be passed persistently during a solve.
template<typename T1>
struct casadi_clarabel_data {
  // Pointer to the problem data
  const casadi_clarabel_prob<T1>* prob;
  // The QP data (pointers to cost, bounds, etc.)
  casadi_qp_data<T1>* qp;
  // A return status from Clarabel
  int return_status;
};

// Allocate and initialize memory for Clarabel
template<typename T1>
int clarabel_init_mem(casadi_clarabel_data<T1>* d) {
  return 0;
}

template<typename T1>
void clarabel_free_mem(casadi_clarabel_data<T1>* d) {
}

// The work vector sizes are determined from the QP problem.
template<typename T1>
void casadi_clarabel_work(const casadi_clarabel_prob<T1>* p,
                          casadi_int* sz_arg,
                          casadi_int* sz_res,
                          casadi_int* sz_iw,
                          casadi_int* sz_w) {
  casadi_qp_work(p->qp, sz_arg, sz_res, sz_iw, sz_w);
}

// (Optionally, initialize pointer arrays here; we leave this empty.)
template<typename T1>
void casadi_clarabel_init(casadi_clarabel_data<T1>* d,
                          const T1*** arg,
                          T1*** res,
                          casadi_int** iw,
                          T1** w) {
  // No special initialization needed in this example.
}

// The main solve routine: builds the Clarabel problem, calls the solver,
// and retrieves the solution.
template<typename T1>
int casadi_clarabel_solve(casadi_clarabel_data<T1>* d,
                          const double** arg,
                          double** res,
                          casadi_int* iw,
                          double* w) {
  /*const casadi_clarabel_prob<T1>* p = d->prob;
  const casadi_qp_prob<T1>* p_qp = p->qp;
  casadi_qp_data<T1>* d_qp = d->qp;

  // Build the CSC matrix for the quadratic cost P.
  ClarabelCscMatrix P;
  clarabel_CscMatrix_init(&P,
      p_qp->nx,          // number of rows (variables)
      p_qp->nx,          // number of columns
      reinterpret_cast<uintptr_t*>(const_cast<int*>(p->colindp)),
      reinterpret_cast<uintptr_t*>(const_cast<int*>(p->rowp)),
      d_qp->h           // pointer to the nonzero Hessian entries
  );

  // Build the CSC matrix for the constraints A.
  ClarabelCscMatrix A;
  clarabel_CscMatrix_init(&A,
      p_qp->na,         // number of constraints (rows)
      p_qp->nx,         // number of variables (columns)
      reinterpret_cast<uintptr_t*>(const_cast<int*>(p->colinda)),
      reinterpret_cast<uintptr_t*>(const_cast<int*>(p->rowa)),
      d_qp->a          // pointer to the nonzero entries of A
  );

  // For the cost, use q stored in d_qp->g.
  const ClarabelFloat* q = d_qp->g;
  // For the constraints, we assume that the right-hand side is stored in d_qp->b.
  const ClarabelFloat* b = d_qp->b;

  // For the cones, we simply assume that all constraints are equalities.
  // (Clarabel represents an equality as a zero cone.)
  ClarabelSupportedConeT cones[1];
  cones[0] = ClarabelZeroConeT(p_qp->na);

  // Get default settings for Clarabel.
  ClarabelDefaultSettings settings = clarabel_DefaultSettings_default();

  // Create a Clarabel solver instance with the problem data.
  ClarabelDefaultSolver* solver = clarabel_DefaultSolver_new(&P, q, &A, b, 1, cones, &settings);

  // Solve the problem.
  clarabel_DefaultSolver_solve(solver);

  // Retrieve the solution.
  ClarabelDefaultSolution solution = clarabel_DefaultSolver_solution(solver);
  for (int i = 0; i < p_qp->nx; ++i) {
    d_qp->x[i] = solution.x[i];
  }
  if (d_qp->f) {
    *d_qp->f = solution.obj_val;
  }
  d_qp->success = (solution.status == ClarabelSolved);
  d->return_status = solution.status;

  clarabel_DefaultSolver_free(solver);*/

// QP Example

    /* From dense matrix:
     * [[6., 0.],
     *  [0., 4.]]
     */
    ClarabelCscMatrix P;

    uintptr_t colptr[3] = { 0, 1, 2 };
    uintptr_t rowval[2] = { 0, 1 };
    ClarabelFloat nzval[2] = { 6., 4. };    
    clarabel_CscMatrix_init(
        &P,
        2,                          // row
        2,                          // col
        colptr,   // colptr
        rowval,      // rowval
        nzval// nzval
    );

    ClarabelFloat q[2] = { -1., -4. };

    /* From dense matrix:
     * [[ 1., -2.], // <-- LHS of equality constraint (lower bound)
     *  [ 1.,  0.], // <-- LHS of inequality constraint (upper bound)
     *  [ 0.,  1.], // <-- LHS of inequality constraint (upper bound)
     *  [-1.,  0.], // <-- LHS of inequality constraint (lower bound)
     *  [ 0., -1.]] // <-- LHS of inequality constraint (lower bound)
     */
    ClarabelCscMatrix A;

    uintptr_t colptr_A[3] = { 0, 2, 5 };
    uintptr_t rowval_A[5] = { 0, 1, 0, 2, 3 };
    ClarabelFloat nzval_A[5] = { 1., 1., -2., 1., -1. };
    clarabel_CscMatrix_init(
        &A,
        5,                                             // row
        2,                                             // col
        colptr_A,                      // colptr
        rowval_A,             // rowval
        nzval_A // nzval
    );

    ClarabelFloat b[5] = { 0., 1., 1., 1., 1. };

    ClarabelSupportedConeT cones[2];

    cones[0].tag = ClarabelZeroConeT_Tag;
    cones[0].zero_cone_t = 1;
    cones[1].tag = ClarabelNonnegativeConeT_Tag;
    cones[1].nonnegative_cone_t = 4;

    // Settings
    ClarabelDefaultSettings settings = clarabel_DefaultSettings_default();

    // Build solver
    ClarabelDefaultSolver *solver = clarabel_DefaultSolver_new(
        &P, // P
        q,  // q
        &A, // A
        b,  // b
        2,  // n_cones
        cones, &settings
    );

    // Solve
    clarabel_DefaultSolver_solve(solver);

    // Get solution
    ClarabelDefaultSolution solution = clarabel_DefaultSolver_solution(solver);
    printf("sol.x[0]: %e", solution.x[0]);

    // Free the matrices and the solver
    clarabel_DefaultSolver_free(solver);

  return 0;
}