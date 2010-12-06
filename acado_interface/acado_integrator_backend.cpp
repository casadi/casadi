#include <acado/utils/acado_utils.hpp>
#include <acado/matrix_vector/matrix_vector.hpp>
#include <acado/symbolic_expression/symbolic_expression.hpp>
#include <acado/function/function_.hpp>
#include <acado/integrator/integrator.hpp>
#include "acado_integrator_backend.hpp"
#include <cassert>
namespace ACADO{

Integrator* AcadoIntegratorBackend::create(void *user_data){
  return new AcadoIntegratorBackend(user_data);
}

AcadoIntegratorBackend::AcadoIntegratorBackend(void *user_data)
              :Integrator( ){

   initializeVariables();
}


AcadoIntegratorBackend::AcadoIntegratorBackend( const DifferentialEquation& rhs_ )
              :Integrator( ){
                assert(0);

    init( rhs_ );
}

AcadoIntegratorBackend::AcadoIntegratorBackend( const AcadoIntegratorBackend& arg )
              :Integrator( arg ){

    constructAll( arg );
}



AcadoIntegratorBackend::~AcadoIntegratorBackend( ){

    deleteAll();
}


AcadoIntegratorBackend& AcadoIntegratorBackend::operator=( const AcadoIntegratorBackend& arg ){
                assert(0);

    if ( this != &arg ){
        deleteAll();
        Integrator::operator=( arg );
        constructAll(arg);
    }
    return *this;
}



returnValue AcadoIntegratorBackend::init( const DifferentialEquation &rhs_ ){

  
    // RHS:
    // ---------
    rhs = new DifferentialEquation( rhs_ );
    m   = rhs->getDim ();
    ma  = rhs->getNXA ();
    mdx = rhs->getNDX ();
    mn  = rhs->getN   ();
    mu  = rhs->getNU  ();
    mui = rhs->getNUI ();
    mp  = rhs->getNP  ();
    mpi = rhs->getNPI ();
    mw  = rhs->getNW  ();
    md  = rhs->getNumDynamicEquations();

    rhs->makeImplicit();
    allocateMemory();

    return SUCCESSFUL_RETURN;
}


void AcadoIntegratorBackend::allocateMemory( ){


    int run1, run2, run3;

    if( m < 1 ){
        ACADOERROR(RET_TRIVIAL_RHS);
        return;
    }

    if( mdx > m-ma ){
        ACADOERROR(RET_TO_MANY_DIFFERENTIAL_STATE_DERIVATIVES);
        return;
    }

    if( m < md+ma ) md = m - ma;


    initializeVariables();

    initialAlgebraicResiduum = new double[m];
    relaxationConstant       = 5.0;


    // BUTCHER TABLEAU:
    // ----------------
    dim = 7;
    A  = new double*[dim];
    b4 = new double [dim];
    b5 = new double [dim];
    c  = new double [dim];

    for( run1 = 0; run1 < dim; run1++ ){
        A[run1] = new double[dim];
    }

    initializeButcherTableau();

    // RK-STARTER:
    // -----------
    ndir = rhs->getNumberOfVariables() + 1 + 2*md;

    eta4  = new double [m];
    eta5  = new double [m];
    eta4_ = new double [m];
    eta5_ = new double [m];

    for( run1 = 0; run1 < m; run1++ ){

        eta4 [run1] = 0.0;
        eta5 [run1] = 0.0;
        eta4_[run1] = 0.0;
        eta5_[run1] = 0.0;
    }


    k     = new double**[4];
    k2    = new double**[4];
    l     = new double**[4];
    l2    = new double**[4];

    for( run3 = 0; run3 < 4; run3++ ){

        k   [run3]  = new double*[dim];
        k2  [run3]  = new double*[dim];
        l   [run3]  = new double*[dim];
        l2  [run3]  = new double*[dim];

        for( run1 = 0; run1 < dim; run1++ ){
            k    [run3][run1] = new double[m];
            k2   [run3][run1] = new double[m];

            for( run2 = 0; run2 < m; run2++ ){
                k [run3][run1][run2] = 0.0;
                k2[run3][run1][run2] = 0.0;
            }

            l    [run3][run1] = new double[ndir];
            l2   [run3][run1] = new double[ndir];

            for( run2 = 0; run2 < ndir; run2++ ){
                l [run3][run1][run2] = 0.0;
                l2[run3][run1][run2] = 0.0;
            }
        }
    }


    iseed = new double[ndir];

    x     = new double [ndir];

    for( run1 = 0; run1 < ndir; run1++ ){
        x[run1]     = 0.0;
        iseed[run1] = 0.0;
    }

    t = 0.0;


    // BDF-METHOD:
    // -----------
    nstep     = 5;
    psi       = (double**)calloc(1,sizeof(double*));
    psi_      = new double[nstep];
    gamma     = (double**)calloc(1,sizeof(double*));

    c2.init(m);
    c2.setZero();

    nOfNewtonSteps = (int*)calloc(1,sizeof(int));
    eta            = new double*[4];
    eta2           = new double*[4];

    nOfNewtonSteps[0] = 0;

    psi   [0] = (double*)calloc(nstep,sizeof(double));
    gamma [0] = (double*)calloc(nstep,sizeof(double));

    nablaY.init( nstep, m );
    nablaY.setZero();

    nablaY_.init( nstep, m );
    nablaY_.setZero();

    phi.init( nstep, m );
    phi.setZero();

    delta.init( nstep, m );
    delta.setZero();


    for( run1 = 0; run1 < nstep; run1++ ){
         psi_    [run1] = 0.0;
         psi  [0][run1] = 0.0;
         gamma[0][run1] = 0.0;
    }

    maxNM = 1;
    M     = (Matrix**)calloc(maxNM,sizeof(Matrix*));
    M_index   = (int*)calloc(maxNM,sizeof(int));

    M_index[0] = 0;
    M      [0] = 0;
    nOfM       = 0;


    for( run1 = 0; run1 < 4; run1++ ){
        eta [run1] = new double[m];
        eta2[run1] = new double[m];
        for( run2 = 0; run2 < m; run2++ ){
            eta [run1][run2] = 0.0;
            eta2[run1][run2] = 0.0;
        }
    }


    F  = new double[m];
    F2 = new double[m];


    // INTERNAL INDEX LISTS:
    // ---------------------
    diff_index = new int[m];

    for( run1 = 0; run1 < md; run1++ ){
        diff_index[run1] = rhs->getStateEnumerationIndex( run1 );
        if( diff_index[run1] == rhs->getNumberOfVariables() ){
            diff_index[run1] = diff_index[run1] + 1 + run1;
        }
    }

    ddiff_index = new int[md];

    for( run1 = 0; run1 < md; run1++ ){
        ddiff_index[run1] = rhs->index( VT_DDIFFERENTIAL_STATE, run1 );
        if( ddiff_index[run1] == rhs->getNumberOfVariables() ){
            ddiff_index[run1] = ddiff_index[run1] + 1 + md + run1;
        }
    }


    alg_index = new int[ma];

    for( run1 = 0; run1 < ma; run1++ ){
        alg_index [run1]    = rhs->index( VT_ALGEBRAIC_STATE, run1 );
        diff_index[md+run1] = alg_index [run1];
    }

    control_index       = new int[mu ];

    for( run1 = 0; run1 < mu; run1++ ){
        control_index[run1] = rhs->index( VT_CONTROL, run1 );
    }

    parameter_index     = new int[mp ];

    for( run1 = 0; run1 < mp; run1++ ){
        parameter_index[run1] = rhs->index( VT_PARAMETER, run1 );
    }

    int_control_index   = new int[mui];

    for( run1 = 0; run1 < mui; run1++ ){
        int_control_index[run1] = rhs->index( VT_INTEGER_CONTROL, run1 );
    }

    int_parameter_index = new int[mpi];

    for( run1 = 0; run1 < mpi; run1++ ){
        int_parameter_index[run1] = rhs->index( VT_INTEGER_PARAMETER, run1 );
    }

    disturbance_index   = new int[mw ];

    for( run1 = 0; run1 < mw; run1++ ){
        disturbance_index[run1] = rhs->index( VT_DISTURBANCE, run1 );
    }

    time_index = rhs->index( VT_TIME, 0 );


    // OTHERS:
    // -------
    diff_scale.init(m);

    for( run1 = 0; run1 < md; run1++ ){
        diff_scale(run1) = rhs->scale( VT_DIFFERENTIAL_STATE, run1 );
    }
    for( run1 = 0; run1 < ma; run1++ ){
        diff_scale(md+run1) = rhs->scale( VT_ALGEBRAIC_STATE, run1 );
    }

    initial_guess    = new double[m];

    for( run1 = 0; run1 < m; run1++ ){

        initial_guess[run1] = 0.0;
    }


    // SENSITIVITIES:
    // --------------
    c2G .init(md);
    c2G2.init(md);
    c2G3.init(md);

    phiG.init(nstep,m);
    phiG2.init(nstep,m);
    phiG3.init(nstep,m);

    deltaG.init(nstep,m);
    deltaG2.init(nstep,m);
    deltaG3.init(nstep,m);

    kH         = new double**[dim];
    for( run1 = 0; run1 < dim; run1++ ){
         kH[run1] = new double*[nstep];
         for( run2 = 0; run2 < nstep; run2++ ){
             kH[run1][run2] = new double[ndir];
         }
    }
    zH         = new double*[6];
    for( run1 = 0; run1 < 6; run1++ ){
         zH[run1] = new double[ndir];
    }

    kH2        = new double**[dim];
    for( run1 = 0; run1 < dim; run1++ ){
         kH2[run1] = new double*[nstep];
         for( run2 = 0; run2 < nstep; run2++ ){
             kH2[run1][run2] = new double[ndir];
         }
    }
    zH2        = new double*[6];
    for( run1 = 0; run1 < 6; run1++ ){
         zH2[run1] = new double[ndir];
    }
    kH3        = new double**[dim];
    for( run1 = 0; run1 < dim; run1++ ){
         kH3[run1] = new double*[nstep];
         for( run2 = 0; run2 < nstep; run2++ ){
             kH3[run1][run2] = new double[ndir];
         }
    }
    zH3        = new double*[6];
    for( run1 = 0; run1 < 6; run1++ ){
         zH3[run1] = new double[ndir];
    }


    // STORAGE:
    // --------
    maxAlloc = 1;
}


void AcadoIntegratorBackend::initializeVariables(){

    initialAlgebraicResiduum = 0  ;
    relaxationConstant       = 0.0;

    dim = 0; A = 0; b4 = 0; b5 = 0; c = 0;

    eta4  = 0; eta5  = 0;
    eta4_ = 0; eta5_ = 0;

    k = 0; k2 = 0; l = 0; l2 = 0;

    iseed = 0; nstep = 0; psi   = 0;
    psi_  = 0; gamma = 0; eta   = 0;
    eta2  = 0; x     = 0; t     = 0;

    nOfNewtonSteps = 0;
    maxNM = 0; M = 0; M_index = 0; nOfM = 0;

    F  = 0; F2 = 0;

    initial_guess = 0;

    ndir = 0; G = 0; etaG = 0;

    G2  = 0; G3   = 0; etaG2 = 0; etaG3 = 0;
    H   = 0; etaH = 0; kH    = 0; zH    = 0;
    H2  = 0; H3   = 0; etaH2 = 0; etaH3 = 0;
    kH2 = 0; zH2  = 0; kH3   = 0; zH3   = 0;

    maxAlloc = 0;

    nFcnEvaluations = 0;
    nJacEvaluations = 0;
}


void AcadoIntegratorBackend::constructAll( const AcadoIntegratorBackend& arg ){

    int run1, run2, run3;

    rhs = new DifferentialEquation( *arg.rhs );

    m   = arg.m;
    ma  = arg.ma;
    mdx = arg.mdx;
    mn  = arg.mn;
    mu  = arg.mu;
    mui = arg.mui;
    mp  = arg.mp;
    mpi = arg.mpi;
    mw  = arg.mw;
    md  = m-ma;

    initialAlgebraicResiduum = new double[m]         ;
    relaxationConstant       = arg.relaxationConstant;


    // BUTCHER TABLEAU:
    // ----------------
    dim = arg.dim;
    A  = new double*[dim];
    b4 = new double [dim];
    b5 = new double [dim];
    c  = new double [dim];

    for( run1 = 0; run1 < dim; run1++ ){
        A[run1] = new double[dim];
    }

    initializeButcherTableau();


    // RK-STARTER:
    // -----------
    ndir = rhs->getNumberOfVariables() + 1 + 2*md;

    eta4  = new double [m];
    eta5  = new double [m];
    eta4_ = new double [m];
    eta5_ = new double [m];

    for( run1 = 0; run1 < m; run1++ ){

        eta4 [run1] = 0.0;
        eta5 [run1] = 0.0;
        eta4_[run1] = 0.0;
        eta5_[run1] = 0.0;
    }


    k     = new double**[4];
    k2    = new double**[4];
    l     = new double**[4];
    l2    = new double**[4];

    for( run3 = 0; run3 < 4; run3++ ){

        k   [run3]  = new double*[dim];
        k2  [run3]  = new double*[dim];
        l   [run3]  = new double*[dim];
        l2  [run3]  = new double*[dim];

        for( run1 = 0; run1 < dim; run1++ ){
            k    [run3][run1] = new double[m];
            k2   [run3][run1] = new double[m];

            for( run2 = 0; run2 < m; run2++ ){
                k [run3][run1][run2] = 0.0;
                k2[run3][run1][run2] = 0.0;
            }

            l    [run3][run1] = new double[ndir];
            l2   [run3][run1] = new double[ndir];

            for( run2 = 0; run2 < ndir; run2++ ){
                l [run3][run1][run2] = 0.0;
                l2[run3][run1][run2] = 0.0;
            }
        }
    }

    iseed = new double[ndir];

    x     = new double [ndir];

    for( run1 = 0; run1 < ndir; run1++ ){
        x[run1]     = 0.0;
        iseed[run1] = 0.0;
    }

    t     = 0.0;


    // BDF-METHOD:
    // -----------
    nstep     = 5;
    psi       = (double**)calloc(1,sizeof(double*));
    psi_      = new double[nstep];
    gamma     = (double**)calloc(1,sizeof(double*));

    c2.init(m);
    c2.setZero();

    nOfNewtonSteps = (int*)calloc(1,sizeof(int));
    eta            = new double*[4];
    eta2           = new double*[4];

    nOfNewtonSteps[0] = 0;

    psi   [0] = (double*)calloc(nstep,sizeof(double));
    gamma [0] = (double*)calloc(nstep,sizeof(double));

    nablaY.init( nstep, m );
    nablaY.setZero();

    nablaY_.init( nstep, m );
    nablaY_.setZero();

    phi.init( nstep, m );
    phi.setZero();

    delta.init( nstep, m );
    delta.setZero();

    for( run1 = 0; run1 < nstep; run1++ ){

         psi_    [run1] = 0.0;
         psi  [0][run1] = 0.0;
         gamma[0][run1] = 0.0;
    }


    maxNM = 1;
    M     = (Matrix**)calloc(maxNM,sizeof(Matrix*));
    M_index   = (int*)calloc(maxNM,sizeof(int));

    M_index[0] = 0;
    M      [0] = 0;
    nOfM       = 0;

    las = arg.las;

    for( run1 = 0; run1 < 4; run1++ ){
        eta [run1] = new double[m];
        eta2[run1] = new double[m];
        for( run2 = 0; run2 < m; run2++ ){
            eta [run1][run2] = 0.0;
            eta2[run1][run2] = 0.0;
        }
    }

    F  = new double[m];
    F2 = new double[m];


    // SETTINGS:
    // ---------
    h    = (double*)calloc(1,sizeof(double));
    h[0]  = 0.001    ;
    hini  = 0.001    ;
    hmin = 0.000001  ;
    hmax = 1.0e10    ;

    tune  = 0.5      ;
    TOL   = 0.000001 ;


    // INTERNAL INDEX LISTS:
    // ---------------------
    diff_index = new int[m];

    for( run1 = 0; run1 < md; run1++ ){
        diff_index[run1] = rhs->index( VT_DIFFERENTIAL_STATE, run1 );
        if( diff_index[run1] == rhs->getNumberOfVariables() ){
            diff_index[run1] = diff_index[run1] + 1 + run1;
        }
    }

    ddiff_index = new int[md];

    for( run1 = 0; run1 < md; run1++ ){
        ddiff_index[run1] = rhs->index( VT_DDIFFERENTIAL_STATE, run1 );
        if( ddiff_index[run1] == rhs->getNumberOfVariables() ){
            ddiff_index[run1] = ddiff_index[run1] + 1 + md + run1;
        }
    }

    alg_index = new int[ma];

    for( run1 = 0; run1 < ma; run1++ ){
        alg_index [run1]    = rhs->index( VT_ALGEBRAIC_STATE, run1 );
        diff_index[md+run1] = alg_index [run1];
    }

    control_index       = new int[mu ];

    for( run1 = 0; run1 < mu; run1++ ){
        control_index[run1] = rhs->index( VT_CONTROL, run1 );
    }

    parameter_index     = new int[mp ];

    for( run1 = 0; run1 < mp; run1++ ){
        parameter_index[run1] = rhs->index( VT_PARAMETER, run1 );
    }

    int_control_index   = new int[mui];

    for( run1 = 0; run1 < mui; run1++ ){
        int_control_index[run1] = rhs->index( VT_INTEGER_CONTROL, run1 );
    }

    int_parameter_index = new int[mpi];

    for( run1 = 0; run1 < mpi; run1++ ){
        int_parameter_index[run1] = rhs->index( VT_INTEGER_PARAMETER, run1 );
    }

    disturbance_index   = new int[mw ];

    for( run1 = 0; run1 < mw; run1++ ){
        disturbance_index[run1] = rhs->index( VT_DISTURBANCE, run1 );
    }

    time_index = rhs->index( VT_TIME, 0 );


    // OTHERS:
    // -------
    maxNumberOfSteps = 1000;
    count            = 0   ;
    count2           = 0   ;
    count3           = 0   ;

    diff_scale.init(m);

    for( run1 = 0; run1 < md; run1++ ){
        diff_scale(run1) = rhs->scale( VT_DIFFERENTIAL_STATE, run1 );
    }
    for( run1 = 0; run1 < ma; run1++ ){
        diff_scale(md+run1) = rhs->scale( VT_ALGEBRAIC_STATE, run1 );
    }

    initial_guess    = new double[m];

    for( run1 = 0; run1 < m; run1++ ){

        initial_guess[run1] = 0.0;
    }


    // PRINT-LEVEL:
    // ------------
    PrintLevel = LOW;


    // SENSITIVITIES:
    // --------------
    nFDirs     = 0   ;
    nBDirs     = 0   ;

    nFDirs2    = 0   ;
    nBDirs2    = 0   ;

    G          = NULL;
    etaG       = NULL;

    G2         = NULL;
    G3         = NULL;
    etaG2      = NULL;
    etaG3      = NULL;

    c2G.init(md);
    c2G2.init(md);
    c2G3.init(md);

    phiG.init(nstep,m);
    phiG2.init(nstep,m);
    phiG3.init(nstep,m);

    deltaG.init(nstep,m);
    deltaG2.init(nstep,m);
    deltaG3.init(nstep,m);


    H          = NULL;
    etaH       = NULL;

    kH         = new double**[dim];
    for( run1 = 0; run1 < dim; run1++ ){
         kH[run1] = new double*[nstep];
         for( run2 = 0; run2 < nstep; run2++ ){
             kH[run1][run2] = new double[ndir];
         }
    }
    zH         = new double*[6];
    for( run1 = 0; run1 < 6; run1++ ){
         zH[run1] = new double[ndir];
    }

    H2         = NULL;
    H3         = NULL;
    etaH2      = NULL;
    etaH3      = NULL;

    kH2        = new double**[dim];
    for( run1 = 0; run1 < dim; run1++ ){
         kH2[run1] = new double*[nstep];
         for( run2 = 0; run2 < nstep; run2++ ){
             kH2[run1][run2] = new double[ndir];
         }
    }
    zH2        = new double*[6];
    for( run1 = 0; run1 < 6; run1++ ){
         zH2[run1] = new double[ndir];
    }
    kH3        = new double**[dim];
    for( run1 = 0; run1 < dim; run1++ ){
         kH3[run1] = new double*[nstep];
         for( run2 = 0; run2 < nstep; run2++ ){
             kH3[run1][run2] = new double[ndir];
         }
    }
    zH3        = new double*[6];
    for( run1 = 0; run1 < 6; run1++ ){
         zH3[run1] = new double[ndir];
    }

    // THE STATE OF AGGREGATION:
    // -------------------------
    soa        = SOA_UNFROZEN;


    // STORAGE:
    // --------
    maxAlloc = 1;

    nFcnEvaluations = 0;
    nJacEvaluations = 0;
}


Integrator* AcadoIntegratorBackend::clone() const{

    return new AcadoIntegratorBackend(*this);
}


void AcadoIntegratorBackend::deleteAll(){

    int run1, run2;

    if( initialAlgebraicResiduum != 0 )
        delete[] initialAlgebraicResiduum;


    // BUTCHER-
    // TABLEAU:
    // ----------
    for( run1 = 0; run1 < dim; run1++ ){
        delete[] A[run1];
    }

    delete[] A;
    delete[] b4;
    delete[] b5;
    delete[] c ;


    // RK-ALGORITHM:
    // -------------
    if( eta4 != NULL ){
        delete[] eta4;
    }
    if( eta4_ != NULL ){
        delete[] eta4_;
    }
    if( eta5 != NULL ){
        delete[] eta5;
    }
    if( eta5_ != NULL ){
        delete[] eta5_;
    }



    for( run2 = 0; run2 < 4; run2++ ){

         for( run1 = 0; run1 < dim; run1++ ){
            if( k[run2]  != NULL )
                delete[] k[run2][run1] ;
            if( k2[run2] != NULL )
                delete[] k2[run2][run1];
            if( l[run2]  != NULL )
                delete[] l[run2][run1] ;
            if( l2[run2] != NULL )
                delete[] l2[run2][run1];
         }

        if( k != NULL )
            delete[] k[run2] ;

        if( k2 != NULL )
            delete[] k2[run2];

        if( l != NULL )
            delete[] l[run2] ;

        if( l2 != NULL )
            delete[] l2[run2];
   }


    if( k != NULL )
        delete[] k;

    if( k2!= NULL )
        delete[] k2;

    if( l != NULL )
        delete[] l ;

    if( l2!= NULL )
        delete[] l2;

    if( iseed != NULL ){
        delete[] iseed;
    }

    if( x != NULL )
        delete[] x;



    for( run1 = 0; run1 < maxAlloc; run1++ ){

        if( psi[run1] != NULL ){
            free(psi[run1]);
        }
        if( gamma[run1] != NULL ){
            free(gamma[run1]);
        }
    }

    if( psi != NULL ){
        free(psi);
    }

    if( psi_ != NULL ){
        delete[] psi_;
    }

    if( gamma != NULL )
        free(gamma);

    if( nOfNewtonSteps != NULL )
        free(nOfNewtonSteps);

    for( run1 = 0; run1 < 4; run1++ ){

        if( eta != NULL )
            delete[] eta[run1];
        if( eta2 != NULL )
            delete[] eta2[run1];
    }

    delete[] eta    ;
    delete[] eta2   ;

    for( run1 = 0; run1 < maxNM; run1++ ){
         if( M[run1] != 0  )
             delete M[run1];
    }

    if( M != NULL ){
        free(M);
        free(M_index);
    }

    if( F != NULL )
        delete[] F;
    if( F2 != NULL )
        delete[] F2;


    // OTHERS:
    // -------
    if( initial_guess != NULL ){
        delete[] initial_guess;
    }

    // SENSITIVITIES:
    // --------------

    if( G  != NULL )
        delete[] G;

    if( etaG  != NULL )
        delete[] etaG;

    if( G2  != NULL )
        delete[] G2;

    if( G3  != NULL )
        delete[] G3;

    if( etaG2  != NULL )
        delete[] etaG2;

    if( etaG3  != NULL )
        delete[] etaG3;


    // ----------------------------------------

    if( H  != NULL )
        delete[] H;

    // ----------------------------------------

    if( H2  != NULL )
        delete[] H2;
    if( H3  != NULL )
        delete[] H3;

    for( run1 = 0; run1 < nstep; run1++ ){
        if( etaH != NULL )
            delete[] etaH[run1];
        if( etaH2 != NULL )
            delete[] etaH2[run1];
        if( etaH3 != NULL )
            delete[] etaH3[run1];
    }
    if( etaH != NULL )
        delete[] etaH;
    if( etaH2 != NULL )
        delete[] etaH2;
    if( etaH3 != NULL )
        delete[] etaH3;


    for( run1 = 0; run1 < dim; run1++ ){
        for( run2 = 0; run2 < nstep; run2++ ){
            if( kH != NULL )
                delete[] kH[run1][run2];
            if( kH2 != NULL )
                delete[] kH2[run1][run2];
            if( kH3 != NULL )
                delete[] kH3[run1][run2];
        }
        if( kH != NULL )
            delete[] kH[run1];
        if( kH2 != NULL )
            delete[] kH2[run1];
        if( kH3 != NULL )
            delete[] kH3[run1];
    }
    if( kH != NULL )
        delete[] kH;
    if( kH2 != NULL )
        delete[] kH2;
    if( kH3 != NULL )
        delete[] kH3;


    for( run1 = 0; run1 < 6; run1++ ){
        if( zH != NULL )
            delete[] zH[run1];
        if( zH2 != NULL )
            delete[] zH2[run1];
        if( zH3 != NULL )
            delete[] zH3[run1];
    }
    if( zH != NULL )
        delete[] zH;
    if( zH2 != NULL )
        delete[] zH2;
    if( zH3 != NULL )
        delete[] zH3;
}


returnValue AcadoIntegratorBackend::freezeMesh(){


    if( soa != SOA_UNFROZEN ){
       if( PrintLevel != NONE ){
           return ACADOWARNING(RET_ALREADY_FROZEN);
       }
       return RET_ALREADY_FROZEN;
    }

    soa = SOA_FREEZING_MESH;
    return SUCCESSFUL_RETURN;
}


returnValue AcadoIntegratorBackend::freezeAll(){

    if( soa != SOA_UNFROZEN ){
       if( PrintLevel != NONE ){
           return ACADOWARNING(RET_ALREADY_FROZEN);
       }
       return RET_ALREADY_FROZEN;
    }

    soa = SOA_FREEZING_ALL;
    return SUCCESSFUL_RETURN;
}


returnValue AcadoIntegratorBackend::unfreeze(){

    int run1, run2;

    for( run1 = 1; run1 < maxAlloc; run1++ ){

        if( psi[run1] != NULL ){
            free(psi[run1]);
        }
        if( gamma[run1] != NULL ){
            free(gamma[run1]);
        }
    }

    maxAlloc = 1;

    psi            = (double**)realloc(psi,maxAlloc*sizeof(double*));
    gamma          = (double**)realloc(gamma,maxAlloc*sizeof(double*));
    nOfNewtonSteps = (int*)realloc(nOfNewtonSteps,(maxAlloc)*sizeof(int));

    for( run1 = 0; run1 < maxAlloc; run1++ ){
        nOfNewtonSteps[run1] = 0;
    }

    for( run1 = 0; run1 < maxAlloc; run1++ ){

        for( run2 = 0; run2 < 5; run2++ ){
             psi  [run1][run2] = 0.0;
             gamma[run1][run2] = 0.0;
        }
    }

    for( run1 = 0; run1 < maxNM; run1++ ){
         if( M[run1] != 0  )
             delete M[run1];
         M[run1] = 0;
    }

    maxNM = 1;
    nOfM  = 0;
    M       = (Matrix**)realloc(M,maxNM*sizeof(Matrix*));
    M_index = (int*)realloc(M_index,maxAlloc*sizeof(int));

    h = (double*)realloc(h,maxAlloc*sizeof(double));

    soa = SOA_UNFROZEN;

    return SUCCESSFUL_RETURN;
}


returnValue AcadoIntegratorBackend::evaluate( const Vector &x0  ,
                                     const Vector &xa  ,
                                     const Vector &p   ,
                                     const Vector &u   ,
                                     const Vector &w   ,
                                     const Grid   &t_    ){
    printf("Evaluate %p\n",this);
  
  

    int         run1;
    returnValue returnvalue;

    if( rhs == NULL ){
        return ACADOERROR(RET_TRIVIAL_RHS);
    }

    Integrator::initializeOptions();

    timeInterval = t_;

    xStore.init( md+ma, timeInterval );
    iStore.init( mn   , timeInterval );

    t             = timeInterval.getFirstTime();
    x[time_index] = timeInterval.getFirstTime();

    if( x0.isEmpty() == BT_TRUE ) return ACADOERROR(RET_MISSING_INPUTS);

    if( (int) x0.getDim() < md )
        return ACADOERROR(RET_INPUT_HAS_WRONG_DIMENSION);

    for( run1 = 0; run1 < md; run1++ ){
        eta4[run1]          = x0(run1);
        eta5[run1]          = x0(run1);
        xStore(0,run1)      = x0(run1);
        x[diff_index[run1]] = x0(run1);
    }

    if( soa != SOA_MESH_FROZEN && soa != SOA_EVERYTHING_FROZEN  ){
       h[0] = hini;

       if( timeInterval.getIntervalLength() - h[0] < EPS ){
           h[0] = timeInterval.getIntervalLength();
       }

       if( h[0] < 10.0*EPS )
           return ACADOERROR(RET_TO_SMALL_OR_NEGATIVE_TIME_INTERVAL);
    }

    if( ma > 0 ){

        if( (int) xa.getDim() < ma )
            return ACADOERROR(RET_INPUT_HAS_WRONG_DIMENSION);

        for( run1 = 0; run1 < ma; run1++ ){
            eta4[md+run1]          = xa(run1);
            eta5[md+run1]          = xa(run1);
            xStore(0,md+run1)      = xa(run1);
            x[diff_index[md+run1]] = xa(run1);
            initial_guess[md+run1] = xa(run1);
        }
    }

    if( nFDirs != 0 )
        for( run1 = 0; run1 < md; run1++ )
            etaG[run1] = fseed(diff_index[run1]);

    if( mp > 0 ){

        if( (int) p.getDim() < mp )
            return ACADOERROR(RET_INPUT_HAS_WRONG_DIMENSION);

        for( run1 = 0; run1 < mp; run1++ ){
            x[parameter_index[run1]] = p(run1);
        }
    }

    if( mu > 0 ){

        if( (int) u.getDim() < mu )
            return ACADOERROR(RET_INPUT_HAS_WRONG_DIMENSION);

        for( run1 = 0; run1 < mu; run1++ ){
            x[control_index[run1]] = u(run1);
        }
    }


    if( mw > 0 ){

        if( (int) w.getDim() < mw )
            return ACADOERROR(RET_INPUT_HAS_WRONG_DIMENSION);

        for( run1 = 0; run1 < mw; run1++ ){
            x[disturbance_index[run1]] = w(run1);
        }
    }

	// log initial step
	logCurrentIntegratorStep( x0,xa );

     // start the time measurement:
     // ---------------------------

        totalTime.start();
        nFcnEvaluations = 0;
        nJacEvaluations = 0;


     // initialize the scaling based on the initial states:
     // ---------------------------------------------------

        double atol;
        get( ABSOLUTE_TOLERANCE, atol );

        for( run1 = 0; run1 < m; run1++ )
            diff_scale(run1) = fabs(eta4[run1]) + atol/TOL;


    returnvalue = rhs[0].evaluate( 0, x, initialAlgebraicResiduum );

    if( returnvalue != SUCCESSFUL_RETURN ){

        totalTime.stop();
        return ACADOWARNING(returnvalue);
    }

     // PRINTING:
     // ---------
        if( PrintLevel == MEDIUM ){
            acadoPrintCopyrightNotice( "AcadoIntegratorBackend -- A BDF integrator (order 4)." );
        }
        if( PrintLevel == HIGH ){
            acadoPrintf("BDF: t = %.16e                          ", t );
            for( run1 = 0; run1 < md; run1++ ){
                acadoPrintf("x[%d] = %.16e  ", run1, eta4[run1] );
            }
            for( run1 = 0; run1 < ma; run1++ ){
                acadoPrintf("xa[%d] = %.16e  ", run1, eta4[md+run1] );
            }
            acadoPrintf("\n");
        }


    count3 = 0;

    returnvalue = rk_start();

    if( returnvalue != SUCCESSFUL_RETURN ){

        totalTime.stop();

        if( PrintLevel != NONE )
            return ACADOERROR(returnvalue);

        return returnvalue;
    }

    if( soa != SOA_MESH_FROZEN && soa != SOA_EVERYTHING_FROZEN  ){

       if( timeInterval.getLastTime() - t - h[0] < EPS ){
           h[0] = timeInterval.getLastTime() - t;
       }
    }

    returnvalue = RET_FINAL_STEP_NOT_PERFORMED_YET;

    count = 7;

    while( returnvalue == RET_FINAL_STEP_NOT_PERFORMED_YET && count <= maxNumberOfSteps ){

        returnvalue = step(count);
        count++;
    }

	// log final step
	logCurrentIntegratorStep( );

    count2 = count-1;


    for( run1 = 0; run1 < mn; run1++ )
        iStore( 0, run1 ) = iStore( 1, run1 );


    // stop the measurement of total time:
    // -----------------------------------
       totalTime.stop();


    // SET THE LOGGING INFORMATION:
    // ----------------------------------------------------------------------------------------

       setLast( LOG_TIME_INTEGRATOR                              , totalTime.getTime()           );
       setLast( LOG_NUMBER_OF_INTEGRATOR_STEPS                   , count-1                       );
       setLast( LOG_NUMBER_OF_INTEGRATOR_REJECTED_STEPS          , getNumberOfRejectedSteps()    );
       setLast( LOG_NUMBER_OF_INTEGRATOR_FUNCTION_EVALUATIONS    , nFcnEvaluations               );
       setLast( LOG_NUMBER_OF_BDF_INTEGRATOR_JACOBIAN_EVALUATIONS, nJacEvaluations               );
       setLast( LOG_TIME_INTEGRATOR_FUNCTION_EVALUATIONS         , functionEvaluation.getTime()  );
       setLast( LOG_TIME_BDF_INTEGRATOR_JACOBIAN_EVALUATION      , jacComputation.getTime()      );
       setLast( LOG_TIME_BDF_INTEGRATOR_JACOBIAN_DECOMPOSITION   , jacDecomposition.getTime()    );

    // ----------------------------------------------------------------------------------------



    if( count > maxNumberOfSteps ){

        if( PrintLevel != NONE )
            return ACADOERROR(RET_MAX_NUMBER_OF_STEPS_EXCEEDED);
        return RET_MAX_NUMBER_OF_STEPS_EXCEEDED;
    }


     // PRINTING:
     // ---------
        if( PrintLevel == MEDIUM ){

            acadoPrintf("\n Results at  t =  %.16e   : \n\n", t );
            for( run1 = 0; run1 < md; run1++ ){
                acadoPrintf("x[%d] = %.16e  ", run1, nablaY(0,run1) );
            }
            for( run1 = 0; run1 < ma; run1++ ){
                acadoPrintf("xa[%d] = %.16e  ", run1, nablaY(0,md+run1) );
            }
            acadoPrintf("\n");
            printBDFfinalResults();
        }
		
	int printIntegratorProfile = 0;
	get( PRINT_INTEGRATOR_PROFILE,printIntegratorProfile );
	
	if ( (BooleanType)printIntegratorProfile == BT_TRUE )
	{
		printRunTimeProfile( );
	}
	else
	{
		if( PrintLevel == MEDIUM  || PrintLevel == HIGH )
			acadoPrintf("BDF: number of steps:  %d\n", count-1 );
	}

    return returnvalue;
}


returnValue AcadoIntegratorBackend::setProtectedForwardSeed( const Vector &xSeed ,
                                                    const Vector &pSeed ,
                                                    const Vector &uSeed ,
                                                    const Vector &wSeed ,
                                                    const int    &order   ){

  
    if( order == 2 ){
        return setForwardSeed2( xSeed, pSeed, uSeed, wSeed );
    }
    if( order < 1 || order > 2 ){
        return ACADOERROR(RET_INPUT_OUT_OF_RANGE);
    }

    if( nBDirs > 0 ){
        return ACADOERROR(RET_INPUT_OUT_OF_RANGE);
    }

    int run2;

    if( G  != NULL ){
        delete[] G;
        G = NULL;
    }

    if( etaG  != NULL ){
        delete[] etaG;
        etaG = NULL;
    }

    nFDirs = 1;

    fseed.init(ndir);
    fseed.setZero();

    G = new double[ndir];

    for( run2 = 0; run2 < (ndir); run2++ )
         G[run2] = 0.0;

    etaG = new double[m];

    nablaG.init( nstep, m );


    if( xSeed.getDim() != 0 )
        for( run2 = 0; run2 < md; run2++ )
            fseed(diff_index[run2]) = xSeed(run2);

    if( pSeed.getDim() != 0 )
        for( run2 = 0; run2 < mp; run2++ ){
            fseed(parameter_index[run2]) = pSeed(run2);
            G    [parameter_index[run2]] = pSeed(run2);
        }

    if( uSeed.getDim() != 0 )
        for( run2 = 0; run2 < mu; run2++ ){
            fseed(control_index[run2]) = uSeed(run2);
            G    [control_index[run2]] = uSeed(run2);
        }

    if( wSeed.getDim() != 0 )
        for( run2 = 0; run2 < mw; run2++ ){
            fseed(disturbance_index[run2]) = wSeed(run2);
            G    [disturbance_index[run2]] = wSeed(run2);
        }

    return SUCCESSFUL_RETURN;
}


returnValue AcadoIntegratorBackend::setForwardSeed2( const Vector &xSeed ,
                                            const Vector &pSeed ,
                                            const Vector &uSeed ,
                                            const Vector &wSeed   ){

    int run2;


    if( G2  != NULL ){
        delete[] G2;
        G2 = NULL;
    }

    if( G3  != NULL ){
        delete[] G3;
        G3 = NULL;
    }

    if( etaG2  != NULL ){
        delete[] etaG2;
        etaG2 = NULL;
    }

    if( etaG3  != NULL ){
        delete[] etaG3;
        etaG3 = NULL;
    }


    nFDirs2 = 1;

    fseed2.init(ndir);
    fseed2.setZero();

    G2 = new double[ndir];
    G3 = new double[ndir];

    for( run2 = 0; run2 < (ndir); run2++ ){
         G2[run2] = 0.0;
         G3[run2] = 0.0;
    }

    etaG2 = new double[m];
    etaG3 = new double[m];

    nablaG2.init( nstep, m );
    nablaG3.init( nstep, m );

    if( xSeed.getDim() != 0 )
        for( run2 = 0; run2 < md; run2++ )
            fseed2(diff_index[run2]) = xSeed(run2);

    if( pSeed.getDim() != 0 )
        for( run2 = 0; run2 < mp; run2++ ){
            fseed2(parameter_index[run2]) = pSeed(run2);
            G2    [parameter_index[run2]] = pSeed(run2);
        }

    if( uSeed.getDim() != 0 )
        for( run2 = 0; run2 < mu; run2++ ){
            fseed2(control_index[run2]) = uSeed(run2);
            G2    [control_index[run2]] = uSeed(run2);
        }

    if( wSeed.getDim() != 0 )
        for( run2 = 0; run2 < mw; run2++ ){
            fseed2(disturbance_index[run2]) = wSeed(run2);
            G2    [disturbance_index[run2]] = wSeed(run2);
        }

    return SUCCESSFUL_RETURN;
}


returnValue AcadoIntegratorBackend::setProtectedBackwardSeed( const Vector &seed, const int &order ){

    if( order == 2 ){
        return setBackwardSeed2(seed);
    }
    if( order < 1 || order > 2 ){
        return ACADOERROR(RET_INPUT_OUT_OF_RANGE);
    }

    if( nFDirs > 0 ){
        return ACADOERROR(RET_INPUT_OUT_OF_RANGE);
    }

    int run1, run2, run3;

    if( H  != NULL ){
        delete[] H;
        H = NULL;
    }

    for( run1 = 0; run1 < nstep; run1++ ){
        if( etaH != NULL )
            delete[] etaH[run1];
    }
    if( etaH != NULL ){
        delete[] etaH;
        etaH = 0;
    }

    nBDirs = 1;

    bseed.init(m);
    bseed.setZero();
    H = new double[m];

    if( seed.getDim() != 0 )
        for( run2 = 0; run2 < md; run2++ )
            bseed(run2) = seed(run2);

    etaH = new double*[nstep];
    for( run3 = 0; run3 < nstep; run3++ ){
        etaH[run3] = new double[ndir];
        for( run2 = 0; run2 < ndir; run2++ ){
            etaH[run3][run2] = 0.0;
        }
    }

    nablaH.init( nstep, ndir );
    nablaH.setZero();

    nablaH_.init( nstep, ndir );
    nablaH_.setZero();

    deltaH.init( nstep, ndir );
    deltaH.setZero();

    c2H.init(ndir);
    c2H.setZero();

    return SUCCESSFUL_RETURN;
}


returnValue AcadoIntegratorBackend::setBackwardSeed2( const Vector &seed ){
assert(0);
  
    int run1, run2, run3;

    if( H2  != NULL ){
        delete[] H2;
        H2 = NULL;
    }
    if( H3  != NULL ){
        delete[] H3;
        H3 = NULL;
    }

    for( run1 = 0; run1 < nstep; run1++ ){
        if( etaH2 != NULL )
            delete[] etaH2[run1];
        if( etaH3 != NULL )
            delete[] etaH3[run1];
    }
    if( etaH2 != NULL ){
        delete[] etaH2;
        etaH2 = 0;
    }
    if( etaH3 != NULL ){
        delete[] etaH3;
        etaH3 = 0;
    }

    nBDirs2 = 1;

    bseed2.init(m);

    H2 = new double[m];
    H3 = new double[m];

    if( seed.getDim() != 0 ){
        for( run2 = 0; run2 < md; run2++ ){
            bseed2(run2) = seed(run2);
        }
    }

    etaH2 = new double*[nstep];
    for( run3 = 0; run3 < nstep; run3++ ){
        etaH2[run3] = new double[ndir];
        for( run2 = 0; run2 < ndir; run2++ ){
             etaH2[run3][run2] = 0.0;
        }
    }
    etaH3 = new double*[nstep];
    for( run3 = 0; run3 < nstep; run3++ ){
        etaH3[run3] = new double[ndir];
        for( run2 = 0; run2 < ndir; run2++ ){
             etaH3[run3][run2] = 0.0;
        }
    }

    nablaH2.init( nstep, ndir );
    nablaH2.setZero();

    nablaH3.init( nstep, ndir );
    nablaH3.setZero();

    nablaH2_.init( nstep, ndir );
    nablaH2_.setZero();

    nablaH3_.init( nstep, ndir );
    nablaH3_.setZero();

    deltaH2.init( nstep, ndir );
    deltaH2.setZero();

    deltaH3.init( nstep, ndir );
    deltaH3.setZero();

    c2H2.init( ndir );
    c2H2.setZero();

    c2H3.init( ndir );
    c2H3.setZero();

    return SUCCESSFUL_RETURN;
}


returnValue AcadoIntegratorBackend::evaluateSensitivities(){
    printf("eval_sens %p\n",this);

    int         run1;
    returnValue returnvalue;

    if( rhs == NULL ){
        return ACADOERROR(RET_TRIVIAL_RHS);
    }

    if( soa != SOA_EVERYTHING_FROZEN ){
        return ACADOERROR(RET_NOT_FROZEN);
    }

    if( nBDirs2 == 0 && nFDirs != 0 ){
        t = timeInterval.getFirstTime();
        dxStore.init ( md+ma, timeInterval );
        for( run1 = 0; run1 < md; run1++ ){
             etaG[run1] = fseed(diff_index[run1]);
        }
    }

    nablaH.setZero();

    if( nBDirs != 0 ){
        for( run1 = 0; run1 < md; run1++ ){
            nablaH(0,diff_index[run1]) = bseed(run1);
        }
    }

    if( nFDirs2 != 0 ){
        t = timeInterval.getFirstTime();
        ddxStore.init( md+ma, timeInterval );
        for( run1 = 0; run1 < md; run1++ ){
            etaG2[run1] = fseed2(diff_index[run1]);
            etaG3[run1] = 0.0;
        }
    }

    nablaH2.setZero();
    nablaH3.setZero();

    if( nBDirs2 != 0 ){
        for( run1 = 0; run1 < md; run1++ ){
            nablaH2(0,diff_index[run1]) = bseed2(run1);
        }
    }


    returnvalue = RET_FINAL_STEP_NOT_PERFORMED_YET;

    if( nBDirs > 0 || nBDirs2 > 0 ){

        t = timeInterval.getLastTime();
        h[0] = 1.0;

        int oldCount = count;

        count--;
        while( returnvalue == RET_FINAL_STEP_NOT_PERFORMED_YET && count >= 7 ){

            returnvalue = step( count );
            count--;
        }

        if( returnvalue != SUCCESSFUL_RETURN &&
            returnvalue != RET_FINAL_STEP_NOT_PERFORMED_YET ){
            count = oldCount;
            if( PrintLevel != NONE )
                return ACADOERROR(returnvalue);
            return returnvalue;
        }

        returnvalue = rk_start();

        if( returnvalue == SUCCESSFUL_RETURN ){

           // PRINTING:
           // ---------
              if( PrintLevel == MEDIUM )
                  printBDFfinalResults();

            count = oldCount;

            return SUCCESSFUL_RETURN;
        }
        count = oldCount;
    }
    else{

        t = timeInterval.getFirstTime();

        returnvalue = rk_start();

        if( returnvalue != SUCCESSFUL_RETURN ){
            if( PrintLevel != NONE )
                return ACADOERROR(returnvalue);
            return returnvalue;
        }

        returnvalue = RET_FINAL_STEP_NOT_PERFORMED_YET;
        count = 7;
        while( returnvalue == RET_FINAL_STEP_NOT_PERFORMED_YET &&
               count <= maxNumberOfSteps ){

            returnvalue = step(count);
            count++;
        }

        if( nBDirs2 == 0 && nFDirs != 0 )
            for( run1 = 0; run1 < m; run1++ )
                dxStore( 0, run1 ) = dxStore( 1, run1 );

        if( nFDirs2 != 0 )
            for( run1 = 0; run1 < m; run1++ )
                ddxStore( 0, run1 ) = ddxStore( 1, run1 );

        if( count > maxNumberOfSteps ){
            if( PrintLevel != NONE )
                return ACADOERROR(RET_MAX_NUMBER_OF_STEPS_EXCEEDED);
            return RET_MAX_NUMBER_OF_STEPS_EXCEEDED;
        }

        // PRINTING:
        // ---------
           if( PrintLevel == MEDIUM )
               printBDFfinalResults();

        return SUCCESSFUL_RETURN;
    }

    if( PrintLevel != NONE )
        return ACADOERROR(returnvalue);

    return returnvalue;
}


returnValue AcadoIntegratorBackend::step(int number_){
  printf(".");
    int run1;
    double E = EPS;

     if( soa != SOA_EVERYTHING_FROZEN && soa != SOA_MESH_FROZEN ){

         int oldAlloc = maxAlloc;

         if( number_ >= maxAlloc ){
             maxAlloc = 2*maxAlloc;

             psi            = (double**)realloc(psi,maxAlloc*sizeof(double*));
             gamma          = (double**)realloc(gamma,maxAlloc*sizeof(double*));
             nOfNewtonSteps = (int*)realloc(nOfNewtonSteps,(maxAlloc)*sizeof(int));

             for( run1 = oldAlloc; run1 < maxAlloc; run1++ ){
                 nOfNewtonSteps[run1] = 0;
             }
             for( run1 = oldAlloc; run1 < maxAlloc; run1++ ){
                 psi  [run1] = (double*)calloc(nstep,sizeof(double));
                 gamma[run1] = (double*)calloc(nstep,sizeof(double));
             }

             M_index      = (int*)realloc(M_index,maxAlloc*sizeof(int));
             h = (double*)realloc(h,maxAlloc*sizeof(double));
         }
     }


    if( soa == SOA_EVERYTHING_FROZEN || soa == SOA_MESH_FROZEN ){
        rel_time_scale = h[0];
        h[0] = h[number_];
    }

    if( soa == SOA_FREEZING_MESH ||
        soa == SOA_MESH_FROZEN   ||
        soa == SOA_UNFROZEN      ){

        if( number_ == 7 ){
           E = determinePredictor(7, BT_TRUE );
        }
        else{
           E = determinePredictor(7, BT_FALSE );
        }
    }
    if( soa == SOA_FREEZING_ALL ){

        if( number_ == 7 ){
           E = determinePredictor(number_, BT_TRUE );
        }
        else{
           E = determinePredictor(number_, BT_FALSE );
        }
    }


    if( soa != SOA_EVERYTHING_FROZEN && soa != SOA_MESH_FROZEN ){

        int number_of_rejected_steps = 0;

        if( E < 0.0 ){
            return ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
        }

        // REJECT THE STEP IF GIVEN TOLERANCE IS NOT ACHIEVED:
        // ---------------------------------------------------
        while( E >= TOL ){

            if( PrintLevel == HIGH ){

                acadoPrintf("STEP REJECTED: error estimate           = %.16e \n", E   );
                acadoPrintf("               required local tolerance = %.16e \n", TOL );
            }

            number_of_rejected_steps++;

            if( soa == SOA_FREEZING_MESH ||
                soa == SOA_UNFROZEN      ){

                psi[7][0] = psi_[0];
                psi[7][1] = psi_[1];
                psi[7][2] = psi_[2];
                psi[7][3] = psi_[3];
            }
            if( soa == SOA_FREEZING_ALL ){

                psi[number_][0] = psi_[0];
                psi[number_][1] = psi_[1];
                psi[number_][2] = psi_[2];
                psi[number_][3] = psi_[3];
            }


            if( h[0] <= hmin + EPS ){
                return ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }
            h[0] = 0.5*h[0];
            if( E > 0.9*INFTY ) h[0] = 0.2*h[0];

            if( h[0] < hmin ){
                h[0] = hmin;
            }

            if( soa == SOA_FREEZING_MESH ||
                soa == SOA_UNFROZEN      ){

                E = determinePredictor( 7, BT_FALSE );
            }
            if( soa == SOA_FREEZING_ALL ){

                E = determinePredictor( number_, BT_FALSE );
            }

            if( E < 0.0 ){
                return ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }
        }

        count3 += number_of_rejected_steps;
    }

    // PROCEED IF THE STEP IS ACCEPTED:
    // --------------------------------

    if( soa != SOA_EVERYTHING_FROZEN ) nablaY = nablaY_;

	// log current step
	logCurrentIntegratorStep( );


     // compute forward derivatives if requested:
     // ------------------------------------------

     if( nFDirs > 0 && nBDirs2 == 0 && nFDirs2 == 0 ){

         if( nBDirs != 0 ){
             return ACADOERROR(RET_WRONG_DEFINITION_OF_SEEDS);
         }

         if( soa == SOA_FREEZING_ALL || soa == SOA_EVERYTHING_FROZEN ){
             determineBDFEtaGForward(number_);
         }
         else{
             determineBDFEtaGForward(7);
         }
     }
     if( nBDirs > 0 ){

         if( soa != SOA_EVERYTHING_FROZEN ){
             return ACADOERROR(RET_NOT_FROZEN);
         }
         if( nFDirs != 0 || nBDirs2 != 0 || nFDirs2 != 0 ){
             return ACADOERROR(RET_WRONG_DEFINITION_OF_SEEDS);
         }
         determineBDFEtaHBackward(number_);
     }
     if( nFDirs2 > 0 ){

         if( soa != SOA_EVERYTHING_FROZEN ){
             return ACADOERROR(RET_NOT_FROZEN);
         }
         if( nBDirs != 0 || nBDirs2 != 0 || nFDirs != 1 ){
             return ACADOERROR(RET_WRONG_DEFINITION_OF_SEEDS);
         }
         determineBDFEtaGForward2(number_);
     }
     if( nBDirs2 > 0 ){

         if( soa != SOA_EVERYTHING_FROZEN ){
             return ACADOERROR(RET_NOT_FROZEN);
         }
         if( nBDirs != 0 || nFDirs2 != 0 || nFDirs != 1 ){
             return ACADOERROR(RET_WRONG_DEFINITION_OF_SEEDS);
         }

         determineBDFEtaHBackward2(number_);
     }


     // increase the time:
     // ----------------------------------------------

     if( nBDirs > 0 || nBDirs2 > 0 ){

         t = t - h[0];
     }
     else{

         t = t + h[0];
     }

     // PRINTING:
     // ---------
        printBDFIntermediateResults();


     // STORAGE:
     // --------

     if( soa == SOA_FREEZING_MESH || soa == SOA_FREEZING_ALL ){

         if( number_ >= maxAlloc){

             maxAlloc = 2*maxAlloc;
             h = (double*)realloc(h,maxAlloc*sizeof(double));
         }

         h[number_] = h[0];
     }

     if( nFDirs == 0 && nBDirs  == 0 && nFDirs2 == 0 && nBDirs == 0 ){

         interpolate( number_, nablaY, xStore );

         int i1 = timeInterval.getFloorIndex( t             );
         int i2 = timeInterval.getFloorIndex( t + c[6]*h[0] );
         int jj;

         for( jj = i1+1; jj <= i2; jj++ )
             for( run1 = 0; run1 < mn; run1++ )
                 iStore( jj, run1 ) = x[rhs->index( VT_INTERMEDIATE_STATE, run1 )];
     }
     if( nFDirs  > 0 && nBDirs2 == 0 && nFDirs2 == 0 ) interpolate( number_, nablaG , dxStore  );
     if( nFDirs2 > 0                                 ) interpolate( number_, nablaG3, ddxStore );


     if( nBDirs == 0 || nBDirs2 == 0 ){

     // Stop the algorithm if  t >= te:
     // ----------------------------------------------
        if( t >= timeInterval.getLastTime() - EPS ){
            x[time_index] = timeInterval.getLastTime();
            for( run1 = 0; run1 < m; run1++ ){
                x[diff_index[run1]] = nablaY(0,run1);
            }

            if( soa == SOA_FREEZING_MESH ){
                soa = SOA_MESH_FROZEN;
            }
            if( soa == SOA_FREEZING_ALL ){
                soa = SOA_EVERYTHING_FROZEN;
            }

            return SUCCESSFUL_RETURN;
        }
     }


     if( soa != SOA_EVERYTHING_FROZEN && soa != SOA_MESH_FROZEN ){


     // recompute the scaling based on the actual states:
     // -------------------------------------------------

        double atol;
        get( ABSOLUTE_TOLERANCE, atol );

        for( run1 = 0; run1 < m; run1++ )
            diff_scale(run1) = fabs(nablaY(0,run1)) + atol/TOL;


     // apply a numeric stabilization of the step size control:
     // -------------------------------------------------------
        double Emin = 1e-3*sqrt(TOL)*pow( hini, 2.5 );

        if( E < Emin     ) E = Emin    ;
        if( E < 10.0*EPS ) E = 10.0*EPS;



     // determine the new step size:
     // ----------------------------------------------
        rel_time_scale = h[0];
        h[0] = h[0]*pow( tune*(TOL/E) , 0.20 );

        if( h[0] / rel_time_scale > 2.0 ) h[0] = 2.0*rel_time_scale;

        if( h[0] > hmax ){
          h[0] = hmax;
        }
        if( h[0] < hmin ){
          h[0] = hmin;
        }

        if( t + h[0] >= timeInterval.getLastTime() ){
          h[0] = timeInterval.getLastTime()-t;
        }
    }

    return RET_FINAL_STEP_NOT_PERFORMED_YET;
}



returnValue AcadoIntegratorBackend::stop(){
assert(0);

    return ACADOERROR(RET_NOT_IMPLEMENTED_YET);
}



returnValue AcadoIntegratorBackend::getProtectedX( Vector *xEnd ) const{
  printf("\ngetresults\n");

    int run1;

    if( (int) xEnd[0].getDim() != m )
        return RET_INPUT_HAS_WRONG_DIMENSION;

    for( run1 = 0; run1 < m; run1++ )
        xEnd[0](run1) = nablaY(0,run1);

    return SUCCESSFUL_RETURN;
}


returnValue AcadoIntegratorBackend::getProtectedForwardSensitivities( Matrix *Dx, int order ) const{
    printf("\ngetfwd\n");

  
    int run1;

    if( Dx == NULL ){
        return SUCCESSFUL_RETURN;
    }

    if( order == 1 && nFDirs2 == 0 ){

        if( (int) Dx[0].getNumCols() != nFDirs )
            return RET_INPUT_HAS_WRONG_DIMENSION;
        if( (int) Dx[0].getNumRows() != m )
            return RET_INPUT_HAS_WRONG_DIMENSION;

        for( run1 = 0; run1 < m; run1++ ){
            Dx[0](run1,0) = nablaG(0,run1);
        }
    }

    if( order == 2 ){

        if( (int) Dx[0].getNumCols() != nFDirs2 )
            return RET_INPUT_HAS_WRONG_DIMENSION;
        if( (int) Dx[0].getNumRows() != m )
            return RET_INPUT_HAS_WRONG_DIMENSION;

        for( run1 = 0; run1 < m; run1++ ){
            Dx[0](run1,0) = nablaG3(0,run1);
        }
    }

    if( order == 1 && nFDirs2 > 0 ){

        if( (int) Dx[0].getNumCols() != nFDirs2 )
            return RET_INPUT_HAS_WRONG_DIMENSION;
        if( (int) Dx[0].getNumRows() != m )
            return RET_INPUT_HAS_WRONG_DIMENSION;

        for( run1 = 0; run1 < m; run1++ ){
            Dx[0](run1,0) = nablaG2(0,run1);
        }
    }

    if( order < 1 || order > 2 ){
        return ACADOERROR(RET_INPUT_OUT_OF_RANGE);
    }

    return SUCCESSFUL_RETURN;
}


returnValue AcadoIntegratorBackend::getProtectedBackwardSensitivities(Vector &Dx_x0,
                                                             Vector &Dx_p ,
                                                             Vector &Dx_u ,
                                                             Vector &Dx_w ,
                                                             int order      ) const{
    printf("\ngetbwd\n");

                                                               
    if( order == 1 && nBDirs2 == 0 )
        copyBackward( Dx_x0, Dx_p, Dx_u, Dx_w, nablaH );

    if( order == 1 && nBDirs2 > 0 )
        copyBackward( Dx_x0, Dx_p, Dx_u, Dx_w, nablaH2 );

    if( order == 2 )
        copyBackward( Dx_x0, Dx_p, Dx_u, Dx_w, nablaH3 );

    if( order < 1 || order > 2 )
        return ACADOERROR(RET_INPUT_OUT_OF_RANGE);

    return SUCCESSFUL_RETURN;
}



returnValue AcadoIntegratorBackend::setDxInitialization( double *dx0 ){

    int run1;

    if( dx0 != NULL ){

        for( run1 = 0; run1 < md; run1++ ){
            initial_guess[run1] = dx0[run1];
        }
    }
    else{

        for( run1 = 0; run1 < md; run1++ ){
            initial_guess[run1] = 0.0;
        }
    }

    for( run1 = 0; run1 < md; run1++ ){
        initial_guess[md+run1] = 0.0;
    }

    return SUCCESSFUL_RETURN;
}


int AcadoIntegratorBackend::getNumberOfSteps() const{

    return count2-2;
}

int AcadoIntegratorBackend::getNumberOfRejectedSteps() const{

    return count3;
}




double AcadoIntegratorBackend::getStepSize() const{

    return h[0];
}



//
// PROTECTED MEMBER FUNCTIONS:
//


double AcadoIntegratorBackend::determinePredictor( int number_, BooleanType ini ){

    int run1;
    double pp;

    if( soa != SOA_EVERYTHING_FROZEN && soa != SOA_MESH_FROZEN ){

        if( soa != SOA_FREEZING_ALL ){

            psi_[0] = psi[number_][0];
            psi_[1] = psi[number_][1];
            psi_[2] = psi[number_][2];
            psi_[3] = psi[number_][3];

            psi[number_][3] = psi[number_][2] + h[0];
            psi[number_][2] = psi[number_][1] + h[0];
            psi[number_][1] = psi[number_][0] + h[0];
            psi[number_][0] = h[0];

            gamma[number_][4] =   1.0/psi[number_][0] + 1.0/psi[number_][1]
                                + 1.0/psi[number_][2] + 1.0/psi[number_][3];
        }
        else{

            psi[number_][3] = psi[number_-1][2] + h[0];
            psi[number_][2] = psi[number_-1][1] + h[0];
            psi[number_][1] = psi[number_-1][0] + h[0];
            psi[number_][0] = h[0];

            gamma[number_][4] =   1.0/psi[number_][0] + 1.0/psi[number_][1]
                                + 1.0/psi[number_][2] + 1.0/psi[number_][3];
        }
    }

    const double scalT = h[0]/rel_time_scale;

    pp = psi[number_][0]*scalT/h[0];
    for( run1 = 0; run1 < m; run1++ )
        phi(1,run1) = pp*nablaY(1,run1);

    pp *= psi[number_][1]*scalT/h[0];
    for( run1 = 0; run1 < m; run1++ )
        phi(2,run1) = pp*nablaY(2,run1);

    pp *= psi[number_][2]*scalT/h[0];
    for( run1 = 0; run1 < m; run1++ )
        phi(3,run1) = pp*nablaY(3,run1);

    pp *= psi[number_][3]*scalT/h[0];
    for( run1 = 0; run1 < m; run1++ )
        phi(4,run1) = pp*nablaY(4,run1);


    for( run1 = 0; run1 < m; run1++ )
        delta(1,run1) = nablaY(0,run1) + phi(1,run1);

    for( run1 = 0; run1 < m; run1++ )
        delta(2,run1) = delta(1,run1) + phi(2,run1);

    for( run1 = 0; run1 < m; run1++ )
        delta(3,run1) = delta(2,run1) + phi(3,run1);

    for( run1 = 0; run1 < m; run1++ )
        eta[0][run1] = delta(3,run1) + phi(4,run1);


    for( run1 = 0; run1 < md; run1++ ){

        c2(run1) = - nablaY(0,run1)/psi[number_][0]
                   - delta (1,run1)/psi[number_][1]
                   - delta (2,run1)/psi[number_][2]
                   - delta (3,run1)/psi[number_][3];
    }

    returnValue returnvalue;

    correctorTime.start();

    returnvalue = determineCorrector(number_, ini);

    correctorTime.stop();

    if( returnvalue == RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF ){
        return -1;
    }
    if( returnvalue != SUCCESSFUL_RETURN ){
        return INFTY;
    }

    for( run1 = 0; run1 < m; run1++ ){

        nablaY_(0,run1) = eta[nOfNewtonSteps[number_]][run1];

        nablaY_(1,run1) = nablaY_(0,run1) - nablaY(0,run1);
        nablaY_(2,run1) = (nablaY_(1,run1) - nablaY(1,run1)*scalT)*(h[0]/psi[number_][1]);
        nablaY_(3,run1) = (nablaY_(2,run1) - nablaY(2,run1)*scalT*scalT)*(h[0]/psi[number_][2]);
        nablaY_(4,run1) = (nablaY_(3,run1) - nablaY(3,run1)*scalT*scalT*scalT)*(h[0]/psi[number_][3]);
    }

    double EE = EPS;
    for( run1 = 0; run1 < md; run1++ ){
        if( (eta[0][run1]-nablaY_(0,run1))/diff_scale(run1) >= EE  ){
            EE = (eta[0][run1]-nablaY_(0,run1))/diff_scale(run1);
        }
        if( (eta[0][run1]-nablaY_(0,run1))/diff_scale(run1) <= -EE ){
            EE = -(eta[0][run1]-nablaY_(0,run1))/diff_scale(run1);
        }
    }

    return 0.05*EE;
}


returnValue AcadoIntegratorBackend::determineCorrector( int stepnumber, BooleanType ini ){

    int run1, run2;
    int newtonsteps;

    double norm1 = 0.0;
    double norm2 = 0.0;

    BooleanType COMPUTE_JACOBIAN  = BT_FALSE;
    BooleanType JACOBIAN_COMPUTED = BT_FALSE;

    if( soa != SOA_MESH_FROZEN && soa != SOA_EVERYTHING_FROZEN ){
        COMPUTE_JACOBIAN = ini;
    }

    if( soa != SOA_MESH_FROZEN && soa != SOA_EVERYTHING_FROZEN ){
       nOfNewtonSteps[stepnumber] = 0;
    }

    if( stepnumber > 0 ){
        M_index[stepnumber] = M_index[stepnumber-1];
    }

    newtonsteps = 0;

    while( newtonsteps < 3 ){

       // evaluate F:
       // -----------
       x[time_index] = t + h[0];
       for( run2 = 0; run2 < m; run2++ ){
            x[diff_index[run2]] = eta[newtonsteps][run2];
       }
       for( run2 = 0; run2 < md; run2++ ){
           x[ddiff_index[run2]] = gamma[stepnumber][4]*eta[newtonsteps][run2]+c2(run2);
       }

       functionEvaluation.start();

       if( rhs[0].evaluate( 3*stepnumber+newtonsteps, x, F ) != SUCCESSFUL_RETURN ){
           return ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
       }

       relaxAlgebraic( F, x[time_index] );

       nFcnEvaluations++;
       functionEvaluation.stop();

       if( COMPUTE_JACOBIAN == BT_TRUE ){

           jacComputation.start();

           if( PrintLevel == HIGH ){
               acadoPrintf("(RE-) COMPUTE-JACOBIAN... \n");
           }

           if( soa == SOA_FREEZING_MESH || soa == SOA_FREEZING_ALL ){

               if( nOfM >= maxNM ){
                   int oldMN = maxNM;
                   maxNM += maxNM;
                   M = (Matrix**)realloc(M, maxNM*sizeof(Matrix*));
                   for( run1 = oldMN; run1 < maxNM; run1++ )
                       M[run1] = 0;
               }
               M_index[stepnumber] = nOfM;
               M[nOfM] = new Matrix(m,m);
               nOfM++;
           }
           else{
               if( M[0] == 0 ) M[0] = new Matrix(m,m);
               M_index[stepnumber] = 0;
               M[0]->init(m,m);
           }

           for( run1 = 0; run1 < md; run1++ ){
               iseed[ddiff_index[run1]] = gamma[stepnumber][4];
               iseed[ diff_index[run1]] = 1.0;
               if( rhs[0].AD_forward( 3*stepnumber+newtonsteps, iseed,
                                      k2[0][0] ) != SUCCESSFUL_RETURN ){
                  return ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
               }
               for( run2 = 0; run2 < m; run2++ )
                   M[M_index[stepnumber]]->operator()(run2,run1) = k2[0][0][run2];

               iseed[ddiff_index[run1]] = 0.0;
               iseed[ diff_index[run1]] = 0.0;
           }

           for( run1 = 0; run1 < ma; run1++ ){
               iseed[ diff_index[md+run1]] = 1.0;
               if( rhs[0].AD_forward( 3*stepnumber+newtonsteps, iseed,
                                        k2[0][0] ) != SUCCESSFUL_RETURN ){
                  return ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
               }
               for( run2 = 0; run2 < m; run2++ )
                   M[M_index[stepnumber]]->operator()(run2,md+run1) = k2[0][0][run2];

               iseed[diff_index[md+run1]] = 0.0;
           }

           nJacEvaluations++;
           jacComputation.stop();
           jacDecomposition.start();

           if( decomposeJacobian( *M[M_index[stepnumber]] ) != SUCCESSFUL_RETURN )
               return ACADOERROR(RET_THE_DAE_INDEX_IS_TOO_LARGE);

           jacDecomposition.stop();

           JACOBIAN_COMPUTED = BT_TRUE;
           COMPUTE_JACOBIAN  = BT_FALSE;
       }

       norm1 = applyNewtonStep( eta[newtonsteps+1],
                                eta[newtonsteps],
                               *M[M_index[stepnumber]],
                                F );

       if( soa == SOA_MESH_FROZEN || soa == SOA_EVERYTHING_FROZEN ){
           if( newtonsteps == nOfNewtonSteps[stepnumber] ){
               return SUCCESSFUL_RETURN;
           }
       }


       double subTOL = 1e-3*TOL;

       if( norm1 < subTOL ){
           nOfNewtonSteps[stepnumber] = newtonsteps+1;
           return SUCCESSFUL_RETURN;
       }

       if( newtonsteps == 0 ){
           norm2 = norm1;
       }

       if( newtonsteps == 1 && (norm1/norm2 > 0.33 || norm1/norm2 > sqrt(0.33*TOL/norm2)) ){

           if( JACOBIAN_COMPUTED == BT_FALSE ){
               COMPUTE_JACOBIAN = BT_TRUE;
               newtonsteps      = -1     ;
           }
       }

       if( newtonsteps == 1 ){
           norm2 = norm1;
       }

       if( newtonsteps == 2 ){

           if( norm1 > 0.33*TOL ){

               if( JACOBIAN_COMPUTED == BT_FALSE ){
                   COMPUTE_JACOBIAN = BT_TRUE;
                   newtonsteps      = -1     ;
               }
               else{
                   return RET_INFO_UNDEFINED;
               }
           }
           else{
               nOfNewtonSteps[stepnumber] = newtonsteps+1;
               return SUCCESSFUL_RETURN;
           }
       }

       newtonsteps++;
       nOfNewtonSteps[stepnumber] = newtonsteps;
    }

    return RET_INFO_UNDEFINED;
}


returnValue AcadoIntegratorBackend::rk_start(){

    int run1, run2, run3;
    double E = EPS;

     if( soa == SOA_EVERYTHING_FROZEN || soa == SOA_MESH_FROZEN ){
         rel_time_scale = h[0];
         h[0] = h[1];
     }

     // MEMORY REALLOCATION:
     // --------------------
     if( soa != SOA_EVERYTHING_FROZEN && soa != SOA_MESH_FROZEN ){

         int oldAlloc = maxAlloc;

         if( maxAlloc < 8 ){
             maxAlloc = 8;

             psi            = (double**)realloc(psi,maxAlloc*sizeof(double*));
             gamma          = (double**)realloc(gamma,maxAlloc*sizeof(double*));
             nOfNewtonSteps = (int*)realloc(nOfNewtonSteps,(maxAlloc)*sizeof(int));

             for( run1 = oldAlloc; run1 < maxAlloc; run1++ ){
                 nOfNewtonSteps[run1] = 0;
             }
             for( run1 = oldAlloc; run1 < maxAlloc; run1++ ){

                 psi  [run1] = (double*)calloc(nstep,sizeof(double));
                 gamma[run1] = (double*)calloc(nstep,sizeof(double));
             }
         }

         M_index      = (int*)realloc(M_index,maxAlloc*sizeof(int));
         h = (double*)realloc(h,maxAlloc*sizeof(double));
     }


    if( soa != SOA_EVERYTHING_FROZEN ){

    returnValue returnvalue;
    int step_rejections = 0;

    while( step_rejections < 1 ){

    // determine k:
    // ------------
    run1 = 0;

    for( run2 = 0; run2 < m; run2++ ){
        k[0][run1][run2] = initial_guess[run2];
    }

    while( run1 < dim ){

        if( run1 > 0 ){
            for( run2 = 0; run2 < m; run2++ ){
                k[0][run1][run2] = k[nOfNewtonSteps[run1-1]][run1-1][run2];
            }
        }

        returnvalue = rk_start_solve(run1);

        if( returnvalue != SUCCESSFUL_RETURN ){

            if( PrintLevel == HIGH ){
                acadoPrintf("RUNGE-KUTTA STARTER: \n");
                acadoPrintf("STEP REJECTED: error estimate           = %.16e \n", E        );
                acadoPrintf("               required local tolerance = %.16e \n", TOL      );
                count3++;
            }

            if( soa == SOA_EVERYTHING_FROZEN || soa == SOA_MESH_FROZEN ){
                return ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }

            if( h[0] <= hmin + EPS ){
                return ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }
            h[0] = 0.2*h[0];
            if( h[0] < hmin ){
                h[0] = hmin;
            }
            run1 = -1;
        }

        run1++;
    }


    // save previous eta4:
    // ----------------------------------------------

       for( run1 = 0; run1 < m; run1++ ){
           eta4_[run1]  = eta4[run1];
           eta5 [run1]  = eta4[run1];
       }

    // determine eta4 and eta5:
    // ----------------------------------------------

       for( run1 = 0; run1 < dim; run1++ ){
           for( run2 = 0; run2 < md; run2++ ){
               eta4[run2] = eta4[run2] + b4[run1]*h[0]*k[nOfNewtonSteps[run1]][run1][run2];
               eta5[run2] = eta5[run2] + b5[run1]*h[0]*k[nOfNewtonSteps[run1]][run1][run2];
           }
       }
       for( run1 = 0; run1 < ma; run1++ ){
           eta4[md+run1] = k[nOfNewtonSteps[dim-1]][dim-1][md+run1];
       }

    // determine local error estimate E:
    // ----------------------------------------------

       E = EPS;
       for( run2 = 0; run2 < md; run2++ ){
           if( (eta4[run2]-eta5[run2])/diff_scale(run2) >= E  ){
               E = (eta4[run2]-eta5[run2])/diff_scale(run2);
           }
           if( (eta4[run2]-eta5[run2])/diff_scale(run2) <= -E ){
               E = (-eta4[run2]+eta5[run2])/diff_scale(run2);
           }
       }

       if( E < TOL*h[0] || soa == SOA_EVERYTHING_FROZEN || soa == SOA_MESH_FROZEN  ){
          break;
       }
       else{

            if( PrintLevel == HIGH ){
                acadoPrintf("RUNGE-KUTTA STARTER: \n");
                acadoPrintf("STEP REJECTED: error estimate           = %.16e \n", E        );
                acadoPrintf("               required local tolerance = %.16e \n", TOL*h[0] );
                count3++;
            }

            for( run1 = 0; run1 < m; run1++ ){
                eta4[run1] = eta4_[run1];
            }
            if( h[0] <= hmin + EPS ){
                return ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }
            h[0] = 0.2*h[0];
            if( h[0] < hmin ){
                h[0] = hmin;
            }
       }

       step_rejections++;
    }

    }


    // PROCEED IF THE STEP IS ACCEPTED:
    // --------------------------------

    // plan the first step of the BDF method:
    // --------------------------------------

    if( soa != SOA_EVERYTHING_FROZEN ){

    for( run1 = 0; run1 < nstep-1; run1++ ){
        for( run2 = 0; run2 < md; run2++ ){
            nablaY(nstep-2-run1,run2) = eta4_[run2];
            for( run3 = 0; run3 <= run1+3; run3++ ){
                nablaY(nstep-2-run1,run2) = nablaY(nstep-2-run1,run2) +
                  h[0]*A[run1+3][run3]*k[nOfNewtonSteps[run3]][run3][run2];
            }
        }
        for( run2 = 0; run2 < ma; run2++ ){
            nablaY(nstep-2-run1,md+run2) = k[nOfNewtonSteps[run1+3]][run1+3][md+run2];
        }
    }

    }


     // compute forward derivatives if requested:
     // ------------------------------------------

     if( nFDirs > 0 && nBDirs2 == 0 && nFDirs2 == 0 ){

         if( nBDirs != 0 ){
             return ACADOERROR(RET_WRONG_DEFINITION_OF_SEEDS);
         }

         determineRKEtaGForward();
     }

     if( nBDirs > 0 ){

         if( soa != SOA_EVERYTHING_FROZEN ){
             return ACADOERROR(RET_NOT_FROZEN);
         }
         if( nFDirs != 0 || nBDirs2 != 0 || nFDirs2 != 0 ){
             return ACADOERROR(RET_WRONG_DEFINITION_OF_SEEDS);
         }
         determineRKEtaHBackward();
     }
     if( nFDirs2 > 0 ){

         if( soa != SOA_EVERYTHING_FROZEN ){
             return ACADOERROR(RET_NOT_FROZEN);
         }
         if( nBDirs != 0 || nBDirs2 != 0 || nFDirs != 1 ){
             return ACADOERROR(RET_WRONG_DEFINITION_OF_SEEDS);
         }
         determineRKEtaGForward2();
     }
     if( nBDirs2 > 0 ){

         if( soa != SOA_EVERYTHING_FROZEN ){
             return ACADOERROR(RET_NOT_FROZEN);
         }
         if( nBDirs != 0 || nFDirs2 != 0 || nFDirs != 1 ){
             return ACADOERROR(RET_WRONG_DEFINITION_OF_SEEDS);
         }

         determineRKEtaHBackward2();
     }

    // Printing:
    // ---------
    printRKIntermediateResults();

    const double hhh = c[3]*h[0];


    // STORAGE:
    // --------
    if( soa == SOA_FREEZING_MESH || soa == SOA_FREEZING_ALL ){

        h[1] = h[0];
        for( run3 = 3; run3 >= 0; run3-- )
            h[2+run3] = hhh;
    }


    int i1 = timeInterval.getFloorIndex( t             );
    int i2 = timeInterval.getFloorIndex( t + c[6]*h[0] );
    int jj;

    for( jj = i1+1; jj <= i2; jj++ ){

        for( run1 = 0; run1 < m; run1++ ){
            if( nFDirs == 0 && nBDirs  == 0 && nFDirs2 == 0 && nBDirs == 0 )   xStore( jj, run1 ) = nablaY (3,run1);
            if( nFDirs  > 0 && nBDirs2 == 0 && nFDirs2 == 0                )  dxStore( jj, run1 ) = nablaG (3,run1);
            if( nFDirs2 > 0                                                ) ddxStore( jj, run1 ) = nablaG3(3,run1);
        }
        for( run1 = 0; run1 < mn; run1++ )
            iStore( jj, run1 ) = x[rhs->index( VT_INTERMEDIATE_STATE, run1 )];
    }


    // PREPARE MODIFIED DIVIDED DIFFERENCES:
    // -------------------------------------
    if( soa != SOA_EVERYTHING_FROZEN )
        prepareDividedDifferences( nablaY );

    if( nFDirs > 0 && nBDirs2 == 0 && nFDirs2 == 0 )
        prepareDividedDifferences( nablaG );

    if( nFDirs2 > 0 ){
        prepareDividedDifferences( nablaG2 );
        prepareDividedDifferences( nablaG3 );
    }


    if( soa != SOA_EVERYTHING_FROZEN ){

        if( soa != SOA_FREEZING_ALL ){

            psi[7][0] = hhh    ;
            psi[7][1] = 2.0*hhh;
            psi[7][2] = 3.0*hhh;
            psi[7][3] = 4.0*hhh;
        }
        else{

            psi[6][0] = hhh    ;
            psi[6][1] = 2.0*hhh;
            psi[6][2] = 3.0*hhh;
            psi[6][3] = 4.0*hhh;
        }
    }

    // increase the time:
    // ----------------------------------------------

    t = t + c[6]*h[0];

    if( soa != SOA_EVERYTHING_FROZEN && soa != SOA_MESH_FROZEN ){

    // determine the new step size:
    // ----------------------------------------------

        rel_time_scale = 0.2*h[0];
        h[0] = 0.2*h[0]*pow( tune*(TOL*h[0]/E) , 1.0/3.0 );

        if( h[0] > hmax ){
            h[0] = hmax;
        }
        if( h[0] < hmin ){
            h[0] = hmin;
        }

        if( t + h[0] >= timeInterval.getLastTime() ){
            h[0] = timeInterval.getLastTime()-t;
        }
    }
    else{
        h[0] = 0.2*h[0];
    }

    return SUCCESSFUL_RETURN;
}


returnValue AcadoIntegratorBackend::rk_start_solve( int stepnumber ){

    int run1, run2, run3;
    int newtonsteps;

    double norm1 = 0.0;
    double norm2 = 0.0;

    BooleanType COMPUTE_JACOBIAN;
    BooleanType JACOBIAN_COMPUTED = BT_FALSE;

    if( stepnumber == 0 && soa != SOA_MESH_FROZEN && soa != SOA_EVERYTHING_FROZEN ){
        COMPUTE_JACOBIAN = BT_TRUE;
    }
    else{
        COMPUTE_JACOBIAN = BT_FALSE;
    }

    if( soa != SOA_MESH_FROZEN && soa != SOA_EVERYTHING_FROZEN ){
       nOfNewtonSteps[stepnumber] = 0;
    }

    if( stepnumber > 0 ){
        M_index[stepnumber] = M_index[stepnumber-1];
    }

    newtonsteps = 0;

    while( newtonsteps < 3 ){

       // evaluate F:
       // -----------
       x[time_index] = t + c[stepnumber]*h[0];
       for( run2 = 0; run2 < md; run2++ ){
            x[diff_index[run2]] = eta4[run2];
            for( run3 = 0; run3 < stepnumber+1; run3++ ){
                 x[diff_index[run2]] = x[diff_index[run2]] +
                          A[stepnumber][run3]*h[0]*k[nOfNewtonSteps[run3]][run3][run2];
            }
       }
       for( run2 = 0; run2 < ma; run2++ ){
           x[diff_index[md+run2]] = k[newtonsteps][stepnumber][md+run2];
       }
       for( run2 = 0; run2 < md; run2++ ){
           x[ddiff_index[run2]] = k[newtonsteps][stepnumber][run2];
       }

       functionEvaluation.start();

       if( rhs[0].evaluate( 3*stepnumber+newtonsteps, x, F ) != SUCCESSFUL_RETURN ){
           return ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
       }

       relaxAlgebraic( F, x[time_index] );

       nFcnEvaluations++;
       functionEvaluation.stop();


       if( COMPUTE_JACOBIAN == BT_TRUE ){

           jacComputation.start();

           if( PrintLevel == HIGH ){
               acadoPrintf("(RE-) COMPUTE-JACOBIAN... \n");
           }

           const double ise = h[0]*A[stepnumber][stepnumber];

           if( soa == SOA_FREEZING_MESH || soa == SOA_FREEZING_ALL ){

               if( nOfM >= maxNM ){
                   int oldMN = maxNM;
                   maxNM += maxNM;
                   M = (Matrix**)realloc(M, maxNM*sizeof(Matrix*));
                   for( run1 = oldMN; run1 < maxNM; run1++ )
                       M[run1] = 0;
               }
               M_index[stepnumber] = nOfM;
               M[nOfM] = new Matrix(m,m);
               nOfM++;
           }
           else{
               if( M[0] == 0 ) M[0] = new Matrix(m,m);
               M_index[stepnumber] = 0;
               M[0]->init(m,m);
           }

           for( run1 = 0; run1 < md; run1++ ){
               iseed[ddiff_index[run1]] = 1.0;
               iseed[ diff_index[run1]] = ise;
               if( rhs[0].AD_forward( 3*stepnumber+newtonsteps, iseed,
                                      k2[0][0] ) != SUCCESSFUL_RETURN ){
                  return ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
               }

               for( run2 = 0; run2 < m; run2++ )
                   M[M_index[stepnumber]]->operator()(run2,run1) = k2[0][0][run2];

               iseed[ddiff_index[run1]] = 0.0;
               iseed[ diff_index[run1]] = 0.0;
           }
           for( run1 = 0; run1 < ma; run1++ ){
               iseed[ diff_index[md+run1]] = 1.0;
               if( rhs[0].AD_forward( 3*stepnumber+newtonsteps, iseed,
                                      k2[0][0] ) != SUCCESSFUL_RETURN ){
                  return ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
               }

               for( run2 = 0; run2 < m; run2++ )
                   M[M_index[stepnumber]]->operator()(run2,md+run1) = k2[0][0][run2];

               iseed[diff_index[md+run1]] = 0.0;
           }

           nJacEvaluations++;
           jacComputation.stop();
           jacDecomposition.start();

           if( decomposeJacobian( *M[M_index[stepnumber]] ) != SUCCESSFUL_RETURN )
                   return ACADOERROR(RET_THE_DAE_INDEX_IS_TOO_LARGE);

           jacDecomposition.stop();

           JACOBIAN_COMPUTED = BT_TRUE;
           COMPUTE_JACOBIAN  = BT_FALSE;
       }

       norm1 = applyNewtonStep( k[newtonsteps+1][stepnumber],
                                k[newtonsteps][stepnumber]  ,
                               *M[M_index[stepnumber]],
                                F );

       if( soa == SOA_MESH_FROZEN || soa == SOA_EVERYTHING_FROZEN ){
           if( newtonsteps == nOfNewtonSteps[stepnumber] ){
               return SUCCESSFUL_RETURN;
           }
       }

       double subTOL = 1e-3*TOL;

       if( norm1 < subTOL ){

           nOfNewtonSteps[stepnumber] = newtonsteps+1;
           return SUCCESSFUL_RETURN;
       }

       if( newtonsteps == 0 ){
           norm2 = norm1;
       }

       if( newtonsteps == 1 && (norm1/norm2 > 0.33 || norm1/norm2 > sqrt(0.33*TOL/norm2)) ){

           if( JACOBIAN_COMPUTED == BT_FALSE ){
               COMPUTE_JACOBIAN = BT_TRUE;
               newtonsteps      = -1     ;
           }
       }

       if( newtonsteps == 1 ){
           norm2 = norm1;
       }

       if( newtonsteps == 2 ){

           if( JACOBIAN_COMPUTED == BT_FALSE ){
               COMPUTE_JACOBIAN = BT_TRUE;
               newtonsteps      = -1     ;
           }
           else{

               if( norm1 > 0.33*TOL ){
                   // printf("RK - REJECTION !!!!!!  norm1 = %.16e \n", norm1 );
                   return RET_INFO_UNDEFINED;
               }

               nOfNewtonSteps[stepnumber] = newtonsteps+1;
               return SUCCESSFUL_RETURN;
           }
       }

       newtonsteps++;
       nOfNewtonSteps[stepnumber] = newtonsteps;
    }

    return RET_INFO_UNDEFINED;
}


void AcadoIntegratorBackend::determineRKEtaGForward(){

    int run1, run2, run3, newtonsteps;


    for( run2 = 0; run2 < m; run2++ ){
        k[0][0][run2] = 0.0;
    }


    // determine k:
    // -----------------------------------------------
       for( run1 = 0; run1 < dim; run1++ ){

           if( run1 > 0 ){
               for( run2 = 0; run2 < m; run2++ ){
                   k[0][run1][run2] = k[nOfNewtonSteps[run1-1]][run1-1][run2];
               }
           }

           newtonsteps = 0;
           while( newtonsteps < nOfNewtonSteps[run1] ){

               for( run2 = 0; run2 < md; run2++ ){
                   G[diff_index[run2]] = etaG[run2];
                   for( run3 = 0; run3 < run1; run3++ ){
                       G[diff_index[run2]] = G[diff_index[run2]] +
                                   A[run1][run3]*h[0]*k[nOfNewtonSteps[run3]][run3][run2];
                   }
                   G[diff_index[run2]] = G[diff_index[run2]] +
                               A[run1][run1]*h[0]*k[newtonsteps][run1][run2];
               }
               for( run2 = 0; run2 < ma; run2++ ){
                   G[diff_index[md+run2]] = k[newtonsteps][run1][md+run2];
               }
               for( run2 = 0; run2 < md; run2++ ){
                   G[ddiff_index[run2]] = k[newtonsteps][run1][run2];
               }

               if( rhs[0].AD_forward( 3*run1+newtonsteps, G, F )
                                      != SUCCESSFUL_RETURN ){
                   ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
                   return;
               }


               applyNewtonStep( k[newtonsteps+1][run1],
                                k[newtonsteps][run1]  ,
                               *M[M_index[run1]],
                                F );

               newtonsteps++;
           }
       }

    // determine nablaG:
    // ----------------------------------------------

       for( run1 = 0; run1 < nstep-1; run1++ ){
           for( run2 = 0; run2 < md; run2++ ){
               nablaG(nstep-2-run1,run2) = etaG[run2];
               for( run3 = 0; run3 <= run1+3; run3++ ){
                   nablaG(nstep-2-run1,run2) = nablaG(nstep-2-run1,run2) +
                     h[0]*A[run1+3][run3]*k[nOfNewtonSteps[run3]][run3][run2];
               }
           }
           for( run2 = 0; run2 < ma; run2++ ){
               nablaG(nstep-2-run1,md+run2) = k[nOfNewtonSteps[run1+3]][run1+3][md+run2];
           }
        }
}


void AcadoIntegratorBackend::determineRKEtaHBackward(){

    int run1, run2, run3, newtonsteps;
    const double hhh = c[3]*h[0];
    int number_;
    double *Hstore = new double[ndir];

    const double scalT = rel_time_scale/hhh;

        for( run3 = 0; run3 < 4; run3++ ){
            for( run1 = 0; run1 < dim; run1++ ){
                for( run2 = 0; run2 < ndir; run2++ ){
                    l[run3][run1][run2] = 0.0;
                }
            }
        }

        for( run1 = 0; run1 < ndir; run1++ ){
            Hstore[run1] = nablaH(0,run1);
        }

        for( run1 = 0; run1 < ndir; run1++ ){

            zH[5][run1] = nablaH(3,run1)*scalT*scalT;
            zH[4][run1] = nablaH(2,run1)*scalT + zH[5][run1]/(3.0);
            zH[3][run1] = - zH[5][run1]/(3.0);
            zH[2][run1] = nablaH(1,run1) + zH[4][run1]/(2.0);
            zH[1][run1] = (zH[3][run1]-zH[4][run1])/(2.0);
            zH[0][run1] = - zH[3][run1]/(2.0);

            nablaH(3,run1) = nablaH(0,run1)*h[0] + zH[2][run1]*(h[0]/hhh);
            nablaH(2,run1) = (zH[1][run1]-zH[2][run1])*(h[0]/hhh);
            nablaH(1,run1) = (zH[0][run1]-zH[1][run1])*(h[0]/hhh);
            nablaH(0,run1) = -zH[0][run1]*(h[0]/hhh);
        }

        for( run1 = 0; run1 < ndir; run1++ ){
            kH[6][nOfNewtonSteps[6]][run1] = nablaH(3,run1)/h[0];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH[6][nOfNewtonSteps[6]][diff_index[run1]] =
            kH[6][nOfNewtonSteps[6]][diff_index[run1]]*h[0]*A[6][6];
        }

        newtonsteps = nOfNewtonSteps[6];
        newtonsteps--;

        number_ = 6;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H );

            if( rhs[0].AD_backward( 3*number_+newtonsteps, H, l[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH[number_][newtonsteps][run1] =   kH[number_][newtonsteps+1][run1]
                                                - l[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH[number_][newtonsteps][diff_index[run1]] = 
                   kH[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }


        for( run1 = 0; run1 < ndir; run1++ ){
            kH[5][nOfNewtonSteps[5]][run1] = kH[6][0][run1] + nablaH(2,run1)/h[0];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH[5][nOfNewtonSteps[5]][diff_index[run1]] = kH[6][0][diff_index[run1]]
            +  nablaH(3,diff_index[run1])*A[6][5]
            +  nablaH(2,diff_index[run1])*A[5][5];
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH[5][nOfNewtonSteps[5]][diff_index[run1]] =
               kH[5][nOfNewtonSteps[5]][diff_index[run1]]
               -l[run3][number_][diff_index[run1]]*h[0]*A[6][5];
            }
        }


        newtonsteps = nOfNewtonSteps[5];
        newtonsteps--;

        number_ = 5;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H );

            if( rhs[0].AD_backward( 3*number_+newtonsteps, H, l[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH[number_][newtonsteps][run1] =   kH[number_][newtonsteps+1][run1]
                                                - l[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH[number_][newtonsteps][diff_index[run1]] = 
                   kH[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }


        for( run1 = 0; run1 < ndir; run1++ ){
            kH[4][nOfNewtonSteps[4]][run1] = kH[5][0][run1] + nablaH(1,run1)/h[0];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH[4][nOfNewtonSteps[4]][diff_index[run1]] = kH[5][0][diff_index[run1]]
            +  nablaH(3,diff_index[run1])*A[6][4]
            +  nablaH(2,diff_index[run1])*A[5][4]
            +  nablaH(1,diff_index[run1])*A[4][4];
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH[4][nOfNewtonSteps[4]][diff_index[run1]] =
               kH[4][nOfNewtonSteps[4]][diff_index[run1]]
               -l[run3][5][diff_index[run1]]*h[0]*A[5][4];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH[4][nOfNewtonSteps[4]][diff_index[run1]] =
               kH[4][nOfNewtonSteps[4]][diff_index[run1]]
               -l[run3][6][diff_index[run1]]*h[0]*A[6][4];
            }
        }


        newtonsteps = nOfNewtonSteps[4];
        newtonsteps--;

        number_ = 4;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H );

            if( rhs[0].AD_backward( 3*number_+newtonsteps, H, l[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH[number_][newtonsteps][run1] =   kH[number_][newtonsteps+1][run1]
                                                - l[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH[number_][newtonsteps][diff_index[run1]] = 
                   kH[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }


        for( run1 = 0; run1 < ndir; run1++ ){
            kH[3][nOfNewtonSteps[3]][run1] = kH[4][0][run1] + nablaH(0,run1)/h[0];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH[3][nOfNewtonSteps[3]][diff_index[run1]] = kH[4][0][diff_index[run1]]
            +  nablaH(3,diff_index[run1])*A[6][3]
            +  nablaH(2,diff_index[run1])*A[5][3]
            +  nablaH(1,diff_index[run1])*A[4][3]
            +  nablaH(0,diff_index[run1])*A[3][3];
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++){
               kH[3][nOfNewtonSteps[3]][diff_index[run1]] =
               kH[3][nOfNewtonSteps[3]][diff_index[run1]]
               -l[run3][4][diff_index[run1]]*h[0]*A[4][3];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH[3][nOfNewtonSteps[3]][diff_index[run1]] =
               kH[3][nOfNewtonSteps[3]][diff_index[run1]]
               -l[run3][5][diff_index[run1]]*h[0]*A[5][3];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH[3][nOfNewtonSteps[3]][diff_index[run1]] =
               kH[3][nOfNewtonSteps[3]][diff_index[run1]]
               -l[run3][6][diff_index[run1]]*h[0]*A[6][3];
            }
        }


        newtonsteps = nOfNewtonSteps[3];
        newtonsteps--;

        number_ = 3;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H );

            if( rhs[0].AD_backward( 3*number_+newtonsteps, H, l[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH[number_][newtonsteps][run1] =   kH[number_][newtonsteps+1][run1]
                                                - l[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH[number_][newtonsteps][diff_index[run1]] = 
                   kH[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }


        for( run1 = 0; run1 < ndir; run1++ ){
            kH[2][nOfNewtonSteps[2]][run1] = kH[3][0][run1];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH[2][nOfNewtonSteps[2]][diff_index[run1]] = kH[3][0][diff_index[run1]]
            +  nablaH(3,diff_index[run1])*A[6][2]
            +  nablaH(2,diff_index[run1])*A[5][2]
            +  nablaH(1,diff_index[run1])*A[4][2]
            +  nablaH(0,diff_index[run1])*A[3][2];
            for( run3 = 0; run3 < nOfNewtonSteps[3]; run3++){
               kH[2][nOfNewtonSteps[2]][diff_index[run1]] =
               kH[2][nOfNewtonSteps[2]][diff_index[run1]]
               -l[run3][3][diff_index[run1]]*h[0]*A[3][2];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++){
               kH[2][nOfNewtonSteps[2]][diff_index[run1]] =
               kH[2][nOfNewtonSteps[2]][diff_index[run1]]
               -l[run3][4][diff_index[run1]]*h[0]*A[4][2];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH[2][nOfNewtonSteps[2]][diff_index[run1]] =
               kH[2][nOfNewtonSteps[2]][diff_index[run1]]
               -l[run3][5][diff_index[run1]]*h[0]*A[5][2];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH[2][nOfNewtonSteps[2]][diff_index[run1]] =
               kH[2][nOfNewtonSteps[2]][diff_index[run1]]
               -l[run3][6][diff_index[run1]]*h[0]*A[6][2];
            }
        }


        newtonsteps = nOfNewtonSteps[2];
        newtonsteps--;

        number_ = 2;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H );

            if( rhs[0].AD_backward( 3*number_+newtonsteps, H, l[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH[number_][newtonsteps][run1] =   kH[number_][newtonsteps+1][run1]
                                                - l[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH[number_][newtonsteps][diff_index[run1]] = 
                   kH[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }


        for( run1 = 0; run1 < ndir; run1++ ){
            kH[1][nOfNewtonSteps[1]][run1] = kH[2][0][run1];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH[1][nOfNewtonSteps[1]][diff_index[run1]] = kH[2][0][diff_index[run1]]
            +  nablaH(3,diff_index[run1])*A[6][1]
            +  nablaH(2,diff_index[run1])*A[5][1]
            +  nablaH(1,diff_index[run1])*A[4][1]
            +  nablaH(0,diff_index[run1])*A[3][1];
            for( run3 = 0; run3 < nOfNewtonSteps[2]; run3++){
               kH[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l[run3][2][diff_index[run1]]*h[0]*A[2][1];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[3]; run3++){
               kH[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l[run3][3][diff_index[run1]]*h[0]*A[3][1];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++){
               kH[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l[run3][4][diff_index[run1]]*h[0]*A[4][1];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l[run3][5][diff_index[run1]]*h[0]*A[5][1];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l[run3][6][diff_index[run1]]*h[0]*A[6][1];
            }
        }


        newtonsteps = nOfNewtonSteps[1];
        newtonsteps--;

        number_ = 1;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H );

            if( rhs[0].AD_backward( 3*number_+newtonsteps, H, l[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH[number_][newtonsteps][run1] =   kH[number_][newtonsteps+1][run1]
                                                - l[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH[number_][newtonsteps][diff_index[run1]] = 
                   kH[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }



        for( run1 = 0; run1 < ndir; run1++ ){
            kH[0][nOfNewtonSteps[0]][run1] = kH[1][0][run1];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH[0][nOfNewtonSteps[0]][diff_index[run1]] = kH[1][0][diff_index[run1]]
            +  nablaH(3,diff_index[run1])*A[6][0]
            +  nablaH(2,diff_index[run1])*A[5][0]
            +  nablaH(1,diff_index[run1])*A[4][0]
            +  nablaH(0,diff_index[run1])*A[3][0];
            for( run3 = 0; run3 < nOfNewtonSteps[1]; run3++){
               kH[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l[run3][1][diff_index[run1]]*h[0]*A[1][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[2]; run3++){
               kH[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l[run3][2][diff_index[run1]]*h[0]*A[2][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[3]; run3++){
               kH[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l[run3][3][diff_index[run1]]*h[0]*A[3][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++){
               kH[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l[run3][4][diff_index[run1]]*h[0]*A[4][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l[run3][5][diff_index[run1]]*h[0]*A[5][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l[run3][6][diff_index[run1]]*h[0]*A[6][0];
            }
        }

        newtonsteps = nOfNewtonSteps[0];
        newtonsteps--;

        number_ = 0;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H );

            if( rhs[0].AD_backward( 3*number_+newtonsteps, H, l[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }

            for( run1 = 0; run1 < ndir; run1++ ){
                kH[number_][newtonsteps][run1] =   kH[number_][newtonsteps+1][run1]
                                                - l[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH[number_][newtonsteps][diff_index[run1]] = 
                   kH[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }

        for( run1 = 0; run1 < ndir; run1++ ){
            nablaH(0,run1) = Hstore[run1];
        }
        for( run1 = 0; run1 < ndir; run1++ ){
            for( run3 = 0; run3 < nOfNewtonSteps[0]; run3++)
               nablaH(0,run1) -= l[run3][0][run1];
            for( run3 = 0; run3 < nOfNewtonSteps[1]; run3++)
               nablaH(0,run1) -= l[run3][1][run1];
            for( run3 = 0; run3 < nOfNewtonSteps[2]; run3++)
               nablaH(0,run1) -= l[run3][2][run1];
            for( run3 = 0; run3 < nOfNewtonSteps[3]; run3++)
               nablaH(0,run1) -= l[run3][3][run1];
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++)
               nablaH(0,run1) -=  l[run3][4][run1];
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++)
               nablaH(0,run1) -= l[run3][5][run1];
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++)
               nablaH(0,run1) -= l[run3][6][run1];
        }

    delete[] Hstore;
}


void AcadoIntegratorBackend::determineRKEtaGForward2(){

    int run1, run2, run3, newtonsteps;


    for( run2 = 0; run2 < m; run2++ ){
        k [0][0][run2] = 0.0;
        k2[0][0][run2] = 0.0;
    }


    // determine k:
    // -----------------------------------------------
       for( run1 = 0; run1 < dim; run1++ ){

           if( run1 > 0 ){
               for( run2 = 0; run2 < m; run2++ ){
                   k [0][run1][run2] = k [nOfNewtonSteps[run1-1]][run1-1][run2];
                   k2[0][run1][run2] = k2[nOfNewtonSteps[run1-1]][run1-1][run2];
               }
           }

           newtonsteps = 0;
           while( newtonsteps < nOfNewtonSteps[run1] ){

               for( run2 = 0; run2 < md; run2++ ){
                   G2[diff_index[run2]] = etaG2[run2];
                   G3[diff_index[run2]] = etaG3[run2];
                   for( run3 = 0; run3 < run1; run3++ ){
                       G2[diff_index[run2]] = G2[diff_index[run2]] +
                                   A[run1][run3]*h[0]*k [nOfNewtonSteps[run3]][run3][run2];
                       G3[diff_index[run2]] = G3[diff_index[run2]] +
                                   A[run1][run3]*h[0]*k2[nOfNewtonSteps[run3]][run3][run2];
                   }
                   G2[diff_index[run2]] = G2[diff_index[run2]] +
                               A[run1][run1]*h[0]*k [newtonsteps][run1][run2];
                   G3[diff_index[run2]] = G3[diff_index[run2]] +
                               A[run1][run1]*h[0]*k2[newtonsteps][run1][run2];
               }
               for( run2 = 0; run2 < ma; run2++ ){
                   G2[diff_index[md+run2]] = k [newtonsteps][run1][md+run2];
                   G3[diff_index[md+run2]] = k2[newtonsteps][run1][md+run2];
               }
               for( run2 = 0; run2 < md; run2++ ){
                   G2[ddiff_index[run2]] = k [newtonsteps][run1][run2];
                   G3[ddiff_index[run2]] = k2[newtonsteps][run1][run2];
               }

               if( rhs[0].AD_forward2( 3*run1+newtonsteps, G2, G3, F, F2 )
                                      != SUCCESSFUL_RETURN ){
                   ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
                   return;
               }


               applyNewtonStep( k[newtonsteps+1][run1],
                                k[newtonsteps][run1]  ,
                               *M[M_index[run1]],
                                F );
               applyNewtonStep( k2[newtonsteps+1][run1],
                                k2[newtonsteps][run1]  ,
                               *M[M_index[run1]],
                                F2 );

               newtonsteps++;
           }
       }

    // determine nablaG 2,3:
    // ----------------------------------------------

       for( run1 = 0; run1 < nstep-1; run1++ ){
           for( run2 = 0; run2 < md; run2++ ){
               nablaG2(nstep-2-run1,run2) = etaG2[run2];
               nablaG3(nstep-2-run1,run2) = etaG3[run2];
               for( run3 = 0; run3 <= run1+3; run3++ ){
                   nablaG2(nstep-2-run1,run2) = nablaG2(nstep-2-run1,run2) +
                     h[0]*A[run1+3][run3]*k [nOfNewtonSteps[run3]][run3][run2];
                   nablaG3(nstep-2-run1,run2) = nablaG3(nstep-2-run1,run2) +
                     h[0]*A[run1+3][run3]*k2[nOfNewtonSteps[run3]][run3][run2];
               }
           }
           for( run2 = 0; run2 < ma; run2++ ){
               nablaG2(nstep-2-run1,md+run2) = k [nOfNewtonSteps[run1+3]][run1+3][md+run2];
               nablaG3(nstep-2-run1,md+run2) = k2[nOfNewtonSteps[run1+3]][run1+3][md+run2];
           }
        }
}


void AcadoIntegratorBackend::determineRKEtaHBackward2(){

    int run1, run2, run3, run4, newtonsteps;
    const double hhh = c[3]*h[0];
    int number_;

    double *H1store = new double[ndir];
    double *H2store = new double[ndir];

    const double scalT = rel_time_scale/hhh;

    for( run4 = 0; run4 < nBDirs2; run4++ ){

        for( run3 = 0; run3 < 4; run3++ ){
            for( run1 = 0; run1 < dim; run1++ ){
                for( run2 = 0; run2 < ndir; run2++ ){
                    l [run3][run1][run2] = 0.0;
                    l2[run3][run1][run2] = 0.0;
                }
            }
        }

        for( run1 = 0; run1 < ndir; run1++ ){
            H1store[run1] = nablaH2(0,run1);
            H2store[run1] = nablaH3(0,run1);
        }

        for( run1 = 0; run1 < ndir; run1++ ){

            zH2[5][run1] = nablaH2(3,run1)*scalT*scalT;
            zH2[4][run1] = nablaH2(2,run1)*scalT + zH2[5][run1]/(3.0);
            zH2[3][run1] = - zH2[5][run1]/(3.0);
            zH2[2][run1] = nablaH2(1,run1) + zH2[4][run1]/(2.0);
            zH2[1][run1] = (zH2[3][run1]-zH2[4][run1])/(2.0);
            zH2[0][run1] = - zH2[3][run1]/(2.0);

            zH3[5][run1] = nablaH3(3,run1)*scalT*scalT;
            zH3[4][run1] = nablaH3(2,run1)*scalT + zH3[5][run1]/(3.0);
            zH3[3][run1] = - zH3[5][run1]/(3.0);
            zH3[2][run1] = nablaH3(1,run1) + zH3[4][run1]/(2.0);
            zH3[1][run1] = (zH3[3][run1]-zH3[4][run1])/(2.0);
            zH3[0][run1] = - zH3[3][run1]/(2.0);

            nablaH2(3,run1) = nablaH2(0,run1)*h[0] + zH2[2][run1]*(h[0]/hhh);
            nablaH2(2,run1) = (zH2[1][run1]-zH2[2][run1])*(h[0]/hhh);
            nablaH2(1,run1) = (zH2[0][run1]-zH2[1][run1])*(h[0]/hhh);
            nablaH2(0,run1) = -zH2[0][run1]*(h[0]/hhh);

            nablaH3(3,run1) = nablaH3(0,run1)*h[0] + zH3[2][run1]*(h[0]/hhh);
            nablaH3(2,run1) = (zH3[1][run1]-zH3[2][run1])*(h[0]/hhh);
            nablaH3(1,run1) = (zH3[0][run1]-zH3[1][run1])*(h[0]/hhh);
            nablaH3(0,run1) = -zH3[0][run1]*(h[0]/hhh);
        }


        for( run1 = 0; run1 < ndir; run1++ ){
            kH2[6][nOfNewtonSteps[6]][run1] = nablaH2(3,run1)/h[0];
            kH3[6][nOfNewtonSteps[6]][run1] = nablaH3(3,run1)/h[0];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH2[6][nOfNewtonSteps[6]][diff_index[run1]] =
            kH2[6][nOfNewtonSteps[6]][diff_index[run1]]*h[0]*A[6][6];
            kH3[6][nOfNewtonSteps[6]][diff_index[run1]] =
            kH3[6][nOfNewtonSteps[6]][diff_index[run1]]*h[0]*A[6][6];
        }

        newtonsteps = nOfNewtonSteps[6];
        newtonsteps--;

        number_ = 6;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH2[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H2 );

            applyMTranspose( kH3[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H3 );


            if( rhs[0].AD_backward2( 3*number_+newtonsteps, H2, H3,
                                    l[newtonsteps][number_], l2[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH2[number_][newtonsteps][run1] =   kH2[number_][newtonsteps+1][run1]
                                                 - l[newtonsteps][number_][run1];
                kH3[number_][newtonsteps][run1] =   kH3[number_][newtonsteps+1][run1]
                                                 - l2[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH2[number_][newtonsteps][diff_index[run1]] = 
                   kH2[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
                kH3[number_][newtonsteps][diff_index[run1]] = 
                   kH3[number_][newtonsteps+1][diff_index[run1]]
                 - l2[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l2[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }

        for( run1 = 0; run1 < ndir; run1++ ){
            kH2[5][nOfNewtonSteps[5]][run1] = kH2[6][0][run1] + nablaH2(2,run1)/h[0];
            kH3[5][nOfNewtonSteps[5]][run1] = kH3[6][0][run1] + nablaH3(2,run1)/h[0];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH2[5][nOfNewtonSteps[5]][diff_index[run1]] = kH2[6][0][diff_index[run1]]
            +  nablaH2(3,diff_index[run1])*A[6][5]
            +  nablaH2(2,diff_index[run1])*A[5][5];
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH2[5][nOfNewtonSteps[5]][diff_index[run1]] =
               kH2[5][nOfNewtonSteps[5]][diff_index[run1]]
               -l[run3][number_][diff_index[run1]]*h[0]*A[6][5];
            }
            kH3[5][nOfNewtonSteps[5]][diff_index[run1]] = kH3[6][0][diff_index[run1]]
            +  nablaH3(3,diff_index[run1])*A[6][5]
            +  nablaH3(2,diff_index[run1])*A[5][5];
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH3[5][nOfNewtonSteps[5]][diff_index[run1]] =
               kH3[5][nOfNewtonSteps[5]][diff_index[run1]]
               -l2[run3][number_][diff_index[run1]]*h[0]*A[6][5];
            }
        }


        newtonsteps = nOfNewtonSteps[5];
        newtonsteps--;

        number_ = 5;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH2[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H2 );
            applyMTranspose( kH3[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H3 );

            if( rhs[0].AD_backward2( 3*number_+newtonsteps, H2, H3,
                                     l[newtonsteps][number_], l2[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH2[number_][newtonsteps][run1] =   kH2[number_][newtonsteps+1][run1]
                                                 - l[newtonsteps][number_][run1];
                kH3[number_][newtonsteps][run1] =   kH3[number_][newtonsteps+1][run1]
                                                 - l2[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH2[number_][newtonsteps][diff_index[run1]] = 
                   kH2[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
                kH3[number_][newtonsteps][diff_index[run1]] = 
                   kH3[number_][newtonsteps+1][diff_index[run1]]
                 - l2[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l2[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }


        for( run1 = 0; run1 < ndir; run1++ ){
            kH2[4][nOfNewtonSteps[4]][run1] = kH2[5][0][run1] + nablaH2(1,run1)/h[0];
            kH3[4][nOfNewtonSteps[4]][run1] = kH3[5][0][run1] + nablaH3(1,run1)/h[0];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH2[4][nOfNewtonSteps[4]][diff_index[run1]] = kH2[5][0][diff_index[run1]]
            +  nablaH2(3,diff_index[run1])*A[6][4]
            +  nablaH2(2,diff_index[run1])*A[5][4]
            +  nablaH2(1,diff_index[run1])*A[4][4];
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH2[4][nOfNewtonSteps[4]][diff_index[run1]] =
               kH2[4][nOfNewtonSteps[4]][diff_index[run1]]
               -l[run3][5][diff_index[run1]]*h[0]*A[5][4];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH2[4][nOfNewtonSteps[4]][diff_index[run1]] =
               kH2[4][nOfNewtonSteps[4]][diff_index[run1]]
               -l[run3][6][diff_index[run1]]*h[0]*A[6][4];
            }
            kH3[4][nOfNewtonSteps[4]][diff_index[run1]] = kH3[5][0][diff_index[run1]]
            +  nablaH3(3,diff_index[run1])*A[6][4]
            +  nablaH3(2,diff_index[run1])*A[5][4]
            +  nablaH3(1,diff_index[run1])*A[4][4];
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH3[4][nOfNewtonSteps[4]][diff_index[run1]] =
               kH3[4][nOfNewtonSteps[4]][diff_index[run1]]
               -l2[run3][5][diff_index[run1]]*h[0]*A[5][4];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH3[4][nOfNewtonSteps[4]][diff_index[run1]] =
               kH3[4][nOfNewtonSteps[4]][diff_index[run1]]
               -l2[run3][6][diff_index[run1]]*h[0]*A[6][4];
            }
        }


        newtonsteps = nOfNewtonSteps[4];
        newtonsteps--;

        number_ = 4;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH2[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H2 );
            applyMTranspose( kH3[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H3 );

            if( rhs[0].AD_backward2( 3*number_+newtonsteps, H2, H3,
                                     l[newtonsteps][number_], l2[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH2[number_][newtonsteps][run1] =   kH2[number_][newtonsteps+1][run1]
                                                 - l[newtonsteps][number_][run1];
                kH3[number_][newtonsteps][run1] =   kH3[number_][newtonsteps+1][run1]
                                                 - l2[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH2[number_][newtonsteps][diff_index[run1]] = 
                   kH2[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
                kH3[number_][newtonsteps][diff_index[run1]] = 
                   kH3[number_][newtonsteps+1][diff_index[run1]]
                 - l2[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l2[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }


        for( run1 = 0; run1 < ndir; run1++ ){
            kH2[3][nOfNewtonSteps[3]][run1] = kH2[4][0][run1] + nablaH2(0,run1)/h[0];
            kH3[3][nOfNewtonSteps[3]][run1] = kH3[4][0][run1] + nablaH3(0,run1)/h[0];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH2[3][nOfNewtonSteps[3]][diff_index[run1]] = kH2[4][0][diff_index[run1]]
            +  nablaH2(3,diff_index[run1])*A[6][3]
            +  nablaH2(2,diff_index[run1])*A[5][3]
            +  nablaH2(1,diff_index[run1])*A[4][3]
            +  nablaH2(0,diff_index[run1])*A[3][3];
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++){
               kH2[3][nOfNewtonSteps[3]][diff_index[run1]] =
               kH2[3][nOfNewtonSteps[3]][diff_index[run1]]
               -l[run3][4][diff_index[run1]]*h[0]*A[4][3];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH2[3][nOfNewtonSteps[3]][diff_index[run1]] =
               kH2[3][nOfNewtonSteps[3]][diff_index[run1]]
               -l[run3][5][diff_index[run1]]*h[0]*A[5][3];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH2[3][nOfNewtonSteps[3]][diff_index[run1]] =
               kH2[3][nOfNewtonSteps[3]][diff_index[run1]]
               -l[run3][6][diff_index[run1]]*h[0]*A[6][3];
            }
            kH3[3][nOfNewtonSteps[3]][diff_index[run1]] = kH3[4][0][diff_index[run1]]
            +  nablaH3(3,diff_index[run1])*A[6][3]
            +  nablaH3(2,diff_index[run1])*A[5][3]
            +  nablaH3(1,diff_index[run1])*A[4][3]
            +  nablaH3(0,diff_index[run1])*A[3][3];
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++){
               kH3[3][nOfNewtonSteps[3]][diff_index[run1]] =
               kH3[3][nOfNewtonSteps[3]][diff_index[run1]]
               -l2[run3][4][diff_index[run1]]*h[0]*A[4][3];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH3[3][nOfNewtonSteps[3]][diff_index[run1]] =
               kH3[3][nOfNewtonSteps[3]][diff_index[run1]]
               -l2[run3][5][diff_index[run1]]*h[0]*A[5][3];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH3[3][nOfNewtonSteps[3]][diff_index[run1]] =
               kH3[3][nOfNewtonSteps[3]][diff_index[run1]]
               -l2[run3][6][diff_index[run1]]*h[0]*A[6][3];
            }
        }


        newtonsteps = nOfNewtonSteps[3];
        newtonsteps--;

        number_ = 3;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH2[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H2 );
            applyMTranspose( kH3[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H3 );

            if( rhs[0].AD_backward2( 3*number_+newtonsteps, H2, H3,
                                     l[newtonsteps][number_], l2[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH2[number_][newtonsteps][run1] =   kH2[number_][newtonsteps+1][run1]
                                                 - l[newtonsteps][number_][run1];
                kH3[number_][newtonsteps][run1] =   kH3[number_][newtonsteps+1][run1]
                                                 - l2[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH2[number_][newtonsteps][diff_index[run1]] = 
                   kH2[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
                kH3[number_][newtonsteps][diff_index[run1]] = 
                   kH3[number_][newtonsteps+1][diff_index[run1]]
                 - l2[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l2[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }


        for( run1 = 0; run1 < ndir; run1++ ){
            kH2[2][nOfNewtonSteps[2]][run1] = kH2[3][0][run1];
            kH3[2][nOfNewtonSteps[2]][run1] = kH3[3][0][run1];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH2[2][nOfNewtonSteps[2]][diff_index[run1]] = kH2[3][0][diff_index[run1]]
            +  nablaH2(3,diff_index[run1])*A[6][2]
            +  nablaH2(2,diff_index[run1])*A[5][2]
            +  nablaH2(1,diff_index[run1])*A[4][2]
            +  nablaH2(0,diff_index[run1])*A[3][2];
            for( run3 = 0; run3 < nOfNewtonSteps[3]; run3++){
               kH2[2][nOfNewtonSteps[2]][diff_index[run1]] =
               kH2[2][nOfNewtonSteps[2]][diff_index[run1]]
               -l[run3][3][diff_index[run1]]*h[0]*A[3][2];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++){
               kH2[2][nOfNewtonSteps[2]][diff_index[run1]] =
               kH2[2][nOfNewtonSteps[2]][diff_index[run1]]
               -l[run3][4][diff_index[run1]]*h[0]*A[4][2];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH2[2][nOfNewtonSteps[2]][diff_index[run1]] =
               kH2[2][nOfNewtonSteps[2]][diff_index[run1]]
               -l[run3][5][diff_index[run1]]*h[0]*A[5][2];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH2[2][nOfNewtonSteps[2]][diff_index[run1]] =
               kH2[2][nOfNewtonSteps[2]][diff_index[run1]]
               -l[run3][6][diff_index[run1]]*h[0]*A[6][2];
            }
            kH3[2][nOfNewtonSteps[2]][diff_index[run1]] = kH3[3][0][diff_index[run1]]
            +  nablaH3(3,diff_index[run1])*A[6][2]
            +  nablaH3(2,diff_index[run1])*A[5][2]
            +  nablaH3(1,diff_index[run1])*A[4][2]
            +  nablaH3(0,diff_index[run1])*A[3][2];
            for( run3 = 0; run3 < nOfNewtonSteps[3]; run3++){
               kH3[2][nOfNewtonSteps[2]][diff_index[run1]] =
               kH3[2][nOfNewtonSteps[2]][diff_index[run1]]
               -l2[run3][3][diff_index[run1]]*h[0]*A[3][2];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++){
               kH3[2][nOfNewtonSteps[2]][diff_index[run1]] =
               kH3[2][nOfNewtonSteps[2]][diff_index[run1]]
               -l2[run3][4][diff_index[run1]]*h[0]*A[4][2];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH3[2][nOfNewtonSteps[2]][diff_index[run1]] =
               kH3[2][nOfNewtonSteps[2]][diff_index[run1]]
               -l2[run3][5][diff_index[run1]]*h[0]*A[5][2];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH3[2][nOfNewtonSteps[2]][diff_index[run1]] =
               kH3[2][nOfNewtonSteps[2]][diff_index[run1]]
               -l2[run3][6][diff_index[run1]]*h[0]*A[6][2];
            }
        }


        newtonsteps = nOfNewtonSteps[2];
        newtonsteps--;

        number_ = 2;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH2[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H2 );
            applyMTranspose( kH3[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H3 );

            if( rhs[0].AD_backward2( 3*number_+newtonsteps, H2, H3,
                                     l[newtonsteps][number_], l2[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH2[number_][newtonsteps][run1] =   kH2[number_][newtonsteps+1][run1]
                                                 - l[newtonsteps][number_][run1];
                kH3[number_][newtonsteps][run1] =   kH3[number_][newtonsteps+1][run1]
                                                 - l2[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH2[number_][newtonsteps][diff_index[run1]] = 
                   kH2[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
                kH3[number_][newtonsteps][diff_index[run1]] = 
                   kH3[number_][newtonsteps+1][diff_index[run1]]
                 - l2[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l2[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }


        for( run1 = 0; run1 < ndir; run1++ ){
            kH2[1][nOfNewtonSteps[1]][run1] = kH2[2][0][run1];
            kH3[1][nOfNewtonSteps[1]][run1] = kH3[2][0][run1];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH2[1][nOfNewtonSteps[1]][diff_index[run1]] = kH2[2][0][diff_index[run1]]
            +  nablaH2(3,diff_index[run1])*A[6][1]
            +  nablaH2(2,diff_index[run1])*A[5][1]
            +  nablaH2(1,diff_index[run1])*A[4][1]
            +  nablaH2(0,diff_index[run1])*A[3][1];
            for( run3 = 0; run3 < nOfNewtonSteps[2]; run3++){
               kH2[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH2[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l[run3][2][diff_index[run1]]*h[0]*A[2][1];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[3]; run3++){
               kH2[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH2[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l[run3][3][diff_index[run1]]*h[0]*A[3][1];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++){
               kH2[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH2[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l[run3][4][diff_index[run1]]*h[0]*A[4][1];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH2[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH2[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l[run3][5][diff_index[run1]]*h[0]*A[5][1];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH2[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH2[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l[run3][6][diff_index[run1]]*h[0]*A[6][1];
            }
            kH3[1][nOfNewtonSteps[1]][diff_index[run1]] = kH3[2][0][diff_index[run1]]
            +  nablaH3(3,diff_index[run1])*A[6][1]
            +  nablaH3(2,diff_index[run1])*A[5][1]
            +  nablaH3(1,diff_index[run1])*A[4][1]
            +  nablaH3(0,diff_index[run1])*A[3][1];
            for( run3 = 0; run3 < nOfNewtonSteps[2]; run3++){
               kH3[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH3[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l2[run3][2][diff_index[run1]]*h[0]*A[2][1];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[3]; run3++){
               kH3[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH3[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l2[run3][3][diff_index[run1]]*h[0]*A[3][1];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++){
               kH3[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH3[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l2[run3][4][diff_index[run1]]*h[0]*A[4][1];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH3[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH3[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l2[run3][5][diff_index[run1]]*h[0]*A[5][1];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH3[1][nOfNewtonSteps[1]][diff_index[run1]] =
               kH3[1][nOfNewtonSteps[1]][diff_index[run1]]
               -l2[run3][6][diff_index[run1]]*h[0]*A[6][1];
            }
        }


        newtonsteps = nOfNewtonSteps[1];
        newtonsteps--;

        number_ = 1;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH2[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H2 );
            applyMTranspose( kH3[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H3 );

            if( rhs[0].AD_backward2( 3*number_+newtonsteps, H2, H3,
                                     l[newtonsteps][number_], l2[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH2[number_][newtonsteps][run1] =   kH2[number_][newtonsteps+1][run1]
                                                 - l[newtonsteps][number_][run1];
                kH3[number_][newtonsteps][run1] =   kH3[number_][newtonsteps+1][run1]
                                                 - l2[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH2[number_][newtonsteps][diff_index[run1]] = 
                   kH2[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
                kH3[number_][newtonsteps][diff_index[run1]] = 
                   kH3[number_][newtonsteps+1][diff_index[run1]]
                 - l2[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l2[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }



        for( run1 = 0; run1 < ndir; run1++ ){
            kH2[0][nOfNewtonSteps[0]][run1] = kH2[1][0][run1];
            kH3[0][nOfNewtonSteps[0]][run1] = kH3[1][0][run1];
        }
        for( run1 = 0; run1 < md; run1++ ){
            kH2[0][nOfNewtonSteps[0]][diff_index[run1]] = kH2[1][0][diff_index[run1]]
            +  nablaH2(3,diff_index[run1])*A[6][0]
            +  nablaH2(2,diff_index[run1])*A[5][0]
            +  nablaH2(1,diff_index[run1])*A[4][0]
            +  nablaH2(0,diff_index[run1])*A[3][0];
            for( run3 = 0; run3 < nOfNewtonSteps[1]; run3++){
               kH2[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH2[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l[run3][1][diff_index[run1]]*h[0]*A[1][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[2]; run3++){
               kH2[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH2[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l[run3][2][diff_index[run1]]*h[0]*A[2][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[3]; run3++){
               kH2[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH2[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l[run3][3][diff_index[run1]]*h[0]*A[3][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++){
               kH2[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH2[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l[run3][4][diff_index[run1]]*h[0]*A[4][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH2[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH2[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l[run3][5][diff_index[run1]]*h[0]*A[5][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH2[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH2[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l[run3][6][diff_index[run1]]*h[0]*A[6][0];
            }
            kH3[0][nOfNewtonSteps[0]][diff_index[run1]] = kH3[1][0][diff_index[run1]]
            +  nablaH3(3,diff_index[run1])*A[6][0]
            +  nablaH3(2,diff_index[run1])*A[5][0]
            +  nablaH3(1,diff_index[run1])*A[4][0]
            +  nablaH3(0,diff_index[run1])*A[3][0];
            for( run3 = 0; run3 < nOfNewtonSteps[1]; run3++){
               kH3[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH3[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l2[run3][1][diff_index[run1]]*h[0]*A[1][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[2]; run3++){
               kH3[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH3[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l2[run3][2][diff_index[run1]]*h[0]*A[2][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[3]; run3++){
               kH3[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH3[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l2[run3][3][diff_index[run1]]*h[0]*A[3][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[4]; run3++){
               kH3[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH3[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l2[run3][4][diff_index[run1]]*h[0]*A[4][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[5]; run3++){
               kH3[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH3[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l2[run3][5][diff_index[run1]]*h[0]*A[5][0];
            }
            for( run3 = 0; run3 < nOfNewtonSteps[6]; run3++){
               kH3[0][nOfNewtonSteps[0]][diff_index[run1]] =
               kH3[0][nOfNewtonSteps[0]][diff_index[run1]]
               -l2[run3][6][diff_index[run1]]*h[0]*A[6][0];
            }
        }

        newtonsteps = nOfNewtonSteps[0];
        newtonsteps--;

        number_ = 0;

        while( newtonsteps >= 0 ){

            applyMTranspose( kH2[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H2 );
            applyMTranspose( kH3[number_][newtonsteps+1],
                            *M[M_index[number_]],
                             H3 );

            if( rhs[0].AD_backward2( 3*number_+newtonsteps, H2, H3,
                                     l[newtonsteps][number_], l2[newtonsteps][number_] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }


            for( run1 = 0; run1 < ndir; run1++ ){
                kH2[number_][newtonsteps][run1] =   kH2[number_][newtonsteps+1][run1]
                                                 - l[newtonsteps][number_][run1];
                kH3[number_][newtonsteps][run1] =   kH3[number_][newtonsteps+1][run1]
                                                 - l2[newtonsteps][number_][run1];
            }
            for( run1 = 0; run1 < md; run1++ ){
                kH2[number_][newtonsteps][diff_index[run1]] = 
                   kH2[number_][newtonsteps+1][diff_index[run1]]
                 - l[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l[newtonsteps][number_][ddiff_index[run1]];
                kH3[number_][newtonsteps][diff_index[run1]] = 
                   kH3[number_][newtonsteps+1][diff_index[run1]]
                 - l2[newtonsteps][number_][diff_index[run1]]*h[0]*A[number_][number_]
                 - l2[newtonsteps][number_][ddiff_index[run1]];
            }

            newtonsteps--;
        }


        for( run1 = 0; run1 < ndir; run1++ ){
            nablaH2(0,run1) = H1store[run1];
            nablaH3(0,run1) = H2store[run1];
        }

        for( run1 = 0; run1 < ndir; run1++ ){
            for( run3 = 0; run3 < nOfNewtonSteps[0]; run3++){
               for( run2 = 0; run2 < dim; run2++ ){
                   nablaH2(0,run1) = nablaH2(0,run1)
                   -l[run3][run2][run1];
                   nablaH3(0,run1) = nablaH3(0,run1)
                  -l2[run3][run2][run1];
               }
            }
        }

    }
    delete[] H1store;
    delete[] H2store;
}



void AcadoIntegratorBackend::determineBDFEtaGForward( int number_ ){

    int run1, run2, run4;
    int newtonsteps;
    double pp;

    const double scalT = psi[number_][0]/rel_time_scale;

    for( run4 = 0; run4 < nFDirs; run4++ ){

        for( run1 = 0; run1 < m; run1++ ){
            phiG(1,run1) = nablaG(1,run1)*scalT;
        }

        pp = (psi[number_][1]/psi[number_][0]);
        for( run1 = 0; run1 < m; run1++ ){
            phiG(2,run1) = pp*nablaG(2,run1)*scalT*scalT;
        }

        pp = pp*(psi[number_][2]/psi[number_][0]);
        for( run1 = 0; run1 < m; run1++ ){
            phiG(3,run1) = pp*nablaG(3,run1)*scalT*scalT*scalT;
        }

        pp = pp*(psi[number_][3]/psi[number_][0]);
        for( run1 = 0; run1 < m; run1++ ){
            phiG(4,run1) = pp*nablaG(4,run1)*scalT*scalT*scalT*scalT;
        }

        for( run1 = 0; run1 < m; run1++ ){
            deltaG(1,run1) = nablaG(0,run1) + phiG(1,run1);
        }

        for( run1 = 0; run1 < m; run1++ ){
            deltaG(2,run1) = deltaG(1,run1) + phiG(2,run1);
        }

        for( run1 = 0; run1 < m; run1++ ){
            deltaG(3,run1) = deltaG(2,run1) + phiG(3,run1);
        }

        for( run1 = 0; run1 < m; run1++ ){
            eta[0][run1] = deltaG(3,run1) + phiG(4,run1);
        }


        for( run1 = 0; run1 < md; run1++ ){

            c2G(run1) = - nablaG (0,run1)/psi[number_][0]
                        - deltaG (1,run1)/psi[number_][1]
                        - deltaG (2,run1)/psi[number_][2]
                        - deltaG (3,run1)/psi[number_][3];
        }


        newtonsteps = 0;
        while( newtonsteps < nOfNewtonSteps[number_] ){
            for( run2 = 0; run2 < m; run2++ ){
                 G[diff_index[run2]] = eta[newtonsteps][run2];
            }
            for( run2 = 0; run2 < md; run2++ ){
                 G[ddiff_index[run2]] =
                 gamma[number_][4]*eta[newtonsteps][run2]+c2G(run2);
            }

            if( rhs[0].AD_forward( 3*number_+newtonsteps, G, F )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }

            applyNewtonStep( eta[newtonsteps+1],
                             eta[newtonsteps],
                            *M[M_index[number_]],
                             F );

            newtonsteps++;
        }

        for( run1 = 0; run1 < m; run1++ ){
            nablaY_(0,run1) = eta[nOfNewtonSteps[number_]][run1];
            nablaY_(1,run1) = nablaY_(0,run1) - nablaG(0,run1);
            nablaY_(2,run1) = (nablaY_(1,run1) - nablaG(1,run1)*scalT)*(psi[number_][0]/psi[number_][1]);
            nablaY_(3,run1) = (nablaY_(2,run1) - nablaG(2,run1)*scalT*scalT)*(psi[number_][0]/psi[number_][2]);
            nablaY_(4,run1) = (nablaY_(3,run1) - nablaG(3,run1)*scalT*scalT*scalT)*(psi[number_][0]/psi[number_][3]);
        }

        nablaG = nablaY_;
    }
}


void AcadoIntegratorBackend::determineBDFEtaHBackward( int number_ ){

    int run1, run2;
    int newtonsteps;
    double pp;

    const double scalT = rel_time_scale/psi[number_][0];

    newtonsteps = nOfNewtonSteps[number_];

    for( run1 = 0; run1 < newtonsteps; run1++ ){
        for( run2 = 0; run2 < ndir; run2++ ){
            l[run1][0][run2] = 0.0;
        }
    }

    for( run1 = 0; run1 < ndir; run1++ ){

        nablaH_(4,run1) = nablaH(4,run1)*scalT*scalT*scalT;
        nablaH_(3,run1) = nablaH(3,run1)*scalT*scalT + nablaH_(4,run1)*(psi[number_][0]/psi[number_][3]);
        nablaH_(2,run1) = nablaH(2,run1)*scalT + nablaH_(3,run1)*(psi[number_][0]/psi[number_][2]);
        nablaH_(1,run1) = nablaH(1,run1) + nablaH_(2,run1)*(psi[number_][0]/psi[number_][1]);
        nablaH_(0,run1) = nablaH(0,run1) + nablaH_(1,run1)/psi[number_][0];

        etaH[newtonsteps][run1] = nablaH_(0,run1);
        c2H(run1) = 0.0;
    }

    newtonsteps--;
    while( newtonsteps >= 0 ){

        applyMTranspose( etaH[newtonsteps+1], M[M_index[number_]][0], H );

        if( rhs[0].AD_backward( 3*number_+newtonsteps, H, l[newtonsteps][0] ) != SUCCESSFUL_RETURN )
            ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);

        for( run1 = 0; run1 < ndir; run1++ ){
            etaH[newtonsteps][run1] = etaH[newtonsteps+1][run1] - l[newtonsteps][0][run1];
        }

        for( run1 = 0; run1 < md; run1++ )
            etaH[newtonsteps][diff_index[run1]] -= l[newtonsteps][0][ddiff_index[run1]]*gamma[number_][4];

        for( run1 = 0; run1 < md; run1++ )
            c2H(diff_index[run1]) -= l[newtonsteps][0][ddiff_index[run1]];

        newtonsteps--;
    }


    for( run1 = 0; run1 < ndir; run1++ ){

        deltaH(3,run1) = etaH[0][run1]   - c2H(run1)/psi[number_][3];
        deltaH(2,run1) = deltaH(3,run1) - c2H(run1)/psi[number_][2];
        deltaH(1,run1) = deltaH(2,run1) - c2H(run1)/psi[number_][1];
    }

    for( run1 = 0; run1 < ndir; run1++ ){

        nablaH(0,run1) =  - nablaH_(1,run1) /psi[number_][0]
                                 - c2H(run1) /psi[number_][0]
                                 + deltaH(1,run1);
    }
    pp =  psi[number_][0];
    for( run1 = 0; run1 < ndir; run1++ ){
        nablaH(1,run1) =  - nablaH_(2,run1)*(psi[number_][0]/psi[number_][1])
                                     + pp*deltaH(1,run1);
    }
    pp = pp*(psi[number_][1]/psi[number_][0]);
    for( run1 = 0; run1 < ndir; run1++ ){
        nablaH(2,run1) =  - nablaH_(3,run1)*(psi[number_][0]/psi[number_][2])
                                 + pp*deltaH(2,run1);
    }
    pp = pp*(psi[number_][2]/psi[number_][0]);

    for( run1 = 0; run1 < ndir; run1++ ){
        nablaH(3,run1) =  - nablaH_(4,run1)*(psi[number_][0]/psi[number_][3])
                                 + pp*deltaH(3,run1);
    }
    pp = pp*(psi[number_][3]/psi[number_][0]);

    for( run1 = 0; run1 < ndir; run1++ ){
        nablaH(4,run1) =  pp*etaH[0][run1];
    }
}



void AcadoIntegratorBackend::determineBDFEtaGForward2( int number_ ){

    int run1, run2, run4;
    int newtonsteps;
    double pp;

    const double scalT = psi[number_][0]/rel_time_scale;

    for( run4 = 0; run4 < nFDirs2; run4++ ){

        for( run1 = 0; run1 < m; run1++ ){
            phiG2(1,run1) = nablaG2(1,run1)*scalT;
        }
        for( run1 = 0; run1 < m; run1++ ){
            phiG3(1,run1) = nablaG3(1,run1)*scalT;
        }

        pp = (psi[number_][1]/psi[number_][0]);
        for( run1 = 0; run1 < m; run1++ ){
            phiG2(2,run1) = pp*nablaG2(2,run1)*scalT*scalT;
        }
        for( run1 = 0; run1 < m; run1++ ){
            phiG3(2,run1) = pp*nablaG3(2,run1)*scalT*scalT;
        }

        pp = pp*(psi[number_][2]/psi[number_][0]);
        for( run1 = 0; run1 < m; run1++ ){
            phiG2(3,run1) = pp*nablaG2(3,run1)*scalT*scalT*scalT;
        }
        for( run1 = 0; run1 < m; run1++ ){
            phiG3(3,run1) = pp*nablaG3(3,run1)*scalT*scalT*scalT;
        }

        pp = pp*(psi[number_][3]/psi[number_][0]);
        for( run1 = 0; run1 < m; run1++ ){
            phiG2(4,run1) = pp*nablaG2(4,run1)*scalT*scalT*scalT*scalT;
        }
        for( run1 = 0; run1 < m; run1++ ){
            phiG3(4,run1) = pp*nablaG3(4,run1)*scalT*scalT*scalT*scalT;
        }

        for( run1 = 0; run1 < m; run1++ ){
            deltaG2(1,run1) = nablaG2(0,run1) + phiG2(1,run1);
        }
        for( run1 = 0; run1 < m; run1++ ){
            deltaG3(1,run1) = nablaG3(0,run1) + phiG3(1,run1);
        }

        for( run1 = 0; run1 < m; run1++ ){
            deltaG2(2,run1) = deltaG2(1,run1) + phiG2(2,run1);
        }
        for( run1 = 0; run1 < m; run1++ ){
            deltaG3(2,run1) = deltaG3(1,run1) + phiG3(2,run1);
        }

        for( run1 = 0; run1 < m; run1++ ){
            deltaG2(3,run1) = deltaG2(2,run1) + phiG2(3,run1);
        }
        for( run1 = 0; run1 < m; run1++ ){
            deltaG3(3,run1) = deltaG3(2,run1) + phiG3(3,run1);
        }

        for( run1 = 0; run1 < m; run1++ ){
            eta [0][run1] = deltaG2(3,run1) + phiG2(4,run1);
        }
        for( run1 = 0; run1 < m; run1++ ){
            eta2[0][run1] = deltaG3(3,run1) + phiG3(4,run1);
        }

        for( run1 = 0; run1 < md; run1++ ){

            c2G2(run1) = - nablaG2 (0,run1)/psi[number_][0]
                         - deltaG2 (1,run1)/psi[number_][1]
                         - deltaG2 (2,run1)/psi[number_][2]
                         - deltaG2 (3,run1)/psi[number_][3];
        }
        for( run1 = 0; run1 < md; run1++ ){

            c2G3(run1) = - nablaG3 (0,run1)/psi[number_][0]
                         - deltaG3 (1,run1)/psi[number_][1]
                         - deltaG3 (2,run1)/psi[number_][2]
                         - deltaG3 (3,run1)/psi[number_][3];
        }


        newtonsteps = 0;
        while( newtonsteps < nOfNewtonSteps[number_] ){
            for( run2 = 0; run2 < m; run2++ ){
                 G2[diff_index[run2]] = eta [newtonsteps][run2];
                 G3[diff_index[run2]] = eta2[newtonsteps][run2];
            }
            for( run2 = 0; run2 < md; run2++ ){
                 G2[ddiff_index[run2]] =
                 gamma[number_][4]*eta [newtonsteps][run2]+c2G2(run2);
                 G3[ddiff_index[run2]] =
                 gamma[number_][4]*eta2[newtonsteps][run2]+c2G3(run2);
            }

            if( rhs[0].AD_forward2( 3*number_+newtonsteps, G2, G3, F, F2 )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }

            applyNewtonStep( eta[newtonsteps+1],
                             eta[newtonsteps],
                            *M[M_index[number_]],
                             F );

            applyNewtonStep( eta2[newtonsteps+1],
                             eta2[newtonsteps],
                            *M[M_index[number_]],
                             F2 );

            newtonsteps++;
        }

        for( run1 = 0; run1 < m; run1++ ){
            nablaY_(0,run1) = eta[nOfNewtonSteps[number_]][run1];
            nablaY_(1,run1) = nablaY_(0,run1) - nablaG2(0,run1);
            nablaY_(2,run1) = (nablaY_(1,run1) - nablaG2(1,run1)*scalT)*(psi[number_][0]/psi[number_][1]);
            nablaY_(3,run1) = (nablaY_(2,run1) - nablaG2(2,run1)*scalT*scalT)*(psi[number_][0]/psi[number_][2]);
            nablaY_(4,run1) = (nablaY_(3,run1) - nablaG2(3,run1)*scalT*scalT*scalT)*(psi[number_][0]/psi[number_][3]);
        }
        nablaG2 = nablaY_;


        for( run1 = 0; run1 < m; run1++ ){
            nablaY_(0,run1) = eta2[nOfNewtonSteps[number_]][run1];
            nablaY_(1,run1) = nablaY_(0,run1) - nablaG3(0,run1);
            nablaY_(2,run1) = (nablaY_(1,run1) - nablaG3(1,run1)*scalT)*(psi[number_][0]/psi[number_][1]);
            nablaY_(3,run1) = (nablaY_(2,run1) - nablaG3(2,run1)*scalT*scalT)*(psi[number_][0]/psi[number_][2]);
            nablaY_(4,run1) = (nablaY_(3,run1) - nablaG3(3,run1)*scalT*scalT*scalT)*(psi[number_][0]/psi[number_][3]);
        }
        nablaG3 = nablaY_;
    }
}



void AcadoIntegratorBackend::determineBDFEtaHBackward2( int number_ ){

    int run1, run2, run4;
    int newtonsteps;
    double pp;

    const double scalT = rel_time_scale/psi[number_][0];

    for( run4 = 0; run4 < nBDirs2; run4++ ){

        newtonsteps = nOfNewtonSteps[number_];

        for( run1 = 0; run1 < newtonsteps; run1++ ){
            for( run2 = 0; run2 < ndir; run2++ ){
                l [run1][0][run2] = 0.0;
                l2[run1][0][run2] = 0.0;
            }
        }

        for( run1 = 0; run1 < ndir; run1++ ){

            nablaH2_(4,run1) = nablaH2(4,run1)*scalT*scalT*scalT;
            nablaH2_(3,run1) = nablaH2(3,run1)*scalT*scalT + nablaH2_(4,run1)*(psi[number_][0]/psi[number_][3]);
            nablaH2_(2,run1) = nablaH2(2,run1)*scalT + nablaH2_(3,run1)*(psi[number_][0]/psi[number_][2]);
            nablaH2_(1,run1) = nablaH2(1,run1) + nablaH2_(2,run1)*(psi[number_][0]/psi[number_][1]);
            nablaH2_(0,run1) = nablaH2(0,run1) + nablaH2_(1,run1)/psi[number_][0];

            nablaH3_(4,run1) = nablaH3(4,run1)*scalT*scalT*scalT;
            nablaH3_(3,run1) = nablaH3(3,run1)*scalT*scalT + nablaH3_(4,run1)*(psi[number_][0]/psi[number_][3]);
            nablaH3_(2,run1) = nablaH3(2,run1)*scalT + nablaH3_(3,run1)*(psi[number_][0]/psi[number_][2]);
            nablaH3_(1,run1) = nablaH3(1,run1) + nablaH3_(2,run1)*(psi[number_][0]/psi[number_][1]);
            nablaH3_(0,run1) = nablaH3(0,run1) + nablaH3_(1,run1)/psi[number_][0];

            etaH2[newtonsteps][run1] = nablaH2_(0,run1);
            etaH3[newtonsteps][run1] = nablaH3_(0,run1);
            c2H2              (run1) = 0.0;
            c2H3              (run1) = 0.0;
        }

        newtonsteps--;
        while( newtonsteps >= 0 ){

            applyMTranspose( etaH2[newtonsteps+1],
                            *M[M_index[number_]],
                             H2 );

            applyMTranspose( etaH3[newtonsteps+1],
                            *M[M_index[number_]],
                             H3 );

            if( rhs[0].AD_backward2( 3*number_+newtonsteps, H2, H3,
                                     l[newtonsteps][0], l2[newtonsteps][0] )
                != SUCCESSFUL_RETURN ){
                ACADOERROR(RET_UNSUCCESSFUL_RETURN_FROM_INTEGRATOR_BDF);
            }

            for( run1 = 0; run1 < ndir; run1++ ){
                etaH2[newtonsteps][run1] = etaH2[newtonsteps+1][run1] -
                                                 l [newtonsteps][0][run1];
                etaH3[newtonsteps][run1] = etaH3[newtonsteps+1][run1] -
                                                 l2[newtonsteps][0][run1];
            }

            for( run1 = 0; run1 < md; run1++ ){
                etaH2[newtonsteps][diff_index[run1]] =
                etaH2[newtonsteps][diff_index[run1]] -
                l[newtonsteps][0][ddiff_index[run1]]*gamma[number_][4];
                etaH3[newtonsteps][diff_index[run1]] =
                etaH3[newtonsteps][diff_index[run1]] -
                l2[newtonsteps][0][ddiff_index[run1]]*gamma[number_][4];
            }

            for( run1 = 0; run1 < md; run1++ ){
                c2H2(diff_index[run1]) -= l[newtonsteps][0][ddiff_index[run1]];
                c2H3(diff_index[run1]) -= l2[newtonsteps][0][ddiff_index[run1]];
            }

            newtonsteps--;
        }

        for( run1 = 0; run1 < ndir; run1++ ){

            deltaH2(3,run1) = etaH2[0][run1]   - c2H2(run1)/psi[number_][3];
            deltaH2(2,run1) = deltaH2(3,run1) - c2H2(run1)/psi[number_][2];
            deltaH2(1,run1) = deltaH2(2,run1) - c2H2(run1)/psi[number_][1];

            deltaH3(3,run1) = etaH3[0][run1]   - c2H3(run1)/psi[number_][3];
            deltaH3(2,run1) = deltaH3(3,run1) - c2H3(run1)/psi[number_][2];
            deltaH3(1,run1) = deltaH3(2,run1) - c2H3(run1)/psi[number_][1];
        }

        for( run1 = 0; run1 < ndir; run1++ ){
            nablaH2(0,run1) =  - nablaH2_(1,run1) /psi[number_][0]
                                      - c2H2(run1) /psi[number_][0]
                                      + deltaH2(1,run1);
            nablaH3(0,run1) =  - nablaH3_(1,run1) /psi[number_][0]
                                      - c2H3(run1) /psi[number_][0]
                                      + deltaH3(1,run1);
        }

        pp = psi[number_][0];
        for( run1 = 0; run1 < ndir; run1++ ){
            nablaH2(1,run1) =  - nablaH2_(2,run1)*(psi[number_][0]/psi[number_][1])
                                      + pp*deltaH2(1,run1);
            nablaH3(1,run1) =  - nablaH3_(2,run1)*(psi[number_][0]/psi[number_][1])
                                      + pp*deltaH3(1,run1);
        }
        pp = pp*(psi[number_][1]/psi[number_][0]);
        for( run1 = 0; run1 < ndir; run1++ ){
            nablaH2(2,run1) =  - nablaH2_(3,run1)*(psi[number_][0]/psi[number_][2])
                                      + pp*deltaH2(2,run1);
            nablaH3(2,run1) =  - nablaH3_(3,run1)*(psi[number_][0]/psi[number_][2])
                                      + pp*deltaH3(2,run1);
        }
        pp = pp*(psi[number_][2]/psi[number_][0]);
        for( run1 = 0; run1 < ndir; run1++ ){
            nablaH2(3,run1) =  - nablaH2_(4,run1)*(psi[number_][0]/psi[number_][3])
                                      + pp*deltaH2(3,run1);
            nablaH3(3,run1) =  - nablaH3_(4,run1)*(psi[number_][0]/psi[number_][3])
                                      + pp*deltaH3(3,run1);
        }
        pp = pp*(psi[number_][3]/psi[number_][0]);
        for( run1 = 0; run1 < ndir; run1++ ){
            nablaH2(4,run1) =  pp*etaH2[0][run1];
            nablaH3(4,run1) =  pp*etaH3[0][run1];
        }
    }
}



void AcadoIntegratorBackend::initializeButcherTableau(){


    A[0][0] =  0.096   ;
    A[0][1] =  0.0     ;
    A[0][2] =  0.0     ;
    A[0][3] =  0.0     ;
    A[0][4] =  0.0     ;
    A[0][5] =  0.0     ;
    A[0][6] =  0.0     ;

    A[1][0] =  0.104   ;
    A[1][1] =  0.096   ;
    A[1][2] =  0.0     ;
    A[1][3] =  0.0     ;
    A[1][4] =  0.0     ;
    A[1][5] =  0.0     ;
    A[1][6] =  0.0     ;

    A[2][0] =  0.1310450349650284    ;
    A[2][1] =  0.07295496503497158102;
    A[2][2] =  0.096   ;
    A[2][3] =  0.0     ;
    A[2][4] =  0.0     ;
    A[2][5] =  0.0     ;
    A[2][6] =  0.0     ;

    A[3][0] =  0.2199597787833081951;
    A[3][1] = -0.1447179487179487179;
    A[3][2] =  0.02875816993464052288;
    A[3][3] =  0.096   ;
    A[3][4] =  0.0     ;
    A[3][5] =  0.0     ;
    A[3][6] =  0.0     ;

    A[4][0] =  0.1608848667672197084;
    A[4][1] = -0.0512400094214750;
    A[4][2] = -0.02467973856209150327;
    A[4][3] =  0.2190348812163467684;
    A[4][4] =  0.096   ;
    A[4][5] =  0.0     ;
    A[4][6] =  0.0     ;

    A[5][0] =  0.1761704161863276;
    A[5][1] = -0.2376951929172075;
    A[5][2] =  0.1249803878146932;
    A[5][3] =  0.3034259664066430296;
    A[5][4] =  0.1371184225095437181;
    A[5][5] =  0.096   ;
    A[5][6] =  0.0     ;

    A[6][0] =  0.1822523881347410759;
    A[6][1] = -0.3465441147470548;
    A[6][2] =  0.213542483660130719;
    A[6][3] =  0.3547492429521829614;
    A[6][4] =  0.1     ;
    A[6][5] =  0.2     ;
    A[6][6] =  0.096   ;


    b4[0] =  0.0;
    b4[1] =  0.0;
    b4[2] =  0.0;
    b4[3] =  0.4583333333333333333;
    b4[4] =  0.04166666666666666667;
    b4[5] =  0.04166666666666666667;
    b4[6] =  0.4583333333333333333;

    b5[0] = -0.3413924474339141;
    b5[1] =  0.8554210048117954;
    b5[2] = -0.5726796951403962;
    b5[3] =  0.4498353628579520;
    b5[4] =  0.0888157749045627;
    b5[5] =  0.0400000000000000;
    b5[6] =  0.4800000000000000;

    c[0] = 0.096;
    c[1] = 0.2;
    c[2] = 0.3;
    c[3] = 0.2;
    c[4] = 0.4;
    c[5] = 0.6;
    c[6] = 0.8;
}



void AcadoIntegratorBackend::printBDFfinalResults2( Matrix &div ){

    int run2;

    acadoPrintf("w.r.t. the states:\n");
    for( run2 = 0; run2 < md; run2++ ){
        acadoPrintf("%.16e  ", div(0,diff_index[run2]) );
    }
    acadoPrintf("\n");

    if( mu > 0 ){
        acadoPrintf("w.r.t. the controls:\n");
        for( run2 = 0; run2 < mu; run2++ ){
            acadoPrintf("%.16e  ", div(0,control_index[run2]) );
        }
        acadoPrintf("\n");
    }
    if( mp > 0 ){
        acadoPrintf("w.r.t. the parameters:\n");
        for( run2 = 0; run2 < mp; run2++ ){
            acadoPrintf("%.16e  ", div(0,parameter_index[run2]) );
        }
        acadoPrintf("\n");
    }
    if( mw > 0 ){
        acadoPrintf("w.r.t. the diturbances:\n");
        for( run2 = 0; run2 < mw; run2++ ){
            acadoPrintf("%.16e  ", div(0,disturbance_index[run2]) );
        }
        acadoPrintf("\n");
    }
}



void AcadoIntegratorBackend::printBDFfinalResults(){

    int run1;

        // Forward Sensitivities:
        // ----------------------

        if( nFDirs > 0 && nBDirs2 == 0 && nFDirs2 == 0 ){
            acadoPrintf("BDF: Forward Sensitivities:\n");
            for( run1 = 0; run1 < m; run1++ ){
                acadoPrintf("%.16e  ", nablaG(0,run1) );
            }
            acadoPrintf("\n");
        }

        if( nFDirs2 > 0 ){

            acadoPrintf("BDF: First Order Forward Sensitivities:\n");
            for( run1 = 0; run1 < m; run1++ ){
                acadoPrintf("%.16e  ", nablaG2(0,run1) );
            }
            acadoPrintf("\n");

            acadoPrintf("BDF: Second Order Forward Sensitivities:\n");
            for( run1 = 0; run1 < m; run1++ ){
                acadoPrintf("%.16e  ", nablaG3(0,run1) );
            }
            acadoPrintf("\n");
        }

        // Backward Sensitivities:
        // -----------------------

        if( nBDirs > 0 ){
            acadoPrintf("BDF: t = %.16e  h = %.16e\n", t-c[6]*h[0], h[0] );
            acadoPrintf("BDF: Backward Sensitivities:\n");
            printBDFfinalResults2( nablaH );
        }

        // 2nd order Backward Sensitivities:
        // ---------------------------------

        if( nBDirs2 > 0 ){
            acadoPrintf("BDF: t = %.16e  h = %.16e\n", t-c[6]*h[0], h[0] );
            acadoPrintf("BDF: Backward Sensitivities:\n");
            printBDFfinalResults2( nablaH2 );
            acadoPrintf("BDF: 2nd order Backward Sensitivities:\n");
            printBDFfinalResults2( nablaH3 );
        }
}



void AcadoIntegratorBackend::printBDFIntermediateResults(){

    int run1;

     if( PrintLevel == HIGH ){

       if( nBDirs == 0 && nBDirs2 == 0 )
         acadoPrintf("BDF: t = %.16e  h = %.16e  ", t, h[0] );

       if( soa != SOA_EVERYTHING_FROZEN ){

         for( run1 = 0; run1 < md; run1++ ){
             acadoPrintf("x[%d] = %.16e  ", run1, nablaY(0,run1) );
         }
         for( run1 = 0; run1 < ma; run1++ ){
             acadoPrintf("xa[%d] = %.16e  ", run1, nablaY(0,md+run1) );
         }
         acadoPrintf("\n");
       }
       else{
         if( nBDirs == 0 && nBDirs2 == 0 )
            acadoPrintf("\n");
       }

        // Forward Sensitivities:
        // ----------------------

        if( nFDirs > 0 && nBDirs2 == 0 && nFDirs2 == 0 ){
            acadoPrintf("BDF: Forward Sensitivities:\n");
            for( run1 = 0; run1 < m; run1++ ){
                acadoPrintf("%.16e  ", nablaG(0,run1) );
            }
            acadoPrintf("\n");
        }

        if( nFDirs2 > 0 ){

            acadoPrintf("BDF: First Order Forward Sensitivities:\n");
            for( run1 = 0; run1 < m; run1++ ){
                acadoPrintf("%.16e  ", nablaG2(0,run1) );
            }
            acadoPrintf("\n");

            acadoPrintf("BDF: Second Order Forward Sensitivities:\n");
            for( run1 = 0; run1 < m; run1++ ){
                acadoPrintf("%.16e  ", nablaG3(0,run1) );
            }
            acadoPrintf("\n");
        }

     }
}


void AcadoIntegratorBackend::printRKIntermediateResults(){

    int run1, run3;

    if( PrintLevel == HIGH ){

        for( run3 = 0; run3 < 4; run3++ ){
            if( nBDirs == 0 && nBDirs2 == 0 )
            acadoPrintf("BDF: t = %.16e  h = %.16e  ", t+h[0]*c[3+run3], h[0]*c[3] );
            if( soa != SOA_EVERYTHING_FROZEN ){
                for( run1 = 0; run1 < md; run1++ ){
                    acadoPrintf("x[%d] = %.16e  ", run1, nablaY(3-run3,run1) );
                }
                for( run1 = 0; run1 < ma; run1++ ){
                    acadoPrintf("xa[%d] = %.16e  ", run1, nablaY(3-run3,run1+md) );
                }
                acadoPrintf("\n");
            }
            else{
                if( nBDirs == 0 && nBDirs2 == 0 )
                    acadoPrintf("\n");
            }

            // Forward Sensitivities:
            // ----------------------

            if( nFDirs > 0 && nBDirs2 == 0 && nFDirs2 == 0 ){
                acadoPrintf("BDF: Forward Sensitivities:\n");
                for( run1 = 0; run1 < m; run1++ ){
                    acadoPrintf("%.16e  ", nablaG(3-run3,run1) );
                }
                acadoPrintf("\n");
            }

            if( nFDirs2 > 0 ){

                acadoPrintf("BDF: First Order Forward Sensitivities:\n");
                for( run1 = 0; run1 < m; run1++ ){
                    acadoPrintf("%.16e  ", nablaG2(3-run3,run1) );
                }
                acadoPrintf("\n");

                acadoPrintf("BDF: Second Order Forward Sensitivities:\n");
                for( run1 = 0; run1 < m; run1++ ){
                    acadoPrintf("%.16e  ", nablaG3(3-run3,run1) );
                }
                acadoPrintf("\n");
            }
        }

        // Backward Sensitivities:
        // -----------------------

        if( nBDirs > 0 ){
            acadoPrintf("BDF: t = %.16e  h_RK = %.16e\n", t-c[6]*h[0], h[0] );
            acadoPrintf("BDF: Backward Sensitivities:\n");
            printBDFfinalResults2( nablaH );
        }

        // 2nd order Backward Sensitivities:
        // ---------------------------------

        if( nBDirs2 > 0 ){
            acadoPrintf("BDF: t = %.16e  h_RK = %.16e\n", t-c[6]*h[0], h[0] );
            acadoPrintf("BDF: Backward Sensitivities:\n");
            printBDFfinalResults2( nablaH2 );
            acadoPrintf("BDF: 2nd order Backward Sensitivities:\n");
            printBDFfinalResults2( nablaH3 );
        }
    }
}


returnValue AcadoIntegratorBackend::decomposeJacobian( Matrix &J ) const{

    switch( las ){

        case HOUSEHOLDER_METHOD:
             return J.computeQRdecomposition();

        case SPARSE_LU:
             return J.computeSparseLUdecomposition();

        default:
             return ACADOERROR( RET_NOT_IMPLEMENTED_YET );
    }

    return SUCCESSFUL_RETURN;
}


double AcadoIntegratorBackend::applyNewtonStep( double *etakplus1, const double *etak, const Matrix &J, const double *FFF ){

    int run1;
    Vector bb(m,FFF);
    Vector deltaX;

    switch( las ){

        case      HOUSEHOLDER_METHOD:  deltaX = J.solveQR      ( bb ); break;
        case      SPARSE_LU:           deltaX = J.solveSparseLU( bb ); break;
        default:                       deltaX.setZero          (    ); break;        
    }

    for( run1 = 0; run1 < m; run1++ )
        etakplus1[run1] = etak[run1] - deltaX(run1);

    return deltaX.getNorm( VN_LINF, diff_scale );
}


void AcadoIntegratorBackend::applyMTranspose( double *seed1, Matrix &J, double *seed2 ){

    int run1;
    Vector bb(m);

    for( run1 = 0; run1 < m; run1++ )
        bb(run1) = seed1[diff_index[run1]];

    Vector deltaX;

    switch( las ){

        case      HOUSEHOLDER_METHOD:  deltaX = J.solveTransposeQR      ( bb ); break;
        case      SPARSE_LU:           deltaX = J.solveTransposeSparseLU( bb ); break;
        default:                       deltaX.setZero                   (    ); break;
    }

    for( run1 = 0; run1 < m; run1++ )
        seed2[run1] = deltaX(run1);
}



void AcadoIntegratorBackend::relaxAlgebraic( double *residuum, double timePoint ){

    int           relaxationType  ;
    double        relaxationPar   ;

    get( RELAXATION_PARAMETER, relaxationPar );

    const double  a = relaxationPar*(timeInterval.getIntervalLength());
    const double  b = 1.0;
    double        damping = 1.0;
    int           run1   ;
    double        normRES;

    get( ALGEBRAIC_RELAXATION, relaxationType );

    switch( relaxationType ){

        case ART_EXPONENTIAL:
             damping = exp(relaxationConstant*(timeInterval.getFirstTime()-timePoint)/( timeInterval.getIntervalLength() ));
             for( run1 = 0; run1 < ma; run1++ )
                 residuum[md+run1] -= initialAlgebraicResiduum[md+run1]*damping;
             break;

        case ART_ADAPTIVE_POLYNOMIAL:
             normRES = 0.0;
             for( run1 = 0; run1 < ma; run1++ )
                 if( fabs(initialAlgebraicResiduum[md+run1]) > normRES )
                     normRES = fabs(initialAlgebraicResiduum[md+run1]);

             damping = pow( 1.0 - ( timePoint-timeInterval.getFirstTime())/( timeInterval.getIntervalLength() ) ,
                                    b + a/(normRES+100.0*EPS) );

             break;

        case ART_UNKNOWN:
             damping = 1.0;
             break;
    }

    for( run1 = 0; run1 < ma; run1++ )
        residuum[md+run1] -= initialAlgebraicResiduum[md+run1]*damping;

    return;
}



int AcadoIntegratorBackend::getDim() const{

    return md+ma;
}

int AcadoIntegratorBackend::getDimX() const{

    return md;
}


void AcadoIntegratorBackend::prepareDividedDifferences( Matrix &div ){

    int run1;

    for( run1 = 0; run1 < m; run1++ ){

        div(3,run1) = div(2,run1) - div(3,run1);

        div(2,run1) = div(1,run1) - div(2,run1);
        div(3,run1) = ( div(2,run1) - div(3,run1) )/(2.0);

        div(1,run1) = div(0,run1) - div(1,run1);
        div(2,run1) = ( div(1,run1) - div(2,run1) )/(2.0);
        div(3,run1) = ( div(2,run1) - div(3,run1) )/(3.0);

        div(4,run1) = 0.0; // the RK-starter has order 3 only.
    }
    return;
}


void AcadoIntegratorBackend::copyBackward( Vector       &Dx_x0,
                                  Vector       &Dx_p ,
                                  Vector       &Dx_u ,
                                  Vector       &Dx_w ,
                                  const Matrix &div    ) const{

    int run1;

    if( Dx_x0.getDim() != 0 )
        for( run1 = 0; run1 < m; run1++ )
            Dx_x0(run1)= div(0,diff_index[run1]);

    if( Dx_u.getDim() != 0 )
        for( run1 = 0; run1 < mu; run1++ )
            Dx_u(run1) = div(0,control_index[run1]);

    if( Dx_p.getDim() != 0 )
        for( run1 = 0; run1 < mp; run1++ )
            Dx_p(run1) = div(0,parameter_index[run1]);

    if( Dx_w.getDim() != 0 )
        for( run1 = 0; run1 < mw; run1++ )
            Dx_w(run1) = div(0,disturbance_index[run1]);

    return;
}


void AcadoIntegratorBackend::interpolate( int number_, Matrix &div, VariablesGrid &poly ){

    int i1 = timeInterval.getFloorIndex( t-h[0] );
    int i2 = timeInterval.getFloorIndex( t      );
    int jj, run1;

    for( jj = i1+1; jj <= i2; jj++ ){

        for( run1 = 0; run1 < m; run1++ )
            poly( jj, run1 ) = div(0,run1);

        double pp1 = (timeInterval.getTime(jj) - t)/h[0];
        double pp2 = pp1*(pp1*h[0] + h[0])/h[0];
        double pp3 = pp2*(pp1*h[0] + psi[number_][1])/h[0];
        double pp4 = pp3*(pp1*h[0] + psi[number_][2])/h[0];

        for( run1 = 0; run1 < m; run1++ )
            poly( jj, run1 ) += pp1*div(1,run1);

        for( run1 = 0; run1 < m; run1++ )
            poly( jj, run1 ) += pp2*div(2,run1);

        for( run1 = 0; run1 < m; run1++ )
            poly( jj, run1 ) += pp3*div(3,run1);

        for( run1 = 0; run1 < m; run1++ )
            poly( jj, run1 ) += pp4*div(4,run1);
    }
}


void AcadoIntegratorBackend::logCurrentIntegratorStep(	const Vector& currentX,
												const Vector& currentXA
												)
{
	return;
	
	int run1;

	// log differential states
	if ( currentX.isEmpty( ) == BT_TRUE )
	{
		Vector currentDiffStates(md);
		for( run1 = 0; run1 < md; run1++ )
			currentDiffStates( run1 ) = nablaY( 0,run1 );
		setLast( LOG_DIFFERENTIAL_STATES,currentDiffStates,t );
	}
	else
	{
		setLast( LOG_DIFFERENTIAL_STATES,currentX,t );
	}

	// log algebraic states
	if ( currentX.isEmpty( ) == BT_TRUE )
	{
		Vector currentAlgStates(ma);
		for( run1 = 0; run1 < ma; run1++ )
			currentAlgStates( run1 ) = nablaY( 0,md+run1 );
		setLast( LOG_ALGEBRAIC_STATES,currentAlgStates,t );
	}
	else
	{
		setLast( LOG_ALGEBRAIC_STATES,currentXA,t );
	}
}

returnValue AcadoIntegratorBackend::init( const DifferentialEquation &rhs_,
                                        const Transition           &trs_ ){

  assert(0);
    return Integrator::init( rhs_, trs_ );
}


} // namespace ACADO

// end of file.
