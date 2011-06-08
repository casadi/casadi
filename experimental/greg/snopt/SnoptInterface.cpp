// SnoptInterface.cpp
// Greg Horn
// Casadi 2011

#include <stdio.h>
#include <string.h>
#include <iostream>

#include <cstdlib>

#include "SnoptInterface.hpp"


static SnoptInterface * si;

using namespace std;
using namespace CasADi;

SnoptInterface::~SnoptInterface()
{
	delete []iGfun;
	delete []jGvar;

	delete []iAfun;
	delete []jAvar;
	delete []A;

	delete []x;
	delete []xlow;
	delete []xupp;
	delete []xmul;
	delete []xstate;

	delete []F;
	delete []Flow;
	delete []Fupp;
	delete []Fmul;
	delete []Fstate;
}

SnoptInterface::SnoptInterface(const SXFunction& user_F) : F_(user_F)
{
	F_.init();

	G_ = F_.jacobian();
	G_.init();

	si = this;

	init();
}


void
SnoptInterface::init()
{
	// design variables
	n = F_.input().size();
	x = new doublereal[n];
	xlow = new doublereal[n];
	xupp = new doublereal[n];
	xmul = new doublereal[n];
	xstate = new integer[n];
	for (int k=0; k<n; k++){
		x[k] = 0;
		xlow[k] = -SNOPT_INFINITY;
		xupp[k] = SNOPT_INFINITY;
		xmul[k] = 0;
		xstate[k] = 0;
	}

	// objective/constraint functions
	neF = F_.output().size();
	objAdd = 0;
	objRow = FIRST_FORTRAN_INDEX;
	F = new doublereal[neF];
	Flow = new doublereal[neF];
	Fupp = new doublereal[neF];
	Fmul = new doublereal[neF];
	Fstate = new integer[neF];
	for (int k=0; k<neF; k++){
		F[k] = 0;
		Flow[k] = -SNOPT_INFINITY;
		Fupp[k] = 0;
		Fmul[k] = 0;
		Fstate[k] = 0;
	}
	Fupp[0] = SNOPT_INFINITY;


	// linear jacobian
	lenA = 0; // implement linearity auto-detection
	if (lenA == 0)
		lenA = 1;
	iAfun = new integer[lenA];
	jAvar = new integer[lenA];
	A = new doublereal[lenA];
	neA = 0;

	// nonlinear jacobian
	vector<int> rowind,col;
	G_.output().sparsity().getSparsityCRS(rowind,col);

	lenG = col.size();
	iGfun = new integer[lenG];
	jGvar = new integer[lenG];

	// populate nonlinear jacobian sparsity pattern
	neG = 0;
	for(int r=0; r<rowind.size()-1; ++r)
        for(int el=rowind[r]; el<rowind[r+1]; ++el){
			iGfun[neG] = r + FIRST_FORTRAN_INDEX;
			jGvar[neG] = col[el] + FIRST_FORTRAN_INDEX;
            neG++;
        }
}

void
SnoptInterface::run()
{
	cerr << "Initializing SnoptInterface" << endl;

	fprintf(stderr,"Initializing SnoptInterface...\n");
	fflush(stderr);
 #define LENRW 20000
 #define LENIW 10000
 #define LENCW 500

//#define LENRW 2000000
//#define LENIW 100000
//#define LENCW 5000

	integer    minrw, miniw, mincw;
	integer    lenrw = LENRW, leniw = LENIW, lencw = LENCW;
	doublereal rw[LENRW];
	integer    iw[LENIW];
	char       cw[8*LENCW];

	integer    Cold = 0, Basis = 1, Warm = 2;

	integer    INFO;

	integer    nxname = 1, nFname = 1, npname;
	char       xnames[1*8], Fnames[1*8];
	char       Prob[200];

	integer    iSpecs = 4,  spec_len;
	integer    iSumm  = 6;
	integer    iPrint = 9,  prnt_len;

	char       printname[200];
	char       specname[200];

	integer    nS, nInf;
	doublereal sInf;
	integer    DerOpt, Major, iSum, iPrt, strOpt_len;
	char       strOpt[200];

	/* open output files using snfilewrappers.[ch] */
	sprintf(specname ,   "%s", "sntoya.spc");   spec_len = strlen(specname);
	sprintf(printname,   "%s", "sntoya.out");   prnt_len = strlen(printname);

	/* Open the print file, fortran style */
	snopenappend_
		( &iPrint, printname,   &INFO, prnt_len );

	/*     ================================================================== */
	/*     First,  sninit_ MUST be called to initialize optional parameters   */
	/*     to their default values.                                           */
	/*     ================================================================== */

	sninit_
		( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );



	/* Read in specs file (optional) */
	/* snfilewrapper_ will open the specs file, fortran style, */
	/* then call snspec_ to read in specs.                        */

	// snfilewrapper_
	// 	( specname, &iSpecs, &INFO, cw, &lencw,
	// 	  iw, &leniw, rw, &lenrw, spec_len, 8*lencw);

	// if( INFO != 101 )
    // {
	// 	printf("Warning: trouble reading specs file %s \n", specname);
    // }



	/*     ------------------------------------------------------------------ */
	/*     Tell SnoptA that userfg computes derivatives.                      */
	/*     The parameters iPrt and iSum may refer to the Print and Summary    */
	/*     file respectively.  Setting them to 0 suppresses printing.         */
	/*     ------------------------------------------------------------------ */

	DerOpt = 1;
	iPrt   = 0;
	iSum   = 0;
	sprintf(strOpt,"%s","Derivative option");
	strOpt_len = strlen(strOpt);
	snseti_
		( strOpt, &DerOpt, &iPrt, &iSum, &INFO,
		  cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8*500 );

	strcpy(Prob,"mysnopta");
	INFO = 0;

	Major = 250;
	strcpy( strOpt,"Major Iteration limit");
	strOpt_len = strlen(strOpt);
	snseti_
		( strOpt, &Major, &iPrint, &iSumm, &INFO,
		  cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8*500 );

	/*     ------------------------------------------------------------------ */
	/*     Solve the problem                                                  */
	/*     ------------------------------------------------------------------ */
	snopta_
		( &Cold, &neF, &n, &nxname, &nFname,
		  &objAdd, &objRow, Prob, (U_fp)toyusrfg_,
		  iAfun, jAvar, &lenA, &neA, A,
		  iGfun, jGvar, &lenG, &neG,
		  xlow, xupp, xnames, Flow, Fupp, Fnames,
		  x, xstate, xmul, F, Fstate, Fmul,
		  &INFO, &mincw, &miniw, &minrw,
		  &nS, &nInf, &sInf,
		  cw, &lencw, iw, &leniw, rw, &lenrw,
		  cw, &lencw, iw, &leniw, rw, &lenrw,
		  npname, 8*nxname, 8*nFname,
		  8*500, 8*500);

	snclose_( &iPrint );
	snclose_( &iSpecs );
}


int SnoptInterface::toyusrfg_
( integer    *Status, integer *n,    doublereal x[],
  integer    *needF,  integer *neF,  doublereal F[],
  integer    *needG,  integer *neG,  doublereal G[],
  char       *cu,     integer *lencu,
  integer    iu[],    integer *leniu,
  doublereal ru[],    integer *lenru )
{
	if( *needF > 0 ) {
		si->F_.setInput(x);
		si->F_.evaluate();
		si->F_.getOutput(F);
	}

	if( *needG > 0 ){
		si->G_.setInput(x);
		si->G_.evaluate();
		si->G_.getOutput(G);
	}
	return 0;
}
