#ifndef SNOPTPROBLEM_H
#define SNOPTPROBLEM_H

#include "snopt.h"

/* File snoptProblem.hpp
 *   C++ interface for SNOPT
 *
 * 10 Jul 2014: First version (based on previous work by Josh Griffin).
 */

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void summaryOff();
void summaryOn();

class snoptProblem {
private:
  int probID;
  int iprint, isumm;
  void init2zero();

protected:
  snoptProblem();
  snoptProblem(const char*name);
  snoptProblem(const char*name, const char *prtfile);
  ~snoptProblem();

  char    Prob[30];

  int     inform;
  int     initCalled, memCalled, allocA, allocG;

  int     leniw, lenrw;
  double *rw;
  int    *iw;

  int     lenru, leniu;
  double *ru;
  int    *iu;

  isnLog  snLog;
  isnLog2 snLog2;
  isqLog  sqLog;
  isnSTOP snSTOP;

  int  errMsgExit      ( const char *var );
  void allocI          ( int leniw );
  void allocR          ( int lenrw );
  void reallocI        ( int leniw );
  void reallocR        ( int lenrw );

  virtual void setWorkspace () = 0;

public:
  virtual int solve   ( int starttype ) = 0;

  int getParameter    ( const char *stroptin, char *stroptout );
  int getIntParameter ( const char *stropt,   int    &opt );
  int getRealParameter( const char *stropt,   double &opt );
  int setParameter    ( const char *stroptin );
  int setIntParameter ( const char *stropt,   int     opt );
  int setRealParameter( const char *stropt,   double  opt );

  void setProbName    ( const char *Prob );

  int  setSpecsFile   ( const char *specname );
  void setPrintFile   ( const char *prtname );

  void setUserI       ( int    *iu, int leniu );
  void setUserR       ( double *ru, int lenru );

  void setUserspace   ( int    *iu, int leniu,
			double *ru, int lenru );

  void setLog         ( isnLog snLog, isnLog2 snLog2, isqLog sqLog );
  void setSTOP        ( isnSTOP snSTOP );
};

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class snoptProblemA : public snoptProblem {
private:
  void init2zero();

protected:
  int     jacComputed;
  int     n, neF;
  int     ObjRow;
  double  ObjAdd;

  double *x, *xlow, *xupp, *xmul;
  double *F, *Flow, *Fupp, *Fmul;

  int    *xstate, *Fstate;

  int     lenA, lenG, neA, neG;
  int    *iAfun, *jAvar, *iGfun, *jGvar;
  double *A;

  snFunA  usrfunA;

  int     userDataSet();

public:
  snoptProblemA();
  snoptProblemA( const char *name );
  snoptProblemA( const char *name, const char *prtfile );
  ~snoptProblemA();

  int  computeJac    ( int &neA, int &neG );
  int  solve         ( int starttype );
  void setWorkspace  ();

  void setProblemSize( int n, int neF );
  void setObjective  ( int ObjRow, double ObjAdd );

  void setA          ( int lenA, int neA, int *iAfun, int *jAvar, double *A );
  void setG          ( int lenG, int neG, int *iGfun, int *jGvar );

  void setX          ( double *x, double *xlow, double *xupp,
                       double *xmul, int *xstate );
  void setF          ( double *F, double *Flow, double *Fupp,
                       double *Fmul, int *Fstate );
  void setUserFun    ( snFunA usrfun );
};

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class snoptProblemC : public snoptProblem {
private:
  void init2zero();

  snFunC  usrfunC;

protected:
  int     m, n, ne, nnCon, nnObj, nnJac, negCon, iObj;
  double  ObjAdd;

  int    *hs;
  double *x, *bl, *bu, *pi, *rc;

  int    *indJ, *locJ;
  double *Jval;

  int     userDataSet();

public:
  snoptProblemC();
  snoptProblemC( const char*name );
  snoptProblemC( const char*name, const char *prtfile );
  ~snoptProblemC();

  int  solve         ( int starttype );
  void setWorkspace  ();

  void setProblemSize( int m, int n, int nnCon, int nnJac, int nnObj );
  void setObjective  ( int iObj, double ObjAdd );

  void setJ          ( int ne, double *Jval, int *indJ, int *locJ );

  void setX          ( double *bl, double *bu, double *x,
		       double *pi, double *rc, int *hs );
  void setUserFun    ( snFunC usrfun );
};

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class snoptProblemB : public snoptProblemC {
private:
  void init2zero();
  int  userDataSet();

protected:
  snObjB  funobj;
  snConB  funcon;

public:
  snoptProblemB();
  snoptProblemB( const char*name );
  snoptProblemB( const char*name, const char *prtfile );
  ~snoptProblemB();

  int  solve         ( int starttype );
  void setFuncon     ( snConB funcon );
  void setFunobj     ( snObjB funobj );
};

#endif /* SNOPTPROBLEM_H */
