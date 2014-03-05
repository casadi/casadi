 
typedef void (*UserFun)(int * mode, int* nnObj, int * nnCon, int *nJac, int *nnL, int * neJac, double *x, double *fObj, double *gObj, double * fCon, double* gCon, int* nState, char* cu, int* lencu, int* iu, int* leniu, double* ru, int *lenru);

extern "C" {

 
  extern void snopt_init(const int * iPrint,const int * iSumm, char * cw, const int * lencw, int *iw, const  int *leniw, double * rw, const int * lenrw);
  
 typedef int integer;
 typedef double doublereal;
 typedef long ftnlen;
 typedef UserFun U_fp;
 
 extern int snoptc_(
const char *start
,integer *m
,integer *n
,integer *ne
,integer *nName
,integer *nnCon
,integer *nnObj
,integer *nnJac
,integer *iObj
,doublereal *ObjAdd
,const char *prob
,U_fp userfun
,const doublereal *jcol
,const integer* indJ
,const integer* locJ
,doublereal* bl
,doublereal* bu
,const char* names
,integer* hs
,doublereal* x
,doublereal* pi
,doublereal* rc
,integer *info
,integer *mincw
,integer *miniw
,integer *minrw
,integer *ns
,integer *ninf
,doublereal *sinf
,doublereal *obj

,char* cu
,integer *lencu
,integer *iu
,integer *leniu
,doublereal* ru
,integer *lenru
,char* cw
,integer *lencw
,integer *iw
,integer *leniw
,doublereal* rw
,integer *lenrw
,ftnlen start_len
,ftnlen prob_len
,ftnlen names_len
,ftnlen cu_len
,ftnlen cw_len
 );
 
 
extern int sngeti_
( char *buffer, integer *ivalue, integer *inform,
 char *cw, integer *lencw, integer *iw,
 integer *leniw, doublereal *rw, integer *lenrw,
 ftnlen buffer_len, ftnlen cw_len);

extern int snseti_
( const char *buffer, integer *ivalue, integer *iprint,
 integer *isumm, integer *inform, char *cw,
 integer *lencw, integer *iw, integer *leniw,
 doublereal *rw, integer *lenrw, ftnlen buffer_len,
 ftnlen cw_len );

extern int snsetr_
( const char *buffer, doublereal *rvalue, integer * iprint,
 integer *isumm, integer *inform, char *cw,
 integer *lencw, integer *iw, integer *leniw,
 doublereal *rw, integer *lenrw, ftnlen buffer_len,
 ftnlen cw_len );
 
 extern int sninit_
( integer *iPrint, integer *iSumm, char *cw,
  integer *lencw, integer *iw, integer *leniw,
  doublereal *rw, integer *lenrw, ftnlen cw_len );
  

  extern void snopt_c(
    const char * Start, const int * lenstart,const int * m, const int * n, const int * neA, const int * nName, const int *nnCon, const int *nnObj, const int *nnJac, const int *iObj, const double *ObjAdd, const char* Prob , UserFun userfun,
    
    const double* Acol,const int* indA,const int *locA, double* bl, double* u,
    
    char* Names,
    
    // Initial values
    int* hs, double* x, double* pi, double * rc,
    
    // Outputs
    int *INFO, int* mincw, int* miniw, int* minrw, int * nS, int* nInf,double* sInf, double* Obj,
    
    // Working spaces for usrfun
    char* cu, const int* lencu, int* iu, const int* leniu, double* ru, const int* lenru,
    // Working spaces for SNOPT
    char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);
    
  // cvalue is always length 8
  extern void snopt_getc (const char *buffer, const int* lenbuffer, char *cvalue, int* Errors, char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);
  extern void snopt_geti (const char *buffer, const int* lenbuffer, int *ivalue, int* Errors, char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);
  
  extern void snopt_getr (const char *buffer, const int* lenbuffer, double *rvalue, int* Errors, char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);

  extern void snopt_set (const char *buffer, const int* lenbuffer, int * iPrint, int * iSumm, int* Errors, char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);
  extern void snopt_seti (const char *buffer, const int* lenbuffer, const int *ivalue, int * iPrint, int * iSumm, int* Errors, char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);
 
  extern void snopt_setr (const char *buffer, const int* lenbuffer, const double *rvalue, int * iPrint, int * iSumm, int* Errors, char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);
  
  
  extern void snopt_spec (const int *iSpecs, int* INFO, char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);
  
  extern void snopt_memb (int *INFO, const int* m,const int* n,const int* neA,const int*  negCon, const int* nnCon, const int* nnJac,const int*  nnObj, int* mincw, int* miniw, int* minrw, char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);

}



