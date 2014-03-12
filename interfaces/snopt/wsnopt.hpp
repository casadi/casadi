 
typedef void (*UserFun)(int * mode, int* nnObj, int * nnCon, int *nJac, int *nnL, int * neJac, double *x, double *fObj, double *gObj, double * fCon, double* gCon, int* nState, char* cu, int* lencu, int* iu, int* leniu, double* ru, int *lenru);

extern "C" {
 
  extern void snopt_init(const int * iPrint,const int * iSumm, char * cw, const int * lencw, int *iw, const  int *leniw, double * rw, const int * lenrw);

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

    // direct calls to the fortran library
  extern void sninit_(const int * iPrint, const int * iSumm,
                      char* cw, const int* lencw,
                      int* iw, const int* leniw,
                      double* rw, const int* lenrw,
                      const long cw_len8);

  extern void snseti_(const char *buffer, const int *ivalue, int * iPrint, int * iSumm, int* Errors,
                      char* cw, const int* lencw,
                      int* iw, const int* leniw,
                      double* rw, const int* lenrw,
                      const long buffer_ftn_len, const long cw_len8);

  extern void snsetr_(const char *buffer, const double *ivalue, int * iPrint, int * iSumm, int* Errors,
                      char* cw, const int* lencw,
                      int* iw, const int* leniw,
                      double* rw, const int* lenrw,
                      const long buffer_ftn_len, const long cw_len8);

  extern void snset_(const char *buffer, int * iPrint, int * iSumm, int* Errors,
                     char* cw, const int* lencw,
                     int* iw, const int* leniw,
                     double* rw, const int* lenrw,
                     const long buffer_ftn_len, const long cw_len8);


}
