#include "muscod_internal.hpp"


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <csignal>
#include <cmath>
extern "C"{
#include "lac_util.h"
#include "model.h"
#include "def_qps.h"
#include "mssqp.h"
#include "minlp.h"
#include "util.h"
#include "def_usrmod.h"
#include "def_ind.h"
}

#define PSFILE "mtest.ps"

using namespace std;
namespace CasADi{

int MuscodInternal::instance_counter = 0;

MuscodInternal::MuscodInternal(muscodSetupFcn setupFcn) : setupFcn_(setupFcn){
  // Increse the counter
  if(instance_counter>0) throw CasadiException("Only one instance of MuscodInternal is allowed, delete this before creating a new instance");
  instance_counter++;

  addOption("datfile",OT_STRING);
  addOption("resfile",OT_STRING,"muscod_results.txt");
  addOption("logfile",OT_STRING,"log.txt");
  addOption("restartfile",OT_STRING,"restart.bin");
  addOption("backupfile",OT_STRING,"restart.bin");
  addOption("ppath", OT_STRING, "PAR");
  addOption("acc",OT_REAL,1e-6); // accuracy of NLP solution
  addOption("tol",OT_REAL,1e-7); // integration tolerance
  addOption("itol",OT_REAL,1e-7); // initial integration tolerance
  addOption("nhtopy",OT_INTEGER,0); // # of homotopy steps
  addOption("sf", OT_INTEGER,0); // index of first SQP iteration using fixed discretization
  addOption("mf", OT_INTEGER,0); // maximum # of successive ``fixed'' iterations
  addOption("mi",OT_INTEGER,100); // maximum total # of SQP iterations
  addOption("levmar",OT_REAL,0.0); // Levenberg Marquardt regularization of hessian
  addOption("eflag",OT_INTEGER,0); // use two gradient evaluations per SQP iteration
  addOption("sflag",OT_INTEGER,0); // stop after each SQP iteration
  addOption("wflag",OT_INTEGER,0); // warm/cool start using final data of previous run
  addOption("rfac",OT_INTEGER,0.0);
  addOption("bflag",OT_INTEGER,-1);
  addOption("cflag",OT_INTEGER,0);
  addOption("name", OT_STRING, "muscod_problem");
  
}

MuscodInternal::~MuscodInternal(){
  instance_counter--;
}
    
void MuscodInternal::print(std::ostream& stream) const{
  stream << "Muscod interface";
}
    
void MuscodInternal::solve(){
  const double acc = getOption("acc").toDouble(); // accuracy of NLP solution
  const double tol = getOption("tol").toDouble();  // integration tolerance
  const double itol = getOption("itol").toDouble(); // initial integration tolerance
  const long nhtopy = getOption("nhtopy").toInt(); // # of homotopy steps
  const long frstart = getOption("sf").toInt(); // index of first SQP iteration using fixed discretization
  const long frmax   = getOption("mf").toInt(); // maximum # of successive ``fixed'' iterations
  const long itmax   = getOption("mi").toInt(); // maximum total # of SQP iterations
  const long plevel  = 0;  //
  const double levmar  = getOption("levmar").toDouble(); // Levenberg Marquardt regularization of hessian
  const long eflag   = getOption("eflag").toInt(); // use two gradient evaluations per SQP iteration
  const long sflag   = getOption("sflag").toInt(); // stop after each SQP iteration
  const long wflag   = getOption("wflag").toInt(); // warm/cool start using final data of previous run
  const double rfac    = getOption("rfac").toDouble();
  const long bflag   = getOption("bflag").toInt();
  const long cflag   = getOption("cflag").toInt();
  
  MInfo *modptr;
  signal(SIGINT, sigproc);
  set_lformat(24, 16, 10, 10);
  
  char* pname = const_cast<char*>(getOption("name").toString().c_str());
  set_pname(pname);

  // Log file
  char *logname = const_cast<char*>(getOption("logfile").toString().c_str());
  FILE *log = fopen(logname,"w");

  // What is this?
  if (bflag >= 0) { 
    minlp ( setupFcn_, pname, plevel, def_QPSOL, wflag, cflag, acc, rfac,
            nhtopy, itol, tol, frmax, frstart, itmax, eflag, sflag, bflag );
    return;
  }

  // Read dat-file
  char *datname = const_cast<char*>(getOption("datfile").toString().c_str());
  FILE *data = fopen(datname,"r");
  char *ppath = const_cast<char*>(getOption("ppath").toString().c_str());
  ini_MODEL(ppath, setupFcn_, data, log, plevel, &modptr );
  fclose(data);

  // Initialize the SQP algorithm
  char* restartname = const_cast<char*>(getOption("restartfile").toString().c_str());
  FILE *restart = NULL;
  if (wflag || cflag)
    restart = fopen(restartname,"rb");
  ini_MSSQP(ppath, restart, log, plevel, modptr, def_QPSOL);
  if (wflag || cflag)
    fclose(restart);

  // Solve
  char* backup_name = const_cast<char*>(getOption("backupfile").toString().c_str());
  char* results_name = const_cast<char*>(getOption("resfile").toString().c_str());
  FILE *backup  = fopen(backup_name,"wb");
  FILE *results = fopen(results_name,"w");
  char psfile[]= PSFILE;
  mssqp(
    (wflag) ? MS_WARM : MS_COLD, 
    acc, 
    rfac,
    levmar,
    nhtopy, 
    itol,
    tol,
    frmax,
    frstart,
    itmax, 
    eflag,
    sflag, 
#ifdef PSFILE
    psfile,
#endif
    log, 
    plevel, 
    backup, 
    10, 
    results
  );
  fclose(results);
  fclose(backup);
  
}



} // namespace CasADi

