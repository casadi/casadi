//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/*! \brief Algorithm for solving (nonconvex) quadratic programs
 *
 * \author Ashutosh Mahajan, Jeff Linderoth, Mahdi Namazifar, MINOTAUR Team
 */

#include <iomanip>
#include <iostream>

#include <MinotaurConfig.h>

#include "AMPLInterface.h"
#include "BranchAndBound.h"
#include "Constraint.h"
#include "EngineFactory.h"
#include "Environment.h"
#include "Handler.h"
#include "IntVarHandler.h"
#include "LinearHandler.h"
#include "LPEngine.h"
#include "LPProcessor.h"
#include "Logger.h"
#include "MaxVioBrancher.h"
#include "NLPEngine.h"
#include "NLPMultiStart.h"
#include "NlPresHandler.h"
#include "NodeIncRelaxer.h"
#include "Option.h"
#include "Problem.h"
#include "ProblemSize.h"
#include "Presolver.h"
#include "Relaxation.h"
#include "ReliabilityBrancher.h"
#include "Solution.h"
#include "Timer.h"
#include "TransPoly.h"
#include "TreeManager.h"

#undef DEBUG_POLLYMAIN

using namespace Minotaur;

void show_help()
{
  std::cout << "Usage:" << std::endl
            << "To show version: polly -v (or --show_version yes) " 
            << std::endl
            << "To show all options: polly -= (or --show_options yes)" 
            << std::endl
            << "To solve an instance: polly --option1 [value] "
            << "--option2 [value] ... " << " .nl-file" << std::endl;
}


// Add options specific to handling of multilinear
void add_options(OptionDBPtr options)
{
  DoubleOptionPtr d_option;
  IntOptionPtr i_option;  
  StringOptionPtr s_option;
  std::string str;

  // polly works only with AMPL
  options->findString("interface_type")->setValue("AMPL");

  // Now add your options
  str = "Grouping Strategy for Multlinear: ";
  str += "TT (Term by term), CONV (Convex Hull), TC (Term Cover)";
  s_option = (StringOptionPtr) new Option<std::string>("ml_group_strategy",
     str, true, "TC");
  options->insert(s_option);

  i_option = (IntOptionPtr) new Option<int>("ml_max_group_size",
     "Maximum size of individual element in grouping: >= 1, <= 20", true, 6);
  options->insert(i_option);

  d_option = (DoubleOptionPtr) new Option<double>(
             "ml_cover_augmentation_factor", 
             "Covering augmentation factor for ml grouping: >= 1", true, 2.0);
  options->insert(d_option);

  d_option = (DoubleOptionPtr) new Option<double>("ml_feastol", 
     "Feasibility Tolerance for handling multlinear terms: > 0", true, 1.0e-6);
  options->insert(d_option);
}


PresolverPtr createPres(EnvPtr env, ProblemPtr p, size_t ndefs, 
                        HandlerVector &handlers)
{
  // create handlers for presolve
  PresolverPtr pres = PresolverPtr(); // NULL
  const std::string me("glob: ");

  p->calculateSize();
  if (env->getOptions()->findBool("presolve")->getValue() == true) {
    LinearHandlerPtr lhandler = (LinearHandlerPtr) new LinearHandler(env, p);
    handlers.push_back(lhandler);
    if (p->isQP() || p->isQuadratic() || p->isLinear() ||
        true==env->getOptions()->findBool("use_native_cgraph")->getValue()) {
      lhandler->setPreOptPurgeVars(true);
      lhandler->setPreOptPurgeCons(true);
      lhandler->setPreOptCoeffImp(true);
    } else {
      lhandler->setPreOptPurgeVars(false);
      lhandler->setPreOptPurgeCons(false);
      lhandler->setPreOptCoeffImp(false);
    }
    if (ndefs>0) {
      lhandler->setPreOptDualFix(false);
    } else {
      lhandler->setPreOptDualFix(true);
    }

    if (!p->isLinear() && 
         true==env->getOptions()->findBool("use_native_cgraph")->getValue() && 
         true==env->getOptions()->findBool("nl_presolve")->getValue() 
         ) {
      NlPresHandlerPtr nlhand = (NlPresHandlerPtr) new NlPresHandler(env, p);
      handlers.push_back(nlhand);
    }

    // write the names.
    std::cout << me << "handlers used in presolve:" << std::endl;
    for (HandlerIterator h = handlers.begin(); h != handlers.end(); 
        ++h) {
      std::cout<<(*h)->getName()<<std::endl;
    }
  }

  pres = (PresolverPtr) new Presolver(p, env, handlers);
  return pres;
}


int main(int argc, char* argv[])
{
  EnvPtr      env         = (EnvPtr) new Environment();
  Timer       *timer      = env->getNewTimer();
  OptionDBPtr options     = env->getOptions();
  TransPoly*  transformer = 0;
  int         err         = 0;
  const std::string me("polly main: ");
  ProblemPtr orig_prob, ref_prob;
  HandlerVector handlers; 
  //VarVector *orig_v=0;
  PresolverPtr pres;


  MINOTAUR_AMPL::AMPLInterfacePtr iface = (MINOTAUR_AMPL::AMPLInterfacePtr) 
    new MINOTAUR_AMPL::AMPLInterface(env);

  timer->start();

  options->findBool("presolve")->setValue(true);
  options->findBool("nl_presolve")->setValue(true);
  //  JTL -- I don't want to presolve linear constraints.  Makes my lambdas weird
  //options->findBool("lin_presolve")->setValue(false);
  options->findBool("lin_presolve")->setValue(true);

  env->readOptions(argc, argv);
  
  // Serdar added all the ML options?
  //add_options(options);
  options->findBool("use_native_cgraph")->setValue(true);

  // Change it to best first
  options->findString("tree_search")->setValue("bfs");

  // Change it to group size
  //options->findInt("ml_max_group_size")->setValue(6);


  // check if user needs help.
  if (options->findBool("show_options")->getValue() ||
      options->findFlag("=")->getValue()) {
    options->write(std::cout);
    exit(0);
  }

  if (options->findBool("show_help")->getValue() ||
      options->findFlag("?")->getValue()) {
    show_help();
    exit(0);
  }

  if (options->findBool("show_version")->getValue() ||
      options->findFlag("v")->getValue()) {
    std::cout << "Minotaur version " << env->getVersion() << std::endl;
#if DEBUG
    std::cout << me; 
    env->writeFullVersion(std::cout);
    std::cout << std::endl;
#endif
    exit(0);
  }

  if (options->findString("problem_file")->getValue()=="") {
    show_help();
    exit(0);
  }

  // read MINLP from AMPL & create Hessian/Jacobian for NLP solves
  orig_prob = iface->readInstance(options->findString("problem_file")
                                      ->getValue());
  std::cout << "time used in reading instance = " << std::fixed 
    << std::setprecision(2) << timer->query() << std::endl;

  // get presolver.
  handlers.clear();
  //orig_v = new VarVector(orig_prob->varsBegin(), orig_prob->varsEnd());
  pres = createPres(env, orig_prob, iface->getNumDefs(), handlers);
  if (options->findBool("presolve")->getValue() == true) {
    std::cout << me << "Presolving ... " << std::endl;
    pres->solve();
    std::cout << me << "Finished presolving." << std::endl;
    for (HandlerVector::iterator it=handlers.begin(); it!=handlers.end(); ++it) {
      (*it)->writeStats(env->getLogger()->msgStream(LogExtraInfo));
    }
  }
  handlers.clear();

  // Transform the problem to the 'standard' form.
  transformer = new TransPoly(env, orig_prob);
  transformer->reformulate(ref_prob, handlers, err);
  assert(0==err);

#if defined(DEBUG_POLLYMAIN)
  std::cout << me << "reformulated problem: " << std::endl;
  ref_prob->write(std::cout);
#endif

  //JTL.  Just getting feet wet
  // (1) Print which handler is handline which constraints

#if defined(DEBUG_POLLYMAIN)
  for (HandlerIterator h = handlers.begin(); h != handlers.end(); 
       ++h) {
    std::cout<<(*h)->getName()<< " will handle constraints with Id: ";
    for (ConstraintVector::const_iterator it = (*h)->consBegin(); it != (*h)->consEnd(); ++it) {
      std::cout << (*it)->getId() << " ";
    }
    std::cout << std::endl;
  }
#endif

  //XXX Maybe transpoly should do this?.  Give it an int var handler if necessary
  ref_prob->calculateSize();
  UInt nbin = ref_prob->getSize()->bins;
  UInt nint = ref_prob->getSize()->ints;
  
  if (nbin + nint > 0) {
    IntVarHandlerPtr v_hand = (IntVarHandlerPtr) new IntVarHandler(env, ref_prob);
    handlers.push_back(v_hand);
  }
  
  //XXX Should set max group size based on max degree of polynomial
  

  // Make your Branch and Bound
  BranchAndBound *bab = new BranchAndBound(env, ref_prob);
  
  //Need engine
  EngineFactory *efac = new EngineFactory(env);
  EnginePtr e = efac->getLPEngine();

  // Need a node processor
  LPProcessorPtr nproc;
  nproc = (LPProcessorPtr) new LPProcessor(env, e, handlers);
  bab->setNodeProcessor(nproc);
  
  // Need a (for now) simple brancher.

  BrancherPtr br;

  if (env->getOptions()->findString("brancher")->getValue() == "maxvio") {  
    MaxVioBrancherPtr maxviol_br;
    maxviol_br    = (MaxVioBrancherPtr) new MaxVioBrancher(env, handlers);
    br = maxviol_br;
  }
  else if (env->getOptions()->findString("brancher")->getValue() == "rel") {
    ReliabilityBrancherPtr rel_br;
    rel_br    = (ReliabilityBrancherPtr) new ReliabilityBrancher(env, handlers);
    rel_br->setEngine(e);
    rel_br->setThresh(4);
    rel_br->setMaxDepth(100);
    rel_br->setIterLim(10000);
    std::cout << me << "reliability branching iteration limit = " 
               << rel_br->getIterLim() << std::endl;
    br = rel_br;
  }

  nproc->setBrancher(br);

  // Need a relaxer
  NodeIncRelaxerPtr nr = (NodeIncRelaxerPtr) new NodeIncRelaxer(env, handlers);
  nr->setEngine(e);
  RelaxationPtr rel;
  nr->setRelaxation(rel);
  bab->setNodeRelaxer(nr);
  bab->shouldCreateRoot(true);

  if (env->getOptions()->findBool("msheur")->getValue() == true) {
    std::cout << "performing NLP Multistart heuristic" << std::endl;

    // Trying to give it a heuristic
    EnginePtr nlp_e = efac->getNLPEngine();
    orig_prob->setNativeDer();
    NLPMSPtr ms_heur = (NLPMSPtr) new NLPMultiStart(env, orig_prob, nlp_e);
    bab->addPreRootHeur(ms_heur); 
  }

  bab->solve();

  std::cout << me << "status of branch-and-bound: " 
            << getSolveStatusString(bab->getStatus()) << std::endl;

  SolutionPtr sol, sol2;
  
  // write solution
  sol = bab->getSolution(); // presolved solution needs translation.
  if (sol) {

#if 0
    //Don't write to ampl now
    if (options->findFlag("AMPL")->getValue()) {
      iface->writeSolution(sol);
    }
#endif
    //sol->write(std::cout);    
  }


  std::cout << me << "nodes created in branch-and-bound = " << 
    bab->getTreeManager()->getSize() << std::endl;
  std::cout << me << std::fixed << std::setprecision(4) << 
    "best bound estimate from remaining nodes = " <<  bab->getLb() 
            << std::endl;
  std::cout << me << std::fixed << std::setprecision(4) << 
    "best solution value = " <<  bab->getUb() << std::endl;

  std::cout << me << "time used = " << std::fixed << std::setprecision(2) 
            << timer->query() << std::endl;
  
  bab->writeStats(env->getLogger()->msgStream(LogExtraInfo));
  e->writeStats(env->getLogger()->msgStream(LogExtraInfo));

  //XXX Clean up memory

  if (bab) delete (bab);
  if (efac) delete (efac);
  if (timer) delete (timer);
  if (iface) delete (iface);
}

// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
