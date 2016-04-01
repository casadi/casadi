//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file Environment.cpp
 * \brief Implement the Environment class
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "MinotaurConfig.h"
#include "Environment.h"
#include "Logger.h"
#include "Option.h"
#include "Timer.h"
#include "Version.h"

using namespace Minotaur;

const std::string Environment::me_ = "Environment: ";

Environment::Environment()
{
  logger_     = (LoggerPtr) new Logger();
  options_    = (OptionDBPtr) new OptionDB();
  timerFac_   = new TimerFactory();
  timer_      = timerFac_->getTimer();
  createDefaultOptions_();
}


Environment::~Environment()
{
  delete timer_;
  delete timerFac_;
}


void Environment::createDefaultOptions_()
{
  BoolOptionPtr b_option;
  IntOptionPtr i_option;
  DoubleOptionPtr d_option;
  StringOptionPtr s_option;

  // bool options
  b_option = (BoolOptionPtr) new Option<bool>("show_version", 
      "Write version number of Minotaur: <0/1>",
      true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("show_options", 
      "Write all options and their default values: <0/1>",
      true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("show_help", 
      "Show usage and exit: <0/1>",
      true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("use_internal_quad", 
      "Should quadratic functions be treated natively: <0/1>",
      true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("use_warmstart", 
      "Should warm starts be used to solve relaxations: <0/1>",
      true, true);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("presolve", 
      "Should presolve be used: <0/1>", true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("persp_ref", 
      "Should prespective reformulation be used: <0/1>", true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("quad_cone_ref", 
      "tighten quadratic constraints with linear binary variables <0/1>",
      true, true);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("display_problem", 
      "Should display problem before solving: <0/1>", true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("display_presolved_problem", 
      "Should display problem after presolving: <0/1>", true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("display_size", 
      "Should display problem size before solving: <0/1>", true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("display_presolved_size", 
      "Should display problem size after presolving: <0/1>", true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("guided_dive", 
      "Should perform guided dives after branching: <0/1>", true, true);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("solve", 
      "Should solve the problem: <0/1>", true, true);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("lin_presolve", 
      "Should presolve using linear handler: <0/1>", true, true);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("quad_presolve", 
      "Should presolve using quadratic handler: <0/1>", true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("ignore_SOS1", 
      "Ignore all SOS1 constraints: <0/1>", true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("ignore_SOS2", 
      "Ignore all SOS2 constraints: <0/1>", true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("nl_presolve", 
      "Should presolve using nonlinear presolver: <0/1>", true, true);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("lin_show_stats", 
      "Should show statistics of linear handler: <0/1>", true, 
      true);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("MSheur", 
      "Use multi-start heuristic for continuous nonlinear problem: <0/1>", 
      true, false);
  options_->insert(b_option);
  
  b_option = (BoolOptionPtr) new Option<bool>("FPump", 
      "Use feasibility pump heuristic for MINLP: <0/1>", true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("modify_rel_only", 
      "If true, apply all modifications to relaxation only  <0/1>",
      true, true);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("convexity_type",
      "Is the problem convex or quasiconvex: convexity <=> 0 and quasiconvexity <=> 1", 
      true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("qg_accpm",
      "Do analytic center cutting plane cuts in qg algorithm: <0/1>", 
      true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("partial_BB",
      "Fix the bounds partially in QGHandler: <0/1>",
      true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("expand_poly", 
      "Fully expand and save polynomials in objective and constraints: <0/1>",
      true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("expand_quad", 
      "Fully expand and save quadratic objective and constraints: <0/1>",
      true, true);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool> ("use_native_cgraph", 
     "If true, use Minotaur's computational graph to evaluate nonlinear functions and their derivatives. <0/1>", true, false);
  options_->insert(b_option);

  b_option = (BoolOptionPtr) new Option<bool>("msheur", 
      "Enable multi-start initial heuristic: <0/1>", true, false);
  options_->insert(b_option);

  // reset, so that we don't accidently add it again.
  b_option.reset();


  // int options
  i_option = (IntOptionPtr) new Option<int>("bnb_node_limit", 
      "Limit on number of nodes created in branch-and-bound: >0",
      true, 1000000000);
  options_->insert(i_option);

  i_option = (IntOptionPtr) new Option<int>("bnb_sol_limit", 
      "Limit on the number of solutions found: >0", true, 
      1000000000);
  options_->insert(i_option);

  i_option = (IntOptionPtr) new Option<int>("ampl_log_level", 
      "Verbosity of ampl interface: 0-6", true, LogInfo);
  options_->insert(i_option);

  i_option = (IntOptionPtr) new Option<int>("bqpd_log_level", 
      "Verbosity of bqpd engine: 0-6", true, LogInfo);
  options_->insert(i_option);

  i_option = (IntOptionPtr) new Option<int>("bqpd_ws_mode", 
      "Warm starting mode for bqpd: 0-6", true, 6);
  options_->insert(i_option);

  i_option = (IntOptionPtr) new Option<int>("br_log_level", 
      "Verbosity of brancher: 0-6", true, LogInfo);
  options_->insert(i_option);

  i_option = (IntOptionPtr) new Option<int>("log_level", 
      "Verbosity of the main solving process: 0-6", true, LogInfo);
  options_->insert(i_option);

  i_option = (IntOptionPtr) new Option<int>("divheur", 
      "Use diving heuristic for MINLP: <-1/0/1>", 
      true, -1);
  options_->insert(i_option);

  i_option = (IntOptionPtr) new Option<int>("engine_log_level", 
      "Verbosity of engine: 0-6", true, LogInfo);
  options_->insert(i_option);

  i_option = (IntOptionPtr) new Option<int>("handler_log_level", 
      "Verbosity of handlers: 0-6", true, LogInfo);
  options_->insert(i_option);

  i_option = (IntOptionPtr) new Option<int>("heur_log_level", 
      "Verbosity of Multi Start Heuristic: 0-6", true, LogInfo);
  options_->insert(i_option);
  
  // Serdar added these options for MultilinearTermsHandler class
  i_option = (IntOptionPtr) new Option<int>("ml_max_group_size",
       "Maximum size of individual element in grouping: >= 2, <= 20", true, 6);
  options_->insert(i_option);
  // Serdar ended.

  i_option = (IntOptionPtr) new Option<int>("node_processor_log_level", 
      "Verbosity of node processor: 0-6", true, LogInfo);
  options_->insert(i_option);

  i_option = (IntOptionPtr) new Option<int>("presolve_log_level", 
      "Verbosity of presolver: 0-6", true, LogInfo);
  options_->insert(i_option);

  i_option = (IntOptionPtr) new Option<int>("rand_seed", 
      "Seed to random number generator: >=0 (0 = time(NULL))", true, 0);
  options_->insert(i_option);

   i_option = (IntOptionPtr) new Option<int>("strbr_pivot_limit",
      "Limit on number of iterations allowed during strong branching: >0",
      true, 25);
  options_->insert(i_option);
 
  i_option = (IntOptionPtr) new Option<int>("trans_log_level", 
      "Verbosity of Transformer ", true, LogInfo);
  options_->insert(i_option);

  i_option.reset();

  
  // double options
  
  // Serdar added these options for MultilinearTermsHandler class
  d_option = (DoubleOptionPtr) new Option<double>("ml_feastol", 
      "MultilinearTermsHandler feasibility tolerance.", true, 0.00001);
  options_->insert(d_option);

  d_option = (DoubleOptionPtr) new Option<double>("bnb_time_limit", 
      "Limit on time in branch-and-bound in seconds: >0",
      true, 1e20);
  options_->insert(d_option);
  
  d_option = (DoubleOptionPtr) new Option<double>("bnb_log_interval", 
      "Display interval in seconds for branch-and-bound status: >0", true, 
      5.);
  options_->insert(d_option);

  d_option = (DoubleOptionPtr) new Option<double>("int_tol", 
      "Tolerance for checking integrality",
      true, 0.000001);
  options_->insert(d_option);

  // Serdar added these options for MultilinearTermsHandler class
  d_option = (DoubleOptionPtr) new Option<double>("ml_cover_augmentation_factor", 
      "Covering augmentation factor for ml grouping: >= 1", true, 2.0);
  options_->insert(d_option);
  // Serdar ended.

  d_option = (DoubleOptionPtr) new Option<double>("obj_cut_off", 
      "Nodes with objective value above obj_cut_off are assumed infeasible",
      true, INFINITY);
  options_->insert(d_option);

  d_option = (DoubleOptionPtr) new Option<double>("obj_gap_percent", 
      "Stop if the objective gap percent falls below this level",
      true, 0.0);
  options_->insert(d_option);

  d_option.reset();
 
  // string options
  s_option = (StringOptionPtr) new Option<std::string>("brancher", 
      "Name of brancher: rel, maxvio, lex, rand, maxfreq", 
      true, "rel");
  options_->insert(s_option);

  s_option = (StringOptionPtr) new Option<std::string>("config_file", 
      "Name of file that contains parameters or options", true, "");
  options_->insert(s_option);

  s_option = (StringOptionPtr) new Option<std::string>("interface_type", 
      "What interface is this environment being used with: AMPL or C++",
      true, "C++");
  options_->insert(s_option);

  s_option = (StringOptionPtr) new Option<std::string>("lp_engine", 
      "Engine for solving Linear Relxations: Osi, None", 
      true, "OsiClp");
  options_->insert(s_option);

  // Serdar added default options for MultilinearTermsHandler.
  s_option = (StringOptionPtr) new Option<std::string>("ml_group_strategy",
      "Group strategy", true, "TC");
  options_->insert(s_option);
  // Serdar ended

  s_option = (StringOptionPtr) new Option<std::string>("nlp_engine", 
      "Engine for solving nonlinear relaxations: Filter-SQP, IPOPT, None", 
      true, "Filter-SQP");
  options_->insert(s_option);

  s_option = (StringOptionPtr) new Option<std::string>("problem_file", 
      "Name of file that contains the instance to be solved", 
      true, "");
  options_->insert(s_option);

  s_option = (StringOptionPtr) new Option<std::string>("qp_engine", 
      "Engine for solving QP relaxations: bqpd, None", 
      true, "bqpd");
  options_->insert(s_option);

  s_option = (StringOptionPtr) new Option<std::string>("tree_search", 
      "Strategy for tree search: dfs, bfs, BthenD", true, "BthenD");
  options_->insert(s_option);

  s_option = (StringOptionPtr) new Option<std::string>("vbc_file", 
      "File name for storing tree information for Vbctool", true, "");
  options_->insert(s_option);

  s_option.reset();
}


void Environment::convertAndAddOneOption_(BoolOptionPtr &b_option, 
                                          IntOptionPtr &i_option,
                                          DoubleOptionPtr &d_option, 
                                          StringOptionPtr &s_option,
                                          FlagOptionPtr &f_option,
                                          std::string &name,
                                          std::string &value,
                                          std::ostringstream &logstr)
{
  std::stringstream mystream; // to convert string to int, bool, double.
  std::string off = "  ";
  // now add the option.
  if (b_option) {
    bool b_value = getBoolValue_(value);
    b_option->setValue(b_value);
    logstr << off << b_option->getName() << " = " <<
      b_option->getValue() << std::endl; 
  } else if (i_option) {
    int i_value;
    mystream << std::flush;
    mystream << value;
    mystream >> i_value;
    i_option->setValue(i_value);
    logstr << off << i_option->getName() << " = " <<
      i_option->getValue() << std::endl; 
  } else if (d_option) {
    double d_value;
    mystream << std::flush;
    mystream << value;
    mystream >> d_value;
    d_option->setValue(d_value);
    logstr << off << d_option->getName() << " = " <<
      d_option->getValue() << std::endl; 
  } else if (s_option) {
    s_option->setValue(value);
    logstr << off << s_option->getName() << " = " <<
      s_option->getValue() << std::endl; 
  } else {
    // option is unknown. We will add this option with the following
    // assumption. If the string has an '=' in it (i.e. value!=""), 
    // then it is a string
    // option, otherwise it is a 'flag' option. We can not have 'int' or
    // 'double' option here.
    logstr << off << name << " is not a known option. ";
    if (value == "") {
      logstr << "Added as a flag." << std::endl;
      f_option = (FlagOptionPtr) new Option<bool>(name, 
          "Flag added by user", false, true);
      options_->insert(f_option, true);
    } else {
      logstr << "Added as a string option with value \"" 
        << value << "\"" << std::endl;
      s_option = (StringOptionPtr) new Option<std::string>(name, 
          "Option added by user", false, value);
      options_->insert(s_option);
    }
  }
}


void Environment::findOption_(const std::string &name,
                              BoolOptionPtr &b_option,
                              IntOptionPtr &i_option,
                              DoubleOptionPtr &d_option,
                              StringOptionPtr &s_option,
                              FlagOptionPtr &f_option)
{
  b_option = BoolOptionPtr();   // NULL
  i_option = IntOptionPtr();    // NULL
  d_option = DoubleOptionPtr(); // NULL
  s_option = StringOptionPtr(); // NULL
  f_option = FlagOptionPtr();   // NULL

  b_option = options_->findBool(name);
  if (b_option) {
    return;
  } 

  i_option = options_->findInt(name);
  if (i_option) {
    return;
  } 

  d_option = options_->findDouble(name);
  if (d_option) {
    return;
  }

  s_option = options_->findString(name);
  if (s_option) {
    return;
  }

  f_option = options_->findFlag(name);
}


bool Environment::getBoolValue_(const std::string & str)
{
  if (str=="") {
    return false;
  }

  if (str=="Y" || str=="YES" || str == "y" || str=="yes" || str == "T" ||
     str=="TRUE" || str=="t" || str=="true" || str=="1") {
    return true;
  }

  return false;
}


LoggerPtr Environment::getLogger() const
{
  return logger_;
}


LogLevel Environment::getLogLevel() const
{
  return logger_->getMaxLevel();
}


Timer* Environment::getNewTimer() 
{
  return timerFac_->getTimer();
}


OptionDBPtr Environment::getOptions()
{
  return options_;
}


double Environment::getTime(int &err)
{
  if (timer_) {
    err = 0;
    return timer_->query();
  } 
#if SPEW
  logger_->msgStream(LogError) << me_ <<
    "timer queried before it is started." << std::endl;
#endif
  err = 1;
  return 0.0;
}


const Timer* Environment::getTimer()
{
  return timer_;
}


std::string Environment::getVersion() 
{
  std::stringstream name_stream;
  name_stream << MINOTAUR_MAJOR_VERSION <<  '.' <<  MINOTAUR_MINOR_VERSION
              << " patch " <<  MINOTAUR_PATCH_VERSION;
  return name_stream.str();
}


void Environment::readConfigFile_(std::string fname, UInt &num_p)
{
  std::string line, w1, w2;
  std::ifstream istr(fname.c_str());
  const std::string spc = " \t\r\n";
  char comment = '#';
  BoolOptionPtr b_option;
  IntOptionPtr i_option;
  DoubleOptionPtr d_option;
  StringOptionPtr s_option;
  FlagOptionPtr f_option;
  std::ostringstream logstr;

  if (!istr.is_open()) {
    logger_->errStream() << me_ << "cannot open file " << fname << std::endl; 
    return;
  }

  while (istr.good()) {
    getline(istr, line);
    if (line.empty()) {
      continue;
    }
    // remove everything after comment.
    line = line.substr(0, line.find_first_of(comment));

    // find first word. It is parameter name.
    w1 = line.substr(0, line.find_first_of(spc));

    // remove first word from line. Get second word. It is the value.
    line = line.substr(w1.length()); 
    w2 = line.substr(0,line.find_first_not_of(spc)); 
    if (w2.length()==line.length()) {
      w2="";
    } else {
      w2 = line.substr(w2.length());
    }
    w2 = w2.substr(0,w2.find_first_of(spc));

    if (w1.empty()) {
      continue;
    } 
    if (w2.empty()) {
      logger_->errStream() << me_ << "parameter " << w1 << " in configuration file " 
        << fname << " has no value. Ignored this parameter." << std::endl;
      continue;
    }
    if ("problem_file"==w1) {
      ++num_p;
    }
    // add the option.
    findOption_(w1, b_option, i_option, d_option, s_option, f_option);
    convertAndAddOneOption_(b_option, i_option, d_option, s_option,
        f_option, w1, w2, logstr);
  }
  istr.close();
}


void Environment::readOptions(int argc, char **argv)
{
  std::string name, s_value;
  int leading_dashes;
  UInt num_p = 0; // number of filenames provided for solve.
  BoolOptionPtr b_option;
  IntOptionPtr i_option;
  DoubleOptionPtr d_option;
  StringOptionPtr s_option;
  FlagOptionPtr f_option;
  std::string offset = "  ";
  std::ostringstream ostr;

  if (argc<2) {
    logger_->msgStream(LogInfo) << me_ 
      << "User provided no options." << std::endl;
  } else {
    ostr << me_ << "User provided options:" << std::endl;
  }
  for (int i=1; i<argc; ++i) {
    name = argv[i];

    // removing leading dashes
    leading_dashes = 0;
    leading_dashes = removeDashes_(name);
    assert(name.size()>0);

    if (leading_dashes==0) {
      // looks like an name of instance
      s_option = options_->findString("problem_file");
      s_option->setValue(name);
      ++num_p;
      ostr << offset << s_option->getName() << " = " <<
        s_option->getValue() << std::endl; 
    } else {
      // looks like an option. remove leading minotaur dot, if any.
      removeMinotaurPrepend_(name);

      // we have an option for which we need an argument also. First check
      // if the option contains an '=' sign followed by a string.
      s_value = separateEqualToArg_(name);

      // find if it is a recognized option
      findOption_(name, b_option, i_option, d_option, s_option, f_option);

      if (f_option) { // true means we have a known flag.
        if (s_value=="") { 
          f_option->setValue(true);
        } else {
          f_option->setValue(getBoolValue_(s_value));
        }
        ostr << offset << f_option->getName() << " = " <<
          f_option->getValue() << std::endl; 
      } else {
        if (s_value=="" && (b_option || i_option || d_option || s_option)) {
          // the next argv is the argument
          assert(i+1<argc);       // throw exception instead
          ++i;
          s_value = argv[i];
        }
        convertAndAddOneOption_(b_option, i_option, d_option, s_option,
            f_option, name, s_value, ostr);
        if (s_option && "config_file"==s_option->getName()) {
          readConfigFile_(s_option->getValue(), num_p);
        }
      }
    }
  }
  if (argc>1) {
    ostr << me_ << "End of user provided options." << std::endl << std::endl; 
  }


  // update the log level if set by the user
  logger_->setMaxLevel((LogLevel)getOptions()->findInt("log_level")
                       ->getValue());
  // display all the new options set.
  logger_->msgStream(LogInfo) << ostr.str();

  if (num_p>1) {
    logger_->msgStream(LogInfo) << me_ 
      << "more than one filename given as input."
      << std::endl
      << me_ << "Only file \"" 
      << options_->findString("problem_file")->getValue()
      << "\" will be read." << std::endl << std::endl;
  } 
#if SPEW
  if (num_p==0) {
    logger_->msgStream(LogInfo) << me_ 
      << "No filename provided as input." << std::endl;
  }
#endif
}


UInt Environment::removeDashes_(std::string &name)
{
  size_t first_occ;

  // exit if string is NULL.
  if (name.length()<1) {
    return 0;
  }

  // check for '--'
  first_occ = name.find("--");
  if (first_occ==0) {
    name = name.substr(2); // todo catch exception out_of_range.
    return 2;
  } 

  // check for '-'
  first_occ = name.find("-");
  if (first_occ==0) {
    name = name.substr(1); // todo catch exception out_of_range.
    return 1;
  }

  return 0;
}


/// \todo implement me.
void Environment::removeMinotaurPrepend_(std::string & )
{
  // std::cout << fullname << std::endl;
}


std::string Environment::separateEqualToArg_(std::string &name)
{
  if (name.size()>1) {
    // search for the '=' sign from position 1 only (can't be at 0).
    size_t pos = name.find('=', 1);
    if (pos == std::string::npos) {
      // '=' not found
      return "";
    } else {
      std::string value = "";
      value = name.substr(pos+1);
      name = name.substr(0, pos);
      return value;
    }
  }
  return "";
}


void Environment::setLogLevel(LogLevel l) 
{
  logger_->setMaxLevel(l);
}


void Environment::startTimer(int &err)
{
  if (timer_) {
    err = 0;
    timer_->start();
  } else {
    err = 1;
  }
}


void Environment::stopTimer(int &err)
{
  if (timer_) {
    err = 0;
    timer_->stop();
  } else {
    err = 1;
  }
}


void Environment::writeFullVersion(std::ostream &out)
{
  out << "Minotaur full version " << getVersion() << " svn version "
      << MINOTAUR_SVN_VERSION;
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
