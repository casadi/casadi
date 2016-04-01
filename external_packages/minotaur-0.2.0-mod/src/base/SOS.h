//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file SOS.h
 * \brief Declare data structures for an SOS-I and SOS-II constraint
 * \author Ashutosh Mahajan, IIT Bombay
 */

#ifndef MINOTAURSOS_H
#define MINOTAURSOS_H

#include "Types.h"
#include "Variable.h"

namespace Minotaur {


class SOS {
public:
  /// Default constructor
  SOS();

  SOS(int n, SOSType type, const double *weights, const VarVector &vars,
      int priority, int id, std::string name);

  /// Destroy
  virtual ~SOS();

  int getId() const;
  int getNz();
  SOSType getType();
  int getPriority() const;

  std::string getName() const;
  const double* getWeights();

  VariableConstIterator varsBegin() const;
  VariableConstIterator varsEnd() const;

private:
  /// Index within the problem.
  int id_;

  /// Number of elements in the SOS
  int n_;

  /// Priority of this SOS over others.
  int priority_;

  /// \brief Type of SOS
  SOSType type_;
      
  /// Values
  double * weights_;

  /// Vector of variables.
  VarVector vars_;

  /// Name of the sos-constraint
  std::string name_;

};
}
#endif

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
