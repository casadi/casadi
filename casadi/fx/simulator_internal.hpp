#ifndef SIMULATOR_INTERNAL_HPP
#define SIMULATOR_INTERNAL_HPP

#include "simulator.hpp"

namespace CasADi{

/** \brief Simulator data storage classs
  \author Joel Andersson 
  \date 2010
*/
class SimulatorInternal : public FXNode{
public:
  
  /** \brief  Constructor */
  SimulatorInternal(const Integrator& integrator, const FX& output_fcn, const std::vector<double>& grid);
  
  /** \brief  Destructor */
  virtual ~SimulatorInternal();
  
  /** \brief  initialize */
  virtual void init();

  /** \brief  Integrate */
  virtual void evaluate(int fsens_order, int asens_order);

  Integrator integrator_;
  FX output_fcn_;
  std::vector<double> grid_;
};
  
} // namespace CasADi

#endif // SIMULATOR_INTERNAL_HPP
