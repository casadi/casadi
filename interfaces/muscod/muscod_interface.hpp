#ifndef MUSCOD_INTERFACE_HPP
#define MUSCOD_INTERFACE_HPP

#include <casadi/options_functionality.hpp>

namespace CasADi{
  
  // Pointer to a user defined setup function
  typedef void (*muscodSetupFcn)();
  
  // Forward declaration
  class MuscodInternal;
  
  // Smart pointer class
  class MuscodInterface : public OptionsFunctionality{
    public:

      /** \brief  Default constructor */
      MuscodInterface();

      /** \brief  Create instance */
      explicit MuscodInterface(muscodSetupFcn setupFcn);
      
      /** \brief  Access functions and members of the node */
      MuscodInternal* operator->();

      /** \brief  Const access functions and members of the node */
      const MuscodInternal* operator->() const;
      
      /** \brief  Make sure that the pointer points towards a valid object */
      virtual bool checkNode()() const;
      
      /** \brief  Solve the problem  */
      void solve();

};



} // namespace CasADi

#endif //MUSCOD_INTERFACE_HPP
