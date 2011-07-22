#ifndef MUSCOD_INTERNAL_HPP
#define MUSCOD_INTERNAL_HPP

#include "muscod_interface.hpp"
#include "muscod_function.hpp"
namespace CasADi{
  
class MuscodInternal : public OptionsFunctionalityNode{
  friend class MuscodInterface;
  
  /** \brief  Constructor only accessable from the MuscodInterface pointer class */
  explicit MuscodInternal(muscodSetupFcn setupFcn);
  
  public:
    
    /** \brief  Destructor */
    virtual ~MuscodInternal();

    /** \brief  Clone function */
    virtual MuscodInternal* clone() const{ return new MuscodInternal(*this);}

    
    /** \brief  Solve the problem  */
    void solve();

    /** \brief  Solve the problem */
    virtual void print(std::ostream& stream) const;

    
    muscodSetupFcn setupFcn_;
    static int instance_counter;
};

} // namespace CasADi

#endif //MUSCOD_INTERNAL_HPP
