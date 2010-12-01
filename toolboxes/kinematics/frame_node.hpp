#ifndef FRAME_NODE_HPP
#define FRAME_NODE_HPP

#include "frame.hpp"

namespace KINEMATICS{
using namespace CasADi;

/** \brief Internal class to make Frame trees reference-safe

*/
class FrameNode{
  public:
    friend class Frame;
      
/** \brief  Constructor */
    explicit FrameNode(const std::string& name, const SXMatrix & q,const SXMatrix & dq,const SXMatrix & ddq);
    FrameNode(const std::string& name, const Frame &ref,const SXMatrix & T);
    ~FrameNode();

/** \brief  Print */
    void print(std::ostream &stream);

    
  protected:
  std::string name;
  int count;
  Frame ref; 
  SXMatrix q;
  SXMatrix dq;
  SXMatrix ddq;
  SXMatrix R; // the 3x3 rotation matrix
  SXMatrix p; // 3 vectors
};

} // namespace KINEMATICS

#endif //FRAME_NODE_HPP

