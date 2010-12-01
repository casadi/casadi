#ifndef ELEMENT_HPP
#define ELEMENT_HPP

namespace CasADi{


/** Element class for Matrix 
  Element is the return type for operator() of the Matrix class.
  From several alternative solutions, this was the only one which didn't cause ambigousity.
  Suggestions for improvement are welcome.
  \author Joel Andersson 
  \date 2010
*/

template<typename M, typename T>
class Element{
  public:
    /// Constructor
    Element(M& mat, int i, int j);
    
    //@{
    /// Methods that modify a part of the parent obejct (A[i] = ?, A[i] += ?, etc.)
    M& operator=(const T &y);
    M& operator+=(const T &y);
    M& operator-=(const T &y);
    M& operator*=(const T &y);
    M& operator/=(const T &y);
    //@}

    /// Get a reference to the element (? = A[i], ? += A[i], etc.)
    operator T();
  
  private:
    M& mat_;
    int i_, j_;
};


// Implementation

template<typename M, typename T>
M& Element<M,T>::operator=(const T &y){
  mat_.getElementRef(i_,j_) = y;
  return mat_;
}

template<typename M, typename T>
M& Element<M,T>::operator+=(const T &y){
  mat_.getElementRef(i_,j_) += y;
  return mat_;
}

template<typename M, typename T>
M& Element<M,T>::operator-=(const T &y){
  mat_.getElementRef(i_,j_) -= y;
  return mat_;
}

template<typename M, typename T>
M& Element<M,T>::operator*=(const T &y){
  mat_.getElementRef(i_,j_) *= y;
  return mat_;
}

template<typename M, typename T>
M& Element<M,T>::operator/=(const T &y){
  mat_.getElementRef(i_,j_) /= y;
  return mat_;
}

template<typename M, typename T>
Element<M,T>::operator T(){
  return mat_.getElement(i_,j_);
}

template<typename M, typename T>
Element<M,T>::Element(M& mat, int i, int j) : mat_(mat), i_(i), j_(j){
}



} // namespace CasADi


#endif // ELEMENT_HPP
