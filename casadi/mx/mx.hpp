#ifndef MX_HPP
#define MX_HPP

#include "../matrix_size.hpp"
#include "../shared_object.hpp"
#include "../sx/sx.hpp"
#include <vector>

namespace CasADi{

/** \brief  Forward declaration */
class MXNode;

/** \brief The main matrix class.
  Operations on MX are lazy on the matrix level.
  \author Joel Andersson 
  \date 2010
*/
class MX : public SharedObject{
public:
 
/** \brief  Default constructor */
  MX();

/** \brief  Construct a symbolic matrix (matrix variable)
Internally represented by SymbolicMatrix
*/
  explicit MX(const std::string& name, int n=1, int m=1);

/** \brief  Create constant */
  MX(double x);
  MX(const std::vector<double> &x);
  MX(const std::vector<double> &x, int n, int m=1, char order='R');
#ifdef WITH_CPP0X
  MX(std::initializer_list<double> list);
#endif

/** \brief  Destructor */
  ~MX();

/** \brief  Immutable matrix style access
   Get an element of the matrix (without changing this object)
 */  
 const MX operator()(int i, int j) const;
/** \brief  Immutable vector style access
   Get an element of the matrix (without changing this object)
 */  
 const MX operator[](int k) const;          

/** \brief  Element of the matrix which is allowed to change the object (vector style index)
  \author Joel Andersson 
  \date 2010
  Allows to get elements from an MX and changing them.
  Elements are scalars, use other techniques for having sliced access.
  \see MatrixElement	
*/
  class Element : public PrintableObject{
    public:
      /** \brief Vector style access constructor */
      Element(MX& mx, int k);
      
      /** \brief  Print */
      virtual void print(std::ostream &stream=std::cout) const;
      
      /** \brief  Automatic type conversion (? =  A[i], sqrt(A[i]) etc.) */
      operator MX() const;
      //@{
      /** \brief  Objects that modify a part of the parent obejct (A[i] = ?, A[i] += ?, etc.) */
      MX& operator=(const MX &y);
      MX& operator+=(const MX &y);
      MX& operator-=(const MX &y);
      MX& operator*=(const MX &y);
      MX& operator/=(const MX &y);
      //@}
    
    protected:
      MX& mx;
      int k;
  };

/** \brief  Mutable matrix style access
   Get an element of the matrix (possibly changing this object)
 */  
  Element operator()(int i, int j);
/** \brief  Mutable vector style access
   Get an element of the matrix (possibly changing this object)
 */  
  Element operator[](int k);

/** \brief  Get the size */
  const MatrixSize& size() const;
/** \brief  Get the number of elements */
  int numel() const;
/** \brief get the first dimension
For an n-by-m matrix, returns n
  */
  int size1() const;
/** \brief get the second dimension
For an n-by-m matrix, returns m
  */
  int size2() const;

//@{
/** \brief  Operators that changes the object */
  MX& operator+=(const MX &y);
  MX& operator-=(const MX &y);
  MX& operator*=(const MX &y);
  MX& operator/=(const MX &y);
//@}

//@{
/** \brief  Operators */
  friend MX operator+(const MX &x, const MX &y);
  friend MX operator-(const MX &x, const MX &y);
  friend MX operator*(const MX &x, const MX &y);
  friend MX operator/(const MX &x, const MX &y);
  MX operator-() const;
//@}

  /** \brief  Check if the matrix expression is empty */
  bool isEmpty() const;
  
/** \brief  Initialize the tree */
/** \brief    void init(); */

/** \brief  Lowlevel. Get a pointer to the node
A regular user should never use this.
*/
  MXNode* get();
/** \brief  Lowlevel. Get a pointer to the node
A regular user should never use this.
*/
  const MXNode* get() const;

//@{
/** \brief  Quick access a member of the node */
  MXNode* operator->();
  const MXNode* operator->() const;
//@}

//@{
  /** \brief  Create nodes by their ID */
  static MX binary(int op, const MX &x, const MX &y);
  static MX unary(int op, const MX &x);
  static MX scalar_matrix(int op, const MX &x, const MX &y);
  static MX matrix_scalar(int op, const MX &x, const MX &y);
  static MX matrix_matrix(int op, const MX &x, const MX &y);
//@}

  /** \brief  Matrix of all zeros */  
  static MX zeros(int nrow, int ncol);
  
  /** \brief  Matrix of all ones */  
  static MX ones(int nrow, int ncol);
  
  //! \brief delayed setting or getting an element
  MX getElement(int k) const;
  MX& setElement(const MX& el, int k);
  
  /** \brief  Create an all zero matrix with size sz */
  explicit MX(const MatrixSize& sz);

};

//@{
/** \brief  concatenate */
MX vertcat(const std::vector<MX>& comp);
MX horzcat(const std::vector<MX>& comp);
MX vertcat(const MX& a, const MX& b);
MX horzcat(const MX& a, const MX& b);
//@}

/** \brief  Take the 2-norm of a MX
Internally represented by Norm2
*/
MX norm_2(const MX &x);
/** \brief  Take the 1-norm of a MX
Internally represented by Norm1
*/
MX norm_1(const MX &x);
/** \brief  Take the infinity-norm of a MX
Internally represented by NormInf
*/
MX norm_inf(const MX &x);
/** \brief  Take the transpose of a MX 
Internally represented by Transpose
*/
MX trans(const MX &x); // transpose
/** \brief  Take the matrix product of 2 MX objects */
MX prod(const MX &x, const MX &y); // matrix product
/** \brief  Take the inner product of two vectors 
	Equals
	\code
	trans(x)*y
	\endcode
	with x and y vectors
*/
MX inner_prod(const MX &x, const MX &y); // 
/** \brief  Take the outer product of two vectors 
	Equals
	\code
	x*trans(y)
	\endcode
	 with x and y vectors
*/
MX outer_prod(const MX &x, const MX &y); // x*trans(y) with x and y vectors

/** \brief Branching on MX nodes
Ternary operator, "cond ? if_true : if_false"
Internally represented by IfElseNode.
*/
MX if_else(const MX &cond, const MX &if_true, const MX &if_false); 

} // namespace CasADi

namespace std{
//@{
/** \brief  Functions with c equivalents: The implementation and syntax mirrors the standard c functions in math.h */
#define MX CasADi::MX
MX sqrt(const MX &x);
MX sin(const MX &x);
MX cos(const MX &x);
MX tan(const MX &x);
MX atan(const MX &x);
MX asin(const MX &x);
MX acos(const MX &x);
MX exp(const MX &x);
MX log(const MX &x);
MX pow(const MX &x, const MX &n);
MX abs(const MX &x);
MX fabs(const MX &x); // same as abs
MX floor(const MX &x);
MX ceil(const MX &x);
MX erf(const MX &x);
MX fmin(const MX &a, const MX &b);
MX fmax(const MX &a, const MX &b);
#undef MX
//@}
} // namespace std

#endif // MX_HPP
