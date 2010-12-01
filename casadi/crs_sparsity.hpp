#ifndef CRS_SPARSITY_HPP
#define CRS_SPARSITY_HPP

#include "shared_object.hpp"
#include <vector>

namespace CasADi{

// Forward declaration
class CRSSparsityNode;
  
class CRSSparsity : public SharedObject{
  public:
  
    /// \brief Default constructor, null pointer
    CRSSparsity();
    
    /// \brief Construct a sparsity pattern (dense)
    CRSSparsity(int nrow, int ncol);

    /// \brief Construct a sparsity pattern from vectors
    CRSSparsity(int nrow, int ncol, std::vector<int> col, std::vector<int> rowind);
    
    /// \brief Access a member function or object
    CRSSparsityNode* operator->();
    const CRSSparsityNode* operator->() const;
  
    /// \brief Assert that the node is pointing to the right type of object
    void assertNode() const;
    
    /// \brief Scalar expression
    static const CRSSparsity scalar;
};

class CRSSparsityNode : public SharedObjectNode{
  public:
    /// \brief Construct a sparsity pattern from vectors
    CRSSparsityNode(int nrow, int ncol, std::vector<int> col, std::vector<int> rowind);

    //! \brief Number of rows
    int nrow_;
    //! \brief Number of columns
    int ncol_;
    //! vector of length nnz containing the columns for all the indices of the non-zero elements
    std::vector<int> col_;
    //! vector of length n+1 containing the index of the last non-zero element up till each row 
    std::vector<int> rowind_;
};

} // namespace CasADi

#endif // CRS_SPARSITY_HPP
