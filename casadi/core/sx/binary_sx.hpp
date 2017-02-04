/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#ifndef CASADI_BINARY_SX_HPP
#define CASADI_BINARY_SX_HPP

#include "sx_node.hpp"
#include <stack>


/// \cond INTERNAL

#include <unordered_map>

namespace casadi {

// hash the key type of the BinarySX cache
struct BinarySXCacheKeyHash
{
  size_t operator()(const std::tuple<int, size_t, size_t>& v) const
  {                                              
    std::size_t ret=0;
    hash_combine(ret, std::get<0>(v));
    hash_combine(ret, std::get<1>(v));
    hash_combine(ret, std::get<2>(v));
    return ret;
  }                                              
};

/** \brief Represents a basic binary operation on two SXElem nodes
  \author Joel Andersson
  \date 2010
*/
class CASADI_EXPORT BinarySX : public SXNode {
  private:

    /** \brief  Constructor is private, use "create" below */
    BinarySX(unsigned char op, const SXElem& dep0, const SXElem& dep1) :
        op_(op), dep0_(dep0), dep1_(dep1) {}

  public:

    /** \brief  Create a binary expression */
    inline static SXElem create(unsigned char op, const SXElem& dep0, const SXElem& dep1) {
      if (dep0.is_constant() && dep1.is_constant()) {
        // Evaluate constant
        double dep0_val(dep0);
        double dep1_val(dep1);
        double ret_val;
        casadi_math<double>::fun(op, dep0_val, dep1_val, ret_val);
        return ret_val;
      } else {
        // The key of the cache
        auto key = std::make_tuple(op, dep0.__hash__(), dep1.__hash__());
        // Find the key in the cache
        auto it = cached_binarysx_.find(key);
        // If not found and op is commutative try with swapped deps.
        // This way the same cache for e.g. x+y and y+x is used.
        if (it==cached_binarysx_.end() && operation_checker<CommChecker>(op))
          it = cached_binarysx_.find(std::make_tuple(op, dep1.__hash__(), dep0.__hash__()));
        if (it==cached_binarysx_.end()) { // not found -> create new object and return it
          // Allocate a new object
          BinarySX* n = new BinarySX(op, dep0, dep1);

          // Add to hash_table
          cached_binarysx_.emplace_hint(it, key, n);

          // Return it to caller
          return SXElem::create(n);
        } else { // found -> return cached object
          return SXElem::create(it->second);
        }
      }
    }

    /** \brief Destructor
    This is a rather complex destructor which is necessary since the default destructor
    can cause stack overflow due to recursive calling.
    */
    virtual ~BinarySX() {
      assert(count == 0 || count == SXNode::countToBeDeleted);

      // If count == 0 then this destructor was called directly (not via the none recursive ~BinarySX hack).
      // In this case remove this node from the cache. If count is SXNode::countToBeDeleted then this node
      // was already removed from the cache by the call to assignNoDelete in ~BinarySX.
      if(count==0)
        removeFromCache();

      // Start destruction method if any of the dependencies has dependencies
      for (int c1=0; c1<2; ++c1) {
        // Get the node of the dependency and remove it from the smart pointer
        SXNode* n1 = const_cast<SXElem&>(dep(c1)).assignNoDelete(casadi_limits<SXElem>::nan);

        // Check if this was the last reference but node is not deleted already
        if (n1->count==SXNode::countToBeDeleted) {

          // Check if binary
          if (!n1->hasDep()) { // n1 is not binary

            delete n1; // Delete straight away

          } else { // n1 is binary

            // Stack of expressions to be deleted
            std::stack<SXNode*> deletion_stack;

            // Add the node to the deletion stack
            deletion_stack.push(n1);

            // Process stack
            while (!deletion_stack.empty()) {

              // Top element
              SXNode *t = deletion_stack.top();

              // Check if the top element has dependencies with dependencies
              bool added_to_stack = false;
              for (int c2=0; c2<t->ndep(); ++c2) { // for all dependencies of the dependency

                // Get the node of the dependency of the top element
                // and remove it from the smart pointer
                SXNode *n2 = const_cast<SXElem&>(t->dep(c2)).assignNoDelete(casadi_limits<SXElem>::nan);

                // Check if this is the only reference to the element but node is not deleted already
                if (n2->count == SXNode::countToBeDeleted) {

                  // Check if binary
                  if (!n2->hasDep()) {

                    // Delete straight away if not binary
                    delete n2;

                  } else {

                    // Add to deletion stack
                    deletion_stack.push(n2);
                    added_to_stack = true;
                  }
                }
              }

              // Delete and pop from stack if nothing added to the stack
              if (!added_to_stack) {
                delete deletion_stack.top();
                deletion_stack.pop();
              }
            } // while
          }
        }
      }
    }

    void removeFromCache() {
      // remove from cache
      size_t num_erased = cached_binarysx_.erase(std::make_tuple(op_, dep0_.__hash__(), dep1_.__hash__()));
      assert(num_erased==1);
      (void)num_erased;
    }

    virtual bool is_smooth() const { return operation_checker<SmoothChecker>(op_);}

    virtual bool hasDep() const { return true; }

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool is_equal(const SXNode* node, int depth) const {
      const BinarySX* n = dynamic_cast<const BinarySX*>(node);
      if (n==0) return false;
      if (n->op_ != op_) return false;
      if (SXElem::is_equal(n->dep0_, dep0_, depth-1)
          && SXElem::is_equal(n->dep1_, dep1_, depth-1)) return true;
      if (operation_checker<CommChecker>(op_)
          && SXElem::is_equal(n->dep1_, dep0_, depth-1)
          && SXElem::is_equal(n->dep0_, dep1_, depth-1)) return true;
      return false;
    }

    /** \brief  Number of dependencies */
    virtual int ndep() const { return 2;}

    /** \brief  get the reference of a dependency */
    virtual const SXElem& dep(int i) const { return i==0 ? dep0_ : dep1_;}

    /** \brief  Get the operation */
    virtual int op() const { return op_;}

    /** \brief  Print expression */
    virtual std::string print(const std::string& arg1, const std::string& arg2) const {
      std::stringstream ss;

      // Print the prefix
      casadi_math<double>::printPre(op_, ss);

      // Print the first dependency
      ss << arg1;

      //Print the infix
      casadi_math<double>::printSep(op_, ss);

      // Print the second dependency
      ss << arg2;

      // Print the suffix
      casadi_math<double>::printPost(op_, ss);

      return ss.str();
    }

  protected:
    /** \brief  The binary operation as an 1 byte integer (allows 256 values) */
    const unsigned char op_;

    /** \brief  The dependencies of the node */
    const SXElem dep0_, dep1_;

    /** \brief Hash map of all binary SX currently allocated
     * (storage is allocated for it in sx_element.cpp) */
    static std::unordered_map<std::tuple<int, size_t, size_t>, BinarySX*, BinarySXCacheKeyHash> cached_binarysx_;
};

} // namespace casadi
/// \endcond

#endif // CASADI_BINARY_SX_HPP
