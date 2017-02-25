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


#ifndef UNARY_SXElem_HPP
#define UNARY_SXElem_HPP

#include "sx_node.hpp"
#include "casadi/core/global_options.hpp"
#include <stack>

/// \cond INTERNAL

#include <unordered_map>

namespace casadi {

// hash the key type of the UnarySX cache
struct UnarySXCacheKeyHash
{
  size_t operator()(const std::tuple<int, size_t>& v) const
  {                                              
    std::size_t ret=0;
    hash_combine(ret, std::get<0>(v));
    hash_combine(ret, std::get<1>(v));
    return ret;
  }                                              
};

/** \brief Represents a basic unary operation on an SXElem node
  \author Joel Andersson
  \date 2012
*/
class CASADI_EXPORT UnarySX : public SXNode {
  private:

    /** \brief  Constructor is private, use "create" below */
    UnarySX(unsigned char op, const SXElem& dep) : op_(op), dep_(dep) {}

  public:

    /** \brief  Create a unary expression */
    inline static SXElem create(unsigned char op, const SXElem& dep) {
      if (dep.is_constant()) {
        // Evaluate constant
        double dep_val(dep);
        double ret_val;
        casadi_math<double>::fun(op, dep_val, dep_val, ret_val);
        return ret_val;
      } else {
        if(!GlobalOptions::getSXCaching()) {
          // Expression containing free variables
          return SXElem::create(new UnarySX(op, dep));
        }

        // The key of the cache
        auto key = std::make_tuple(op, dep.__hash__());
        // Find the key in the cache
        auto it = cached_unarysx_.find(key);
        if (it==cached_unarysx_.end()) { // not found -> create new object and return it
          // Allocate a new object
          UnarySX* n = new UnarySX(op, dep);

          // Add to hash_table
          cached_unarysx_.emplace_hint(it, key, n);

          // Return it to caller
          return SXElem::create(n);
        } else { // found -> return cached object
          return SXElem::create(it->second);
        }
      }
    }

    /** \brief Destructor */
    virtual ~UnarySX() {
      assert(count == 0 || count == SXNode::countToBeDeleted);

      // If count == 0 then this destructor was called directly (not via the none recursive ~BinarySX hack).
      // In this case remove this node from the cache. If count is SXNode::countToBeDeleted then this node
      // was already removed from the cache by the call to assignNoDelete in ~BinarySX.
      if(count==0)
        removeFromCache();
    }

    void removeFromCache() {
      size_t num_erased = cached_unarysx_.erase(std::make_tuple(op_, dep_.__hash__()));
      // assert(num_erased==1);  only possible if the cache is enabled permanently
      (void)num_erased;
    }

    virtual bool is_smooth() const { return operation_checker<SmoothChecker>(op_);}

    virtual bool hasDep() const { return true; }

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool is_equal(const SXNode* node, int depth) const {
      const UnarySX* n = dynamic_cast<const UnarySX*>(node);
      return n && n->op_ == op_ &&  SXElem::is_equal(n->dep_, dep_, depth-1);
    }

    /** \brief  Number of dependencies */
    virtual int ndep() const { return 1;}

    /** \brief  get the reference of a dependency */
    virtual const SXElem& dep(int i) const { return dep_; }

    /** \brief  Get the operation */
    virtual int op() const { return op_;}

    /** \brief  Print expression */
    virtual std::string print(const std::string& arg1, const std::string& arg2) const {
      std::stringstream ss;

      // Print the prefix
      casadi_math<double>::printPre(op_, ss);

      // Print the dependency
      ss << arg1;

      // Print the suffix
      casadi_math<double>::printPost(op_, ss);

      return ss.str();
    }

  protected:
    /** \brief  The binary operation as an 1 byte integer (allows 256 values) */
    const unsigned char op_;

    /** \brief  The dependencies of the node */
    const SXElem dep_;

    /** \brief Hash map of all binary SX currently allocated
     * (storage is allocated for it in sx_element.cpp) */
    static std::unordered_map<std::tuple<int, size_t>, UnarySX*, UnarySXCacheKeyHash> cached_unarysx_;
};

} // namespace casadi

/// \endcond
#endif // UNARY_SXElem_HPP
