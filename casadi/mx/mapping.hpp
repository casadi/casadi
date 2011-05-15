/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef MAPPING_HPP
#define MAPPING_HPP

#include "mx_node.hpp"
#include <map>

namespace CasADi{
  /** \brief Maps non-zero elements
  \author Joel Andersson
  \date 2011
*/
class Mapping : public MXNode{
  public:

    /// Default constructor
    Mapping(const CRSSparsity& sp);

    /// Clone function
    Mapping* clone() const;
      
    /// Destructor
    virtual ~Mapping(){}
    
    /// Evaluate the function and store the result in the node
    virtual void evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj);

    /// Print
    virtual void print(std::ostream &stream, const std::vector<std::string>& args) const;
    
    /// Is a mapping matrix
    virtual bool isMapping() const{return true;}
    
    /// Add a dependency (index given)
    virtual void addDependency(int depind, const std::vector<int>& nz_d, const std::vector<int>& nz);
    
    /// Add a dependency
    virtual void addDependency(const MX& d, const std::vector<int>& nz_d, const std::vector<int>& nz);
    
    /// Add a dependency
    virtual void addDependency(const MX& d, const std::vector<int>& nz_d);
    
    /// Add a dependency (indices nz to be placed in (i,j))
    void addDepend(const MX& d, std::vector<int> nz, std::vector<int> i, std::vector<int> j);
    
    /// Symbolic forward sensitivities
    virtual MX adFwd(const std::vector<MX>& jx);

  protected:
    
    /// Check if the mapping is ready
    bool isReady() const;
    
    /// Mapping from the output non-zero to the dependency nonzero index
    std::vector<int> nzind_;

    /// Mapping from the output non-zero index of the dependency index
    std::vector<int> depind_;

    /// Map to locate the dependencies
    std::map<const MXNode*, int> depmap_;
};

} // namespace CasADi

#endif // MAPPING_HPP
