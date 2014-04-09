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

#ifndef IO_SCHEME_INTERNAL_HPP
#define IO_SCHEME_INTERNAL_HPP

#include "io_scheme.hpp"
#include "../casadi_exception.hpp"
#include <map>


/// \cond INTERNAL
class SharedObjectNode;

namespace casadi{

class CASADI_SYMBOLIC_EXPORT IOSchemeInternal : public SharedObjectNode{
  public:

    /// Name of the scheme
    virtual std::string name() const=0;

    /// List available entries
    virtual std::string entryNames() const=0;

    /// Get index by entry name
    virtual int index(const std::string &name) const=0;

    /// Number of entries
    virtual int size() const=0;

    /// Get the entry name by index
    virtual std::string entry(int i) const=0;

    /// Get the entry enum name by index
    virtual std::string entryEnum(int i) const=0;

    /// Print a description of the object
    virtual void print(std::ostream &stream=std::cout) const=0;

    /// Print a representation of the object
    virtual void repr(std::ostream &stream=std::cout) const=0;

    /// Describe the index
    virtual std::string describe(int i) const=0;

};

class CASADI_SYMBOLIC_EXPORT IOSchemeBuiltinInternal : public IOSchemeInternal {
  public:

    /// Default constructor
    IOSchemeBuiltinInternal(InputOutputScheme scheme);

    /// Name of the scheme
    virtual std::string name() const;

    /// List available entries
    virtual std::string entryNames() const;

    /// Get index by entry name
    virtual int index(const std::string &name) const;

    /// Number of entries
    virtual int size() const;

    /// Get the entry name by index
    virtual std::string entry(int i) const;

    /// Get the entry enum name by index
    virtual std::string entryEnum(int i) const;

    /// Describe the index
    virtual std::string describe(int i) const;

    /// Print a description of the object
    virtual void print(std::ostream &stream=std::cout) const;

    /// Print a representation of the object
    virtual void repr(std::ostream &stream=std::cout) const;

    /// Clone
    virtual IOSchemeBuiltinInternal* clone() const{ return new IOSchemeBuiltinInternal(scheme_); }

  private:
    InputOutputScheme scheme_;

};

class CASADI_SYMBOLIC_EXPORT IOSchemeCustomInternal : public IOSchemeInternal {
  public:

    /// Default constructor
    IOSchemeCustomInternal(const std::vector<std::string> &entries, const std::vector<std::string> &descriptions=std::vector<std::string>());

    /// Name of the scheme
    virtual std::string name() const;

    /// List available entries
    virtual std::string entryNames() const;

    /// Get index by entry name
    virtual int index(const std::string &name) const;

    /// Number of entries
    virtual int size() const;

    /// Get the entry name by index
    virtual std::string entry(int i) const;

    /// Get the entry enum name by index
    virtual std::string entryEnum(int i) const;

    /// Describe the index
    virtual std::string describe(int i) const;

    /// Print a description of the object
    virtual void print(std::ostream &stream=std::cout) const;

    /// Print a representation of the object
    virtual void repr(std::ostream &stream=std::cout) const;

    /// Clone
    virtual IOSchemeCustomInternal* clone() const{ return new IOSchemeCustomInternal(entries_,descriptions_); }

  private:
    std::vector<std::string> entries_;
    std::vector<std::string> descriptions_;
    std::map<std::string,int> entrymap_;
};



} // namespace casadi
/// \endcond

#endif // IO_SCHEME_INTERNAL_HPP
