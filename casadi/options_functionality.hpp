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

#ifndef OPTIONS_FUNCTIONALITY_HPP
#define OPTIONS_FUNCTIONALITY_HPP

#include "option.hpp"
#include <map>



namespace CasADi{
    
  // Forward declaration
  class OptionsFunctionalityNode;
  
/** \brief Prvides options setting/getting functionality
  Gives a derived class the ability to set and retrieve options in a convenient way.
  It also contains error checking, making sure that the option exists and that the value type is correct.
  
  A derived class should add option names, types and default values to the corresponding vectors.
  \author Joel Andersson 
  \date 2010
  Joel Andersson, K.U. Leuven 2010
  joel.andersson@esat.kuleuven.be
*/
class OptionsFunctionality : public SharedObject{
  public:

    
    /// Default constructor
    OptionsFunctionality();
    
    /// Destructor
    ~OptionsFunctionality();
    
    /// Get a pointer to the node
    const OptionsFunctionalityNode* get() const;
    OptionsFunctionalityNode* get();

    /// Access a member function or object
    OptionsFunctionalityNode* operator->();
    const OptionsFunctionalityNode* operator->() const;
        
    /** \brief  set option */
    void setOption(const std::string &str, const Option& val);

    /** \brief  get an option value */
    Option getOption(const std::string &str) const;

    /** \brief  check if there is an option str */
    bool hasOption(const std::string &str) const;

    /** \brief  check if the user has there is an option str */
    bool hasSetOption(const std::string &str) const;

    /** \brief  Print options to a stream */
    void printOptions(std::ostream &stream=std::cout) const;

  /// Assert that the node is pointing to the right type of object
    virtual bool checkNode() const;

    /** \brief  Copy all options from another object*/
    void copyOptions(const OptionsFunctionality& obj);

};
      
/** \brief Internal class
  \author Joel Andersson 
  \date 2010
*/
class OptionsFunctionalityNode : public SharedObjectNode{

public:
  
/// Constructor, destructor
OptionsFunctionalityNode();
virtual ~OptionsFunctionalityNode();
 
  /** \brief  set option */
  void setOption(const std::string &str, const Option& val);

  /** \brief  check if there is an option str */
  bool hasOption(const std::string &str) const;

  /** \brief  check if the user has there is an option str */
  bool hasSetOption(const std::string &str) const;

  /** \brief  Print options to a stream */
  void printOptions(std::ostream &stream=std::cout) const;
    
  /** \brief  get an option value */
  Option getOption(const std::string &str) const;

  /** \brief  Print */
  virtual void print(std::ostream &stream) const = 0;

  /** \brief  Copy all options from another object*/
  void copyOptions(const OptionsFunctionality& obj);
  
protected:

  void addOption(const std::string &str, const opt_type& type, const Option &def_val=Option());

private:
/** \brief  Allowed options  */
  std::map<std::string, opt_type> allowed_options;

/** \brief  User-set options */
  std::map<std::string, Option> options;

};

} // namespace CasADi


#endif // OPTIONS_FUNCTIONALITY
