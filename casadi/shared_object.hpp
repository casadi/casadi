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

#ifndef SHARED_OBJECT_HPP
#define SHARED_OBJECT_HPP

#include "printable_object.hpp"

namespace CasADi{

/** \brief Forward declaration
  \author Joel Andersson 
  \date 2010
*/
class SharedObjectNode;

/** \brief Shared pointer class
  \author Joel Andersson 
  \date 2010	
*/
class SharedObject : public PrintableObject{
  template<class B> friend B shared_cast(SharedObject& A);
  template<class B> friend const B shared_cast(const SharedObject& A);

  public:
    /// Default constructor
    SharedObject();
    
    /// Copy constructor
    SharedObject(const SharedObject& ref);

    // Destructor
    ~SharedObject();
    
    /// Assignment operator
    SharedObject& operator=(const SharedObject& ref);
    
    /// Get a pointer to the node
    const SharedObjectNode* get() const;
    SharedObjectNode* get();

    /// Access a member function or object
    SharedObjectNode* operator->();
    const SharedObjectNode* operator->() const;

    /// Print a representation of the object
    virtual void repr(std::ostream &stream) const;

    /// Print a destription of the object
    virtual void print(std::ostream &stream=std::cout) const;
    
    //! Initialize the object: more documentation in the node class (SharedObjectNode and derived classes)
    void init();
    
    //! Is a null pointer?
    bool isNull() const; 

    /// Assert that the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Assign the node to something
    void assignNode(SharedObjectNode* node);
    
    /// If there are other references to the object, then make a deep copy of it and point to this new object
    void makeUnique();
    
  private:
    SharedObjectNode *node;
    void count_up(); // increase counter of the node
    void count_down(); // decrease counter of the node
};

// Node class
class SharedObjectNode{
  friend class SharedObject;
  public:
  
  //! Default constructor
  SharedObjectNode();

  //! Copy constructor
  SharedObjectNode(const SharedObjectNode& node);
  
  //! Assignment operator
  SharedObjectNode& operator=(const SharedObjectNode& node);
  
  //! Destructor
  virtual ~SharedObjectNode() = 0;  

  /// Make a deep copy of the instance  
  virtual SharedObjectNode* clone() const;

  //! Initialize the object
  virtual void init();

  /// Print a representation of the object
  virtual void repr(std::ostream &stream) const;

  /// Print a destription of the object
  virtual void print(std::ostream &stream) const;

  private:
    //! Number of references pointing to the object
    unsigned int count;
};

// Typecast a shared object to another one using a dynamic cast
template<class B>
B shared_cast(SharedObject& A){
  
  //! Get a pointer to the node
  SharedObjectNode* ptr = A.get();
  
  //! Create a return object
  B ret;
  
  //! Assign node of B and return
  ret.assignNode(ptr);
      
  //! Null pointer if not pointing towards the right type of object
  if(!ret.checkNode()) ret.assignNode(0);

  return ret;
}

// Typecast a const shared object to another one using a dynamic cast
template<class B>
const B shared_cast(const SharedObject& A){
  SharedObject A_copy = A;
  return shared_cast<B>(A_copy);
}


} // namespace CasADi


#endif // SHARED_OBJECT_HPP

