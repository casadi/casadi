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

#ifndef PRINTABLE_OBJECT_HPP
#define PRINTABLE_OBJECT_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <streambuf>
#include <vector>


/// Enable print limiting with this.
#define LIMIT_PRINTING

/// Maximum amount of characters that gets printed
#define printLimit 20000

/// The symbol that is appended to a stream that reaches the printLimit
#define ellipsis "..."

/// If print limit is reached, return void
#define STREAMLIMITTEST  \
{ \
      Lostream* lostream = dynamic_cast<Lostream*>(&stream); \
      if (lostream) { \
           Lbuffer *s = dynamic_cast<Lbuffer*>(stream.rdbuf()); \
           if (s) { \
             if (s->isfull()) return; \
           } \
      } \
}


#ifdef LIMIT_PRINTING
#ifndef SWIG
#define limited(stream) Lostream(stream).getStream()
#define limitedTo(stream,limit) Lostream(stream,limit).getStream()
#else
#define limited(stream) stream
#define limitedTo(stream,limit) stream
#endif // SWIG
#else
#define limited(stream) stream
#define limitedTo(stream,limit) stream
#endif // LIMIT_PRINTING

#ifndef SWIG
/**
* \brief Limiting buffer
* A buffer that keeps track of how many characters it has processed (counter).
* If counter exceeds limit, the buffer writes an ellipsis '...' to it's sink and ignores any remaining characters.
*
* If limit is negative, the limit is ignored.
*
* This class is based on ideas from http://www.mr-edd.co.uk/blog/beginners_guide_streambuf
*/
class Lbuffer : public std::streambuf
{
    public:
        explicit Lbuffer(std::ostream &sink, int limit = printLimit, std::size_t buff_sz = 0);
        bool isfull() { return counter > limit;}
        
    private:
        bool do_process();
        int_type overflow(int_type ch);
        int sync();

        /// The amount of characters that this buffer has seen so far
        int counter;
        /// The limiting amount of characters
        int limit;

        /// The stream to which all characters are passed
        std::ostream &sink_;
        /// The character buffer
        std::vector<char> buffer_;

};


/* \brief A Limiting version of std::ostream
*
*/
class Lostream : public std::ostream {
  public:
    Lostream(std::ostream &stream_,int limit = printLimit) : buffer(stream_, limit) {
       rdbuf(&buffer);
    };
    
    /** \brief return a reference to self */
    Lostream& getStream() {  return *this; } 
   
  private:
    Lbuffer buffer;
};

#endif // SWIG

namespace CasADi{

/** \brief Base class for objects that have a natural string representation
  \author Joel Andersson 
  \date 2010
*/
class PrintableObject{
  public:
    
#ifndef SWIG
    /// Print a destription of the object
    virtual void print(std::ostream &stream=limited(std::cout)) const;

    /// Print a representation of the object
    virtual void repr(std::ostream &stream=limited(std::cout)) const;

    /// Print a representation of the object to a stream
    friend std::ostream& operator<<(std::ostream &stream, const PrintableObject& obj);
    
    /// Return a string with a representation (for SWIG)
    std::string getRepresentation() const;

    /// Return a string with a destription (for SWIG)
    std::string getDescription() const;
#endif // SWIG    
};

#ifdef SWIG
%extend PrintableObject{
  std::string __str__()  { return $self->getDescription(); }
  std::string __repr__()  { return $self->getRepresentation(); }
}
#endif // SWIG    

} // namespace CasADi


#endif // PRINTABLE_OBJECT_HPP

