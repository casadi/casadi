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

#include "printable_object.hpp"
#include <typeinfo>
#include <sstream>


Lbuffer::Lbuffer(std::ostream &sink, int limit_, std::size_t buff_sz) :
    sink_(sink),  buffer_(buff_sz + 1),
    counter(0), limit(limit_)
{
    sink_.clear();
    char *base = &buffer_.front();
    setp(base, base + buffer_.size() - 1);
}

Lbuffer::int_type Lbuffer::overflow(int_type ch)
{
    if (isfull()) return traits_type::eof();
    
    if (sink_ && ch != traits_type::eof())
    {
        counter++;
        if (isfull()) {
          char p[] = ellipsis;
          sink_.write(p, 3);
          if (!dynamic_cast<std::stringstream*>(&sink_)) {
            char p2[] = "\n";
            sink_.write(p2,1);
          }
          return traits_type::eof();
        }
        *pptr() = char(ch);
        pbump(1);
        if (do_process()) return ch;
    }

    return traits_type::eof();
}


int Lbuffer::sync()
{
	return do_process() ? 0 : -1;
}


bool Lbuffer::do_process()
{
    char *p = pbase();
    std::ptrdiff_t n = pptr() - pbase();
    pbump(-n);

    return sink_.write(pbase(), n);
}


using namespace std;
namespace CasADi{

ostream& operator<<(ostream &stream, const PrintableObject& obj){
  obj.repr(stream);
  return stream;
}

void PrintableObject::repr(std::ostream &stream) const{
  // Print description by default
  print(stream);
}

void PrintableObject::print(std::ostream &stream) const{
  // Print name by default
  stream << typeid(this).name();
}

string PrintableObject::getRepresentation() const{
  stringstream ss;
  repr(limited(ss));
  return ss.str();
}

string PrintableObject::getDescription() const{
  stringstream ss;
  print(limited(ss));
  return ss.str();  
}

    
} // namespace CasADi
    


