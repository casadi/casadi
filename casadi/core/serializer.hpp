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


#ifndef CASADI_SERIALIZER_HPP
#define CASADI_SERIALIZER_HPP


#include "function.hpp"
#include "linsol.hpp"
#include <sstream>
#include <map>
#include <set>

namespace casadi {
  class Slice;
  class Linsol;
  class Sparsity;
  class SparsityInternal;
  class Function;
  class FunctionInternal;
  class MX;
  class MXNode;
  class SXElem;
  class SXNode;
  /** \brief Helper class for Serialization
      \author Joris Gillis
      \date 2018
  */
  class CASADI_EXPORT DeSerializer {
  public:
    DeSerializer(std::istream &in_s);
    void unpack(Sparsity& e);
    void unpack(MX& e);
    void unpack(SXElem& e);
    void unpack(Linsol& e);
    template <class T>
    void unpack(Matrix<T>& e) {
      e = Matrix<T>::deserialize(*this);
    }
    void unpack(Function& e);
    void unpack(Slice& e);
    void unpack(int& e);
    void unpack(bool& e);
    void unpack(casadi_int& e);
    void unpack(std::string& e);
    void unpack(double& e);
    void unpack(char& e);
    template <class T>
    void unpack(std::vector<T>& e) {
      char t;
      unpack(t);
      casadi_assert_dev(t=='V');
      casadi_int s;
      unpack(s);
      e.resize(s);
      for (auto & i : e) unpack(i);
    }


    template <class T>
    void unpack(const std::string& descr, T& e) {
      std::string d;
      unpack(d);
      //uout() << "unpack started: " << descr << std::endl;
      casadi_assert(d==descr, "Mismatch: '" + descr + "' expected, got '" + d + "'.");
      unpack(e);
      //uout() << "unpack: " << descr << ": " << e << std::endl;
    }


    template <class T, class M>
    void shared_unpack(T& e, M& cache) {
      char i;
      unpack("Shared::flag", i);
      switch (i) {
        case 'd': // definition
          e = T::deserialize(*this);
          cache.push_back(e);
          break;
        case 'r': // reference
          {
            casadi_int k;
            unpack("Shared::reference", k);
            e = cache.at(k);
          }
          break;
        default:
          casadi_assert_dev(false);
      }
    }

    void assert_decoration(char e);


  private:
    std::istream& in;
    std::vector<MX> nodes;
    std::vector<Function> functions;
    std::vector<SXElem> sx_nodes;
    std::vector<Sparsity> sparsities;
    std::vector<Linsol> linsols;
  };



  /** \brief Helper class for Serialization


      \author Joris Gillis
      \date 2018
  */
  class CASADI_EXPORT Serializer {
  public:
    /// Constructor
    Serializer(std::ostream& out, const Dict& opts = Dict());

    /// Add a function
    casadi_int add(const Function& f);

    void pack(const Sparsity& e);
    void pack(const MX& e);
    void pack(const SXElem& e);
    void pack(const Linsol& e);
    template <class T>
    void pack(const Matrix<T>& e) {
      e.serialize(*this);
    }
    void pack(const Function& e);
    void pack(const Slice& e);
    void pack(int e);
    void pack(bool e);
    void pack(casadi_int e);
    void pack(double e);
    void pack(const std::string& e);
    void pack(char e);
    template <class T>
    void pack(const std::vector<T>& e) {
      pack('V');
      pack(casadi_int(e.size()));
      for (const auto & i : e) pack(i);
    }

    template <class T>
    void pack(const std::string& descr, const T& e) {
      //uout() << "  pack started: " << descr << std::endl;
      pack(descr);
      pack(e);
      //uout() << "  pack: " << descr << ": " << e << std::endl;
    }

    void decorate(char e);

    template <class T, class M>
    void shared_pack(const T& e, M& map) {
      auto it = map.find(e.get());
      if (it==map.end()) {
        // Not found
        pack("Shared::flag", 'd'); // definition
        e.serialize(*this);
        casadi_int r = map.size();
        map[e.get()] = r;
      } else {
        pack("Shared::flag", 'r'); // reference
        pack("Shared::reference", it->second);
      }
    }

  private:
    std::vector<Function> added_functions_;

    std::map<MXNode*, casadi_int> MX_nodes_;
    std::map<FunctionInternal*, casadi_int> functions_;
    std::map<SXNode*, casadi_int> SX_nodes_;
    std::map<SparsityInternal*, casadi_int> sparsities_;
    std::map<SharedObjectInternal*, casadi_int> linsols_;

    std::ostream& out;

  };


} // namespace casadi

#endif // CASADI_SERIALIZER_HPP
