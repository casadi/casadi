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


#include "kernel_sum_2d_internal.hpp"

namespace casadi {
  using namespace std;

  KernelSum2D::KernelSum2D() {
  }

  KernelSum2D::KernelSum2D(const std::string& name, const Function& f,
                           const std::pair<int, int> & size,
                           double r,
                           int n,
                 const Dict& opts) {
    assignNode(new KernelSum2DInternal(name, f, size, r, n));
    setOption(opts);
    init();
  }

  KernelSum2DInternal* KernelSum2D::operator->() {
    return static_cast<KernelSum2DInternal*>(Function::operator->());
  }

  const KernelSum2DInternal* KernelSum2D::operator->() const {
    return static_cast<const KernelSum2DInternal*>(Function::operator->());
  }

  bool KernelSum2D::test_cast(const SharedObjectNode* ptr) {
    return dynamic_cast<const KernelSum2DInternal*>(ptr)!=0;
  }

} // namespace casadi
