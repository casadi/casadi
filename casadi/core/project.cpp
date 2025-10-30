/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include "project.hpp"
#include "casadi_misc.hpp"
#include <sstream>
#include <vector>

namespace casadi {

  Project::Project(const MX& x, const Sparsity& sp) {
    set_dep(x);
    set_sparsity(Sparsity(sp));
  }

  std::string Project::disp(const std::vector<std::string>& arg) const {
    if (sparsity().is_dense()) {
      return "dense(" + arg.at(0) + ")";
    } else {
      return "project(" + arg.at(0) + ")";
    }
  }

  template<typename T>
  int Project::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    casadi_project(arg[0], dep().sparsity(), res[0], sparsity(), w);
    return 0;
  }

  int Project::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int Project::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  void Project::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res,
      const std::vector<bool>& unique) const {
    bool unique_arg0 = !unique.empty() && unique[0];
    res[0] = arg[0]->get_project(sparsity(), unique_arg0);
  }

  MX Project::get_project(const Sparsity& sp, bool unique) const {
    if (unique && sp.is_subset(dep(0).sparsity())) {
      return dep(0)->get_project(sp);
    } else {
      return MXNode::get_project(sp, unique);
    }
  }

  MX Project::get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz,
      bool unique) const {
    if (unique) {
      // Decay projection into GetNonzeros
      // Defer simplification logic to GetNonzeros::get_nzref

      std::vector<unsigned char> mapping;
      dep(0).sparsity().intersect(sparsity(), mapping);

      // Nonzero indices of the original matrix
      std::vector<casadi_int> dnz;
      casadi_int i=0;
      for (auto e : mapping) {
        if (e & 2) { // entry present in projection
          if (e & 1) {
            dnz.push_back(i); // present in original
          } else {
            dnz.push_back(-1); // not present in original
          }
        }
        if (e & 1) i++; // Only count a nonzero when present
      }

      return dep(0)->get_nzref(sparsity(), dnz)->get_nzref(sp, nz);
    } else {
      return MXNode::get_nzref(sp, nz);
    }
  }

  void Project::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    casadi_int nfwd = fsens.size();
    for (casadi_int d=0; d<nfwd; ++d) {
      fsens[d][0] = project(fseed[d][0], sparsity() * dep().sparsity(), true);
    }
  }

  void Project::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    casadi_int nadj = aseed.size();
    for (casadi_int d=0; d<nadj; ++d) {
      asens[d][0] += project(aseed[d][0], sparsity() * dep().sparsity(), true);
    }
  }

  int Project::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    sparsity().set(res[0], arg[0], dep().sparsity());
    return 0;
  }

  int Project::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    dep().sparsity().bor(arg[0], res[0], sparsity());
    std::fill(res[0], res[0]+nnz(), 0);
    return 0;
  }

  void Project::generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res,
                          const std::vector<bool>& arg_is_ref,
                          std::vector<bool>& res_is_ref,
                          bool prefer_inline) const {
    g << g.project(g.work(arg.front(), dep().nnz(), arg_is_ref.front()), dep(0).sparsity(),
                           g.work(res.front(), nnz(), false), sparsity(), "w") << "\n";
  }

  void Project::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("Project::type", 'n');
  }

  void Densify::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s); // NOLINT
    s.pack("Project::type", 'd');
  }

  void Sparsify::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s); // NOLINT
    s.pack("Project::type", 's');
  }

  MXNode* Project::deserialize(DeserializingStream& s) {
    char t;
    s.unpack("Project::type", t);
    switch (t) {
      case 'n':
        return new Project(s);
      case 'd':
        return new Densify(s);
      case 's':
        return new Sparsify(s);
      default:
        casadi_assert_dev(false);
    }
  }

  void Densify::generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res,
                          const std::vector<bool>& arg_is_ref,
                          std::vector<bool>& res_is_ref,
                          bool prefer_inline) const {
      casadi_int nrow_x, ncol_x, i, el;
      const casadi_int *colind_x, *row_x;
      const casadi_int* sp_x = dep(0).sparsity();
      nrow_x = sp_x[0]; ncol_x = sp_x[1];
      colind_x = sp_x+2; row_x = sp_x+ncol_x+3;
      std::vector<casadi_int> streaks_target;
      std::vector<casadi_int> streaks_len;
      casadi_int streak_len = 0;
      casadi_int streak_target = 0;

      casadi_int iy_prev = -1;
      casadi_int iy;
      casadi_int offset = 0;
      for (i=0; i<ncol_x; ++i) {
        for (el=colind_x[i]; el!=colind_x[i+1]; ++el) {
          iy = offset+row_x[el];
          if (iy!=iy_prev+1) {
            if (streak_len>0) {
              streaks_len.push_back(streak_len);
              streaks_target.push_back(streak_target);
              streak_len = 0;
              streak_target = iy;
            }
            if (streak_len==0) streak_target = iy;
          }
          streak_len++;
          iy_prev = iy;
        }
        offset += nrow_x;
      }
      if (streak_len>0) {
        streaks_len.push_back(streak_len);
        streaks_target.push_back(streak_target);
        streak_target = iy;
      }
      casadi_assert_dev(sum(streaks_len)==dep().sparsity().nnz());

      if (false && size1()==1 && dep().nnz()>1) {
        // Codegen the indices
        std::string col = g.constant(dep().sparsity().T().get_row());
        g.local("cii", "const casadi_int", "*");
        g.local("ss", "const casadi_real", "*");
        std::string src = g.work(arg.front(), dep().nnz(), arg_is_ref.front());
        std::string dst = g.work(res.front(), nnz(), false);

        g << g.clear(dst, nnz()) << "\n";
        g << "for (cii=" << col << ", ss=" << src << "; cii!=" << col << "+" <<
          dep().nnz() << "; ++cii) " << dst << "[*cii] = *ss++;\n";
      } else {
        std::string src = g.work(arg.front(), dep().nnz(), arg_is_ref.front());
        std::string dst = g.work(res.front(), nnz(), false);

        if (streaks_len.size()<4) {
          g.comment("compact_densify" + str(streaks_len)+str(streaks_target));
          casadi_int offset = 0;
          casadi_int target_prev = 0;
          for (casadi_int i=0;i<streaks_len.size();++i) {
            if (streaks_target.at(i)-target_prev)
              g << g.clear(dst + "+" + str(target_prev), streaks_target.at(i)-target_prev) << "\n";
            g << g.copy(src + "+" + str(offset), streaks_len.at(i), dst + "+" + str(streaks_target.at(i))) << "\n";
            offset += streaks_len.at(i);
            target_prev = streaks_target.at(i)+streaks_len.at(i);
          }
          if (nnz()-target_prev)
            g << g.clear(dst + "+" + str(target_prev), nnz()-target_prev) << "\n";
        } else {
          g << g.densify(src, dep(0).sparsity(), dst);
          g << "// " << dep().sparsity().nnz() << " of " << nnz();
          g << "\n";
        }
      }
  }

  void Sparsify::generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res,
                          const std::vector<bool>& arg_is_ref,
                          std::vector<bool>& res_is_ref,
                          bool prefer_inline) const {
    g << g.sparsify(g.work(arg.front(), dep().nnz(), arg_is_ref.front()),
                           g.work(res.front(), nnz(), false), sparsity()) << "\n";
  }

  template<typename T>
  int Densify::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    casadi_densify(arg[0], dep().sparsity(), res[0], false);
    return 0;
  }

  template<typename T>
  int Sparsify::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    casadi_sparsify(arg[0], res[0], sparsity(), false);
    return 0;
  }

  int Densify::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int Densify::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  int Sparsify::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int Sparsify::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

} // namespace casadi
