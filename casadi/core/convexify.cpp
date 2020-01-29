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


#include "convexify.hpp"

using namespace std;
namespace casadi {

  std::string strategy_to_string(Convexify::Strategy s) {
    switch (s) {
      case Convexify::CVX_REGULARIZE: return "regularize";
      case Convexify::CVX_EIGEN_REFLECT: return "eigen-reflect";
      case Convexify::CVX_EIGEN_CLIP: return "eigen-clip";
    }
    return "unknown";
  }

  Convexify::Convexify(const MX& H, const Dict& opts) {

    // Validate and categorize matrix input sparsity
    casadi_assert(H.is_square(), "Convexify ");
    if (H.sparsity().is_symmetric()) {
      type_in_ = SYMM;
    } else if (H.is_tril()) {
      type_in_ = TRIL;
    } else if (H.is_triu()) {
      type_in_ = TRIU;
    } else {
      casadi_error("Convexify operation requires symmetric or triangular input");
    }
    set_dep(H);

    // Read options
    margin_ = 1e-7;
    max_iter_eig_ = 50;
    std::string strategy = "regularize";

    for (auto&& op : opts) {
      if (op.first=="strategy") {
        strategy = op.second.to_string();
      } else if (op.first=="margin") {
        margin_ = op.second;
        casadi_assert(margin_>=0, "Margin must be >=0");
      } else if (op.first=="max_iter_eig") {
        max_iter_eig_ = op.second;
      } else {
        casadi_error("Unknown option '" + op.first + "'.");
      }
    }

    // Interpret strategy
    if (strategy=="regularize") {
      strategy_ = CVX_REGULARIZE;
      casadi_assert(type_in_==SYMM, "Only truly symmetric matrices supported");
    } else if (strategy=="eigen-reflect") {
      strategy_ = CVX_EIGEN_REFLECT;
    } else if (strategy=="eigen-clip") {
      strategy_ = CVX_EIGEN_CLIP;
    } else {
      casadi_error("Invalid convexify strategy. "
        "Choose from regularize|eigen-reflect|eigen-clip");
    }

    Sparsity Hrsp = H.sparsity()+H.sparsity().T();

    Sparsity Hsp;
    if (strategy_==CVX_EIGEN_REFLECT || strategy_==CVX_EIGEN_CLIP) {
      // Uncover strongly connected components
      std::vector<casadi_int> scc_index;
      casadi_int scc_nb = Hrsp.scc(scc_index, scc_offset_);

      // Represent Hessian as block-dense in permuted space
      std::vector<Sparsity> sp;
      for (casadi_int i=0;i<scc_nb;++i) {
        casadi_int block = scc_offset_.at(i+1)-scc_offset_.at(i);
        Sparsity stencil;
        if (type_in_==SYMM) {
          stencil = Sparsity::dense(block, block);
        } else if (type_in_==TRIL) {
          stencil = Sparsity::lower(block);
        } else {
          stencil = Sparsity::upper(block);
        }
        sp.push_back(stencil);
      }

      std::vector<casadi_int> ssc_perm = lookupvector(scc_index);
      std::vector<casadi_int> mapping_dummy;
      Hsp = diagcat(sp).sub(ssc_perm, ssc_perm, mapping_dummy);
      scc_sp_ = Hsp.sub(scc_index, scc_index, scc_mapping_);

      // Find out size of maximum block
      block_size_ = 0;
      for (casadi_int i=0;i<scc_nb;++i) {
        casadi_int block = scc_offset_.at(i+1)-scc_offset_.at(i);
        if (block>block_size_) block_size_ = block;
      }

      scc_transform_ = scc_offset_!=range(H.size1());
    } else if (strategy_==CVX_REGULARIZE) {
      Hsp = Hrsp + Sparsity::diag(H.size1());
    } else {
      Hsp = Hrsp;
    }
    Hsp_project_ = Hsp!=Hrsp;

    if (type_in_==TRIL) {
      Hsp = Sparsity::tril(Hsp);
    } else if (type_in_==TRIU) {
      Hsp = Sparsity::triu(Hsp);
    }
    set_sparsity(Hsp);
  }

  size_t Convexify::sz_iw() const {
    size_t sz_iw_cvx = 0;
    if (strategy_==CVX_EIGEN_REFLECT || strategy_==CVX_EIGEN_CLIP) {
      sz_iw_cvx = 1+3*max_iter_eig_;
    }
    return sz_iw_cvx;
  }

  size_t Convexify::sz_w() const {
    casadi_int sz_w_cvx = 0;
    if (strategy_==CVX_EIGEN_REFLECT || strategy_==CVX_EIGEN_CLIP) {
      sz_w_cvx = max(block_size_, 2*(block_size_-1)*max_iter_eig_);
      if (Hsp_project_) sz_w_cvx = max(sz_w_cvx, size1());
      if (scc_transform_) sz_w_cvx += block_size_*block_size_;
    }
    size_t sz_w = size1()+sz_w_cvx;
    return sz_w;
  }

  std::string Convexify::disp(const std::vector<std::string>& arg) const {
    return "convexify(" + arg.at(0) + ")";
  }

  void Convexify::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    Dict options;
    options["strategy"] = strategy_to_string(strategy_);
    options["margin"] = margin_;
    options["max_iter_eig"] = max_iter_eig_;
    res[0] = convexify(arg[0], options);
  }

  int Convexify::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    double* Bk =  res[0];
    const Sparsity& Hsp = sparsity();

    if (Hsp_project_) {
      casadi_project(arg[0], dep(0).sparsity(), Bk, Hsp, w);
    } else {
      casadi_copy(arg[0], nnz(), Bk);
    }

    if (strategy_==CVX_REGULARIZE) {
      // Determing regularization parameter with Gershgorin theorem
      double reg = margin_-casadi_lb_eig(Hsp, Bk);
      if (reg > 0) casadi_regularize(Hsp, Bk, reg);
    } else if (strategy_==CVX_EIGEN_REFLECT || strategy_==CVX_EIGEN_CLIP) {
      casadi_int offset = 0;

      // Loop over Hessian blocks
      for (casadi_int k=0;k<scc_offset_.size()-1;++k) {
        casadi_int block_size = scc_offset_.at(k+1)-scc_offset_.at(k);

        double *H_block = w;
        double *w_cvx = w;

        // Set w_cvx to dense Hessian block from Bk
        if (scc_transform_) {
          casadi_int kk=0;

          // Loop over columns of block
          for (casadi_int i=0;i<block_size;++i) {
            // Loop over elements in column
            if (type_in_==SYMM) {
              for (casadi_int j=0;j<block_size;++j, ++kk) {
                H_block[kk] = Bk[scc_mapping_[offset+kk]];
              }
            } else if (type_in_==TRIU) {
              for (casadi_int j=0;j<i+1;++j, ++kk) {
                double e = Bk[scc_mapping_[offset+kk]];
                H_block[i*block_size+j] = e;
                H_block[i+block_size*j] = e;
              }
            } else {
              for (casadi_int j=i;j<block_size;++j, ++kk) {
                double e = Bk[scc_mapping_[offset+kk]];
                H_block[i*block_size+j] = e;
                H_block[i+block_size*j] = e;
              }
            }
          }
          w_cvx += block_size*block_size;
        } else {
          H_block = Bk+offset;
        }

        // Perform convexification
        int ret = casadi_cvx(block_size, H_block, margin_, 1e-10,
          strategy_==CVX_EIGEN_REFLECT, max_iter_eig_, w_cvx, iw);
        casadi_assert(!ret, "Failure in convexification.");

        // Fill in upper-rectangular part
        for (casadi_int i=0;i<block_size;++i) {
          for (casadi_int j=0;j<i+1;++j) {
            H_block[block_size*i+j] = H_block[block_size*j+i];
          }
        }

        // Put results back in Bk
        if (scc_transform_) {
          casadi_int kk=0;
          // Loop over columns of block
          for (casadi_int i=0;i<block_size;++i) {
            // Loop over elements in column
            if (type_in_==SYMM) {
              for (casadi_int j=0;j<block_size;++j, ++kk) {
                Bk[scc_mapping_[offset+kk]] = H_block[kk];
              }
            } else if (type_in_==TRIU) {
              for (casadi_int j=0;j<i+1;++j, ++kk) {
                Bk[scc_mapping_[offset+kk]] = H_block[block_size*i+j];
              }
            } else {
              for (casadi_int j=i;j<block_size;++j, ++kk) {
                Bk[scc_mapping_[offset+kk]] = H_block[block_size*i+j];
              }
            }
          }
        }
        if (type_in_==SYMM) {
          offset += block_size*block_size;
        } else {
          offset += block_size*(block_size+1)/2;
        }
      }
    }

    return 0;
  }

  void Convexify::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.version("Convexify", 1);
    s.pack("Convexify::type_in", static_cast<int>(type_in_));
    s.pack("Convexify::strategy", static_cast<int>(strategy_));
    s.pack("Convexify::margin", margin_);
    s.pack("Convexify::max_iter_eig", max_iter_eig_);
    s.pack("Convexify::scc_offset", scc_offset_);
    s.pack("Convexify::scc_mapping", scc_mapping_);
    s.pack("Convexify::block_size", block_size_);
    s.pack("Convexify::scc_sp", scc_sp_);
    s.pack("Convexify::Hsp_project", Hsp_project_);
    s.pack("Convexify::scc_transform", scc_transform_);
  }

  Convexify::Convexify(DeserializingStream& s) : MXNode(s) {
    s.version("Convexify", 1);
    int type_in;
    s.unpack("Convexify::type_in", type_in);
    type_in_ = static_cast<TypeIn>(type_in);
    int strategy;
    s.unpack("Convexify::strategy", strategy);
    strategy_ = static_cast<Strategy>(strategy);
    s.unpack("Convexify::margin", margin_);
    s.unpack("Convexify::max_iter_eig", max_iter_eig_);
    s.unpack("Convexify::scc_offset", scc_offset_);
    s.unpack("Convexify::scc_mapping", scc_mapping_);
    s.unpack("Convexify::block_size", block_size_);
    s.unpack("Convexify::scc_sp", scc_sp_);
    s.unpack("Convexify::Hsp_project", Hsp_project_);
    s.unpack("Convexify::scc_transform", scc_transform_);
  }

  void Convexify::generate(CodeGenerator& g,
                       const std::vector<casadi_int>& arg,
                       const std::vector<casadi_int>& res) const {
    std::string Bproj = "w";
    const Sparsity& Hsp = sparsity();
    std::string Bk = g.work(res[0], nnz());
    std::string arg0 = g.work(arg[0], dep(0).nnz());

    if (Hsp_project_) {
      g << g.project(arg0, dep(0).sparsity(),  Bk, sparsity(), Bproj) << "\n";
    } else {
      g << g.copy(arg0, nnz(), Bk) << "\n";
    }

    if (strategy_==CVX_REGULARIZE) {
      g.comment("Determing regularization parameter with Gershgorin theorem");
      g.local("r", "casadi_real");
      g << "r = " << margin_ << "-" + g.lb_eig(Hsp, Bk) << ";\n";
      g << "if (r>0) " << g.regularize(Hsp, Bk, "r") << "\n";
    } else if (strategy_==CVX_EIGEN_REFLECT || strategy_==CVX_EIGEN_CLIP) {
      g.add_auxiliary(CodeGenerator::AUX_CVX);
      g.local("offset", "casadi_int");
      g << "offset = 0;\n";

      g.comment("Loop over Hessian blocks");
      g.local("k", "casadi_int");
      g << "for (k=0;k<" << scc_offset_.size()-1 << ";++k) {\n";
      g.local("block_size", "casadi_int");
      std::string scc_offset = g.constant(scc_offset_);
      std::string scc_mapping = g.constant(scc_mapping_);
      g << "block_size = " << scc_offset << "[k+1]-" << scc_offset << "[k];\n";

      std::string H_block, w_cvx;

      g.comment("Set w_cvx to dense Hessian block from Bk");
      if (scc_transform_) {
        g.local("kk", "casadi_int");
        g.local("i", "casadi_int");
        g.local("j", "casadi_int");
        g << "for (i=0,kk=0;i<block_size;++i) {\n";
        if (type_in_==SYMM) {
          g << "for (j=0;j<block_size;++j,++kk) {\n";
          g << "w[kk] = " << Bk << "[" << scc_mapping << "[offset+kk]];\n";
          g << "}\n";
        } else {
          if (type_in_==TRIU) g << "for (j=0;j<i+1;++j,++kk) {\n";
          if (type_in_==TRIL) g << "for (j=i;j<block_size;++j,++kk) {\n";
          g.local("e", "casadi_real");
          g << "e = " << Bk << "[" << scc_mapping << "[offset+kk]];\n";
          g << "w[i*block_size+j]=e;\n";
          g << "w[j*block_size+i]=e;\n";
          g << "}\n";
        }
        g << "}\n";
        w_cvx = "w+block_size*block_size";
        H_block = "w";
      } else {
        H_block = Bk + "+offset";
        w_cvx = "w";
      }

      g.comment("Perform convexification");
      g << "if (casadi_cvx(block_size, " + H_block + ", " + str(margin_) +
        ", 1e-10, " + str(strategy_==CVX_EIGEN_REFLECT ? 1 : 0) + ", " +
        str(max_iter_eig_) + ", " + w_cvx + ", iw)) return 1;\n";

      g.comment("Fill in upper-rectangular part");
      g << "for (i=0;i<block_size;++i) {\n";
      g << "for (j=0;j<i+1;++j) {\n";
      g << H_block + "[block_size*i+j] = " + H_block + "[block_size*j+i];\n";
      g << "}\n";
      g << "}\n";

      if (scc_transform_) {
        g.comment("Put results back in Bk");
        g << "for (i=0,kk=0;i<block_size;++i) {\n";
        if (type_in_==SYMM) {
          g << "for (j=0;j<block_size;++j,++kk) {\n";
          g << Bk << "[" << scc_mapping << "[offset+kk]] = " + H_block + "[kk];\n";
          g << "}\n";
        } else {
          if (type_in_==TRIU) g << "for (j=0;j<i+1;++j,++kk) {\n";
          if (type_in_==TRIL) g << "for (j=i;j<block_size;++j,++kk) {\n";
          g << Bk << "[" << scc_mapping << "[offset+kk]] = " + H_block + "[block_size*i+j];\n";
          g << "}\n";
        }
        g << "}\n";
      }

      g << "offset += kk;\n";
      g << "}\n";
    }
  }

} // namespace casadi
