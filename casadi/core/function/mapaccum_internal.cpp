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


#include "mapaccum_internal.hpp"

using namespace std;

namespace casadi {

  MapAccumInternal::MapAccumInternal(const std::string& name, const Function& f, int n,
                                     const std::vector<bool>& input_accum,
                                     const std::vector<int>& output_accum,
                                     bool reverse)
    : FunctionInternal(name), f_(f), n_(n),
      input_accum_(input_accum), output_accum_(output_accum), reverse_(reverse) {

    casadi_assert(n>0);
    casadi_assert(f.n_in()>=1);
    casadi_assert(f.n_out()>=1);

    casadi_assert(input_accum.size()==f_.n_in());
    int n_accum=0;
    for (int i=0;i<input_accum.size();++i) {
      if (input_accum[i]) {
        casadi_assert(n_accum<output_accum.size());
        casadi_assert(output_accum[n_accum]>=0);
        casadi_assert(output_accum[n_accum]<f_.n_out());

        casadi_assert_message(f.sparsity_in(i)==f.sparsity_out(output_accum[n_accum]),
                              "Input #" << i << " and output #" << output_accum[n_accum] <<
                              " must have matching sparsity. " <<
                              "Got " << f.sparsity_in(i).dim() << " and " <<
                              f.sparsity_out(output_accum[n_accum]).dim() << ".");
        n_accum++;
      }
    }
    casadi_assert(output_accum.size()==n_accum);
    casadi_assert(isMonotone(output_accum));
  }

  MapAccumInternal::~MapAccumInternal() {
  }

  void MapAccumInternal::init() {

    int num_in = f_.n_in(), num_out = f_.n_out();

    // Initialize the functions, get input and output sparsities
    // Input and output sparsities

    ibuf_.resize(num_in);
    obuf_.resize(num_out);

    for (int i=0;i<num_in;++i) {
      // Sparsity of the original input
      Sparsity in_sp = f_.input(i).sparsity();

      if (!input_accum_[i]) in_sp = repmat(in_sp, 1, n_);

      // Allocate space for input
      input(i) = DMatrix::zeros(in_sp);
    }

    for (int i=0;i<num_out;++i) {
      // Sparsity of the original output
      Sparsity out_sp = f_.output(i).sparsity();

      out_sp = repmat(out_sp, 1, n_);

      // Allocate space for output
      output(i) = DMatrix::zeros(out_sp);
    }

    step_in_.resize(num_in, 0);
    step_out_.resize(num_out, 0);

    for (int i=0;i<num_in;++i) {
      step_in_[i] = f_.input(i).nnz();
    }

    for (int i=0;i<num_out;++i) {
      step_out_[i] = f_.output(i).nnz();
    }

    // Call the initialization method of the base class
    FunctionInternal::init();

    // We cannot rely on the output buffer to use as accumulator:
    // a null-pointer may be passed as output
    nnz_accum_ = 0;
    for (int i=0;i<num_in;++i) {
      if (input_accum_[i]) nnz_accum_+= step_in_[i];
    }

    // We need double this space: one end to serve as input, one end to dump the output

    alloc_w(f_.sz_w()+2*nnz_accum_);
    alloc_iw(f_.sz_iw());
    alloc_arg(2*f_.sz_arg());
    alloc_res(2*f_.sz_res());
  }

  template<typename T, typename R>
  void MapAccumInternal::evalGen(const T** arg, T** res, int* iw, T* w,
                                 R reduction) {
    int num_in = f_.n_in(), num_out = f_.n_out();

    // Catch: must accomodate scenario where res[j] of the accumulator = 0.

    T* accum = w+f_.sz_w();

    const T** arg1 = arg+f_.sz_arg();
    T** res1 = res+f_.sz_res();

    // Copy the initial values to the accumulators
    for (int j=0; j<num_in; ++j) {
      if (input_accum_[j]) {
        copy(arg[j], arg[j]+step_in_[j], accum);
        accum+= step_in_[j];
      }
    }

    // Set the function accum inputs
    accum = w+f_.sz_w();
    for (int j=0; j<num_in; ++j) {
      if (input_accum_[j]) {
        arg1[j] = accum;
        accum += step_in_[j];
      }
    }

    for (int iter=0; iter<n_; ++iter) {

      int i = reverse_ ? n_-iter-1: iter;

      // Set the function non-accum inputs
      for (int j=0; j<num_in; ++j) {
        accum = w+f_.sz_w();
        if (!input_accum_[j]) {
          arg1[j] = (arg[j]==0) ? 0: arg[j]+i*step_in_[j];
        }
      }

      // Set the function outputs
      for (int j=0; j<num_out; ++j) {
        res1[j] = (res[j]==0) ? 0: res[j]+i*step_out_[j];
      }

      // Point the accumulator outputs to temporary storage
      accum = w+f_.sz_w()+nnz_accum_;
      for (int j=0; j<output_accum_.size(); ++j) {
        int jj = output_accum_[j];
        res1[jj] = accum;
        accum += step_out_[jj];
      }

      // Evaluate the function
      f_(arg1, res1, iw, w);

      // Copy the temporary storage to the accumulator
      copy(w+f_.sz_w()+nnz_accum_, w+f_.sz_w()+nnz_accum_*2, w+f_.sz_w());

      // Copy the accumulator to the global output ...
      accum = w+f_.sz_w();
      for (int j=0; j<output_accum_.size(); ++j) {
        int jj = output_accum_[j];
        // ... but beware of a null pointer
        if (res[jj]!=0) {
          copy(accum, accum+step_out_[jj], res[jj]+i*step_out_[jj]);
          accum += step_out_[jj];
        }
      }
    }
  }

  void MapAccumInternal::evalD(const double** arg, double** res,
                                int* iw, double* w) {
    evalGen(arg, res, iw, w, std::plus<double>());
  }

  void MapAccumInternal::evalSX(const SXElement** arg, SXElement** res,
                                int* iw, SXElement* w) {
    evalGen(arg, res, iw, w, std::plus<SXElement>());
  }

  void MapAccumInternal::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    evalGen(arg, res, iw, w, orop);
  }

  std::vector<MX> bisect(const MX& a, int b) {
    std::vector<int> idx(3, 0);
    idx[1] = b;
    idx[2] = a.size2();
    return horzsplit(a, idx);
  }

  Function MapAccumInternal
  ::getDerForward(const std::string& name, int nfwd, Dict& opts) {

    // Obtain forward mode of the primitive function
    /*  
        x1, y0 = f(x0, u0)
        x2, y1 = f(x1, u1)
        x3, y2 = f(x2, u2)

        F :    x0, [u0 u1 u2] -> [x1 x2 x3], [y0 y1 y2]

        df:    x, u, xp, y, x_dot, u_dot -> xp_dot, y_dot

    */

    Function df = f_.derForward(nfwd);

    /*
      Reverse mode of F looks like:

        F :    x0, [u0 u1 u2], [x1 x2 x3], [y0 y1 y2],
                    x0_dot, [u0_dot u1_dot u2_dot] -> 
                    [x1_dot x2_dot x3_dot], [y0_dot y1_dot y2_dot]

      In terms of the derivate of the primitive:
        

        df (x0, u0, x1, y0, x0_dot, u0_dot) -> x1_dot, y0_dot
        df (x1, u1, x2, y1, x1_dot, u1_dot) -> x2_dot, y1_dot   (A)
        df (x2, u2, x3, y2, x2_dot, u2_dot) -> x3_dot, y2_dot

      This maps to MapAccum.

    */

    // Construct the new MapAccum's input_accum
    /*  
        fb:  x, u, xp, y, x_dot, u_dot -> xp_dot, y_dot
        
        input_accum:
              0, 0, 0, 0,  1,      0
      
    */
    std::vector<bool> input_accum(n_in()+n_out(), false);
    for (int i=0;i<nfwd;++i) {
      input_accum.insert(input_accum.end(), input_accum_.begin(), input_accum_.end());
    }

    // Construct the new MapAccum's output_accum
    /*  
       fb:  x, u, xp, y, x_dot, u_dot -> xp_dot, y_dot
        
        output_accum:
              [0]
      
    */
    std::vector<int> output_accum;
    int offset = 0;
    for (int i=0;i<nfwd;++i) {
      for (int j=0;j<output_accum_.size();++j) {
        output_accum.push_back(offset+output_accum_[j]);
      }
      offset+= n_out();
    }

    // Construct the new MapAccum
    Function ma = df.mapaccum("map", n_, input_accum, output_accum, reverse_);

    /*

      Recall, the function we need to return looks like:

        F :    x0, [u0 u1 u2], [x1 x2 x3], [y0 y1 y2],
                    x0_dot, [u0_dot u1_dot u2_dot] -> 
                    [x1_dot x2_dot x3_dot], [y0_dot y1_dot y2_dot]

      While the just constructed MapAccum has the following signature:

        ma :   [x0 x1 x2], [u0 u1 u2], [x1 x2 x3], [y0 y1 y2],
                    x0_dot, [u0_dot u1_dot u2_dot] -> 
                    [x1_dot x2_dot x3_dot], [y0_dot y1_dot y2_dot]

      This requires an extra wrapper

    */

    std::vector<MX> ins = mx_in();
    std::vector<MX> outs = mx_out();

    std::vector<MX> f_der_ins;
    std::vector<MX> der_ins = ins;
    der_ins.insert(der_ins.end(), outs.begin(), outs.end());

    // For the inputs to the new mapaccum
    // which correspond to accumulators, supply these with the whole history
    /** Takes care of:  [x0 x1 x2]   */
    int i_output_accum=0;
    for (int i=0;i<n_in();++i) {
      if (input_accum_[i]) {
        // [x0 x1 x2 x3] -> all
        // [x0 x1 x2] -> in
        if (reverse_) {
          MX all = horzcat(outs[output_accum_[i_output_accum++]], ins[i]);
          MX in = bisect(all, ins[i].size2())[1];
          f_der_ins.push_back(in);
        } else {
          MX all = horzcat(ins[i], outs[output_accum_[i_output_accum++]]);
          MX in = bisect(all, ins[i].size2()*n_)[0];
          f_der_ins.push_back(in);
        }
      } else {
        f_der_ins.push_back(ins[i]);
      }
    }
    f_der_ins.insert(f_der_ins.end(), outs.begin(), outs.end());

    for (int i=0;i<nfwd;++i) {
      std::vector<MX> ins = mx_in();
      der_ins.insert(der_ins.end(), ins.begin(), ins.end());
      f_der_ins.insert(f_der_ins.end(), ins.begin(), ins.end());
    }

    // Construct the wrapper
    return MX::fun("df", der_ins, ma(f_der_ins));
  }

  Function MapAccumInternal
  ::getDerReverse(const std::string& name, int nadj, Dict& opts) {

    // Obtain Reverse mode of the primitive function
    /*  
        x1, y0 = f(x0, u0)
        x2, y1 = f(x1, u1)
        x3, y2 = f(x2, u2)

        F :    x0, [u0 u1 u2] -> [x1 x2 x3], [y0 y1 y2]

        fb:    x, u, xp, y, xp_bar, y_bar -> x_bar, u_bar

    */
    Function fb = f_.derReverse(nadj);

    /*
      Reverse mode of F looks like:

        F :    x0, [u0 u1 u2], [x1 x2 x3], [y0 y1 y2],
                    [X1_bar X2_bar X3_bar], [Y0_bar Y1_bar Y2_bar] -> 
                    X0_bar, [u0_bar, u1_bar u2_bar]

      In terms of the derivate of the primitive:
        

        fb (x2, u2, x3, y2,          X3_bar, y2_bar) -> x2_bar, u2_bar
        fb (x1, u1, x2, y1, x2_bar + X2_bar, y1_bar) -> x1_bar, u1_bar   (A)
        fb (x0, u0, x1, y0, x1_bar + X1_bar, y0_bar) -> x0_bar, u0_bar

        x0_bar -> X0_bar


      Let's define an extended version of the derivative of the primitive:

      fbX:  x, u, xp, y xp_bar, y_bar, X -> x_bar + X, u_bar

      With this defintion, (A) can be rewritten as:
        
         fbX (x2, u2, x3, y2, X3_bar         , y2_bar, X2_bar) -> x2_bar + X2_bar, u2_bar
         fbX (x1, u1, x2, y1, x2_bar + X2_bar, y1_bar, X1_bar) -> x1_bar + X1_bar, u2_bar
         fbX (x0, u0, x1, y0, x1_bar + X1_bar, y0_bar, 0)      -> x0_bar + 0,      u2_bar

      This maps to MapAccum.

    */

    // Extend the primitive function with extra inputs:
    // for each reverse direction,
    // add one extra input per accumulated input.

    /*  
        fbX:  x, u, xp, y, xp_bar, y_bar, X -> ...
      
    */
    std::vector<MX> ins  = f_.mx_in();

    // Copy nominal inputs
    std::vector<MX> f_der_ins = ins;
    std::vector<MX> der_ins   = ins;

    // Copy nominal outputs
    std::vector<MX> outs = f_.mx_out();
    der_ins.insert(der_ins.end(), outs.begin(), outs.end());
    f_der_ins.insert(f_der_ins.end(), outs.begin(), outs.end());

    for (int i=0;i<nadj;++i) {
      std::vector<MX> outs = f_.mx_out();
      // Copy seeds
      der_ins.insert(der_ins.end(), outs.begin(), outs.end());
      f_der_ins.insert(f_der_ins.end(), outs.begin(), outs.end());

      // Add extra seeds
      for (int j=0;j<output_accum_.size();++j) {
        MX s = MX::sym("x", f_.sparsity_out(output_accum_[j]));
        f_der_ins.push_back(s);
      }
    }

    // Construct the outputs of the extended primitive function.
    // The number of outputs is the same as the primitive function.
    // The outputs corresponding to an accumulator need a summation
    // with the above extra inputs.

    /*  
        fbX:  x, u, xp, y xp_bar, y_bar, X -> x_bar + X, u_bar
      
    */
    // Call the dervative of the primitive
    std::vector<MX> der_outs = fb(der_ins);

    std::vector<MX> f_der_outs;
    for (int i=0;i<nadj;++i) {

      // Obtain the sensitivities of direction i.
      f_der_outs.insert(f_der_outs.end(), der_outs.begin()+i*n_in(), der_outs.begin()+(i+1)*n_in());

      int i_output_accum = 0;
      for (int j=0;j<n_in();++j) {
        if (input_accum_[j]) {
          // Add the extra argument to the output
          int idx = n_in()+n_out()+i*(n_out()+output_accum_.size())+n_out()+i_output_accum;
          f_der_outs[i*n_in()+j] += f_der_ins[idx];
          i_output_accum ++;
        }
      }
    }

    Function fbX = MX::fun("f", f_der_ins, f_der_outs);

    // Construct the new MapAccum's input_accum
    /*  
        fbX:  x, u, xp, y, xp_bar, y_bar, X -> ...
        
        input_accum:
              0, 0, 0, 0,  1,      0,     0
      
    */
    std::vector<bool> input_accum(n_in()+n_out(), false);
    std::vector<bool> input_accum_rev(n_out(), false);
    for (int j=0;j<output_accum_.size();++j) {
      input_accum_rev[output_accum_[j]] = true;
    }

    for (int i=0;i<nadj;++i) {
      input_accum.insert(input_accum.end(), input_accum_rev.begin(), input_accum_rev.end());
      input_accum.insert(input_accum.end(), output_accum_.size(), false);
    }

    // Construct the new MapAccum's output_accum
    /*  
        fbX:  x, u, xp, y xp_bar, y_bar, X -> x_bar + X, u_bar
        
        output_accum:
              [0]
      
    */
    std::vector<int> output_accum_rev;
    for (int j=0;j<input_accum_.size();++j) {
      if (input_accum_[j]) output_accum_rev.push_back(j);
    }

    std::vector<int> output_accum;
    int offset = 0;
    for (int i=0;i<nadj;++i) {
      for (int j=0;j<output_accum_.size();++j) {
        output_accum.push_back(offset+output_accum_rev[j]);
      }
      offset+= n_in();
    }

    // Create the new MapAccum
    Function ma = fbX.mapaccum("map", n_, input_accum, output_accum, !reverse_);

    /*

      Recall, the function we need to return looks like:

                 F :    x0, [u0 u1 u2], [x1 x2 x3], [y0 y1 y2],
                    [X1_bar X2_bar X3_bar], [Y0_bar Y1_bar Y2_bar] -> 
                    X0_bar, [u0_bar, u1_bar u2_bar]

      While the just constructed MapAccum has the following signature:

        ma :   [x0 x1 x2], [u0 u1 u2], [x1 x2 x3], [y0 y1 y2],
                    X3_bar, [Y0_bar Y1_bar Y2_bar], [0 X1_bar X2_bar] -> 
                    [X0_bar X1_bar X2_bar], [u0_bar, u1_bar u2_bar]

      This requires an extra wrapper

    */

    // We need to provision the extra inputs here to the new mapaccum
    der_ins.clear();f_der_ins.clear();
    der_outs.clear();

    ins  = mx_in();
    outs = mx_out();

    // For the inputs to the new mapaccum
    // which correspond to accumulators, supply these with the whole history
    /** Takes care of:  [x0 x1 x2]   */
    der_ins.insert(der_ins.end(), ins.begin(), ins.end());
    der_ins.insert(der_ins.end(), outs.begin(), outs.end());
    int i_output_accum=0;
    for (int i=0;i<n_in();++i) {
      if (input_accum_[i]) {
        // [x0 x1 x2 x3] -> all
        // [x0 x1 x2] -> in
        if (reverse_) {
          MX all = horzcat(outs[output_accum_[i_output_accum++]], ins[i]);
          MX in = bisect(all, ins[i].size2())[1];
          f_der_ins.push_back(in);
        } else {
          MX all = horzcat(ins[i], outs[output_accum_[i_output_accum++]]);
          MX in = bisect(all, ins[i].size2()*n_)[0];
          f_der_ins.push_back(in);
        }
      } else {
        f_der_ins.push_back(ins[i]);
      }
    }
    f_der_ins.insert(f_der_ins.end(), outs.begin(), outs.end());

    // For the seeds to to the new mapaccum
    // which correspond to accumulators, supply
    // these with only the latest seed.
    // Also, supply the extra inputs
    /** Takes care of: X3_bar, [0 X1_bar X2_bar]  */
    for (int i=0;i<nadj;++i) {

      // Seed has the signature of the output
      std::vector<MX> outs = mx_out();

      // Pass seeds unchanged
      der_ins.insert(der_ins.end(), outs.begin(), outs.end());
      f_der_ins.insert(f_der_ins.end(), outs.begin(), outs.end());

      int i_output_accum=0;
      for (int j=0;j<n_in();++j) {
        if (input_accum_[j]) {
          if (reverse_) {
            // 0, [X1_bar X2_bar X3_bar] -> all
            MX all = horzcat(outs[output_accum_[i_output_accum]],
              DMatrix::zeros(ins[j].sparsity()));
            // [0 X1_bar X2_bar], X3_bar -> splits
            std::vector<MX> splits = bisect(all, ins[j].size2());

            // alter the seed corresponding to an accumulator
            f_der_ins[n_in()+n_out()+i*(n_out()+output_accum_.size())+output_accum_[i_output_accum]]
              = splits[0];
            // supply extra inputs
            f_der_ins.push_back(splits[1]);
          } else {
            // 0, [X1_bar X2_bar X3_bar] -> all
            MX all = horzcat(DMatrix::zeros(ins[j].sparsity()),
              outs[output_accum_[i_output_accum]]);
            // [0 X1_bar X2_bar], X3_bar -> splits
            std::vector<MX> splits = bisect(all, ins[j].size2()*n_);

            // alter the seed corresponding to an accumulator
            f_der_ins[n_in()+n_out()+i*(n_out()+output_accum_.size())+output_accum_[i_output_accum]]
              = splits[1];
            // supply extra inputs
            f_der_ins.push_back(splits[0]);
          }

          i_output_accum++;
        }
      }
    }

    // Call the new MapAccum with the newly constructed inputs
    der_outs = ma(f_der_ins);

    // Chop off the sensitivities corresponding to accumulator inputs
    /*
      The output of the new MapAccum delivers [X2_bar X1_bar X0_bar],
      while we only need to output X0_bar.
      
    */
    for (int i=0;i<nadj;++i) {
      for (int j=0;j<n_in();++j) {
        if (input_accum_[j]) {
          MX& x = der_outs[i*n_in()+j];
          if (reverse_) {
            std::vector<MX> r = bisect(x, f_.sparsity_in(j).size2()*(n_-1));
            x = r[1];
          } else {
            std::vector<MX> r = bisect(x, f_.sparsity_in(j).size2());
            x = r[0];
          }
        }
      }
    }

    // Construct the wrapper
    return MX::fun("df", der_ins, der_outs);
  }

  void MapAccumInternal::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void MapAccumInternal::generateBody(CodeGenerator& g) const {
    g.addAuxiliary(CodeGenerator::AUX_COPY_N);
    int num_in = f_.n_in(), num_out = f_.n_out();

    g.body << "  const real_t** arg1 = arg+" << f_.sz_arg() << ";" << endl
           << "  real_t** res1 = res + " << f_.sz_res() << ";" << endl;

    g.body << "  real_t* accum = w+" << f_.sz_w() << ";" << endl;

    // Copy the initial values to the accumulators
    int offset=0;
    for (int j=0; j<num_in; ++j) {
      if (input_accum_[j]) {
        g.body << "  copy_n(arg[" << j << "], " << step_in_[j] <<
          ", accum + " << offset << ");" << std::endl;
        g.body << "  arg1[" << j << "] = accum + " << offset << ";" << std::endl;
        offset+= step_in_[j];
      }
    }

    g.body << "  int i;" << endl;
    if (reverse_) {
      g.body << "  for (i=" << n_-1 << "; i>=0; --i) {" << endl;
    } else {
      g.body << "  for (i=0; i<" << n_ << "; ++i) {" << endl;
    }

      // Set the function non-accum inputs
      for (int j=0; j<num_in; ++j) {
        if (!input_accum_[j]) {
          g.body << "    arg1[" << j << "] = (arg[" << j << "]==0) ? 0: " <<
            "arg[" << j << "]+i*" << step_in_[j] << ";" << endl;
        }
      }

      // Set the function outputs
      for (int j=0; j<num_out; ++j) {
        g.body << "    res1[" << j << "] = (res[" << j << "]==0) ? 0: "
          "res[" << j << "]+i*" << step_out_[j] << ";" << endl;
      }

      // Point the accumulator outputs to temporary storage
      offset = 0;
      for (int j=0; j<output_accum_.size(); ++j) {
        int jj = output_accum_[j];
        g.body << "    res1[" << jj << "] = accum+" << nnz_accum_+offset << ";" << endl;
        offset += step_out_[jj];
      }

      // Evaluate the function
      g.body << "    " << g.call(f_, "arg1", "res1", "iw", "w") << ";" << endl;

      // Copy the temporary storage to the accumulator
      g.body << "    copy_n(accum+" << nnz_accum_ << ", " << nnz_accum_ << ", accum);" << endl;

      // Copy the accumulator to the global output ...
      offset = 0;
      for (int j=0; j<output_accum_.size(); ++j) {
        int jj = output_accum_[j];
        // ... but beware of a null pointer
        g.body << "    if (res[" << jj << "]!=0)" << endl;
        g.body << "      copy_n(accum+ " << offset << ", " << step_out_[jj] <<
          ", res[" << jj << "]+i*" << step_out_[jj] << ");" << endl;
        offset += step_out_[jj];
      }
    g.body << "  }" << endl;
  }

  inline string name(const Function& f) {
    if (f.isNull()) {
      return "NULL";
    } else {
      return f.name();
    }
  }

  void MapAccumInternal::print(ostream &stream) const {
    stream << "MapAccum(" << name(f_) << ", " << n_ << ")";
  }

} // namespace casadi
