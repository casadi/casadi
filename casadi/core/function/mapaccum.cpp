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


#include "mapaccum.hpp"

using namespace std;

namespace casadi {

    Function Mapaccum::create(const std::string& name,
      Function& f, int n, const vector<int>& accum_in, const vector<int>& accum_out,
        const Dict& opts, bool reverse) {

      casadi_assert(inBounds(accum_in, f.n_in()) && isUnique(accum_in));
      casadi_assert(inBounds(accum_out, f.n_out()) && isUnique(accum_out));

      casadi_assert(accum_in.size()==accum_out.size());
      int n_accum=accum_in.size();

      // Check that sparsities match
      for (int i=0;i<n_accum;++i) {
        casadi_assert_message(f.sparsity_in(accum_in[i])==f.sparsity_out(accum_out[i]),
                                "Input #" << accum_in[i] << " and output #" << accum_out[i] <<
                                " must have matching sparsity. " <<
                                "Got " << f.sparsity_in(accum_in[i]).dim() << " and " <<
                                f.sparsity_out(accum_out[i]).dim() << ".");
      }

      if (accum_in==range(n_accum) && accum_out==range(n_accum)) {
        // No need to reorder
        Function ret;
        ret.assignNode(new Mapaccum(name, f, n, n_accum, reverse));
        ret->construct(opts);
        return ret;
      } else {
        // Need to do some reordering
        std::vector<int> temp_in = complement(accum_in, f.n_in());
        std::vector<int> order_in = accum_in;
        order_in.insert(order_in.end(), temp_in.begin(), temp_in.end());

        std::vector<int> temp_out = complement(accum_out, f.n_out());
        std::vector<int> order_out = accum_out;
        order_out.insert(order_out.end(), temp_out.begin(), temp_out.end());

        Function fr = f.slice(order_in, order_out);
        Function ma;
        ma.assignNode(new Mapaccum(name, fr, n, n_accum, reverse));
        ma->construct(opts);

        std::vector<int> order_in_inv = lookupvector(order_in, f.n_in());
        std::vector<int> order_out_inv = lookupvector(order_out, f.n_out());

        return ma.slice(order_in_inv, order_out_inv, opts);
      }


  }

  Mapaccum::Mapaccum(const std::string& name, const Function& f, int n, int n_accum, bool reverse)
    : FunctionInternal(name), f_(f), n_(n),
      n_accum_(n_accum), reverse_(reverse) {

    casadi_assert(n>0);
    casadi_assert(f.n_in()>=1);
    casadi_assert(f.n_out()>=1);

    for (int i=0;i<n_accum_;++i) {
      casadi_assert_message(f.sparsity_in(i)==f.sparsity_out(i),
                            "Input and output #" << i <<
                            " must have matching sparsity. " <<
                            "Got " << f.sparsity_in(i).dim() << " and " <<
                            f.sparsity_out(i).dim() << ".");
    }
  }

  Mapaccum::~Mapaccum() {
  }

  void Mapaccum::init(const Dict& opts) {

    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    int num_in = f_.n_in(), num_out = f_.n_out();
    step_in_.resize(num_in, 0);
    step_out_.resize(num_out, 0);
    for (int i=0; i<num_in; ++i) step_in_[i] = f_.nnz_in(i);
    for (int i=0; i<num_out; ++i) step_out_[i] = f_.nnz_out(i);

    // We cannot rely on the output buffer to use as accumulator:
    // a null-pointer may be passed as output
    nnz_accum_ = 0;
    for (int i=0; i<n_accum_; ++i) {
      nnz_accum_+= step_in_[i];
    }

    // We need double this space: one end to serve as input, one end to dump the output
    alloc_w(f_.sz_w()+2*nnz_accum_);
    alloc_iw(f_.sz_iw());
    alloc_arg(2*f_.sz_arg());
    alloc_res(2*f_.sz_res());
  }

  template<typename T, typename R>
  void Mapaccum::evalGen(const T** arg, T** res, int* iw, T* w, R reduction) const {
    int num_in = f_.n_in(), num_out = f_.n_out();

    // Catch: must accomodate scenario where res[j] of the accumulator = 0.

    T* accum = w+f_.sz_w();

    const T** arg1 = arg+f_.sz_arg();
    T** res1 = res+f_.sz_res();

    // Copy the initial values to the accumulators
    for (int j=0; j<n_accum_; ++j) {
      copy(arg[j], arg[j]+step_in_[j], accum);
      accum+= step_in_[j];
    }

    // Set the function accum inputs
    accum = w+f_.sz_w();
    for (int j=0; j<n_accum_; ++j) {
      arg1[j] = accum;
      accum += step_in_[j];
    }

    for (int iter=0; iter<n_; ++iter) {

      int i = reverse_ ? n_-iter-1: iter;

      // Set the function non-accum inputs
      for (int j=n_accum_; j<num_in; ++j) {
        arg1[j] = (arg[j]==0) ? 0: arg[j]+i*step_in_[j];
      }

      // Set the function outputs
      for (int j=0; j<num_out; ++j) {
        res1[j] = (res[j]==0) ? 0: res[j]+i*step_out_[j];
      }

      // Point the accumulator outputs to temporary storage
      accum = w+f_.sz_w()+nnz_accum_;
      for (int j=0; j<n_accum_; ++j) {
        res1[j] = accum;
        accum += step_out_[j];
      }

      // Evaluate the function
      f_(arg1, res1, iw, w, 0);

      // Copy the temporary storage to the accumulator
      copy(w+f_.sz_w()+nnz_accum_, w+f_.sz_w()+nnz_accum_*2, w+f_.sz_w());

      // Copy the accumulator to the global output ...
      accum = w+f_.sz_w();
      for (int j=0; j<n_accum_; ++j) {
        // ... but beware of a null pointer
        if (res[j]!=0) {
          copy(accum, accum+step_out_[j], res[j]+i*step_out_[j]);
          accum += step_out_[j];
        }
      }
    }
  }

  void Mapaccum::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    evalGen(arg, res, iw, w, std::plus<double>());
  }

  void Mapaccum::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    evalGen(arg, res, iw, w, std::plus<SXElem>());
  }

  void Mapaccum::sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    evalGen(arg, res, iw, w, orop);
  }

  void Mapaccum::sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {

    int num_in = f_.n_in(), num_out = f_.n_out();

    bvec_t** arg1 = arg+f_.sz_arg();
    bvec_t** res1 = res+f_.sz_res();

    bvec_t* accum = w+f_.sz_w();

    /**

      [ | | | ]

      |    _______
      L-> |       | -> [ | | | ]
          |       |
       -> |       | -> temp
      acc |_______| -> [ | | | ]


        */


    // Initialize the accumulator
    // with reverse seeds of accumulated output
    accum = w+f_.sz_w();
    for (int j=0; j<n_accum_; ++j) {
      if (res[j]==0) {
        fill(accum, accum+step_out_[j], 0);
      } else {
        if (reverse_) {
          copy(res[j], res[j]+step_out_[j], accum);
        } else {
          copy(res[j]+(n_-1)*step_out_[j], res[j]+n_*step_out_[j], accum);
        }
      }
      accum += step_out_[j];
    }

    // Point the function accum inputs to the accumulator
    accum = w+f_.sz_w();
    for (int j=0; j<n_accum_; ++j) {
      arg1[j] = accum;
      accum += step_in_[j];
    }

    for (int iter=n_-1; iter>=0; iter--) {

      int i = reverse_ ? n_-iter-1: iter;
      int i_prev = reverse_ ? n_-iter: iter-1;

      // Copy the accumulator to the temporary storage
      copy(w+f_.sz_w(), w+f_.sz_w()+nnz_accum_, w+f_.sz_w()+nnz_accum_);

      accum = w+f_.sz_w();
      // Copy the last seed to the accumulator
      for (int j=0; j<n_accum_; ++j) {
        // ... but beware of a null pointer
        if (res[j]!=0) {
          if (iter!=0) {
            copy(res[j]+i_prev*step_out_[j], res[j]+(i_prev+1)*step_out_[j], accum);
          }
          accum += step_out_[j];
        }

      }

      // For the start of the chain, use the reverse sensitivity
      if (iter==0) {
        accum = w+f_.sz_w();
        for (int j=0; j<n_accum_; ++j) {
          copy(arg[j], arg[j]+step_in_[j], accum);
          accum+= step_in_[j];
        }
      }

      // Set the function non-accum inputs
      for (int j=n_accum_; j<num_in; ++j) {
        accum = w+f_.sz_w();
        // non-accum inputs
        arg1[j] = (arg[j]==0) ? 0: arg[j]+i*step_in_[j];
      }

      // Set the function outputs
      for (int j=0; j<num_out; ++j) {
        res1[j] = (res[j]==0) ? 0: res[j]+i*step_out_[j];
      }

      // Point the accumulator outputs to temporary storage
      accum = w+f_.sz_w()+nnz_accum_;
      for (int j=0; j<n_accum_; ++j) {
        res1[j] = accum;
        accum += step_out_[j];
      }

      // Evaluate the function
      f_->sp_rev(arg1, res1, iw, w, 0);

    }

    accum = w+f_.sz_w();
    // Copy the accumulated result to the reverse sensitivities
    for (int j=0; j<n_accum_; ++j) {
      copy(accum, accum+step_in_[j], arg[j]);
      accum+= step_in_[j];
    }

    // Reset all seeds
    for (int j=0; j<num_out; ++j) {
      if (res[j]!=0) {
        fill(res[j], res[j]+f_.nnz_out(j), bvec_t(0));
      }
    }
  }

  std::vector<MX> bisect(const MX& a, int b) {
    std::vector<int> idx(3, 0);
    idx[1] = b;
    idx[2] = a.size2();
    return horzsplit(a, idx);
  }

  Function Mapaccum
  ::get_forward(const std::string& name, int nfwd, Dict& opts) {
    return get_forward_new(name, nfwd, opts);
  }

  Function Mapaccum
  ::get_forward_old(const std::string& name, int nfwd, Dict& opts) {

    // Obtain forward mode of the primitive function
    /*
        x1, y0 = f(x0, u0)
        x2, y1 = f(x1, u1)
        x3, y2 = f(x2, u2)

        F :    x0, [u0 u1 u2] -> [x1 x2 x3], [y0 y1 y2]

        df:    x, u, xp, y, x_dot, u_dot -> xp_dot, y_dot

    */

    Function df = f_.forward(nfwd);

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
    std::vector<int> accum_in;
    int offset = n_in()+n_out();
    for (int i=0;i<nfwd;++i) {
      for (int j=0;j<n_accum_;++j) accum_in.push_back(offset+j);
      offset+= n_in();
    }

    // Construct the new MapAccum's output_accum
    /*
       fb:  x, u, xp, y, x_dot, u_dot -> xp_dot, y_dot

        output_accum:
              [0]

    */
    std::vector<int> accum_out;
    offset = 0;
    for (int i=0;i<nfwd;++i) {
      for (int j=0;j<n_accum_;++j) accum_out.push_back(offset+j);
      offset+= n_out();
    }

    // Construct the new MapAccum
    Function ma = Mapaccum::create(name, df, n_, accum_in, accum_out,
                                   derived_options(), reverse_);

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
    for (int i=0;i<n_accum_;++i) {
      // [x0 x1 x2 x3] -> all
      // [x0 x1 x2] -> in
      if (reverse_) {
        MX all = horzcat(outs[i], ins[i]);
        MX in = bisect(all, ins[i].size2())[1];
        f_der_ins.push_back(in);
      } else {
        MX all = horzcat(ins[i], outs[i]);
        MX in = bisect(all, ins[i].size2()*n_)[0];
        f_der_ins.push_back(in);
      }
    }
    for (int i=n_accum_;i<n_in();++i) f_der_ins.push_back(ins[i]);
    f_der_ins.insert(f_der_ins.end(), outs.begin(), outs.end());

    for (int i=0;i<nfwd;++i) {
      std::vector<MX> ins = mx_in();
      der_ins.insert(der_ins.end(), ins.begin(), ins.end());
      f_der_ins.insert(f_der_ins.end(), ins.begin(), ins.end());
    }

    // Construct the wrapper
    return Function("df", der_ins, ma(f_der_ins), opts);
  }

  Function Mapaccum
  ::get_reverse(const std::string& name, int nadj, Dict& opts) {
    return get_reverse_new(name, nadj, opts);
  }

  Function Mapaccum
  ::get_reverse_old(const std::string& name, int nadj, Dict& opts) {

    // Obtain Reverse mode of the primitive function
    /*
        x1, y0 = f(x0, u0)
        x2, y1 = f(x1, u1)
        x3, y2 = f(x2, u2)

        F :    x0, [u0 u1 u2] -> [x1 x2 x3], [y0 y1 y2]

        fb:    x, u, xp, y, xp_bar, y_bar -> x_bar, u_bar

    */
    Function fb = f_.reverse(nadj);

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
      for (int j=0;j<n_accum_;++j) {
        MX s = MX::sym("x", f_.sparsity_out(j));
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
      for (int j=0;j<n_accum_;++j) {
        // Add the extra argument to the output
        int idx = n_in()+n_out()+i*(n_out()+n_accum_)+n_out()+i_output_accum;
        f_der_outs[i*n_in()+j] += f_der_ins[idx];
        i_output_accum ++;
      }
    }

    Function fbX("f", f_der_ins, f_der_outs);

    // Construct the new MapAccum's input_accum
    /*
        fbX:  x, u, xp, y, xp_bar, y_bar, X -> ...

        input_accum:
              0, 0, 0, 0,  1,      0,     0

    */
    std::vector<int> accum_in;
    int offset = n_in()+n_out();
    for (int i=0;i<nadj;++i) {
      for (int j=0;j<n_accum_;++j) accum_in.push_back(offset+j);
      offset+= n_out()+n_accum_;
    }

    // Construct the new MapAccum's output_accum
    /*
        fbX:  x, u, xp, y xp_bar, y_bar, X -> x_bar + X, u_bar

        output_accum:
              [0]

    */
    std::vector<int> accum_out;
    offset = 0;
    for (int i=0;i<nadj;++i) {
      for (int j=0;j<n_accum_;++j) accum_out.push_back(offset+j);
      offset+= n_in();
    }

    // Create the new MapAccum
    Function ma = Mapaccum::create(name, fbX, n_, accum_in, accum_out,
      derived_options(), !reverse_);
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
    for (int i=0;i<n_accum_;++i) {
      // [x0 x1 x2 x3] -> all
      // [x0 x1 x2] -> in
      if (reverse_) {
        MX all = horzcat(outs[i], ins[i]);
        MX in = bisect(all, ins[i].size2())[1];
        f_der_ins.push_back(in);
      } else {
        MX all = horzcat(ins[i], outs[i]);
        MX in = bisect(all, ins[i].size2()*n_)[0];
        f_der_ins.push_back(in);
      }
    }
    for (int i=n_accum_;i<n_in();++i) f_der_ins.push_back(ins[i]);
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

      for (int j=0;j<n_accum_;++j) {
        if (reverse_) {
          // 0, [X1_bar X2_bar X3_bar] -> all
          MX all = horzcat(outs[j], DM::zeros(ins[j].sparsity()));
          // [0 X1_bar X2_bar], X3_bar -> splits
          std::vector<MX> splits = bisect(all, ins[j].size2());

          // alter the seed corresponding to an accumulator
          f_der_ins[n_in()+n_out()+i*(n_out()+n_accum_)+j]
            = splits[0];
          // supply extra inputs
          f_der_ins.push_back(splits[1]);
        } else {
          // 0, [X1_bar X2_bar X3_bar] -> all
          MX all = horzcat(DM::zeros(ins[j].sparsity()), outs[j]);
          // [0 X1_bar X2_bar], X3_bar -> splits
          std::vector<MX> splits = bisect(all, ins[j].size2()*n_);

          // alter the seed corresponding to an accumulator
          f_der_ins[n_in()+n_out()+i*(n_out()+n_accum_)+j]
            = splits[1];
          // supply extra inputs
          f_der_ins.push_back(splits[0]);
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
      for (int j=0;j<n_accum_;++j) {
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

    // Construct the wrapper
    return Function("df", der_ins, der_outs, opts);
  }

  void Mapaccum::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void Mapaccum::generateBody(CodeGenerator& g) const {
    g.addAuxiliary(CodeGenerator::AUX_COPY);
    int num_in = f_.n_in(), num_out = f_.n_out();

    g.body << "  const real_t** arg1 = arg+" << f_.sz_arg() << ";" << endl
           << "  real_t** res1 = res + " << f_.sz_res() << ";" << endl;

    g.body << "  real_t* accum = w+" << f_.sz_w() << ";" << endl;

    // Copy the initial values to the accumulators
    int offset=0;
    for (int j=0; j<n_accum_; ++j) {
      g.body << "  copy(arg[" << j << "], " << step_in_[j] <<
        ", accum + " << offset << ");" << std::endl;
      g.body << "  arg1[" << j << "] = accum + " << offset << ";" << std::endl;
      offset+= step_in_[j];
    }

    g.body << "  int i;" << endl;
    if (reverse_) {
      g.body << "  for (i=" << n_-1 << "; i>=0; --i) {" << endl;
    } else {
      g.body << "  for (i=0; i<" << n_ << "; ++i) {" << endl;
    }

      // Set the function non-accum inputs
      for (int j=n_accum_; j<num_in; ++j) {
        g.body << "    arg1[" << j << "] = (arg[" << j << "]==0) ? 0: " <<
          "arg[" << j << "]+i*" << step_in_[j] << ";" << endl;
      }

      // Set the function outputs
      for (int j=0; j<num_out; ++j) {
        g.body << "    res1[" << j << "] = (res[" << j << "]==0) ? 0: "
          "res[" << j << "]+i*" << step_out_[j] << ";" << endl;
      }

      // Point the accumulator outputs to temporary storage
      offset = 0;
      for (int j=0; j<n_accum_; ++j) {
        g.body << "    res1[" << j << "] = accum+" << nnz_accum_+offset << ";" << endl;
        offset += step_out_[j];
      }

      // Evaluate the function
      g.body << "    " << g(f_, "arg1", "res1", "iw", "w") << ";" << endl;

      // Copy the temporary storage to the accumulator
      g.body << "    copy(accum+" << nnz_accum_ << ", " << nnz_accum_ << ", accum);" << endl;

      // Copy the accumulator to the global output ...
      offset = 0;
      for (int j=0; j<n_accum_; ++j) {
        // ... but beware of a null pointer
        g.body << "    if (res[" << j << "]!=0)" << endl;
        g.body << "      copy(accum+ " << offset << ", " << step_out_[j] <<
          ", res[" << j << "]+i*" << step_out_[j] << ");" << endl;
        offset += step_out_[j];
      }
    g.body << "  }" << endl;
  }

  void Mapaccum::print(ostream &stream) const {
    stream << "MapAccum(" << f_.name() << ", " << n_ << ")";
  }

} // namespace casadi
