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


#include "kernel_sum.hpp"

using namespace std;

namespace casadi {

  KernelSum::KernelSum(const std::string& name, const Function& f,
                                           const std::pair<int, int> & size,
                                           double r,
                                           int n)
    : FunctionInternal(name), f_(f), size_(size), r_(r), n_(n) {

    casadi_assert(n>=1);

    casadi_assert_message(n==1, "Vectorized form of KernelSum2D not yet implemented.");

    casadi_assert(f.n_in()>=2);
    casadi_assert(f.sparsity_in(0)==Sparsity::dense(2, 1));
    casadi_assert(f.sparsity_in(1)==Sparsity::dense(1, 1));
    casadi_assert(f.sparsity_in(2)==Sparsity::dense(2, 1));
  }

  KernelSum::~KernelSum() {
  }

  void KernelSum::init() {

    int num_in = f_.n_in(), num_out = f_.n_out();

    // Call the initialization method of the base class
    FunctionInternal::init();

    step_out_.resize(num_out, 0);

    for (int i=0; i<num_out; ++i) {
      step_out_[i] = f_.nnz_out(i);
    }

    // Allocate some space to evaluate each function to.
    nnz_out_ = 0;
    for (int i=0;i<num_out;++i) {
      nnz_out_+= step_out_[i];
    }

    alloc_w(f_.sz_w() + nnz_out_+3);
    alloc_iw(f_.sz_iw());
    alloc_arg(2*f_.sz_arg());
    alloc_res(2*f_.sz_res());

  }

  static bvec_t Orring(bvec_t x, bvec_t y) { return x | y; }

  void KernelSum::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem) {
    // First input is non-differentiable

    int num_in = f_.n_in(), num_out = f_.n_out();

    // Clear the accumulators
    bvec_t** sum = res;
    for (int k=0;k<num_out;++k) {
      if (sum[k]!=0) std::fill(sum[k], sum[k]+step_out_[k], 0);
    }

    const bvec_t** arg1 = arg+f_.sz_arg();
    bvec_t** res1 = res+f_.sz_res();

    // Everything except the first argument can be passed as-is to f
    std::copy(arg+1, arg+num_in, arg1+2);

    // The first argument will be the pixel coordinates p_i
    bvec_t* coord = w+f_.sz_w()+nnz_out_;
    arg1[0] = coord;

    // The second argument will be the pixel value v_i
    bvec_t* value = w+f_.sz_w()+nnz_out_+2;
    arg1[1] = value;

    bvec_t* temp_res = w+f_.sz_w();
    // Clear the temp_res storage space
    if (temp_res!=0) std::fill(temp_res, temp_res+nnz_out_, 0);

    // Set the function outputs
    for (int j=0; j<num_out; ++j) {
      // Make the function outputs end up in temp_res
      res1[j] = (res[j]==0)? 0: temp_res;
      temp_res+= step_out_[j];
    }

    // Set the coordinate input
    coord[0] = 0;
    coord[1] = 0;

    // Set the pixel value input
    value[0] = 0;

    // Evaluate the function
    f_(arg1, res1, iw, w, 0);

    // Sum results from temporary storage to accumulator
    for (int k=0;k<num_out;++k) {
      if (res1[k] && sum[k])
        std::transform(res1[k], res1[k]+step_out_[k], sum[k], sum[k], &Orring);
    }
  }

  void KernelSum::
  eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const {
    int num_in = f_.n_in(), num_out = f_.n_out();

    const double* V = arg[0];
    const double* X = arg[1];

    // Clear the accumulators
    double** sum = res;
    for (int k=0;k<num_out;++k) {
      if (sum[k]!=0) std::fill(sum[k], sum[k]+step_out_[k], 0);
    }

    const double** arg1 = arg+f_.sz_arg();
    double** res1 = res+f_.sz_res();

    // Everything except the first argument can be passed as-is to f
    std::copy(arg+1, arg+num_in, arg1+2);

    // The first argument will be the pixel coordinates p_i
    double* coord = w+f_.sz_w()+nnz_out_;
    arg1[0] = coord;

    // The second argument will be the pixel value v_i

    double* value = w+f_.sz_w()+nnz_out_+2;
    arg1[1] = value;

    double* temp_res = w+f_.sz_w();

    // Clear the temp_res storage space
    if (temp_res!=0) std::fill(temp_res, temp_res+nnz_out_, 0);

    // Set the function outputs
    for (int j=0; j<num_out; ++j) {
      // Make the function outputs end up in temp_res
      res1[j] = (res[j]==0)? 0: temp_res;
      temp_res+= step_out_[j];
    }

    //     ---> j,v
    //   |
    //   v  i,u
    int u = round(X[0]);
    int v = round(X[1]);

    for (int i = -round(r_);i<= round(r_);++i) {
      for (int j = -round(r_);j<= round(r_);++j) {

        // Set the coordinate input
        coord[0] = u+i;
        coord[1] = v+j;

        int I = min(max(u+i, 0), size_.first-1);
        int J = min(max(v+j, 0), size_.second-1);

        // Set the pixel value input
        value[0] = V[I+J*size_.first];

        // Evaluate the function
        f_(arg1, res1, iw, w, 0);

        // Sum results from temporary storage to accumulator
        for (int k=0;k<num_out;++k) {
          if (res1[k] && sum[k])
            std::transform(res1[k], res1[k]+step_out_[k], sum[k], sum[k], std::plus<double>());
        }
      }
    }
  }

  Function KernelSum
  ::get_forward(const std::string& name, int nfwd, Dict& opts) {

    /* Write KernelSum2D in linear form:
    *  
    *    S = F(V, X)  = sum_i  f ( P_i, v_i, X)
    *
    *  With a slight abuse of notation, we have:
    *
    *    S_dot = sum_i  f_forward ( P_i, v_i, X, X_dot)
    *
    *  The forward mode of KernelSum is another KernelSum.
    *  There is a bit of houskeeping in selecting the correct inputs/outputs.
    *
    */

    /* More exactly, the forward mode of the primitive is
    * fd( P_i, v_i, X, S, P_i_dot, v_i_dot, X_dot)
    *
    * we need to bring this in the form 
    *
    *     f_forward ( P_i, v_i, X, X_dot)
    *
    */
    Function fd = f_.forward(nfwd);

    vector<MX> f_inputs(f_.n_in());
    for (int i=0; i<f_inputs.size(); ++i) {
      f_inputs[i] = f_->FunctionInternal::mx_in(i); // NOTE(@jaeandersson): Flawed design
    }
    std::vector<MX> f_outputs  = f_.mx_out();

    std::vector<MX> fd_inputs   = f_inputs;

    // Create nodes for the nominal output (S)
    std::vector<MX> f_call_out  = f_(f_inputs);

    fd_inputs.insert(fd_inputs.end(), f_call_out.begin(), f_call_out.end());
    for (int i=0;i<nfwd;++i) {
      // Pad with blanks: we don't consider P_i_dot and v_i_dot
      fd_inputs.push_back(MX());
      fd_inputs.push_back(MX());

      vector<MX> inputs(f_.n_in());
      for (int i=0; i<inputs.size(); ++i) {
        inputs[i] = f_->FunctionInternal::mx_in(i); // NOTE(@jaeandersson): Flawed design
      }

      fd_inputs.insert(fd_inputs.end(), inputs.begin()+2, inputs.end());
      f_inputs.insert(f_inputs.end(), inputs.begin()+2, inputs.end());
    }

    Function f_forward("f", f_inputs, fd(fd_inputs));

    Function ret = f_forward.kernel_sum(name, size_, r_, n_);

    /* Furthermore, we need to return something of signature
    *  der(V,X,S,V_dot,X_dot)
    *
    */
    std::vector<MX> der_inputs = mx_in();
    std::vector<MX> ret_inputs = der_inputs;

    std::vector<MX> outputs = mx_out();
    der_inputs.insert(der_inputs.end(), outputs.begin(), outputs.end());

    for (int i=0;i<nfwd;++i) {
      // Construct dummy matrix for capturing the V_dot argument.
      der_inputs.push_back(MX::sym("x", Sparsity(size_) ));
      std::vector<MX> inputs = mx_in();
      der_inputs.insert(der_inputs.end(), inputs.begin()+1, inputs.end());
      ret_inputs.insert(ret_inputs.end(), inputs.begin()+1, inputs.end());
    }

    Function der("f", der_inputs, ret(ret_inputs), opts);

    return der;
  }

  Function KernelSum
  ::get_reverse(const std::string& name, int nadj, Dict& opts) {
    /* Write KernelSum2D in linear form:
    *  
    *    S = F(V, X)  = sum_i  f ( P_i, v_i, X)
    *
    *  With a slight abuse of notation, we have:
    *
    *    X_bar = sum_i  f_reverse ( P_i, v_i, X, S_bar)
    *
    *  The reverse mode of KernelSum is another KernelSum.
    *  There is a bit of houskeeping in selecting the correct inputs/outputs.
    *
    */

    int num_in = f_.n_in(), num_out = f_.n_out();

    /* More exactly, the reverse mode of the primitive is
    * fd( P_i, v_i, X, S, S_bar) -> P_i_bar, v_i_bar, X_bar
    *
    * we need to bring this in the form 
    *
    *     f_reverse ( P_i, v_i, X, S_bar) -> X_bar
    *
    *
    *
    */
    Function fd = f_.reverse(nadj);

    // NOTE(@jaeandersson): Flawed design
    std::vector<MX> f_inputs   = f_->FunctionInternal::mx_in();
    std::vector<MX> f_outputs  = f_.mx_out();

    std::vector<MX> fd_inputs   = f_inputs;
    // Create nodes for the nominal output (S)
    std::vector<MX> f_call_out  = f_(f_inputs);

    fd_inputs.insert(fd_inputs.end(), f_call_out.begin(), f_call_out.end());

    for (int i=0;i<nadj;++i) {
      std::vector<MX> outputs   = f_.mx_out();
      fd_inputs.insert(fd_inputs.end(), outputs.begin(), outputs.end());
      f_inputs.insert(f_inputs.end(), outputs.begin(), outputs.end());
    }

    std::vector<MX> fd_outputs = fd(fd_inputs);

    // Drop the outputs we do not need: P_i_bar, v_i_bar
    f_outputs.clear();
    int offset = 2;
    for (int i=0;i<nadj;++i) {
      f_outputs.insert(f_outputs.end(), fd_outputs.begin()+offset,
        fd_outputs.begin()+offset+f_.n_in()-2);
      offset+= f_.n_in();
    }

    Function f_reverse("f", f_inputs, f_outputs);

    Function kn = f_reverse.kernel_sum(name, size_, r_, n_);

    /* Furthermore, we need to return something of signature
    *  der(V,X,S,S_bar) -> V_bar, X_bar
    *
    */

    std::vector<MX> ret_inputs = mx_in();
    std::vector<MX> kn_inputs = ret_inputs;
    for (int i=0;i<num_out;++i) {
      // Dummy symbols for the nominal outputs (S)
      MX output = MX::sym("x", Sparsity(f_.sparsity_out(i).size()));
      ret_inputs.push_back(output);
    }
    for (int i=0;i<nadj;++i) {
      std::vector<MX> outputs = mx_out();
      ret_inputs.insert(ret_inputs.end(), outputs.begin(), outputs.end());
      kn_inputs.insert(kn_inputs.end(), outputs.begin(), outputs.end());
    }

    std::vector<MX> kn_outputs = kn(kn_inputs);

    std::vector<MX> ret_outputs;
    offset = 0;
    for (int i=0;i<nadj;++i) {
      // Use sparse zero as V_bar output
      ret_outputs.push_back(MX::zeros(Sparsity(size_)));
      ret_outputs.insert(ret_outputs.end(), kn_outputs.begin()+offset,
        kn_outputs.begin()+offset+num_in-2);
      offset+= num_in-2;
    }

    Function ret(name, ret_inputs, ret_outputs, opts);

    return ret;
  }

  void KernelSum::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void KernelSum::generateBody(CodeGenerator& g) const {
    g.addAuxiliary(CodeGenerator::AUX_COPY);
    g.addAuxiliary(CodeGenerator::AUX_FILL);
    g.addAuxiliary(CodeGenerator::AUX_AXPY);
    int num_in = f_.n_in(), num_out = f_.n_out();

    g.body << "  const real_t* V = arg[0];" << std::endl;
    g.body << "  const real_t* X = arg[1];" << std::endl;

    // Clear the accumulators
    g.body << "  real_t** sum = res;" << std::endl;
    for (int k=0;k<num_out;++k) {
      g.body << "  if (sum[" << k << "]!=0) fill(sum[" << k << "], "
        << step_out_[k] << ", 0);" << std::endl;
    }

    g.body << "  const real_t** arg1 = arg+" << f_.sz_arg() << ";" << std::endl;
    g.body << "  real_t** res1 = res+" << f_.sz_res() << ";" << std::endl;

    // Everything except the first argument can be passed as-is to f
    g.body << "  int ii;" << std::endl;
    g.body << "  for(ii=0;ii<" << num_in -1<< ";++ii) arg1[2+ii] = arg[1+ii];" << std::endl;

    // The first argument will be the pixel coordinates p_i
    g.body << "  real_t* coord = w+" << f_.sz_w()+nnz_out_ << ";" << std::endl;
    g.body << "  arg1[0] = coord;" << std::endl;

    // The second argument will be the pixel value v_i

    g.body << "  real_t* value = w+" << f_.sz_w()+nnz_out_+2 << ";" << std::endl;
    g.body << "  arg1[1] = value;" << std::endl;

    g.body << "  real_t* temp_res = w+" << f_.sz_w() << ";" << std::endl;

    // Clear the temp_res storage space
    g.body << "  if (temp_res!=0) fill(temp_res, " << nnz_out_ << ", 0);" << std::endl;

    // Set the function outputs
    int offset = 0;
    for (int j=0; j<num_out; ++j) {
      // Make the function outputs end up in temp_res
      g.body << "  res1[" << j << "] = (res[" << j << "]==0)? 0: temp_res + "
        << offset << ";" << std::endl;
      offset+= step_out_[j];
    }

    //     ---> j,v
    //   |
    //   v  i,u
    g.body << "  int u = round(X[0]);" << std::endl;
    g.body << "  int v = round(X[1]);" << std::endl;

    g.body << "  int i,j;" << std::endl;
    g.body << "  for (i = -" << round(r_) <<";i<=" <<  round(r_) << ";++i) {" << std::endl;
    g.body << "    for (j = -" << round(r_) <<";j<= " << round(r_) <<";++j) {" << std::endl;

    // Set the coordinate input
    g.body << "      int I = u+i;" << std::endl;
    g.body << "      int J = v+j;" << std::endl;

    g.body << "      coord[0] = I;" << std::endl;
    g.body << "      coord[1] = J;" << std::endl;

    g.body << "      if (I<0) I=0;" << std::endl;
    g.body << "      if (J<0) J=0;" << std::endl;

    g.body << "      if (I>" << size_.first-1 << ") I=" << size_.first-1 << ";" << std::endl;
    g.body << "      if (J>" << size_.second-1 << ") J=" << size_.second-1 << ";" << std::endl;

    // Set the pixel value input
    g.body << "      value[0] = V[I+J*" << size_.first << "];" << std::endl;

    // Evaluate the function
    g.body << "      if (" << g(f_, "arg1", "res1", "iw", "w", "0") << ") return 1;" << endl;

    // Sum results from temporary storage to accumulator
    for (int k=0; k<num_out; ++k) {
      g.body << "      if (res1[" << k << "] && sum[" << k << "])" << endl
             << "       axpy(" << step_out_[k] << ",1," <<
                          "res1["<< k << "],sum[" << k << "]);" << endl;
    }
    g.body << "    }" << std::endl;
    g.body << "  }" << std::endl;

  }

  void KernelSum::print(ostream &stream) const {
    stream << "KernelSum2D(" << f_.name() << ", " << n_ << ")";
  }

} // namespace casadi
