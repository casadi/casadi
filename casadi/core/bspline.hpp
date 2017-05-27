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


#ifndef CASADI_BSPLINE_HPP
#define CASADI_BSPLINE_HPP

#include "function_internal.hpp"

namespace casadi {


    /** Base class for BSpline evaluators
    *
    *
    *
    */
    class BSplineCommon : public FunctionInternal {
      public:

        BSplineCommon(const std::string &name, const std::vector<double>& knots,
          const std::vector<int>& offset, const std::vector<int>& degree, int m);

        /** \brief  Initialize */
        void init(const Dict& opts) override;

        ///@{
        /** \brief Number of function inputs and outputs */
        size_t get_n_in() override { return 1; }
        size_t get_n_out() override { return 1; }
        ///@}

        ///@{
        /** \brief Options */
        static Options options_;
        const Options& get_options() const override { return options_;}
        ///@}

        static void from_knots(const std::vector< std::vector<double> >& knots,
          std::vector<int>& offset, std::vector<double>& stacked);

        std::vector<int> lookup_mode_;
        std::vector<double> knots_;
        std::vector<int> offset_;
        std::vector<int> degree_;
        std::vector<int> strides_;

        std::vector<int> coeffs_dims_;
        int coeffs_size_;
        int m_;
    };

    class BSpline : public BSplineCommon {
    public:
      static Function create(const std::string &name,
        const std::vector< std::vector<double> >& knots,
        const std::vector<double>& coeffs, const std::vector<int>& degree, int m=1,
        const Dict& opts=Dict());

      BSpline(const std::string &name, const std::vector<double>& knots,
        const std::vector<int>& offset, const std::vector<double>& coeffs,
        const std::vector<int>& degree, int m);

      /** \brief  Destructor */
      ~BSpline() override {}

      /// @{
      /** \brief Sparsities of function inputs and outputs */
      Sparsity get_sparsity_in(int i) override { return Sparsity::dense(offset_.size()-1); }
      Sparsity get_sparsity_out(int i) override { return Sparsity::dense(m_, 1); }
      /// @}

      /** \brief  Initialize */
      void init(const Dict& opts) override;

      /** \brief  Evaluate numerically, work vectors given */
      void eval(void* mem, const double** arg, double** res, int* iw, double* w) const override;

      ///@{
      /** \brief Generate a function that calculates \a nfwd forward derivatives */
      Function get_forward(const std::string& name, int nfwd,
                                   const std::vector<std::string>& i_names,
                                   const std::vector<std::string>& o_names,
                                   const Dict& opts) const override;
      int get_n_forward() const override { return 64;}
      ///@}

      ///@{
      /** \brief Generate a function that calculates \a nadj adjoint derivatives */
      Function get_reverse(const std::string& name, int nadj,
                                   const std::vector<std::string>& i_names,
                                   const std::vector<std::string>& o_names,
                                   const Dict& opts) const override;
      int get_n_reverse() const override { return 64;}
      ///@}

      bool has_jacobian() const override { return true;}
      Function get_jacobian(const std::string& name,
            const std::vector<std::string>& i_names,
            const std::vector<std::string>& o_names, const Dict& opts) const override;

      /** \brief  Print description */
      void print(std::ostream &stream) const override { stream << "BSpline"; }

      /** \brief Is codegen supported? */
      bool has_codegen() const override { return true;}

      /** \brief Generate code for the body of the C function */
      void generateBody(CodeGenerator& g) const override;
      void generateDeclarations(CodeGenerator& g) const override {};

      std::string type_name() const override { return "BSpline"; }

      std::vector<double> coeffs_;

    private:
      std::vector<double> derivative_coeff(int i) const;
      MX jac(const MX& x) const;
    };


    class BSplineDual : public BSplineCommon {
    public:
      static Function create(const std::string &name,
        const std::vector< std::vector<double> >& knots,
        const std::vector<double>& x, const std::vector<int>& degree, int m=1, bool reverse=false,
        const Dict& opts=Dict());

      BSplineDual(const std::string &name, const std::vector<double>& knots,
        const std::vector<int>& offset, const std::vector<double>& x,
        const std::vector<int>& degree, int m, bool reverse);

      /** \brief  Destructor */
      ~BSplineDual() override {}

      ///@{
      /** \brief Number of function inputs and outputs */
      size_t get_n_in() override { return 1; }
      size_t get_n_out() override { return 1; }
      ///@}

      /// @{
      /** \brief Sparsities of function inputs and outputs */
      Sparsity get_sparsity_in(int i) override;
      Sparsity get_sparsity_out(int i) override;
      /// @}

      /** \brief  Initialize */
      void init(const Dict& opts) override;

      /** \brief  Evaluate numerically, work vectors given */
      void eval(void* mem, const double** arg, double** res, int* iw, double* w) const override;

      ///@{
      /** \brief Generate a function that calculates \a nfwd forward derivatives */
      Function get_forward(const std::string& name, int nfwd,
                             const std::vector<std::string>& i_names,
                             const std::vector<std::string>& o_names,
                             const Dict& opts) const override;
      int get_n_forward() const override { return 64;}
      ///@}

      ///@{
      /** \brief Generate a function that calculates \a nadj adjoint derivatives */
      Function get_reverse(const std::string& name, int nadj,
                                   const std::vector<std::string>& i_names,
                                   const std::vector<std::string>& o_names,
                                   const Dict& opts) const override;
      int get_n_reverse() const override { return 64;}
      ///@}


      /** \brief  Propagate sparsity forward */
      void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const override;

      /** \brief  Propagate sparsity backwards */
      void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const override;

      ///@{
      /// Is the class able to propagate seeds through the algorithm?
      bool has_spfwd() const override { return true;}
      bool has_sprev() const override { return true;}
      ///@}

      /** \brief  Print description */
      void print(std::ostream &stream) const override { stream << "BSplineDual"; }

      /** \brief Is codegen supported? */
      bool has_codegen() const override { return true;}

      /** \brief Generate code for the body of the C function */
      void generateBody(CodeGenerator& g) const override;
      void generateDeclarations(CodeGenerator& g) const override {};

      std::string type_name() const override { return "BSplineDual"; }

      std::vector<double> x_;
      bool reverse_;
      int N_;
    };


} // namespace casadi

#endif // CASADI_BSPLINE_HPP
