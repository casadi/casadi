
// File: index.xml

// File: classcasadi_1_1Blocksqp.xml
%feature("docstring") casadi::Blocksqp "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1BonMinMessageHandler.xml
%feature("docstring") casadi::BonMinMessageHandler "

[INTERNAL] 
Helper class to direct messages to  uout()

IPOPT has the concept of a Journal/Journalist BONMIN and CBC do not.

";


// File: classcasadi_1_1BSplineInterpolant.xml
%feature("docstring") casadi::BSplineInterpolant "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Callback.xml
%feature("docstring") casadi::Callback "

[INTERNAL] 
 Callback function functionality.

This class provides a public API to the FunctionInternal class that 
can be 
subclassed by the user, who is then able to implement the 
different virtual
 method. Note that the  Function class also provides a public API to 
FunctionInternal, but only allows
 calling, not being called.

The user is responsible for not deleting this class for the lifetime 
of the
 internal function object.

Joris Gillis, Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_o0

C++ includes: callback.hpp
";

%feature("docstring") casadi::Callback::jit "

[INTERNAL] 
Create a just-in-time compiled function from a C language 
string.

The names and sparsity patterns of all the inputs and outputs must be 

provided. If sparsities are not provided, all inputs and outputs are 

assumed to be scalar. Only specify the function body, assuming that 
input 
and output nonzeros are stored in arrays with the specified 
naming 
convension. The data type used is 'casadi_real', which is 
typically equal 
to 'double or another data type with the same API as 'double.

Inputs may be null pointers. This means that the all entries are zero.
 
Outputs may be null points. This means that the corresponding result 
can be
 ignored.

If an error occurs in the evaluation, issue \"return 1;\";

The final generated function will have a structure similar to:

casadi_int fname(const casadi_real** arg, casadi_real** res, 
casadi_int* 
iw, casadi_real* w, void* mem) { const casadi_real *x1, 
*x2; casadi_real 
*r1, *r2; x1 = *arg++; x2 = *arg++; r1 = *res++; r2 =
 *res++; 
<FUNCTION_BODY> return 0; }

Extra doc: https://github.com/casadi/casadi/wiki/L_1v3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L262-L274

>  Function casadi::Function::jit(const std::string &name, const std::string &body, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const std::vector< Sparsity > &sparsity_in, const std::vector< Sparsity > &sparsity_out, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Create a just-in-time compiled function from a C language string.

The names and sparsity patterns of all the inputs and outputs must be 

provided. If sparsities are not provided, all inputs and outputs are 

assumed to be scalar. Only specify the function body, assuming that 
input 
and output nonzeros are stored in arrays with the specified 
naming 
convension. The data type used is 'casadi_real', which is 
typically equal 
to 'double or another data type with the same API as 'double.

Inputs may be null pointers. This means that the all entries are zero.
 
Outputs may be null points. This means that the corresponding result 
can be
 ignored.

If an error occurs in the evaluation, issue \"return 1;\";

The final generated function will have a structure similar to:

casadi_int fname(const casadi_real** arg, casadi_real** res, 
casadi_int* 
iw, casadi_real* w, void* mem) { const casadi_real *x1, 
*x2; casadi_real 
*r1, *r2; x1 = *arg++; x2 = *arg++; r1 = *res++; r2 =
 *res++; 
<FUNCTION_BODY> return 0; }

Extra doc: https://github.com/casadi/casadi/wiki/L_1v3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L262-L274

";

";

%feature("docstring") casadi::Callback::has_jacobian "

[INTERNAL] 
Return Jacobian of all input elements with respect to all output
 
elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_oh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L175

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L92-L94

";

%feature("docstring") casadi::Callback::get_jacobian "

[INTERNAL] 
Return Jacobian of all input elements with respect to all output
 
elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_oh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L176

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L97-L102

";

%feature("docstring") casadi::Callback::has_forward "

[INTERNAL] 
Return function that calculates forward derivatives.

forward(nfwd) returns a cached instance if available, and calls   Function 
get_forward(casadi_int nfwd) if no cached version is available.

Extra doc: https://github.com/casadi/casadi/wiki/L_oi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L190

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L112-L114

";

%feature("docstring") casadi::Callback::get_forward "

[INTERNAL] 
Return function that calculates forward derivatives.

forward(nfwd) returns a cached instance if available, and calls   Function 
get_forward(casadi_int nfwd) if no cached version is available.

Extra doc: https://github.com/casadi/casadi/wiki/L_oi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L191

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L105-L110

";

%feature("docstring") casadi::Callback::has_reverse "

[INTERNAL] 
Return function that calculates adjoint derivatives.

reverse(nadj) returns a cached instance if available, and calls   Function 
get_reverse(casadi_int nadj) if no cached version is available.

Extra doc: https://github.com/casadi/casadi/wiki/L_oj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L205

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L124-L126

";

%feature("docstring") casadi::Callback::get_reverse "

[INTERNAL] 
Return function that calculates adjoint derivatives.

reverse(nadj) returns a cached instance if available, and calls   Function 
get_reverse(casadi_int nadj) if no cached version is available.

Extra doc: https://github.com/casadi/casadi/wiki/L_oj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L206

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L117-L122

";

%feature("docstring") casadi::Callback::has_jac_sparsity "

[INTERNAL] 
Return sparsity of Jacobian of all input elements.

with respect to all output elements

Extra doc: https://github.com/casadi/casadi/wiki/L_ok

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L218

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L218-L218

";

%feature("docstring") casadi::Callback::get_jac_sparsity "

[INTERNAL] 
Return sparsity of Jacobian of all input elements.

with respect to all output elements

Extra doc: https://github.com/casadi/casadi/wiki/L_ok

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L219

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L219-L220

";

%feature("docstring") casadi::Callback::expand "

[INTERNAL] 
Expand a function to SX.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L207

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L286-L302

>  Function casadi::Function::expand(const std::string &name, const Dict &opts=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Expand a function to SX.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L207

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L286-L302

";

";

%feature("docstring") casadi::Callback::size1_in "

[INTERNAL] 
Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240-L240

>  casadi_int casadi::Function::size1_in(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240-L240

";

";

%feature("docstring") casadi::Callback::size2_in "

[INTERNAL] 
Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242-L242

>  casadi_int casadi::Function::size2_in(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242-L242

";

";

%feature("docstring") casadi::Callback::size_in "

[INTERNAL] 
Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244-L246

>  std::pair<casadi_int, casadi_int> casadi::Function::size_in(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244-L246

";

";

%feature("docstring") casadi::Callback::size1_out "

[INTERNAL] 
Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254-L254

>  casadi_int casadi::Function::size1_out(const std::string &oname) const
------------------------------------------------------------------------
[INTERNAL] 
Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254-L254

";

";

%feature("docstring") casadi::Callback::size2_out "

[INTERNAL] 
Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256-L256

>  casadi_int casadi::Function::size2_out(const std::string &oname) const
------------------------------------------------------------------------
[INTERNAL] 
Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256-L256

";

";

%feature("docstring") casadi::Callback::size_out "

[INTERNAL] 
Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258-L260

>  std::pair<casadi_int, casadi_int> casadi::Function::size_out(const std::string &oname) const
------------------------------------------------------------------------
[INTERNAL] 
Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258-L260

";

";

%feature("docstring") casadi::Callback::nnz_in "

[INTERNAL] 
Get number of input nonzeros.

For a particular input or for all of the inputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L271

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L271-L271

>  casadi_int casadi::Function::nnz_in(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get number of input nonzeros.

For a particular input or for all of the inputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L271

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L271-L271

";

";

%feature("docstring") casadi::Callback::nnz_out "

[INTERNAL] 
Get number of output nonzeros.

For a particular output or for all of the outputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L282

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L282-L282

>  casadi_int casadi::Function::nnz_out(const std::string &oname) const
------------------------------------------------------------------------
[INTERNAL] 
Get number of output nonzeros.

For a particular output or for all of the outputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L282

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L282-L282

";

";

%feature("docstring") casadi::Callback::numel_in "

[INTERNAL] 
Get number of input elements.

For a particular input or for all of the inputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1ve

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L293-L293

>  casadi_int casadi::Function::numel_in(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get number of input elements.

For a particular input or for all of the inputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1ve

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L293-L293

";

";

%feature("docstring") casadi::Callback::numel_out "

[INTERNAL] 
Get number of output elements.

For a particular output or for all of the outputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L304

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L304-L304

>  casadi_int casadi::Function::numel_out(const std::string &oname) const
------------------------------------------------------------------------
[INTERNAL] 
Get number of output elements.

For a particular output or for all of the outputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L304

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L304-L304

";

";

%feature("docstring") casadi::Callback::sparsity_in "

[INTERNAL] 
Get sparsity of a given input.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L373

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L967-L973

>  const Sparsity & casadi::Function::sparsity_in(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get sparsity of a given input.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L373

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L967-L973

";

";

%feature("docstring") casadi::Callback::sparsity_out "

[INTERNAL] 
Get sparsity of a given output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L381

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L983-L989

>  const Sparsity & casadi::Function::sparsity_out(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get sparsity of a given output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L381

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L983-L989

";

";

%feature("docstring") casadi::Callback::is_diff_in "

[INTERNAL] 
Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L390

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1007-L1013

>  std::vector< bool > casadi::Function::is_diff_in() const
------------------------------------------------------------------------
[INTERNAL] 
Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L390

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1007-L1013

";

";

%feature("docstring") casadi::Callback::is_diff_out "

[INTERNAL] 
Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L391

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1015-L1021

>  std::vector< bool > casadi::Function::is_diff_out() const
------------------------------------------------------------------------
[INTERNAL] 
Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L391

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1015-L1021

";

";

%feature("docstring") casadi::Callback::sparsity_jac "

[DEPRECATED] Get, if necessary generate, the sparsity of a Jacobian 
block

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L485

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L485-L488

>  const Sparsity casadi::Function::sparsity_jac(const std::string &iind, const std::string &oind, bool compact=false, bool symmetric=false) const
------------------------------------------------------------------------

[DEPRECATED] Get, if necessary generate, the sparsity of a Jacobian 
block

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L485

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L485-L488

";

";

%feature("docstring") casadi::Callback::call "

[INTERNAL] 
Evaluate the function symbolically or numerically.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L524

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1414-L1421

>  void casadi::Function::call(const MXDict &arg, MXDict &res, bool always_inline=false, bool never_inline=false) const
------------------------------------------------------------------------
[INTERNAL] 
Evaluate the function symbolically or numerically.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L524

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1414-L1421

";

";

%feature("docstring") casadi::Callback::call_gen "

[INTERNAL] 
Call using a map.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1129

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1362-L1376

>  void casadi::Function::call_gen(const std::map< std::string, M > &arg, std::map< std::string, M > &res, bool always_inline, bool never_inline) const
------------------------------------------------------------------------
[INTERNAL] 
Call using a map.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1129

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1362-L1376

";

";

%feature("docstring") casadi::Callback::buf_in "

[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L564

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L386-L398

>  std::vector< const double * > casadi::Function::buf_in(MapArg arg) const
------------------------------------------------------------------------
[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L564

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L386-L398

";

";

%feature("docstring") casadi::Callback::buf_out "

[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L568

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L414-L427

>  std::vector< double * > casadi::Function::buf_out(MPrRes res) const
------------------------------------------------------------------------
[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L568

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L414-L427

";

";

%feature("docstring") casadi::Callback::mapaccum "

[INTERNAL] 
Create a mapaccumulated version of this function.

Suppose the function has a signature of:

::

     f: (x, u) -> (x_next , y )
  



The the mapaccumulated version has the signature:

::

     F: (x0, U) -> (X , Y )
  
      with
          U: horzcat([u0, u1, ..., u_(N-1)])
          X: horzcat([x1, x2, ..., x_N])
          Y: horzcat([y0, y1, ..., y_(N-1)])
  
      and
          x1, y0 <- f(x0, u0)
          x2, y1 <- f(x1, u1)
          ...
          x_N, y_(N-1) <- f(x_(N-1), u_(N-1))
  



Mapaccum has the following benefits over writing an equivalent for-
loop:

much faster at construction time

potentially much faster compilation times (for codegen)

offers a trade-off between memory and evaluation time

The base (settable through the options dictionary, default 10), is 
used to 
create a tower of function calls, containing unrolled for-
loops of length 
maximum base.

This technique is much more scalable in terms of memory-usage, but 
slightly
 slower at evaluation, than a plain for-loop. The effect is 
similar to that
 of a for-loop with a check-pointing instruction after 
each chunk of 
iterations with size base.

Set base to -1 to unroll all the way; no gains in memory efficiency 
here.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L697

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L484-L486

>  Function casadi::Function::mapaccum(casadi_int N, const Dict &opts=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Create a mapaccumulated version of this function.

Suppose the function has a signature of:

::

     f: (x, u) -> (x_next , y )
  



The the mapaccumulated version has the signature:

::

     F: (x0, U) -> (X , Y )
  
      with
          U: horzcat([u0, u1, ..., u_(N-1)])
          X: horzcat([x1, x2, ..., x_N])
          Y: horzcat([y0, y1, ..., y_(N-1)])
  
      and
          x1, y0 <- f(x0, u0)
          x2, y1 <- f(x1, u1)
          ...
          x_N, y_(N-1) <- f(x_(N-1), u_(N-1))
  



Mapaccum has the following benefits over writing an equivalent for-
loop:

much faster at construction time

potentially much faster compilation times (for codegen)

offers a trade-off between memory and evaluation time

The base (settable through the options dictionary, default 10), is 
used to 
create a tower of function calls, containing unrolled for-
loops of length 
maximum base.

This technique is much more scalable in terms of memory-usage, but 
slightly
 slower at evaluation, than a plain for-loop. The effect is 
similar to that
 of a for-loop with a check-pointing instruction after 
each chunk of 
iterations with size base.

Set base to -1 to unroll all the way; no gains in memory efficiency 
here.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L697

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L484-L486

";

";

%feature("docstring") casadi::Callback::fold "

[INTERNAL] 
Create a mapaccumulated version of this function.

Suppose the function has a signature of:

::

     f: (x, u) -> (x_next , y )
  



The the mapaccumulated version has the signature:

::

     F: (x0, U) -> (X , Y )
  
      with
          U: horzcat([u0, u1, ..., u_(N-1)])
          X: horzcat([x1, x2, ..., x_N])
          Y: horzcat([y0, y1, ..., y_(N-1)])
  
      and
          x1, y0 <- f(x0, u0)
          x2, y1 <- f(x1, u1)
          ...
          x_N, y_(N-1) <- f(x_(N-1), u_(N-1))
  



Mapaccum has the following benefits over writing an equivalent for-
loop:

much faster at construction time

potentially much faster compilation times (for codegen)

offers a trade-off between memory and evaluation time

The base (settable through the options dictionary, default 10), is 
used to 
create a tower of function calls, containing unrolled for-
loops of length 
maximum base.

This technique is much more scalable in terms of memory-usage, but 
slightly
 slower at evaluation, than a plain for-loop. The effect is 
similar to that
 of a for-loop with a check-pointing instruction after 
each chunk of 
iterations with size base.

Set base to -1 to unroll all the way; no gains in memory efficiency 
here.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L698

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L477-L483

";

%feature("docstring") casadi::Callback::map "

[INTERNAL]

>  Function casadi::Function::map(casadi_int n, const std::string &parallelization, casadi_int max_num_threads) const
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::Callback::generate_in "

[INTERNAL] 
Export an input file that can be passed to generate C code with 
a 
main.

See: 
 generate_out

See: 
 convert_in to convert between dict/map and vector

Extra doc: https://github.com/casadi/casadi/wiki/L_1wx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L855

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1176-L1186

>  std::vector< DM > casadi::Function::generate_in(const std::string &fname)
------------------------------------------------------------------------
[INTERNAL] 
Export an input file that can be passed to generate C code with a 
main.

See: 
 generate_out

See: 
 convert_in to convert between dict/map and vector

Extra doc: https://github.com/casadi/casadi/wiki/L_1wx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L855

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1176-L1186

";

";

%feature("docstring") casadi::Callback::generate_out "

[INTERNAL] 
Export an output file that can be checked with generated C code
 
output.

See: 
 generate_in

See: 
 convert_out to convert between dict/map and vector

Extra doc: https://github.com/casadi/casadi/wiki/L_1wy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L866

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1188-L1198

>  std::vector< DM > casadi::Function::generate_out(const std::string &fname)
------------------------------------------------------------------------
[INTERNAL] 
Export an output file that can be checked with generated C code 
output.

See: 
 generate_in

See: 
 convert_out to convert between dict/map and vector

Extra doc: https://github.com/casadi/casadi/wiki/L_1wy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L866

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1188-L1198

";

";

%feature("docstring") casadi::Callback::export_code "

[INTERNAL] 
Export function in specific language.

Only allowed for (a subset of) SX/MX Functions

Extra doc: https://github.com/casadi/casadi/wiki/L_1wz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L904

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1200-L1203

>  void casadi::Function::export_code(const std::string &lang, std::ostream &stream, const Dict &options=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Export function in specific language.

Only allowed for (a subset of) SX/MX Functions

Extra doc: https://github.com/casadi/casadi/wiki/L_1wz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L904

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1200-L1203

";

";

%feature("docstring") casadi::Callback::serialize "

[INTERNAL] 
Serialize.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L893

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1217-L1221

>  std::string casadi::Function::serialize(const Dict &opts=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Serialize.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L893

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1217-L1221

";

";

%feature("docstring") casadi::Callback::save "

[INTERNAL] 
Save  Function to a file.

See: 
 load

Extra doc: https://github.com/casadi/casadi/wiki/L_240

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L900

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1212-L1215

";

%feature("docstring") casadi::Callback::sx_in "

[INTERNAL] 
Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L944

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1511-L1517

>  const std::vector< SX > casadi::Function::sx_in() const
------------------------------------------------------------------------
[INTERNAL] 
Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L944

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1511-L1517

";

";

%feature("docstring") casadi::Callback::mx_in "

[INTERNAL] 
Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L949

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1535-L1537

>  const std::vector< MX > casadi::Function::mx_in() const
------------------------------------------------------------------------
[INTERNAL] 
Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L949

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1535-L1537

";

";

%feature("docstring") casadi::Callback::sx_out "

[INTERNAL] 
Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L962

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1519-L1525

>  const std::vector< SX > casadi::Function::sx_out() const
------------------------------------------------------------------------
[INTERNAL] 
Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L962

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1519-L1525

";

";

%feature("docstring") casadi::Callback::mx_out "

[INTERNAL] 
Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L967

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1539-L1541

>  const std::vector< MX > casadi::Function::mx_out() const
------------------------------------------------------------------------
[INTERNAL] 
Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L967

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1539-L1541

";

";

%feature("docstring") casadi::Callback::nz_from_in "

[INTERNAL] 
Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L974

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1543-L1545

";

%feature("docstring") casadi::Callback::nz_from_out "

[INTERNAL] 
Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L975

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1547-L1549

";

%feature("docstring") casadi::Callback::nz_to_in "

[INTERNAL] 
Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L976

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1551-L1553

";

%feature("docstring") casadi::Callback::nz_to_out "

[INTERNAL] 
Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L977

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1555-L1557

";

%feature("docstring") casadi::Callback::convert_in "

[INTERNAL] 
Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L996

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1595-L1597

>  std::vector< MX > casadi::Function::convert_in(const MXDict &arg) const
------------------------------------------------------------------------
[INTERNAL] 
Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L996

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1595-L1597

";

";

%feature("docstring") casadi::Callback::convert_out "

[INTERNAL] 
Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L998

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1603-L1605

>  std::vector< MX > casadi::Function::convert_out(const MXDict &arg) const
------------------------------------------------------------------------
[INTERNAL] 
Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L998

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1603-L1605

";

";

%feature("docstring") casadi::Callback::has_spfwd "

[INTERNAL] 
Is the class able to propagate seeds through the algorithm?

Extra doc: https://github.com/casadi/casadi/wiki/L_1xl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1078

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1627-L1629

";

%feature("docstring") casadi::Callback::has_sprev "

[INTERNAL] 
Is the class able to propagate seeds through the algorithm?

Extra doc: https://github.com/casadi/casadi/wiki/L_1xl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1079

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1631-L1633

";

%feature("docstring") casadi::Callback::Callback "

[INTERNAL] 
Copy constructor (throws an error)

Extra doc: https://github.com/casadi/casadi/wiki/L_o3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L64

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L32-L34

>  casadi::Callback::Callback(const Callback &obj)
------------------------------------------------------------------------
[INTERNAL] 
Copy constructor (throws an error)

Extra doc: https://github.com/casadi/casadi/wiki/L_o3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L64

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L32-L34

";

";

%feature("docstring") casadi::Callback::~Callback "

[INTERNAL] 
Destructor.

Extra doc: https://github.com/casadi/casadi/wiki/L_o4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L69

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L44-L51

";

%feature("docstring") casadi::Callback::construct "

[INTERNAL] 
Construct internal object.

This is the step that actually construct the internal object, as the 
class 
constructor only creates a null pointer. It should be called 
from the user 
constructor.

Extra doc: https://github.com/casadi/casadi/wiki/L_o5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L78

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L36-L42

";

%feature("docstring") casadi::Callback::init "

[INTERNAL] 
Initialize the object.

This function is called after the object construction (for the whole 
class 
hierarchy) is complete, but before the finalization step. It is 
called 
recursively for the whole class hierarchy, starting with the 
lowest level.

Extra doc: https://github.com/casadi/casadi/wiki/L_o6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L88

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L88-L88

";

%feature("docstring") casadi::Callback::finalize "

[INTERNAL] 
Finalize the object.

This function is called after the construction and init steps are 

completed, but before user functions are called. It is called 
recursively 
for the whole class hierarchy, starting with the highest 
level.

Extra doc: https://github.com/casadi/casadi/wiki/L_o7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L98

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L98-L98

";

%feature("docstring") casadi::Callback::eval "

[INTERNAL] 
Evaluate numerically, using temporary matrices and work vectors.

This signature is not thread-safe. For guaranteed thread-safety, use  
eval_buffer

Extra doc: https://github.com/casadi/casadi/wiki/L_o8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L106

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L60-L62

";

%feature("docstring") casadi::Callback::eval_buffer "

[INTERNAL] 
A copy-free low level interface.

In Python, you will be passed two tuples of memoryview objects Note 
that 
only the structural nonzeros are present in the memoryview 
objects/buffers.

Make sure to override  has_eval_buffer() to indicate support for this 
method.

Extra doc: https://github.com/casadi/casadi/wiki/L_o9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L116

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L53-L56

";

%feature("docstring") casadi::Callback::has_eval_buffer "

[INTERNAL] 
Does the  Callback class support a copy-free low level interface
 ?

Extra doc: https://github.com/casadi/casadi/wiki/L_265

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L122

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L57-L59

";

%feature("docstring") casadi::Callback::get_n_in "

[INTERNAL] 
Get the number of inputs.

This function is called during construction.

Extra doc: https://github.com/casadi/casadi/wiki/L_oa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L129

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L64-L66

";

%feature("docstring") casadi::Callback::get_n_out "

[INTERNAL] 
Get the number of outputs.

This function is called during construction.

Extra doc: https://github.com/casadi/casadi/wiki/L_ob

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L136

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L68-L70

";

%feature("docstring") casadi::Callback::get_sparsity_in "

[INTERNAL] 
Get the sparsity of an input.

This function is called during construction.

Extra doc: https://github.com/casadi/casadi/wiki/L_oc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L143

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L72-L74

";

%feature("docstring") casadi::Callback::get_sparsity_out "

[INTERNAL] 
Get the sparsity of an output.

This function is called during construction.

Extra doc: https://github.com/casadi/casadi/wiki/L_od

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L76-L78

";

%feature("docstring") casadi::Callback::get_name_in "

[INTERNAL] 
Get the name of an input.

This function is called during construction.

Extra doc: https://github.com/casadi/casadi/wiki/L_oe

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L80-L82

";

%feature("docstring") casadi::Callback::get_name_out "

[INTERNAL] 
Get the name of an output.

This function is called during construction.

Extra doc: https://github.com/casadi/casadi/wiki/L_of

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L164

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L84-L86

";

%feature("docstring") casadi::Callback::uses_output "

[INTERNAL] 
Do the derivative functions need nondifferentiated outputs?

Extra doc: https://github.com/casadi/casadi/wiki/L_og

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L169

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L88-L90

";

%feature("docstring") casadi::Callback::n_in "

[INTERNAL] 
Get the number of function inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L228

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L787-L789

";

%feature("docstring") casadi::Callback::n_out "

[INTERNAL] 
Get the number of function outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L791-L793

";

%feature("docstring") casadi::Callback::name_in "

[INTERNAL] 
Get input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L320

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L943-L949

>  const std::string & casadi::Function::name_in(casadi_int ind) const
------------------------------------------------------------------------
[INTERNAL] 
Get input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L320

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L943-L949

";

";

%feature("docstring") casadi::Callback::name_out "

[INTERNAL] 
Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L325

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L951-L957

>  const std::string & casadi::Function::name_out(casadi_int ind) const
------------------------------------------------------------------------
[INTERNAL] 
Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L325

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L951-L957

";

";

%feature("docstring") casadi::Callback::index_in "

[INTERNAL] 
Find the index for a string describing a particular entry of an 
input 
scheme.

example: schemeEntry(\"x_opt\") -> returns NLPSOL_X if 
FunctionInternal 
adheres to SCHEME_NLPINput

Extra doc: https://github.com/casadi/casadi/wiki/L_1vk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L333

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L927-L933

";

%feature("docstring") casadi::Callback::index_out "

[INTERNAL] 
Find the index for a string describing a particular entry of an 
output
 scheme.

example: schemeEntry(\"x_opt\") -> returns NLPSOL_X if 
FunctionInternal 
adheres to SCHEME_NLPINput

Extra doc: https://github.com/casadi/casadi/wiki/L_1vl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L341

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L935-L941

";

%feature("docstring") casadi::Callback::default_in "

[INTERNAL] 
Get default input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L346

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1423-L1425

";

%feature("docstring") casadi::Callback::max_in "

[INTERNAL] 
Get largest input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L351

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1427-L1429

";

%feature("docstring") casadi::Callback::min_in "

[INTERNAL] 
Get smallest input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L356

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1431-L1433

";

%feature("docstring") casadi::Callback::nominal_in "

[INTERNAL] 
Get nominal input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L361

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1435-L1437

";

%feature("docstring") casadi::Callback::nominal_out "

[INTERNAL] 
Get nominal output value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L366

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1439-L1441

";

%feature("docstring") casadi::Callback::factory "

[INTERNAL] ";

%feature("docstring") casadi::Callback::oracle "

[INTERNAL] 
Get oracle.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L407

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1846-L1852

";

%feature("docstring") casadi::Callback::wrap "

[INTERNAL] 
Wrap in an  Function instance consisting of only one  MX call.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L412

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1854-L1856

";

%feature("docstring") casadi::Callback::wrap_as_needed "

[INTERNAL] 
Wrap in a  Function with options.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L417

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1858-L1860

";

%feature("docstring") casadi::Callback::which_depends "

[INTERNAL] 
Which variables enter with some order.

Parameters:
-----------

order: 
Only 1 (linear) and 2 (nonlinear) allowed

tr: 
Flip the relationship. Return which expressions contain the variables

Extra doc: https://github.com/casadi/casadi/wiki/L_1vx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L425

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1764-L1771

";

%feature("docstring") casadi::Callback::print_dimensions "

[INTERNAL] 
Print dimensions of inputs and outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L432

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1095-L1097

";

%feature("docstring") casadi::Callback::print_options "

[INTERNAL] 
Print options to a stream.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L437

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1099-L1101

";

%feature("docstring") casadi::Callback::print_option "

[INTERNAL] 
Print all information there is to know about a certain option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L442

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1103-L1105

";

%feature("docstring") casadi::Callback::has_option "

[INTERNAL] 
Does a particular option exist.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L447

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1107-L1114

";

%feature("docstring") casadi::Callback::change_option "

[INTERNAL] 
Change option after object creation for debugging.

This is only possible for a selected number of options that do not 
change 
the numerical results of the computation, e.g. to enable a more
 verbose 
output or saving to file.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L455

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1116-L1126

";

%feature("docstring") casadi::Callback::jacobian_old "

[DEPRECATED] Replaced by  Function::factory.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L466

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L856-L862

";

%feature("docstring") casadi::Callback::hessian_old "

[DEPRECATED] Replaced by  Function::factory.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L471

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L864-L872

";

%feature("docstring") casadi::Callback::jacobian "

[INTERNAL] 
Calculate all Jacobian blocks.

Generates a function that takes all non-differentiated inputs and 
outputs 
and calculates all Jacobian blocks. Inputs that are not needed
 by the 
routine are all-zero sparse matrices with the correct 
dimensions.  Output 
blocks that are not calculated, e.g. if the corresponding input or 
output 
is marked non-differentiated are also all-zero sparse. The 
Jacobian blocks 
are sorted starting by all the blocks for the first 
output, then all the 
blocks for the second output and so on. E.g. f : 
(x, y) -> (r, s) results 
in the function jac_f : (x, y, out_r, out_s) 
-> (jac_r_x, jac_r_y, jac_s_x,
 jac_s_y)

This function is cached.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L508

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L884-L890

";

%feature("docstring") casadi::Callback::rev "

[INTERNAL] 
Propagate sparsity backward with temporary memory allocation.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L633

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L460-L475

>  int casadi::Function::rev(std::vector< bvec_t * > arg, std::vector< bvec_t * > res) const
------------------------------------------------------------------------
[INTERNAL] 
Propagate sparsity backward with temporary memory allocation.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L633

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L460-L475

";

";

%feature("docstring") casadi::Callback::mapsum "

[INTERNAL] 
Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization: 
Type of parallelization used: unroll|serial|openmp

Extra doc: https://github.com/casadi/casadi/wiki/L_1wh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L642

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L724-L731

";

%feature("docstring") casadi::Callback::slice "

[INTERNAL] 
returns a new function with a selection of inputs/outputs of the
 
original

Extra doc: https://github.com/casadi/casadi/wiki/L_1wl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L754

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L715-L722

";

%feature("docstring") casadi::Callback::forward "

[INTERNAL] 
Get a function that calculates  nfwd forward derivatives.



::

     Returns a function with <tt>n_in + n_out + n_in</tt> inputs
     and <tt>nfwd</tt> outputs.
     The first <tt>n_in</tt> inputs correspond to nondifferentiated inputs.
     The next <tt>n_out</tt> inputs correspond to nondifferentiated outputs.
     and the last <tt>n_in</tt> inputs correspond to forward seeds,
     stacked horizontally
     The  <tt>n_out</tt> outputs correspond to forward sensitivities,
     stacked horizontally.     *
     <tt>(n_in = n_in(), n_out = n_out())</tt>
  
    The functions returned are cached, meaning that if called multiple timed
    with the same value, then multiple references to the same function will be returned.
  

Extra doc: https://github.com/casadi/casadi/wiki/L_1wq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L800

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1079-L1085

";

%feature("docstring") casadi::Callback::reverse "

[INTERNAL] 
Get a function that calculates  nadj adjoint derivatives.



::

     Returns a function with <tt>n_in + n_out + n_out</tt> inputs
     and <tt>n_in</tt> outputs.
     The first <tt>n_in</tt> inputs correspond to nondifferentiated inputs.
     The next <tt>n_out</tt> inputs correspond to nondifferentiated outputs.
     and the last <tt>n_out</tt> inputs correspond to adjoint seeds,
     stacked horizontally
     The  <tt>n_in</tt> outputs correspond to adjoint sensitivities,
     stacked horizontally.     *
     <tt>(n_in = n_in(), n_out = n_out())</tt>
  
     <tt>(n_in = n_in(), n_out = n_out())</tt>
  
    The functions returned are cached, meaning that if called multiple timed
    with the same value, then multiple references to the same function will be returned.
  

Extra doc: https://github.com/casadi/casadi/wiki/L_1wr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L820

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1087-L1093

";

%feature("docstring") casadi::Callback::jac_sparsity "

[INTERNAL] 
Get, if necessary generate, the sparsity of a single Jacobian 
block.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L830

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L911-L917

>  Sparsity casadi::Function::jac_sparsity(casadi_int oind, casadi_int iind, bool compact=false) const
------------------------------------------------------------------------
[INTERNAL] 
Get, if necessary generate, the sparsity of a single Jacobian block.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L830

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L911-L917

";

";

%feature("docstring") casadi::Callback::generate "

[INTERNAL] 
Export / Generate C code for the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L840

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1132-L1134

>  std::string casadi::Function::generate(const Dict &opts=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Export / Generate C code for the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L840

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1132-L1134

";

";

%feature("docstring") casadi::Callback::generate_dependencies "

[INTERNAL] 
Export / Generate C code for the dependency function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ww

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L845

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1142-L1144

";

%feature("docstring") casadi::Callback::stats "

[INTERNAL] 
Get all statistics obtained at the end of the last evaluate 
call.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L932

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L896-L898

";

%feature("docstring") casadi::Callback::has_free "

[INTERNAL] 
Does the function have free variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1004

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1635-L1637

";

%feature("docstring") casadi::Callback::get_free "

[INTERNAL] 
Get free variables as a string.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1009

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1128-L1130

";

%feature("docstring") casadi::Callback::free_sx "

[INTERNAL] 
Get all the free variables of the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1014

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1611-L1617

";

%feature("docstring") casadi::Callback::free_mx "

[INTERNAL] 
Get all the free variables of the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1019

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1619-L1625

";

%feature("docstring") casadi::Callback::generate_lifted "

[INTERNAL] 
Extract the functions needed for the Lifted  Newton method.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1024

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1639-L1645

";

%feature("docstring") casadi::Callback::n_nodes "

[INTERNAL] 
Number of nodes in the algorithm.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1030

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1703-L1709

";

%feature("docstring") casadi::Callback::n_instructions "

[INTERNAL] 
Number of instruction in the algorithm (SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xe

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1035

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1647-L1653

";

%feature("docstring") casadi::Callback::instruction_id "

[INTERNAL] 
Identifier index of the instruction (SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1040

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1671-L1677

";

%feature("docstring") casadi::Callback::instruction_input "

[INTERNAL] 
Locations in the work vector for the inputs of the instruction.

(SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1047

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1679-L1685

";

%feature("docstring") casadi::Callback::instruction_constant "

[INTERNAL] 
Get the floating point output argument of an instruction 
(SXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1052

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1687-L1693

";

%feature("docstring") casadi::Callback::instruction_output "

[INTERNAL] 
Location in the work vector for the output of the instruction.

(SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1059

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1695-L1701

";

%feature("docstring") casadi::Callback::instruction_MX "

[INTERNAL] 
Get the  MX node corresponding to an instruction (MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1064

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1655-L1661

";

%feature("docstring") casadi::Callback::instructions_sx "

[INTERNAL] 
Get the SX node corresponding to all instructions (SXFunction)

Note: input and output instructions have no SX representation. This 
method 
returns nan for those instructions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1072

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1663-L1669

";

%feature("docstring") casadi::Callback::sz_arg "

[INTERNAL] 
Get required length of arg field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1085

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1027-L1027

";

%feature("docstring") casadi::Callback::sz_res "

[INTERNAL] 
Get required length of res field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1090

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1029-L1029

";

%feature("docstring") casadi::Callback::sz_iw "

[INTERNAL] 
Get required length of iw field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1095

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1031-L1031

";

%feature("docstring") casadi::Callback::sz_w "

[INTERNAL] 
Get required length of w field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1033-L1033

";

%feature("docstring") casadi::Callback::sz_work "

[INTERNAL] 
Get number of temporary variables needed.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1106

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1023-L1025

";

%feature("docstring") casadi::Callback::set_work "

[INTERNAL] 
Set the (persistent) work vectors.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1111

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1052-L1059

";

%feature("docstring") casadi::Callback::set_temp "

[INTERNAL] 
Set the (temporary) work vectors.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1117

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1061-L1068

";

%feature("docstring") casadi::Callback::setup "

[INTERNAL] 
Set the (persistent and temporary) work vectors.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1123

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1070-L1077

";

%feature("docstring") casadi::Callback::name "

[INTERNAL] 
Name of the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1137

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1250-L1257

";

%feature("docstring") casadi::Callback::is_a "

[INTERNAL] 
Check if the function is of a particular type.

Optionally check if name matches one of the base classes (default 
true)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1144

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1607-L1609

";

%feature("docstring") casadi::Callback::assert_size_in "

[INTERNAL] 
Assert that an input dimension is equal so some given value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1723-L1729

";

%feature("docstring") casadi::Callback::assert_size_out "

[INTERNAL] 
Assert that an output dimension is equal so some given value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1192

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1731-L1736

";

%feature("docstring") casadi::Callback::assert_sparsity_out "

[INTERNAL] 
Assert that an output sparsity is a multiple of some given 
sparsity.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1738-L1747

";

%feature("docstring") casadi::Callback::checkout "

[INTERNAL] 
Checkout a memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1199

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1711-L1713

";

%feature("docstring") casadi::Callback::release "

[INTERNAL] 
Release a memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1202

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1715-L1717

";

%feature("docstring") casadi::Callback::memory "

[INTERNAL] 
Get memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1206

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1719-L1721

";

%feature("docstring") casadi::Callback::cache "

[INTERNAL] 
Get all functions in the cache.

Extra doc: https://github.com/casadi/casadi/wiki/L_26i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1212

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1773-L1780

";

%feature("docstring") casadi::Callback::get_function "

[INTERNAL] 
Get a dependency function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1791-L1797

>  Function casadi::Function::get_function(const std::string &name) const
------------------------------------------------------------------------
[INTERNAL] 
Get a dependency function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1791-L1797

";

";

%feature("docstring") casadi::Callback::has_function "

[INTERNAL] 
Check if a particular dependency exists.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1227

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1799-L1806

";

%feature("docstring") casadi::Callback::find_functions "

[INTERNAL] 
Get all functions embedded in the expression graphs.

Parameters:
-----------

max_depth: 
Maximum depth - a negative number indicates no maximum

Extra doc: https://github.com/casadi/casadi/wiki/L_1y6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1234

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1808-L1824

";

%feature("docstring") casadi::Callback::find_function "

[INTERNAL] 
Get a specific function embedded in the expression graphs.

Parameters:
-----------

name: 
Name of function needed

max_depth: 
Maximum depth - a negative number indicates no maximum

Extra doc: https://github.com/casadi/casadi/wiki/L_1y7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1826-L1843

";

%feature("docstring") casadi::Callback::info "

[INTERNAL] 
Obtain information about function

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1245

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1872-L1874

";

%feature("docstring") casadi::Callback::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L99-L101

";

%feature("docstring") casadi::Callback::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L103-L109

";

%feature("docstring") casadi::Callback::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::Callback::is_null "

[INTERNAL] 
Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L73-L75

";

%feature("docstring") casadi::Callback::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L129-L131

";


// File: classcasadi_1_1casadi__limits.xml
%feature("docstring") casadi::casadi_limits "

[INTERNAL] 
 casadi_limits class

The following class, which acts as a complements to the standard 

std::numeric_limits class, allows specifying certain properties of 
scalar 
objects. The template can be specialized for e.g. symbolic 
scalars 
Joel 
Andersson

C++ includes: casadi_limits.hpp
";


// File: classcasadi_1_1casadi__limits_3_01SXElem_01_4.xml
%feature("docstring") casadi::casadi_limits< SXElem > "

[INTERNAL] C++ includes: sx_elem.hpp
";


// File: classcasadi_1_1CasadiException.xml
%feature("docstring") casadi::CasadiException "

[INTERNAL] 
Casadi exception class.



::

  \\\\author Joel Andersson
  \\\\date 2010
  Example for simple exception throwing:
  \\\\code
          throw CasadiException(\"This is a nasty error\");
  \\\\endcode
  Example for exception chaining:
  \\\\code
          try {
                  throw CasadiException(\"This is a nasty error\");
          catch(CasadiException &e) {
                  throw CasadiException(\"Serious error.\") << e;
          }
  \\\\endcode
  

Extra doc: https://github.com/casadi/casadi/wiki/L_7u

C++ includes: exception.hpp
";

%feature("docstring") casadi::CasadiException::CasadiException "

[INTERNAL] 
Form message string.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L84

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L84-L84

>  casadi::CasadiException::CasadiException(const std::string &msg)
------------------------------------------------------------------------
[INTERNAL] 
Form message string.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L84

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L84-L84

";

";

%feature("docstring") casadi::CasadiException::~CasadiException "

[INTERNAL]  throw ()
Destructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L87

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L87-L87

";

%feature("docstring") casadi::CasadiException::what "

[INTERNAL]  throw ()
Display error.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L90

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L90-L92

";


// File: classcasadi_1_1CasadiHandler.xml
%feature("docstring") casadi::CasadiHandler "

[INTERNAL] ";


// File: classcasadi_1_1CasadiMeta.xml
%feature("docstring") casadi::CasadiMeta "

[INTERNAL] 
Collects global CasADi meta information.

Joris Gillis

Extra doc: https://github.com/casadi/casadi/wiki/L_23k

C++ includes: casadi_meta.hpp
";


// File: classcasadi_1_1ClangCompiler.xml
%feature("docstring") casadi::ClangCompiler "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1CodeGenerator.xml
%feature("docstring") casadi::CodeGenerator "

[INTERNAL] 
Helper class for C code generation.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_ru

C++ includes: code_generator.hpp
";

%feature("docstring") casadi::CodeGenerator::CodeGenerator "

[INTERNAL] 
Constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L46

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L35-L180

";

%feature("docstring") casadi::CodeGenerator::add "

[INTERNAL] 
Add a function (name generated)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L49

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L288-L316

";

%feature("docstring") casadi::CodeGenerator::dump "

[INTERNAL] 
Generate a file, return code as string.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L57

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L318-L322

>  std::string casadi::CodeGenerator::dump()
------------------------------------------------------------------------
[INTERNAL] 
Generate a file, return code as string.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L57

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L318-L322

";

";

%feature("docstring") casadi::CodeGenerator::generate "

[INTERNAL] 
Generate file(s)

The \"prefix\" argument will be prepended to the generated files and 
may be
 a directory or a file prefix. returns the filename

Extra doc: https://github.com/casadi/casadi/wiki/L_rv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L66

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L410-L464

";

%feature("docstring") casadi::CodeGenerator::add_include "

[INTERNAL] 
Add an include file optionally using a relative path \"...\" 
instead 
of an absolute path <...>

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L69

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L942-L962

";

%feature("docstring") casadi::CodeGenerator::add_dependency "

[INTERNAL] 
Add a function dependency.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L74

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L207-L286

";

%feature("docstring") casadi::CodeGenerator::add_external "

[INTERNAL] 
Add an external function declaration.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L77

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L986-L988

";

%feature("docstring") casadi::CodeGenerator::shorthand "

[INTERNAL] 
Add/get a shorthand.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L83

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L995-L1001

>  std::string casadi::CodeGenerator::shorthand(const std::string &name, bool allow_adding=true)
------------------------------------------------------------------------
[INTERNAL] 
Add/get a shorthand.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L83

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L995-L1001

";

";

%feature("docstring") casadi::CodeGenerator::sparsity "

[INTERNAL] ";

%feature("docstring") casadi::CodeGenerator::add_sparsity "

[INTERNAL] ";

%feature("docstring") casadi::CodeGenerator::get_sparsity "

[INTERNAL] 
Get the index of an existing sparsity pattern.

Extra doc: https://github.com/casadi/casadi/wiki/L_rw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L94

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1011-L1013

";

%feature("docstring") casadi::CodeGenerator::get_constant "

[INTERNAL] 
Get or add an integer constant.

Extra doc: https://github.com/casadi/casadi/wiki/L_ry

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L104

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1057-L1078

>  casadi_int casadi::CodeGenerator::get_constant(const std::vector< casadi_int > &v, bool allow_adding=false)
------------------------------------------------------------------------
[INTERNAL] 
Get or add an integer constant.

Extra doc: https://github.com/casadi/casadi/wiki/L_ry

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L104

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1057-L1078

";

";

%feature("docstring") casadi::CodeGenerator::constant "

[INTERNAL]

>  std::string casadi::CodeGenerator::constant(casadi_int v)

>  std::string casadi::CodeGenerator::constant(const std::string &v)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::CodeGenerator::constant_copy "

[INTERNAL] 
Represent an array constant; adding it when new.

Extra doc: https://github.com/casadi/casadi/wiki/L_s0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L121

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1084-L1098

";

%feature("docstring") casadi::CodeGenerator::define_rom_double "

[INTERNAL] 
Allocate file scope double read-only memory.

Extra doc: https://github.com/casadi/casadi/wiki/L_s2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L134

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L727-L732

";

%feature("docstring") casadi::CodeGenerator::rom_double "

[INTERNAL] 
Access file scope double read-only memory.

Extra doc: https://github.com/casadi/casadi/wiki/L_s3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L139

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L734-L739

";

%feature("docstring") casadi::CodeGenerator::define_rom_integer "

[INTERNAL] 
Allocate file scope integer read-only memory.

Extra doc: https://github.com/casadi/casadi/wiki/L_s4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L144

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L741-L746

";

%feature("docstring") casadi::CodeGenerator::rom_integer "

[INTERNAL] 
Access file scope integer read-only memory.

Extra doc: https://github.com/casadi/casadi/wiki/L_s5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L149

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L748-L753

";

%feature("docstring") casadi::CodeGenerator::print_formatted "

[INTERNAL] 
Print without newline characters.

Extra doc: https://github.com/casadi/casadi/wiki/L_s8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L166

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1940-L1964

";

%feature("docstring") casadi::CodeGenerator::flush "

[INTERNAL] 
Flush the buffer to a stream of choice.

Extra doc: https://github.com/casadi/casadi/wiki/L_sa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L181

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1987-L1990

";

%feature("docstring") casadi::CodeGenerator::local "

[INTERNAL] 
Declare a local variable.

Extra doc: https://github.com/casadi/casadi/wiki/L_sb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L186

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1992-L2004

";

%feature("docstring") casadi::CodeGenerator::scope_enter "

[INTERNAL] 
Enter a local scope.

Extra doc: https://github.com/casadi/casadi/wiki/L_sc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L191

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L182-L185

";

%feature("docstring") casadi::CodeGenerator::scope_exit "

[INTERNAL] 
Exit a local scope.

Extra doc: https://github.com/casadi/casadi/wiki/L_sd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L196

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L187-L205

";

%feature("docstring") casadi::CodeGenerator::sx_work "

[INTERNAL] 
Declare a work vector element.

Extra doc: https://github.com/casadi/casadi/wiki/L_se

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2006-L2017

";

%feature("docstring") casadi::CodeGenerator::init_local "

[INTERNAL] 
Specify the default value for a local variable.

Extra doc: https://github.com/casadi/casadi/wiki/L_sf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L206

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2019-L2022

";

%feature("docstring") casadi::CodeGenerator::indent "

[INTERNAL] 
Increase indentation.

Extra doc: https://github.com/casadi/casadi/wiki/L_sg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L211

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L211-L211

";

%feature("docstring") casadi::CodeGenerator::unindent "

[INTERNAL] 
Decrease indentation.

Extra doc: https://github.com/casadi/casadi/wiki/L_sh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L216

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L216-L216

";

%feature("docstring") casadi::CodeGenerator::avoid_stack "

[INTERNAL] 
Avoid stack?

Extra doc: https://github.com/casadi/casadi/wiki/L_si

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L221

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L221-L221

";

%feature("docstring") casadi::CodeGenerator::initializer "

[INTERNAL]

>  std::string casadi::CodeGenerator::initializer(const std::vector< casadi_int > &v)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::CodeGenerator::sanitize_source "

[INTERNAL] 
Sanitize source files for codegen.

Extra doc: https://github.com/casadi/casadi/wiki/L_sl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L239

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2025-L2116

";

%feature("docstring") casadi::CodeGenerator::dot "

[INTERNAL] 
Codegen inner product.

Extra doc: https://github.com/casadi/casadi/wiki/L_sm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L246

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1702-L1708

";

%feature("docstring") casadi::CodeGenerator::mv "

[INTERNAL] 
Codegen dense matrix-vector multiplication.

Extra doc: https://github.com/casadi/casadi/wiki/L_so

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L257

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1902-L1907

>  std::string casadi::CodeGenerator::mv(const std::string &x, casadi_int nrow_x, casadi_int ncol_x, const std::string &y, const std::string &z, bool tr)
------------------------------------------------------------------------
[INTERNAL] 
Codegen dense matrix-vector multiplication.

Extra doc: https://github.com/casadi/casadi/wiki/L_so

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L257

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1902-L1907

";

";

%feature("docstring") casadi::CodeGenerator::axpy "

[INTERNAL] 
Codegen axpy: y += a*x.

Extra doc: https://github.com/casadi/casadi/wiki/L_sp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L263

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1854-L1858

";

%feature("docstring") casadi::CodeGenerator::clip_min "

[INTERNAL] 
Codegen clip_min: Clips the smaller entries in a vector than min
 to 
the min.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L271

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1860-L1864

";

%feature("docstring") casadi::CodeGenerator::clip_max "

[INTERNAL] 
Codegen clip_max: Clips the larger entries in a vector than max 
to the
 max.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L279

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1866-L1870

";

%feature("docstring") casadi::CodeGenerator::vector_fmax "

[INTERNAL] 
Codegen vector_fmax: Takes vectorwise max of a vector and writes
 the 
result to second vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L286

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1872-L1876

";

%feature("docstring") casadi::CodeGenerator::vector_fmin "

[INTERNAL] 
Codegen vector_fmin: Takes vectorwise min of a vector and writes
 the 
result to second vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1878-L1882

";

%feature("docstring") casadi::CodeGenerator::masked_norm_inf "

[INTERNAL] 
codegen masked_norm_inf: The mask tells what entry is used in 
the inf-
norm.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L300

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1884-L1888

";

%feature("docstring") casadi::CodeGenerator::scal "

[INTERNAL] 
What does scal do??

Extra doc: https://github.com/casadi/casadi/wiki/L_sq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L307

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1890-L1893

";

%feature("docstring") casadi::CodeGenerator::mtimes "

[INTERNAL] 
Codegen sparse matrix-matrix multiplication.

Extra doc: https://github.com/casadi/casadi/wiki/L_sr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L312

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1909-L1916

";

%feature("docstring") casadi::CodeGenerator::trilsolve "

[INTERNAL] 
Codegen lower triangular solve.

Extra doc: https://github.com/casadi/casadi/wiki/L_ss

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L320

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1918-L1923

";

%feature("docstring") casadi::CodeGenerator::triusolve "

[INTERNAL] 
Codegen upper triangular solve.

Extra doc: https://github.com/casadi/casadi/wiki/L_st

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L326

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1925-L1930

";

%feature("docstring") casadi::CodeGenerator::bilin "

[INTERNAL] 
Codegen bilinear form.

Extra doc: https://github.com/casadi/casadi/wiki/L_su

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L332

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1710-L1716

";

%feature("docstring") casadi::CodeGenerator::rank1 "

[INTERNAL] 
Rank-1 update.

Extra doc: https://github.com/casadi/casadi/wiki/L_sv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L338

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1718-L1726

";

%feature("docstring") casadi::CodeGenerator::logsumexp "

[INTERNAL] 
\\\\brie LogSumExp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L342

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1933-L1938

";

%feature("docstring") casadi::CodeGenerator::interpn "

[INTERNAL] 
Multilinear interpolation.

Extra doc: https://github.com/casadi/casadi/wiki/L_sw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L347

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1728-L1739

";

%feature("docstring") casadi::CodeGenerator::interpn_grad "

[INTERNAL] 
Multilinear interpolation - calculate gradient.

Extra doc: https://github.com/casadi/casadi/wiki/L_sx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L356

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1741-L1751

";

%feature("docstring") casadi::CodeGenerator::trans "

[INTERNAL] 
Transpose.

Extra doc: https://github.com/casadi/casadi/wiki/L_sy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L366

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1753-L1759

";

%feature("docstring") casadi::CodeGenerator::qr "

[INTERNAL] 
QR factorization.

Extra doc: https://github.com/casadi/casadi/wiki/L_sz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L372

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2151-L2159

";

%feature("docstring") casadi::CodeGenerator::qr_solve "

[INTERNAL] 
QR solve.

Extra doc: https://github.com/casadi/casadi/wiki/L_t0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L381

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2162-L2171

";

%feature("docstring") casadi::CodeGenerator::lsqr_solve "

[INTERNAL] 
\\\\brief LSQR solve

Extra doc: https://github.com/casadi/casadi/wiki/L_t1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L390

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2174-L2179

";

%feature("docstring") casadi::CodeGenerator::ldl "

[INTERNAL] 
LDL factorization.

Extra doc: https://github.com/casadi/casadi/wiki/L_t2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L396

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2182-L2188

";

%feature("docstring") casadi::CodeGenerator::ldl_solve "

[INTERNAL] 
LDL solve.

Extra doc: https://github.com/casadi/casadi/wiki/L_t3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L404

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2191-L2197

";

%feature("docstring") casadi::CodeGenerator::fmax "

[INTERNAL] 
fmax

Extra doc: https://github.com/casadi/casadi/wiki/L_t4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L412

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2200-L2203

";

%feature("docstring") casadi::CodeGenerator::fmin "

[INTERNAL] 
fmin

Extra doc: https://github.com/casadi/casadi/wiki/L_t5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L417

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2206-L2209

";

%feature("docstring") casadi::CodeGenerator::mmax "

[INTERNAL] 
mmax

Extra doc: https://github.com/casadi/casadi/wiki/L_t6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L422

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2248-L2251

";

%feature("docstring") casadi::CodeGenerator::mmin "

[INTERNAL] 
mmin

Extra doc: https://github.com/casadi/casadi/wiki/L_t7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L427

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2254-L2257

";

%feature("docstring") casadi::CodeGenerator::vfmax "

[INTERNAL] 
vfmax

Extra doc: https://github.com/casadi/casadi/wiki/L_ta

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L442

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2224-L2227

>  std::string casadi::CodeGenerator::vfmax(const std::string &x, const std::string &n, const std::string &y)
------------------------------------------------------------------------
[INTERNAL] 
vfmax

Extra doc: https://github.com/casadi/casadi/wiki/L_ta

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L442

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2224-L2227

";

";

%feature("docstring") casadi::CodeGenerator::vfmin "

[INTERNAL] 
vfmin

Extra doc: https://github.com/casadi/casadi/wiki/L_tb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L447

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2230-L2233

>  std::string casadi::CodeGenerator::vfmin(const std::string &x, const std::string &n, const std::string &y)
------------------------------------------------------------------------
[INTERNAL] 
vfmin

Extra doc: https://github.com/casadi/casadi/wiki/L_tb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L447

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2230-L2233

";

";

%feature("docstring") casadi::CodeGenerator::max "

[INTERNAL] 
max

Extra doc: https://github.com/casadi/casadi/wiki/L_tc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L452

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2236-L2239

";

%feature("docstring") casadi::CodeGenerator::min "

[INTERNAL] 
min

Extra doc: https://github.com/casadi/casadi/wiki/L_td

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L457

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2242-L2245

";

%feature("docstring") casadi::CodeGenerator::norm_inf "

[INTERNAL] 
norm_inf

Extra doc: https://github.com/casadi/casadi/wiki/L_te

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L462

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2272-L2275

";

%feature("docstring") casadi::CodeGenerator::norm_2 "

[INTERNAL] 
norm_2



::

       Extra doc: https://github.com/casadi/casadi/wiki/L_256 
  



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L469

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2278-L2281

";

%feature("docstring") casadi::CodeGenerator::max_viol "

[INTERNAL] 
max_viol

Extra doc: https://github.com/casadi/casadi/wiki/L_tf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L474

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2260-L2263

";

%feature("docstring") casadi::CodeGenerator::sum_viol "

[INTERNAL] 
sum_viol

Extra doc: https://github.com/casadi/casadi/wiki/L_tg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L480

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2266-L2269

";

%feature("docstring") casadi::CodeGenerator::bound_consistency "

[INTERNAL] 
bound_consistency

Extra doc: https://github.com/casadi/casadi/wiki/L_th

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L486

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2309-L2314

";

%feature("docstring") casadi::CodeGenerator::lb_eig "

[INTERNAL] 
lb_eig

Extra doc: https://github.com/casadi/casadi/wiki/L_ti

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L492

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2284-L2287

";

%feature("docstring") casadi::CodeGenerator::regularize "

[INTERNAL] 
regularize

Extra doc: https://github.com/casadi/casadi/wiki/L_tj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L497

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2290-L2293

";

%feature("docstring") casadi::CodeGenerator::convexify_eval "

[INTERNAL] 
convexify

Extra doc: https://github.com/casadi/casadi/wiki/L_tk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L502

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2296-L2300

";

%feature("docstring") casadi::CodeGenerator::low "

[INTERNAL] 
low

Extra doc: https://github.com/casadi/casadi/wiki/L_tl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L508

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2303-L2306

";

%feature("docstring") casadi::CodeGenerator::declare "

[INTERNAL] 
Declare a function.

Extra doc: https://github.com/casadi/casadi/wiki/L_tm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L514

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1761-L1772

";

%feature("docstring") casadi::CodeGenerator::comment "

[INTERNAL] 
Write a comment line (ignored if not verbose)

Extra doc: https://github.com/casadi/casadi/wiki/L_tn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L519

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2118-L2122

";

%feature("docstring") casadi::CodeGenerator::add_auxiliary "

[INTERNAL] 
Add a built-in auxiliary function.

Extra doc: https://github.com/casadi/casadi/wiki/L_tp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L612

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1104-L1547

";

%feature("docstring") casadi::CodeGenerator::add_io_sparsities "

[INTERNAL] 
Add io sparsity patterns of a function.

Extra doc: https://github.com/casadi/casadi/wiki/L_tq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L617

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2125-L2148

";

%feature("docstring") casadi::CodeGenerator::work "

[INTERNAL] 
Get work vector name from index

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L622

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L863-L871

";

%feature("docstring") casadi::CodeGenerator::workel "

[INTERNAL] 
Get work vector element from index

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L625

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L873-L879

";

%feature("docstring") casadi::CodeGenerator::print_vector "

[INTERNAL] 
Print real vector to a c file.

Extra doc: https://github.com/casadi/casadi/wiki/L_ts

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L640

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L900-L903

>  void casadi::CodeGenerator::print_vector(std::ostream &s, const std::string &name, const std::vector< double > &v)
------------------------------------------------------------------------
[INTERNAL] 
Print real vector to a c file.

Extra doc: https://github.com/casadi/casadi/wiki/L_ts

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L640

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L900-L903

";

";

%feature("docstring") casadi::CodeGenerator::copy "

[INTERNAL] 
Create a copy operation.

Extra doc: https://github.com/casadi/casadi/wiki/L_tt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L646

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1641-L1648

";

%feature("docstring") casadi::CodeGenerator::copy_check "

[INTERNAL] ";

%feature("docstring") casadi::CodeGenerator::copy_default "

[INTERNAL] ";

%feature("docstring") casadi::CodeGenerator::fill "

[INTERNAL] 
Create a fill operation.

Extra doc: https://github.com/casadi/casadi/wiki/L_tu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L655

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1692-L1700

";

%feature("docstring") casadi::CodeGenerator::clear "

[INTERNAL] 
Create a fill operation.

Extra doc: https://github.com/casadi/casadi/wiki/L_tv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L660

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1670-L1676

";

%feature("docstring") casadi::CodeGenerator::arg "

[INTERNAL] 
Refer to argument.

Extra doc: https://github.com/casadi/casadi/wiki/L_tw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L665

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1678-L1680

";

%feature("docstring") casadi::CodeGenerator::res "

[INTERNAL] 
Refer to resuly.

Extra doc: https://github.com/casadi/casadi/wiki/L_tx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L670

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1682-L1684

";

%feature("docstring") casadi::CodeGenerator::mem "

[INTERNAL] 
Access thread-local memory.

Extra doc: https://github.com/casadi/casadi/wiki/L_ty

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L675

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1686-L1690

";

%feature("docstring") casadi::CodeGenerator::project "

[INTERNAL] 
Sparse assignment.

Extra doc: https://github.com/casadi/casadi/wiki/L_tz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L680

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1775-L1787

";

%feature("docstring") casadi::CodeGenerator::tri_project "

[INTERNAL] 
Project triangular part.

Extra doc: https://github.com/casadi/casadi/wiki/L_u0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L687

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1790-L1798

";

%feature("docstring") casadi::CodeGenerator::densify "

[INTERNAL] 
Densify.

Extra doc: https://github.com/casadi/casadi/wiki/L_u1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L693

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1801-L1809

";

%feature("docstring") casadi::CodeGenerator::sparsify "

[INTERNAL] 
Sparsify.

Extra doc: https://github.com/casadi/casadi/wiki/L_u2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L699

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1812-L1820

";

%feature("docstring") casadi::CodeGenerator::to_mex "

[INTERNAL] 
Create matrix in MATLAB's MEX format.

Extra doc: https://github.com/casadi/casadi/wiki/L_u3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L705

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1549-L1554

";

%feature("docstring") casadi::CodeGenerator::from_mex "

[INTERNAL] 
Get matrix from MATLAB's MEX format.

Extra doc: https://github.com/casadi/casadi/wiki/L_u4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L710

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1556-L1567

";

%feature("docstring") casadi::CodeGenerator::printf "

[INTERNAL]

>  std::string casadi::CodeGenerator::printf(const std::string &str, const std::string &arg1)

>  std::string casadi::CodeGenerator::printf(const std::string &str, const std::string &arg1, const std::string &arg2)

>  std::string casadi::CodeGenerator::printf(const std::string &str, const std::string &arg1, const std::string &arg2, const std::string &arg3)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::CodeGenerator::print_op "

[INTERNAL]

>  std::string casadi::CodeGenerator::print_op(casadi_int op, const std::string &a0, const std::string &a1)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::CodeGenerator::file_slurp "

[INTERNAL] 
Slurp a file.

Extra doc: https://github.com/casadi/casadi/wiki/L_u7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L738

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2317-L2320

";

%feature("docstring") casadi::CodeGenerator::cache_check "

[INTERNAL] 
cache check

Extra doc: https://github.com/casadi/casadi/wiki/L_u8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L743

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2323-L2328

";

%feature("docstring") casadi::CodeGenerator::sz_work "

[INTERNAL] 
Get number of temporary variables needed for all functions.

Extra doc: https://github.com/casadi/casadi/wiki/L_258

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L759

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2330-L2338

";


// File: classcasadi_1_1Collocation.xml
%feature("docstring") casadi::Collocation "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Conic.xml
%feature("docstring") casadi::Conic "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1DaeBuilder.xml


/*
 Variables and equations 
*/

/*
 Symbolic modeling 
*/

/*
Formulate a dynamic system model

*/

/*
 Register an existing variable 
*/

/*
 [DEPRECATED] Specify all variables of a type: Call set_all instead 
*/

/*
 Manipulation 
*/

/*
Reformulate the dynamic optimization problem.

*/

/*
 Functions 
*/

/*
Add or load auxiliary functions

*/

/*
 Import and export 
*/
%feature("docstring") casadi::DaeBuilder "

[INTERNAL] 
A symbolic representation of a differential-algebraic equations 
model.

Variables:
==========





::

  t:      independent variable (usually time)
  c:      constants
  p:      parameters
  d:      dependent parameters (time independent)
  u:      controls
  w:      dependent variables  (time dependent)
  x:      differential states
  z:      algebraic variables
  q:      quadrature states
  y:      outputs
  



Equations:
==========





::

  differential equations: \\\\dot{x} ==  ode(...)
  algebraic equations:          0 ==  alg(...)
  quadrature equations:   \\\\dot{q} == quad(...)
  dependent parameters:         d == ddef(d_prev,p)
  dependent variables:          w == wdef(w_prev,x,z,u,p,t)
  output equations:             y == ydef(...)
  initial equations:     init_lhs == init_rhs(...)
  events:      when when_cond < 0: when_lhs := when_rhs
  



Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_5c

C++ includes: dae_builder.hpp
";

%feature("docstring") casadi::DaeBuilder::t "

[INTERNAL] 
Independent variable (usually time)

Extra doc: https://github.com/casadi/casadi/wiki/L_5e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L93

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L63-L65

";

%feature("docstring") casadi::DaeBuilder::x "

[INTERNAL] 
Differential states.

Extra doc: https://github.com/casadi/casadi/wiki/L_5f

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L98

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L67-L69

";

%feature("docstring") casadi::DaeBuilder::ode "

[INTERNAL] 
Ordinary differential equations (ODE)

Extra doc: https://github.com/casadi/casadi/wiki/L_5g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L103

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L71-L73

";

%feature("docstring") casadi::DaeBuilder::z "

[INTERNAL] 
Algebraic variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_5h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L108

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L75-L77

";

%feature("docstring") casadi::DaeBuilder::alg "

[INTERNAL] 
Algebraic equations.

Extra doc: https://github.com/casadi/casadi/wiki/L_5i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L113

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L79-L81

";

%feature("docstring") casadi::DaeBuilder::q "

[INTERNAL] 
Quadrature states.

Extra doc: https://github.com/casadi/casadi/wiki/L_5j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L118

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L83-L85

";

%feature("docstring") casadi::DaeBuilder::quad "

[INTERNAL] 
Quadrature equations.

Extra doc: https://github.com/casadi/casadi/wiki/L_5k

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L123

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L87-L89

";

%feature("docstring") casadi::DaeBuilder::y "

[INTERNAL] 
 Output variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_5l

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L128

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L91-L93

";

%feature("docstring") casadi::DaeBuilder::ydef "

[INTERNAL] 
Definitions of output variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_5m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L133

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L95-L97

";

%feature("docstring") casadi::DaeBuilder::u "

[INTERNAL] 
Free controls.

Extra doc: https://github.com/casadi/casadi/wiki/L_5n

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L99-L101

";

%feature("docstring") casadi::DaeBuilder::p "

[INTERNAL] 
Parameters.

Extra doc: https://github.com/casadi/casadi/wiki/L_5o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L143

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L103-L105

";

%feature("docstring") casadi::DaeBuilder::c "

[INTERNAL] 
Named constants.

Extra doc: https://github.com/casadi/casadi/wiki/L_5p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L148

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L107-L109

";

%feature("docstring") casadi::DaeBuilder::cdef "

[INTERNAL] 
Definitions of named constants.

Extra doc: https://github.com/casadi/casadi/wiki/L_5q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L111-L113

";

%feature("docstring") casadi::DaeBuilder::d "

[INTERNAL] 
Dependent parameters.

Extra doc: https://github.com/casadi/casadi/wiki/L_5r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L158

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L115-L117

";

%feature("docstring") casadi::DaeBuilder::ddef "

[INTERNAL] 
Definitions of dependent parameters.

Interdependencies are allowed but must be non-cyclic.

Extra doc: https://github.com/casadi/casadi/wiki/L_5s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L165

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L119-L121

";

%feature("docstring") casadi::DaeBuilder::w "

[INTERNAL] 
Dependent variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_5t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L170

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L123-L125

";

%feature("docstring") casadi::DaeBuilder::wdef "

[INTERNAL] 
Dependent variables and corresponding definitions.

Interdependencies are allowed but must be non-cyclic.

Extra doc: https://github.com/casadi/casadi/wiki/L_5u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L177

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L127-L129

";

%feature("docstring") casadi::DaeBuilder::aux "

[INTERNAL] 
Auxiliary variables: Used e.g. to define functions.

Extra doc: https://github.com/casadi/casadi/wiki/L_5v

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L182

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L131-L133

";

%feature("docstring") casadi::DaeBuilder::init_lhs "

[INTERNAL] 
Initial conditions, left-hand-side.

Extra doc: https://github.com/casadi/casadi/wiki/L_5w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L187

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L135-L137

";

%feature("docstring") casadi::DaeBuilder::init_rhs "

[INTERNAL] 
Initial conditions, right-hand-side.

Extra doc: https://github.com/casadi/casadi/wiki/L_5x

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L192

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L139-L141

";

%feature("docstring") casadi::DaeBuilder::when_cond "

[INTERNAL] 
When statement: triggering condition.

Extra doc: https://github.com/casadi/casadi/wiki/L_5y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L197

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L143-L145

";

%feature("docstring") casadi::DaeBuilder::when_lhs "

[INTERNAL] 
When statement: left-hand-side.

Extra doc: https://github.com/casadi/casadi/wiki/L_5z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L202

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L147-L149

";

%feature("docstring") casadi::DaeBuilder::when_rhs "

[INTERNAL] 
When statement: right-hand-side.

Extra doc: https://github.com/casadi/casadi/wiki/L_60

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L207

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L151-L153

";

%feature("docstring") casadi::DaeBuilder::has_t "

[INTERNAL] 
Is there a time variable?

Extra doc: https://github.com/casadi/casadi/wiki/L_64

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L231

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L182-L184

";

%feature("docstring") casadi::DaeBuilder::nx "

[INTERNAL] 
Differential states.

Extra doc: https://github.com/casadi/casadi/wiki/L_65

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L236

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L186-L188

";

%feature("docstring") casadi::DaeBuilder::nz "

[INTERNAL] 
Algebraic variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_66

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L241

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L190-L192

";

%feature("docstring") casadi::DaeBuilder::nq "

[INTERNAL] 
Quadrature states.

Extra doc: https://github.com/casadi/casadi/wiki/L_67

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L246

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L194-L196

";

%feature("docstring") casadi::DaeBuilder::ny "

[INTERNAL] 
 Output variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_68

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L251

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L198-L200

";

%feature("docstring") casadi::DaeBuilder::nu "

[INTERNAL] 
Free controls.

Extra doc: https://github.com/casadi/casadi/wiki/L_69

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L256

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L202-L204

";

%feature("docstring") casadi::DaeBuilder::np "

[INTERNAL] 
Parameters.

Extra doc: https://github.com/casadi/casadi/wiki/L_6a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L261

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L206-L208

";

%feature("docstring") casadi::DaeBuilder::nc "

[INTERNAL] 
Named constants.

Extra doc: https://github.com/casadi/casadi/wiki/L_6b

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L266

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L210-L212

";

%feature("docstring") casadi::DaeBuilder::nd "

[INTERNAL] 
Dependent parameters.

Extra doc: https://github.com/casadi/casadi/wiki/L_6c

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L271

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L214-L216

";

%feature("docstring") casadi::DaeBuilder::nw "

[INTERNAL] 
Dependent variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_6d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L276

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L218-L220

";

%feature("docstring") casadi::DaeBuilder::add_t "

[INTERNAL] 
Add an independent variable (time)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L284

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L430-L435

";

%feature("docstring") casadi::DaeBuilder::add_p "

[INTERNAL] 
Add a new parameter.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L287

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L437-L444

";

%feature("docstring") casadi::DaeBuilder::add_u "

[INTERNAL] 
Add a new control.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L290

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L446-L453

";

%feature("docstring") casadi::DaeBuilder::add_x "

[INTERNAL] 
Add a new differential state.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L455-L462

";

%feature("docstring") casadi::DaeBuilder::add_z "

[INTERNAL] 
Add a new algebraic variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L296

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L464-L471

";

%feature("docstring") casadi::DaeBuilder::add_q "

[INTERNAL] 
Add a new quadrature state.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L299

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L473-L480

";

%feature("docstring") casadi::DaeBuilder::add_c "

[INTERNAL] 
Add a new constant.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L302

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L482-L489

";

%feature("docstring") casadi::DaeBuilder::add_d "

[INTERNAL] 
Add a new dependent parameter.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L305

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L491-L498

";

%feature("docstring") casadi::DaeBuilder::add_w "

[INTERNAL] 
Add a new dependent variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L308

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L500-L507

";

%feature("docstring") casadi::DaeBuilder::add_y "

[INTERNAL] 
Add a new output.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L311

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L509-L516

";

%feature("docstring") casadi::DaeBuilder::set_ode "

[INTERNAL] 
Specify the ordinary differential equation for a state.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L314

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L535-L541

";

%feature("docstring") casadi::DaeBuilder::set_alg "

[INTERNAL] 
Specificy the residual equation for an algebraic variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L317

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L543-L549

";

%feature("docstring") casadi::DaeBuilder::add_aux "

[INTERNAL] 
Add an auxiliary variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L320

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L518-L522

";

%feature("docstring") casadi::DaeBuilder::add_init "

[INTERNAL] 
Add an initial equation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L323

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L524-L527

";

%feature("docstring") casadi::DaeBuilder::add_when "

[INTERNAL] 
Add a when statement.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L326

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L529-L533

";

%feature("docstring") casadi::DaeBuilder::sanity_check "

[INTERNAL] 
Check if dimensions match.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L329

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L551-L557

";

%feature("docstring") casadi::DaeBuilder::register_t "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::register_p "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::register_u "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::register_x "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::register_z "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::register_q "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::register_c "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::register_d "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::register_w "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::register_y "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::set_u "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::set_x "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::set_z "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::set_q "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::set_y "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::clear_in "

[DEPRECATED] Clear input variable: Replaced by clear_all

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L371

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L371-L371

";

%feature("docstring") casadi::DaeBuilder::eliminate_w "

[INTERNAL] 
Eliminate all dependent variables.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L375

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L611-L617

";

%feature("docstring") casadi::DaeBuilder::lift "

[INTERNAL] 
Lift problem formulation by extracting shared subexpressions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L378

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L619-L625

";

%feature("docstring") casadi::DaeBuilder::eliminate_quad "

[INTERNAL] 
Eliminate quadrature states and turn them into ODE states.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L381

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L248-L254

";

%feature("docstring") casadi::DaeBuilder::sort_d "

[INTERNAL] 
Sort dependent parameters.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L384

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L256-L262

";

%feature("docstring") casadi::DaeBuilder::sort_w "

[INTERNAL] 
Sort dependent variables.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L387

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L264-L270

";

%feature("docstring") casadi::DaeBuilder::sort_z "

[INTERNAL] 
Sort algebraic variables.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L390

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L272-L278

";

%feature("docstring") casadi::DaeBuilder::prune "

[INTERNAL] 
Prune unused controls.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L393

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L280-L286

";

%feature("docstring") casadi::DaeBuilder::tear "

[INTERNAL] 
Identify iteration variables and residual equations using naming
 
convention.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L396

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L288-L294

";

%feature("docstring") casadi::DaeBuilder::add_fun "

[INTERNAL] 
Add an external function.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L413

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L770-L773

>  Function casadi::DaeBuilder::add_fun(const std::string &name, const Importer &compiler, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Add an external function.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L413

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L770-L773

";

";

%feature("docstring") casadi::DaeBuilder::has_fun "

[INTERNAL] 
Does a particular function already exist?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L417

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L775-L782

";

%feature("docstring") casadi::DaeBuilder::fun "

[INTERNAL] 
Get all functions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L423

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L814-L816

>  std::vector< Function > casadi::DaeBuilder::fun() const
------------------------------------------------------------------------
[INTERNAL] 
Get all functions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L423

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L814-L816

";

";

%feature("docstring") casadi::DaeBuilder::gather_fun "

[INTERNAL] 
Collect embedded functions from the expression graph.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L426

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L793-L812

";

%feature("docstring") casadi::DaeBuilder::parse_fmi "

[INTERNAL] 
Import existing problem from FMI/XML

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L433

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L433-L433

";

%feature("docstring") casadi::DaeBuilder::provides_directional_derivative "

[INTERNAL] 
Does the FMU provide support for analytic derivatives.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L436

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L230-L238

";

%feature("docstring") casadi::DaeBuilder::load_fmi_description "

[INTERNAL] 
Import problem description from FMI or XML.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L439

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L222-L228

";

%feature("docstring") casadi::DaeBuilder::export_fmu "

[INTERNAL] 
Export instance into an FMU.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L442

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L240-L246

";

%feature("docstring") casadi::DaeBuilder::add_lc "

[INTERNAL] 
Add a named linear combination of output expressions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L445

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L711-L718

";

%feature("docstring") casadi::DaeBuilder::create "

[INTERNAL] 
Load a function from an FMU DLL, standard IO conforming with 

simulator.

Parameters:
-----------

name: 
Name assigned to the resulting function object

opts: 
Optional settings

Extra doc: https://github.com/casadi/casadi/wiki/L_6f

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L472

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L742-L749

>  Function casadi::DaeBuilder::create(const std::string &name, const Dict &opts=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Load a function from an FMU DLL, standard IO conforming with 
simulator.

Parameters:
-----------

name: 
Name assigned to the resulting function object

opts: 
Optional settings

Extra doc: https://github.com/casadi/casadi/wiki/L_6f

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L472

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L742-L749

";

";

%feature("docstring") casadi::DaeBuilder::var "

[INTERNAL] 
Get variable expressions by index.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L700

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L892-L899

>  std::vector< MX > casadi::DaeBuilder::var(const std::vector< size_t > &ind) const
------------------------------------------------------------------------
[INTERNAL] 
Get variable expressions by index.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L700

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L892-L899

";

";

%feature("docstring") casadi::DaeBuilder::beq "

[INTERNAL] 
Get/set the binding equation for a variable

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L490

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L594-L601

";

%feature("docstring") casadi::DaeBuilder::set_beq "

[INTERNAL] 
Get/set the binding equation for a variable

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L491

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L603-L609

";

%feature("docstring") casadi::DaeBuilder::value_reference "

[INTERNAL] 
Get/set value reference

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L496

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L627-L629

";

%feature("docstring") casadi::DaeBuilder::set_value_reference "

[INTERNAL] 
Get/set value reference

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L497

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L631-L633

";

%feature("docstring") casadi::DaeBuilder::description "

[INTERNAL] 
Get/set description

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L502

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L635-L637

";

%feature("docstring") casadi::DaeBuilder::set_description "

[INTERNAL] 
Get/set description

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L503

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L639-L641

";

%feature("docstring") casadi::DaeBuilder::type "

[INTERNAL] 
Get/set the type

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L508

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L643-L652

";

%feature("docstring") casadi::DaeBuilder::set_type "

[INTERNAL] 
Get/set the type

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L509

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L654-L661

";

%feature("docstring") casadi::DaeBuilder::causality "

[INTERNAL] 
Get/set the causality

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L514

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L663-L665

";

%feature("docstring") casadi::DaeBuilder::set_causality "

[INTERNAL] 
Get/set the causality

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L515

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L667-L669

";

%feature("docstring") casadi::DaeBuilder::variability "

[INTERNAL] 
Get/set the variability

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L520

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L671-L673

";

%feature("docstring") casadi::DaeBuilder::set_variability "

[INTERNAL] 
Get/set the variability

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L521

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L675-L677

";

%feature("docstring") casadi::DaeBuilder::initial "

[INTERNAL] 
Get/set the initial property

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L526

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L679-L681

";

%feature("docstring") casadi::DaeBuilder::set_initial "

[INTERNAL] 
Get/set the initial property

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L527

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L683-L685

";

%feature("docstring") casadi::DaeBuilder::unit "

[INTERNAL] 
Get/set the unit

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L532

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L687-L689

";

%feature("docstring") casadi::DaeBuilder::set_unit "

[INTERNAL] 
Get/set the unit

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L533

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L691-L693

";

%feature("docstring") casadi::DaeBuilder::display_unit "

[INTERNAL] 
Get/set the display unit

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L538

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L695-L697

";

%feature("docstring") casadi::DaeBuilder::set_display_unit "

[INTERNAL] 
Get/set the display unit

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L539

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L699-L701

";

%feature("docstring") casadi::DaeBuilder::variable "

[INTERNAL] 
Access a variable by index

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L684

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L862-L868

>  const Variable & casadi::DaeBuilder::variable(size_t ind) const
------------------------------------------------------------------------
[INTERNAL] 
Access a variable by index

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L684

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L862-L868

";

";

%feature("docstring") casadi::DaeBuilder::type_name "

[INTERNAL] 
Readable name of the class.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L74

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L74-L74

";

%feature("docstring") casadi::DaeBuilder::DaeBuilder "

[INTERNAL] 
Construct a  DaeBuilder instance.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L80

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L54-L57

>  casadi::DaeBuilder::DaeBuilder(const std::string &name, const std::string &path=\"\", const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Construct a  DaeBuilder instance.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L80

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L54-L57

";

";

%feature("docstring") casadi::DaeBuilder::name "

[INTERNAL] 
Get variable names by indices.

Extra doc: https://github.com/casadi/casadi/wiki/L_6i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L716

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L929-L938

>  std::vector< std::string > casadi::DaeBuilder::name(const std::vector< size_t > &ind) const
------------------------------------------------------------------------
[INTERNAL] 
Get variable names by indices.

Extra doc: https://github.com/casadi/casadi/wiki/L_6i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L716

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L929-L938

";

";

%feature("docstring") casadi::DaeBuilder::outputs "

[INTERNAL] 
Model structure: outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_61

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L213

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L155-L162

";

%feature("docstring") casadi::DaeBuilder::derivatives "

[INTERNAL] 
Model structure: derivatives.

Extra doc: https://github.com/casadi/casadi/wiki/L_62

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L218

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L164-L171

";

%feature("docstring") casadi::DaeBuilder::initial_unknowns "

[INTERNAL] 
Model structure: initial unknowns.

Extra doc: https://github.com/casadi/casadi/wiki/L_63

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L223

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L173-L180

";

%feature("docstring") casadi::DaeBuilder::clear_all "

[INTERNAL] 
Clear all variables of a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L333

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L414-L420

";

%feature("docstring") casadi::DaeBuilder::set_all "

[INTERNAL] 
Set all variables of a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L336

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L422-L428

";

%feature("docstring") casadi::DaeBuilder::dependent_fun "

[INTERNAL] 
Construct a function for evaluating dependent parameters.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L475

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L827-L836

";

%feature("docstring") casadi::DaeBuilder::der "

[INTERNAL] 
Get the time derivative of an expression, single variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L552

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L568-L581

>  std::string casadi::DaeBuilder::der(const std::string &name) const
------------------------------------------------------------------------
[INTERNAL] 
Get the time derivative of an expression, single variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L552

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L568-L581

";

";

%feature("docstring") casadi::DaeBuilder::numel "

[INTERNAL] 
Get the number of elements of a variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L543

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L703-L705

";

%feature("docstring") casadi::DaeBuilder::dimension "

[INTERNAL] 
Get the dimensions of a variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L546

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L707-L709

";

%feature("docstring") casadi::DaeBuilder::attribute "

[INTERNAL] 
Get an attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L599

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L949-L957

>  std::vector< double > casadi::DaeBuilder::attribute(const std::string &a, const std::vector< std::string > &name) const
------------------------------------------------------------------------
[INTERNAL] 
Get an attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L599

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L949-L957

";

";

%feature("docstring") casadi::DaeBuilder::set_attribute "

[INTERNAL] 
Set an attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L602

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L967-L974

>  void casadi::DaeBuilder::set_attribute(const std::string &a, const std::vector< std::string > &name, const std::vector< double > &val)
------------------------------------------------------------------------
[INTERNAL] 
Set an attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L602

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L967-L974

";

";

%feature("docstring") casadi::DaeBuilder::min "

[INTERNAL] 
Get the lower bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L606

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L985-L992

>  std::vector< double > casadi::DaeBuilder::min(const std::vector< std::string > &name) const
------------------------------------------------------------------------
[INTERNAL] 
Get the lower bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L606

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L985-L992

";

";

%feature("docstring") casadi::DaeBuilder::set_min "

[INTERNAL] 
Set the lower bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L609

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1002-L1008

>  void casadi::DaeBuilder::set_min(const std::vector< std::string > &name, const std::vector< double > &val)
------------------------------------------------------------------------
[INTERNAL] 
Set the lower bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L609

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1002-L1008

";

";

%feature("docstring") casadi::DaeBuilder::max "

[INTERNAL] 
Get the upper bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L612

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1019-L1026

>  std::vector< double > casadi::DaeBuilder::max(const std::vector< std::string > &name) const
------------------------------------------------------------------------
[INTERNAL] 
Get the upper bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L612

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1019-L1026

";

";

%feature("docstring") casadi::DaeBuilder::set_max "

[INTERNAL] 
Set the upper bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L615

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1036-L1042

>  void casadi::DaeBuilder::set_max(const std::vector< std::string > &name, const std::vector< double > &val)
------------------------------------------------------------------------
[INTERNAL] 
Set the upper bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L615

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1036-L1042

";

";

%feature("docstring") casadi::DaeBuilder::nominal "

[INTERNAL] 
Get the nominal value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L618

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1053-L1060

>  std::vector< double > casadi::DaeBuilder::nominal(const std::vector< std::string > &name) const
------------------------------------------------------------------------
[INTERNAL] 
Get the nominal value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L618

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1053-L1060

";

";

%feature("docstring") casadi::DaeBuilder::set_nominal "

[INTERNAL] 
Set the nominal value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L621

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1070-L1076

>  void casadi::DaeBuilder::set_nominal(const std::vector< std::string > &name, const std::vector< double > &val)
------------------------------------------------------------------------
[INTERNAL] 
Set the nominal value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L621

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1070-L1076

";

";

%feature("docstring") casadi::DaeBuilder::start "

[INTERNAL] 
Get the start attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L624

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1088-L1095

>  std::vector< double > casadi::DaeBuilder::start(const std::vector< std::string > &name) const
------------------------------------------------------------------------
[INTERNAL] 
Get the start attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L624

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1088-L1095

";

";

%feature("docstring") casadi::DaeBuilder::set_start "

[INTERNAL] 
Set the start attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L627

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1105-L1111

>  void casadi::DaeBuilder::set_start(const std::vector< std::string > &name, const std::vector< double > &val)
------------------------------------------------------------------------
[INTERNAL] 
Set the start attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L627

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1105-L1111

";

";

%feature("docstring") casadi::DaeBuilder::reset "

[INTERNAL] ";

%feature("docstring") casadi::DaeBuilder::set "

[INTERNAL] 
Set the current value (string)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L633

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1145-L1152

>  void casadi::DaeBuilder::set(const std::vector< std::string > &name, const std::vector< std::string > &val)
------------------------------------------------------------------------
[INTERNAL] 
Set the current value (string)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L633

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1145-L1152

";

";

%feature("docstring") casadi::DaeBuilder::get "

[INTERNAL] 
Evaluate the values for a set of variables at the initial time.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L636

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1158-L1173

>  std::vector< GenericType > casadi::DaeBuilder::get(const std::vector< std::string > &name) const
------------------------------------------------------------------------
[INTERNAL] 
Evaluate the values for a set of variables at the initial time.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L636

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1158-L1173

";

";

%feature("docstring") casadi::DaeBuilder::add_variable "

[INTERNAL] 
Add a new variable from symbolic expressions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L645

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L332-L335

>  void casadi::DaeBuilder::add_variable(const MX &new_v)
------------------------------------------------------------------------
[INTERNAL] 
Add a new variable from symbolic expressions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L645

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L332-L335

";

";

%feature("docstring") casadi::DaeBuilder::add_variable_new "

[INTERNAL] 
Add a new variable from symbolic expressions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L654

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L347-L351

>  size_t casadi::DaeBuilder::add_variable_new(const MX &new_v)
------------------------------------------------------------------------
[INTERNAL] 
Add a new variable from symbolic expressions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L654

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L347-L351

";

";

%feature("docstring") casadi::DaeBuilder::has_variable "

[INTERNAL] 
Check if a particular variable exists.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L657

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L296-L303

";

%feature("docstring") casadi::DaeBuilder::all_variables "

[INTERNAL] 
Get a list of all variables.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L660

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L305-L312

";

%feature("docstring") casadi::DaeBuilder::oracle "

[INTERNAL] 
Get the (cached) oracle, SX or  MX.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L663

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L818-L825

";

%feature("docstring") casadi::DaeBuilder::jac_sparsity "

[INTERNAL] 
Get Jacobian sparsity.

Extra doc: https://github.com/casadi/casadi/wiki/L_6g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L668

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1175-L1183

";

%feature("docstring") casadi::DaeBuilder::new_variable "

[INTERNAL] 
Create a new variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L673

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L314-L320

";

%feature("docstring") casadi::DaeBuilder::find "

[INTERNAL] 
Get indices of variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L706

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L910-L917

>  std::vector< size_t > casadi::DaeBuilder::find(const std::vector< std::string > &name) const
------------------------------------------------------------------------
[INTERNAL] 
Get indices of variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L706

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L910-L917

";

";

%feature("docstring") casadi::DaeBuilder::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L99-L101

";

%feature("docstring") casadi::DaeBuilder::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L103-L109

";

%feature("docstring") casadi::DaeBuilder::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::DaeBuilder::is_null "

[INTERNAL] 
Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L73-L75

";

%feature("docstring") casadi::DaeBuilder::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L129-L131

";


// File: classcasadi_1_1DeserializerBase.xml
%feature("docstring") casadi::DeserializerBase "

[INTERNAL] C++ includes: serializer.hpp
";

%feature("docstring") casadi::DeserializerBase::DeserializerBase "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::~DeserializerBase "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::pop_type "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_sparsity "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_mx "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_dm "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_sx "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_linsol "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_function "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_generictype "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_int "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_double "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_string "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_sparsity_vector
 "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_mx_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_dm_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_sx_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_linsol_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_function_vector
 "

[INTERNAL] ";

%feature("docstring") 
casadi::DeserializerBase::blind_unpack_generictype_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_int_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_double_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::blind_unpack_string_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_sparsity "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_mx "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_dm "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_sx "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_linsol "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_function "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_generictype "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_int "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_double "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_string "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_sparsity_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_mx_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_dm_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_sx_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_linsol_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_function_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_generictype_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_int_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_double_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::unpack_string_vector "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::connect "

[INTERNAL] ";

%feature("docstring") casadi::DeserializerBase::reset "

[INTERNAL] ";


// File: classcasadi_1_1DeserializingStream.xml
%feature("docstring") casadi::DeserializingStream "

[INTERNAL] 
Helper class for Serialization.

Joris Gillis

Extra doc: https://github.com/casadi/casadi/wiki/L_ak

C++ includes: serializing_stream.hpp
";

%feature("docstring") casadi::DeserializingStream::DeserializingStream "

[INTERNAL]

>  casadi::DeserializingStream::DeserializingStream(const DeserializingStream &)=delete
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::DeserializingStream::unpack "

[INTERNAL]

>  void casadi::DeserializingStream::unpack(MX &e)

>  void casadi::DeserializingStream::unpack(SXElem &e)

>  void casadi::DeserializingStream::unpack(Linsol &e)

>  void casadi::DeserializingStream::unpack(Matrix< T > &e)

>  void casadi::DeserializingStream::unpack(Function &e)

>  void casadi::DeserializingStream::unpack(Importer &e)

>  void casadi::DeserializingStream::unpack(GenericType &e)

>  void casadi::DeserializingStream::unpack(std::ostream &s)

>  void casadi::DeserializingStream::unpack(Slice &e)

>  void casadi::DeserializingStream::unpack(int &e)

>  void casadi::DeserializingStream::unpack(bool &e)

>  void casadi::DeserializingStream::unpack(casadi_int &e)

>  void casadi::DeserializingStream::unpack(size_t &e)

>  void casadi::DeserializingStream::unpack(std::string &e)

>  void casadi::DeserializingStream::unpack(double &e)

>  void casadi::DeserializingStream::unpack(char &e)

>  void casadi::DeserializingStream::unpack(std::vector< T > &e)

>  void casadi::DeserializingStream::unpack(std::map< K, V > &e)

>  void casadi::DeserializingStream::unpack(std::pair< A, B > &e)

>  void casadi::DeserializingStream::unpack(const std::string &descr, T &e)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::DeserializingStream::version "

[INTERNAL] ";

%feature("docstring") casadi::DeserializingStream::connect "

[INTERNAL] ";

%feature("docstring") casadi::DeserializingStream::reset "

[INTERNAL] ";


// File: classcasadi_1_1Dple.xml
%feature("docstring") casadi::Dple "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Expm.xml
%feature("docstring") casadi::Expm "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1FastNewton.xml
%feature("docstring") casadi::FastNewton "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Feasiblesqpmethod.xml
%feature("docstring") casadi::Feasiblesqpmethod "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1FileDeserializer.xml
%feature("docstring") casadi::FileDeserializer "

[INTERNAL] C++ includes: serializer.hpp
";

%feature("docstring") casadi::FileDeserializer::FileDeserializer "

[INTERNAL] 
Advanced deserialization of CasADi objects.

See: 
 FileSerializer

Extra doc: https://github.com/casadi/casadi/wiki/L_7t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.hpp#L250

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.cpp#L106-L112

";

%feature("docstring") casadi::FileDeserializer::~FileDeserializer "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::pop_type "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_sparsity "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_mx "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_dm "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_sx "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_linsol "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_function "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_generictype "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_int "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_double "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_string "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_sparsity_vector
 "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_mx_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_dm_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_sx_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_linsol_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_function_vector
 "

[INTERNAL] ";

%feature("docstring") 
casadi::FileDeserializer::blind_unpack_generictype_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_int_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_double_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::blind_unpack_string_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_sparsity "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_mx "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_dm "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_sx "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_linsol "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_function "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_generictype "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_int "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_double "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_string "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_sparsity_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_mx_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_dm_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_sx_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_linsol_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_function_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_generictype_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_int_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_double_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::unpack_string_vector "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::connect "

[INTERNAL] ";

%feature("docstring") casadi::FileDeserializer::reset "

[INTERNAL] ";


// File: classcasadi_1_1FileSerializer.xml
%feature("docstring") casadi::FileSerializer "

[INTERNAL] C++ includes: serializer.hpp
";

%feature("docstring") casadi::FileSerializer::FileSerializer "

[INTERNAL] 
Advanced serialization of CasADi objects.

See: 
 StringSerializer,  FileDeserializer

Extra doc: https://github.com/casadi/casadi/wiki/L_7q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.hpp#L221

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.cpp#L41-L49

";

%feature("docstring") casadi::FileSerializer::~FileSerializer "

[INTERNAL] ";

%feature("docstring") casadi::FileSerializer::pack "

[INTERNAL] ";

%feature("docstring") casadi::FileSerializer::connect "

[INTERNAL] ";

%feature("docstring") casadi::FileSerializer::reset "

[INTERNAL] ";


// File: classcasadi_1_1FixedStepIntegrator.xml
%feature("docstring") casadi::FixedStepIntegrator "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Function.xml
%feature("docstring") casadi::Function "

[INTERNAL] 
 Function object.

A  Function instance is a general multiple-input, multiple-output function 
where 
each input and output can be a sparse matrix.
 For an introduction to
 this class, see the CasADi user guide. Function is a reference counted and 
immutable class; copying a class instance 
is very cheap and its behavior 
(with some exceptions) is not affected 
by calling its member functions.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_1uw

>List of available options

+------------------+-----------------+------------------+------------------+
|        Id        |      Type       |   Description    |     Used in      |
+==================+=================+==================+==================+
| ad_weight        | OT_DOUBLE       | Weighting factor | casadi::Function |
|                  |                 | for derivative   | Internal         |
|                  |                 | calculation.When |                  |
|                  |                 | there is an      |                  |
|                  |                 | option of either |                  |
|                  |                 | using forward or |                  |
|                  |                 | reverse mode     |                  |
|                  |                 | directional      |                  |
|                  |                 | derivatives, the |                  |
|                  |                 | condition ad_wei |                  |
|                  |                 | ght*nf<=(1-      |                  |
|                  |                 | ad_weight)*na is |                  |
|                  |                 | used where nf    |                  |
|                  |                 | and na are       |                  |
|                  |                 | estimates of the |                  |
|                  |                 | number of        |                  |
|                  |                 | forward/reverse  |                  |
|                  |                 | mode directional |                  |
|                  |                 | derivatives      |                  |
|                  |                 | needed. By       |                  |
|                  |                 | default,         |                  |
|                  |                 | ad_weight is     |                  |
|                  |                 | calculated       |                  |
|                  |                 | automatically,   |                  |
|                  |                 | but this can be  |                  |
|                  |                 | overridden by    |                  |
|                  |                 | setting this     |                  |
|                  |                 | option. In       |                  |
|                  |                 | particular, 0    |                  |
|                  |                 | means forcing    |                  |
|                  |                 | forward mode and |                  |
|                  |                 | 1 forcing        |                  |
|                  |                 | reverse mode.    |                  |
|                  |                 | Leave unset for  |                  |
|                  |                 | (class specific) |                  |
|                  |                 | heuristics.      |                  |
+------------------+-----------------+------------------+------------------+
| ad_weight_sp     | OT_DOUBLE       | Weighting factor | casadi::Function |
|                  |                 | for sparsity     | Internal         |
|                  |                 | pattern          |                  |
|                  |                 | calculation calc |                  |
|                  |                 | ulation.Override |                  |
|                  |                 | s default        |                  |
|                  |                 | behavior. Set to |                  |
|                  |                 | 0 and 1 to force |                  |
|                  |                 | forward and      |                  |
|                  |                 | reverse mode     |                  |
|                  |                 | respectively.    |                  |
|                  |                 | Cf. option       |                  |
|                  |                 | \"ad_weight\".     |                  |
|                  |                 | When set to -1,  |                  |
|                  |                 | sparsity is      |                  |
|                  |                 | completely       |                  |
|                  |                 | ignored and      |                  |
|                  |                 | dense matrices   |                  |
|                  |                 | are used.        |                  |
+------------------+-----------------+------------------+------------------+
| always_inline    | OT_BOOL         | Force inlining.  | casadi::Function |
|                  |                 |                  | Internal         |
+------------------+-----------------+------------------+------------------+
| cache            | OT_DICT         | Prepopulate the  | casadi::Function |
|                  |                 | function cache.  | Internal         |
|                  |                 | Default: empty   |                  |
+------------------+-----------------+------------------+------------------+
| compiler         | OT_STRING       | Just-in-time     | casadi::Function |
|                  |                 | compiler plugin  | Internal         |
|                  |                 | to be used.      |                  |
+------------------+-----------------+------------------+------------------+
| custom_jacobian  | OT_FUNCTION     | Override         | casadi::Function |
|                  |                 | CasADi's AD. Use | Internal         |
|                  |                 | together with    |                  |
|                  |                 | 'jac_penalty':   |                  |
|                  |                 | 0. Note: Highly  |                  |
|                  |                 | experimental.    |                  |
|                  |                 | Syntax may break |                  |
|                  |                 | often.           |                  |
+------------------+-----------------+------------------+------------------+
| der_options      | OT_DICT         | Default options  | casadi::Function |
|                  |                 | to be used to    | Internal         |
|                  |                 | populate         |                  |
|                  |                 | forward_options, |                  |
|                  |                 | reverse_options, |                  |
|                  |                 | and              |                  |
|                  |                 | jacobian_options |                  |
|                  |                 | before those     |                  |
|                  |                 | options are      |                  |
|                  |                 | merged in.       |                  |
+------------------+-----------------+------------------+------------------+
| derivative_of    | OT_FUNCTION     | The function is  | casadi::Function |
|                  |                 | a derivative of  | Internal         |
|                  |                 | another          |                  |
|                  |                 | function. The    |                  |
|                  |                 | type of          |                  |
|                  |                 | derivative       |                  |
|                  |                 | (directional     |                  |
|                  |                 | derivative,      |                  |
|                  |                 | Jacobian) is     |                  |
|                  |                 | inferred from    |                  |
|                  |                 | the function     |                  |
|                  |                 | name.            |                  |
+------------------+-----------------+------------------+------------------+
| dump             | OT_BOOL         | Dump function to | casadi::Function |
|                  |                 | file upon first  | Internal         |
|                  |                 | evaluation.      |                  |
|                  |                 | [false]          |                  |
+------------------+-----------------+------------------+------------------+
| dump_dir         | OT_STRING       | Directory to     | casadi::Function |
|                  |                 | dump             | Internal         |
|                  |                 | inputs/outputs   |                  |
|                  |                 | to. Make sure    |                  |
|                  |                 | the directory    |                  |
|                  |                 | exists [.]       |                  |
+------------------+-----------------+------------------+------------------+
| dump_format      | OT_STRING       | Choose file      | casadi::Function |
|                  |                 | format to dump   | Internal         |
|                  |                 | matrices. See    |                  |
|                  |                 | DM.from_file     |                  |
|                  |                 | [mtx]            |                  |
+------------------+-----------------+------------------+------------------+
| dump_in          | OT_BOOL         | Dump numerical   | casadi::Function |
|                  |                 | values of inputs | Internal         |
|                  |                 | to file          |                  |
|                  |                 | (readable with   |                  |
|                  |                 | DM.from_file )   |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| dump_out         | OT_BOOL         | Dump numerical   | casadi::Function |
|                  |                 | values of        | Internal         |
|                  |                 | outputs to file  |                  |
|                  |                 | (readable with   |                  |
|                  |                 | DM.from_file )   |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| enable_fd        | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation by   |                  |
|                  |                 | finite           |                  |
|                  |                 | differencing.    |                  |
|                  |                 | [default:        |                  |
|                  |                 | false]]          |                  |
+------------------+-----------------+------------------+------------------+
| enable_forward   | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | Jacobian-times-  |                  |
|                  |                 | vector products  |                  |
|                  |                 | - typically      |                  |
|                  |                 | using forward    |                  |
|                  |                 | mode AD - if     |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| enable_jacobian  | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | Jacobians of all |                  |
|                  |                 | differentiable   |                  |
|                  |                 | outputs with     |                  |
|                  |                 | respect to all   |                  |
|                  |                 | differentiable   |                  |
|                  |                 | inputs - if      |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| enable_reverse   | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | transposed       |                  |
|                  |                 | Jacobian-times-  |                  |
|                  |                 | vector products  |                  |
|                  |                 | - typically      |                  |
|                  |                 | using reverse    |                  |
|                  |                 | mode AD - if     |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| error_on_fail    | OT_BOOL         | Throw exceptions | casadi::ProtoFun |
|                  |                 | when function    | ction            |
|                  |                 | evaluation fails |                  |
|                  |                 | (default true).  |                  |
+------------------+-----------------+------------------+------------------+
| fd_method        | OT_STRING       | Method for       | casadi::Function |
|                  |                 | finite           | Internal         |
|                  |                 | differencing     |                  |
|                  |                 | [default         |                  |
|                  |                 | 'central']       |                  |
+------------------+-----------------+------------------+------------------+
| fd_options       | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to the    | Internal         |
|                  |                 | finite           |                  |
|                  |                 | difference       |                  |
|                  |                 | instance         |                  |
+------------------+-----------------+------------------+------------------+
| forward_options  | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | forward mode     |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| gather_stats     | OT_BOOL         | Deprecated       | casadi::Function |
|                  |                 | option           | Internal         |
|                  |                 | (ignored):       |                  |
|                  |                 | Statistics are   |                  |
|                  |                 | now always       |                  |
|                  |                 | collected.       |                  |
+------------------+-----------------+------------------+------------------+
| input_scheme     | OT_STRINGVECTOR | Deprecated       | casadi::Function |
|                  |                 | option (ignored) | Internal         |
+------------------+-----------------+------------------+------------------+
| inputs_check     | OT_BOOL         | Throw exceptions | casadi::Function |
|                  |                 | when the         | Internal         |
|                  |                 | numerical values |                  |
|                  |                 | of the inputs    |                  |
|                  |                 | don't make sense |                  |
+------------------+-----------------+------------------+------------------+
| is_diff_in       | OT_BOOLVECTOR   | Indicate for     | casadi::Function |
|                  |                 | each input if it | Internal         |
|                  |                 | should be        |                  |
|                  |                 | differentiable.  |                  |
+------------------+-----------------+------------------+------------------+
| is_diff_out      | OT_BOOLVECTOR   | Indicate for     | casadi::Function |
|                  |                 | each output if   | Internal         |
|                  |                 | it should be     |                  |
|                  |                 | differentiable.  |                  |
+------------------+-----------------+------------------+------------------+
| jac_penalty      | OT_DOUBLE       | When requested   | casadi::Function |
|                  |                 | for a number of  | Internal         |
|                  |                 | forward/reverse  |                  |
|                  |                 | directions, it   |                  |
|                  |                 | may be cheaper   |                  |
|                  |                 | to compute first |                  |
|                  |                 | the full         |                  |
|                  |                 | jacobian and     |                  |
|                  |                 | then multiply    |                  |
|                  |                 | with seeds,      |                  |
|                  |                 | rather than      |                  |
|                  |                 | obtain the       |                  |
|                  |                 | requested        |                  |
|                  |                 | directions in a  |                  |
|                  |                 | straightforward  |                  |
|                  |                 | manner. Casadi   |                  |
|                  |                 | uses a heuristic |                  |
|                  |                 | to decide which  |                  |
|                  |                 | is cheaper. A    |                  |
|                  |                 | high value of    |                  |
|                  |                 | 'jac_penalty'    |                  |
|                  |                 | makes it less    |                  |
|                  |                 | likely for the   |                  |
|                  |                 | heurstic to      |                  |
|                  |                 | chose the full   |                  |
|                  |                 | Jacobian         |                  |
|                  |                 | strategy. The    |                  |
|                  |                 | special value -1 |                  |
|                  |                 | indicates never  |                  |
|                  |                 | to use the full  |                  |
|                  |                 | Jacobian         |                  |
|                  |                 | strategy         |                  |
+------------------+-----------------+------------------+------------------+
| jacobian_options | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | Jacobian         |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| jit              | OT_BOOL         | Use just-in-time | casadi::Function |
|                  |                 | compiler to      | Internal         |
|                  |                 | speed up the     |                  |
|                  |                 | evaluation       |                  |
+------------------+-----------------+------------------+------------------+
| jit_cleanup      | OT_BOOL         | Cleanup up the   | casadi::Function |
|                  |                 | temporary source | Internal         |
|                  |                 | file that jit    |                  |
|                  |                 | creates.         |                  |
|                  |                 | Default: true    |                  |
+------------------+-----------------+------------------+------------------+
| jit_name         | OT_STRING       | The file name    | casadi::Function |
|                  |                 | used to write    | Internal         |
|                  |                 | out code. The    |                  |
|                  |                 | actual file      |                  |
|                  |                 | names used       |                  |
|                  |                 | depend on 'jit_t |                  |
|                  |                 | emp_suffix' and  |                  |
|                  |                 | include          |                  |
|                  |                 | extensions.      |                  |
|                  |                 | Default:         |                  |
|                  |                 | 'jit_tmp'        |                  |
+------------------+-----------------+------------------+------------------+
| jit_options      | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to the    | Internal         |
|                  |                 | jit compiler.    |                  |
+------------------+-----------------+------------------+------------------+
| jit_serialize    | OT_STRING       | Specify          | casadi::Function |
|                  |                 | behaviour when   | Internal         |
|                  |                 | serializing a    |                  |
|                  |                 | jitted function: |                  |
|                  |                 | SOURCE|link|embe |                  |
|                  |                 | d.               |                  |
+------------------+-----------------+------------------+------------------+
| jit_temp_suffix  | OT_BOOL         | Use a temporary  | casadi::Function |
|                  |                 | (seemingly       | Internal         |
|                  |                 | random) filename |                  |
|                  |                 | suffix for       |                  |
|                  |                 | generated code   |                  |
|                  |                 | and libraries.   |                  |
|                  |                 | This is desired  |                  |
|                  |                 | for thread-      |                  |
|                  |                 | safety. This     |                  |
|                  |                 | behaviour may    |                  |
|                  |                 | defeat caching   |                  |
|                  |                 | compiler         |                  |
|                  |                 | wrappers.        |                  |
|                  |                 | Default: true    |                  |
+------------------+-----------------+------------------+------------------+
| max_io           | OT_INT          | Acceptable       | casadi::Function |
|                  |                 | number of inputs | Internal         |
|                  |                 | and outputs.     |                  |
|                  |                 | Warn if          |                  |
|                  |                 | exceeded.        |                  |
+------------------+-----------------+------------------+------------------+
| max_num_dir      | OT_INT          | Specify the      | casadi::Function |
|                  |                 | maximum number   | Internal         |
|                  |                 | of directions    |                  |
|                  |                 | for derivative   |                  |
|                  |                 | functions.       |                  |
|                  |                 | Overrules the    |                  |
|                  |                 | builtin optimize |                  |
|                  |                 | d_num_dir.       |                  |
+------------------+-----------------+------------------+------------------+
| never_inline     | OT_BOOL         | Forbid inlining. | casadi::Function |
|                  |                 |                  | Internal         |
+------------------+-----------------+------------------+------------------+
| output_scheme    | OT_STRINGVECTOR | Deprecated       | casadi::Function |
|                  |                 | option (ignored) | Internal         |
+------------------+-----------------+------------------+------------------+
| post_expand      | OT_BOOL         | After            | casadi::Function |
|                  |                 | construction,    | Internal         |
|                  |                 | expand this      |                  |
|                  |                 | Function .       |                  |
|                  |                 | Default: False   |                  |
+------------------+-----------------+------------------+------------------+
| post_expand_opti | OT_DICT         | Options to be    | casadi::Function |
| ons              |                 | passed to post-  | Internal         |
|                  |                 | construction     |                  |
|                  |                 | expansion.       |                  |
|                  |                 | Default: empty   |                  |
+------------------+-----------------+------------------+------------------+
| print_in         | OT_BOOL         | Print numerical  | casadi::Function |
|                  |                 | values of inputs | Internal         |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| print_out        | OT_BOOL         | Print numerical  | casadi::Function |
|                  |                 | values of        | Internal         |
|                  |                 | outputs          |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| print_time       | OT_BOOL         | print            | casadi::ProtoFun |
|                  |                 | information      | ction            |
|                  |                 | about execution  |                  |
|                  |                 | time. Implies    |                  |
|                  |                 | record_time.     |                  |
+------------------+-----------------+------------------+------------------+
| record_time      | OT_BOOL         | record           | casadi::ProtoFun |
|                  |                 | information      | ction            |
|                  |                 | about execution  |                  |
|                  |                 | time, for        |                  |
|                  |                 | retrieval with   |                  |
|                  |                 | stats() .        |                  |
+------------------+-----------------+------------------+------------------+
| regularity_check | OT_BOOL         | Throw exceptions | casadi::ProtoFun |
|                  |                 | when NaN or Inf  | ction            |
|                  |                 | appears during   |                  |
|                  |                 | evaluation       |                  |
+------------------+-----------------+------------------+------------------+
| reverse_options  | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | reverse mode     |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| user_data        | OT_VOIDPTR      | A user-defined   | casadi::Function |
|                  |                 | field that can   | Internal         |
|                  |                 | be used to       |                  |
|                  |                 | identify the     |                  |
|                  |                 | function or pass |                  |
|                  |                 | additional       |                  |
|                  |                 | information      |                  |
+------------------+-----------------+------------------+------------------+
| verbose          | OT_BOOL         | Verbose          | casadi::ProtoFun |
|                  |                 | evaluation  for  | ction            |
|                  |                 | debugging        |                  |
+------------------+-----------------+------------------+------------------+

C++ includes: function.hpp
";

%feature("docstring") casadi::Function::buf_in "

[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L564

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L386-L398

>  std::vector< const double * > casadi::Function::buf_in(MapArg arg) const
------------------------------------------------------------------------
[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L564

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L386-L398

";

";

%feature("docstring") casadi::Function::buf_out "

[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L568

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L414-L427

>  std::vector< double * > casadi::Function::buf_out(MPrRes res) const
------------------------------------------------------------------------
[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L568

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L414-L427

";

";

%feature("docstring") casadi::Function::jit "

[INTERNAL] 
Create a just-in-time compiled function from a C language 
string.

The names and sparsity patterns of all the inputs and outputs must be 

provided. If sparsities are not provided, all inputs and outputs are 

assumed to be scalar. Only specify the function body, assuming that 
input 
and output nonzeros are stored in arrays with the specified 
naming 
convension. The data type used is 'casadi_real', which is 
typically equal 
to 'double or another data type with the same API as 'double.

Inputs may be null pointers. This means that the all entries are zero.
 
Outputs may be null points. This means that the corresponding result 
can be
 ignored.

If an error occurs in the evaluation, issue \"return 1;\";

The final generated function will have a structure similar to:

casadi_int fname(const casadi_real** arg, casadi_real** res, 
casadi_int* 
iw, casadi_real* w, void* mem) { const casadi_real *x1, 
*x2; casadi_real 
*r1, *r2; x1 = *arg++; x2 = *arg++; r1 = *res++; r2 =
 *res++; 
<FUNCTION_BODY> return 0; }

Extra doc: https://github.com/casadi/casadi/wiki/L_1v3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L262-L274

>  Function casadi::Function::jit(const std::string &name, const std::string &body, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const std::vector< Sparsity > &sparsity_in, const std::vector< Sparsity > &sparsity_out, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Create a just-in-time compiled function from a C language string.

The names and sparsity patterns of all the inputs and outputs must be 

provided. If sparsities are not provided, all inputs and outputs are 

assumed to be scalar. Only specify the function body, assuming that 
input 
and output nonzeros are stored in arrays with the specified 
naming 
convension. The data type used is 'casadi_real', which is 
typically equal 
to 'double or another data type with the same API as 'double.

Inputs may be null pointers. This means that the all entries are zero.
 
Outputs may be null points. This means that the corresponding result 
can be
 ignored.

If an error occurs in the evaluation, issue \"return 1;\";

The final generated function will have a structure similar to:

casadi_int fname(const casadi_real** arg, casadi_real** res, 
casadi_int* 
iw, casadi_real* w, void* mem) { const casadi_real *x1, 
*x2; casadi_real 
*r1, *r2; x1 = *arg++; x2 = *arg++; r1 = *res++; r2 =
 *res++; 
<FUNCTION_BODY> return 0; }

Extra doc: https://github.com/casadi/casadi/wiki/L_1v3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L262-L274

";

";

%feature("docstring") casadi::Function::Function "

[INTERNAL] 
Construct from a file.

Extra doc: https://github.com/casadi/casadi/wiki/L_1uz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L75

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L90-L92

>  casadi::Function::Function(const std::string &fname)
------------------------------------------------------------------------
[INTERNAL] 
Construct from a file.

Extra doc: https://github.com/casadi/casadi/wiki/L_1uz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L75

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L90-L92

";

";

%feature("docstring") casadi::Function::expand "

[INTERNAL] 
Expand a function to SX.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L207

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L286-L302

>  Function casadi::Function::expand(const std::string &name, const Dict &opts=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Expand a function to SX.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L207

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L286-L302

";

";

%feature("docstring") casadi::Function::size1_in "

[INTERNAL] 
Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240-L240

>  casadi_int casadi::Function::size1_in(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240-L240

";

";

%feature("docstring") casadi::Function::size2_in "

[INTERNAL] 
Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242-L242

>  casadi_int casadi::Function::size2_in(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242-L242

";

";

%feature("docstring") casadi::Function::size_in "

[INTERNAL] 
Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244-L246

>  std::pair<casadi_int, casadi_int> casadi::Function::size_in(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244-L246

";

";

%feature("docstring") casadi::Function::size1_out "

[INTERNAL] 
Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254-L254

>  casadi_int casadi::Function::size1_out(const std::string &oname) const
------------------------------------------------------------------------
[INTERNAL] 
Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254-L254

";

";

%feature("docstring") casadi::Function::size2_out "

[INTERNAL] 
Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256-L256

>  casadi_int casadi::Function::size2_out(const std::string &oname) const
------------------------------------------------------------------------
[INTERNAL] 
Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256-L256

";

";

%feature("docstring") casadi::Function::size_out "

[INTERNAL] 
Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258-L260

>  std::pair<casadi_int, casadi_int> casadi::Function::size_out(const std::string &oname) const
------------------------------------------------------------------------
[INTERNAL] 
Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258-L260

";

";

%feature("docstring") casadi::Function::nnz_in "

[INTERNAL] 
Get number of input nonzeros.

For a particular input or for all of the inputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L271

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L271-L271

>  casadi_int casadi::Function::nnz_in(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get number of input nonzeros.

For a particular input or for all of the inputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L271

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L271-L271

";

";

%feature("docstring") casadi::Function::nnz_out "

[INTERNAL] 
Get number of output nonzeros.

For a particular output or for all of the outputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L282

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L282-L282

>  casadi_int casadi::Function::nnz_out(const std::string &oname) const
------------------------------------------------------------------------
[INTERNAL] 
Get number of output nonzeros.

For a particular output or for all of the outputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L282

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L282-L282

";

";

%feature("docstring") casadi::Function::numel_in "

[INTERNAL] 
Get number of input elements.

For a particular input or for all of the inputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1ve

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L293-L293

>  casadi_int casadi::Function::numel_in(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get number of input elements.

For a particular input or for all of the inputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1ve

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L293-L293

";

";

%feature("docstring") casadi::Function::numel_out "

[INTERNAL] 
Get number of output elements.

For a particular output or for all of the outputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L304

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L304-L304

>  casadi_int casadi::Function::numel_out(const std::string &oname) const
------------------------------------------------------------------------
[INTERNAL] 
Get number of output elements.

For a particular output or for all of the outputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L304

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L304-L304

";

";

%feature("docstring") casadi::Function::sparsity_in "

[INTERNAL] 
Get sparsity of a given input.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L373

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L967-L973

>  const Sparsity & casadi::Function::sparsity_in(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get sparsity of a given input.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L373

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L967-L973

";

";

%feature("docstring") casadi::Function::sparsity_out "

[INTERNAL] 
Get sparsity of a given output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L381

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L983-L989

>  const Sparsity & casadi::Function::sparsity_out(const std::string &iname) const
------------------------------------------------------------------------
[INTERNAL] 
Get sparsity of a given output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L381

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L983-L989

";

";

%feature("docstring") casadi::Function::is_diff_in "

[INTERNAL] 
Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L390

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1007-L1013

>  std::vector< bool > casadi::Function::is_diff_in() const
------------------------------------------------------------------------
[INTERNAL] 
Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L390

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1007-L1013

";

";

%feature("docstring") casadi::Function::is_diff_out "

[INTERNAL] 
Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L391

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1015-L1021

>  std::vector< bool > casadi::Function::is_diff_out() const
------------------------------------------------------------------------
[INTERNAL] 
Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L391

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1015-L1021

";

";

%feature("docstring") casadi::Function::sparsity_jac "

[DEPRECATED] Get, if necessary generate, the sparsity of a Jacobian 
block

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L485

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L485-L488

>  const Sparsity casadi::Function::sparsity_jac(const std::string &iind, const std::string &oind, bool compact=false, bool symmetric=false) const
------------------------------------------------------------------------

[DEPRECATED] Get, if necessary generate, the sparsity of a Jacobian 
block

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L485

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L485-L488

";

";

%feature("docstring") casadi::Function::call "

[INTERNAL] 
Evaluate the function symbolically or numerically.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L524

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1414-L1421

>  void casadi::Function::call(const MXDict &arg, MXDict &res, bool always_inline=false, bool never_inline=false) const
------------------------------------------------------------------------
[INTERNAL] 
Evaluate the function symbolically or numerically.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L524

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1414-L1421

";

";

%feature("docstring") casadi::Function::call_gen "

[INTERNAL] 
Call using a map.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1129

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1362-L1376

>  void casadi::Function::call_gen(const std::map< std::string, M > &arg, std::map< std::string, M > &res, bool always_inline, bool never_inline) const
------------------------------------------------------------------------
[INTERNAL] 
Call using a map.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1129

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1362-L1376

";

";

%feature("docstring") casadi::Function::mapaccum "

[INTERNAL] 
Create a mapaccumulated version of this function.

Suppose the function has a signature of:

::

     f: (x, u) -> (x_next , y )
  



The the mapaccumulated version has the signature:

::

     F: (x0, U) -> (X , Y )
  
      with
          U: horzcat([u0, u1, ..., u_(N-1)])
          X: horzcat([x1, x2, ..., x_N])
          Y: horzcat([y0, y1, ..., y_(N-1)])
  
      and
          x1, y0 <- f(x0, u0)
          x2, y1 <- f(x1, u1)
          ...
          x_N, y_(N-1) <- f(x_(N-1), u_(N-1))
  



Mapaccum has the following benefits over writing an equivalent for-
loop:

much faster at construction time

potentially much faster compilation times (for codegen)

offers a trade-off between memory and evaluation time

The base (settable through the options dictionary, default 10), is 
used to 
create a tower of function calls, containing unrolled for-
loops of length 
maximum base.

This technique is much more scalable in terms of memory-usage, but 
slightly
 slower at evaluation, than a plain for-loop. The effect is 
similar to that
 of a for-loop with a check-pointing instruction after 
each chunk of 
iterations with size base.

Set base to -1 to unroll all the way; no gains in memory efficiency 
here.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L697

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L484-L486

>  Function casadi::Function::mapaccum(casadi_int N, const Dict &opts=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Create a mapaccumulated version of this function.

Suppose the function has a signature of:

::

     f: (x, u) -> (x_next , y )
  



The the mapaccumulated version has the signature:

::

     F: (x0, U) -> (X , Y )
  
      with
          U: horzcat([u0, u1, ..., u_(N-1)])
          X: horzcat([x1, x2, ..., x_N])
          Y: horzcat([y0, y1, ..., y_(N-1)])
  
      and
          x1, y0 <- f(x0, u0)
          x2, y1 <- f(x1, u1)
          ...
          x_N, y_(N-1) <- f(x_(N-1), u_(N-1))
  



Mapaccum has the following benefits over writing an equivalent for-
loop:

much faster at construction time

potentially much faster compilation times (for codegen)

offers a trade-off between memory and evaluation time

The base (settable through the options dictionary, default 10), is 
used to 
create a tower of function calls, containing unrolled for-
loops of length 
maximum base.

This technique is much more scalable in terms of memory-usage, but 
slightly
 slower at evaluation, than a plain for-loop. The effect is 
similar to that
 of a for-loop with a check-pointing instruction after 
each chunk of 
iterations with size base.

Set base to -1 to unroll all the way; no gains in memory efficiency 
here.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L697

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L484-L486

";

";

%feature("docstring") casadi::Function::fold "

[INTERNAL] 
Create a mapaccumulated version of this function.

Suppose the function has a signature of:

::

     f: (x, u) -> (x_next , y )
  



The the mapaccumulated version has the signature:

::

     F: (x0, U) -> (X , Y )
  
      with
          U: horzcat([u0, u1, ..., u_(N-1)])
          X: horzcat([x1, x2, ..., x_N])
          Y: horzcat([y0, y1, ..., y_(N-1)])
  
      and
          x1, y0 <- f(x0, u0)
          x2, y1 <- f(x1, u1)
          ...
          x_N, y_(N-1) <- f(x_(N-1), u_(N-1))
  



Mapaccum has the following benefits over writing an equivalent for-
loop:

much faster at construction time

potentially much faster compilation times (for codegen)

offers a trade-off between memory and evaluation time

The base (settable through the options dictionary, default 10), is 
used to 
create a tower of function calls, containing unrolled for-
loops of length 
maximum base.

This technique is much more scalable in terms of memory-usage, but 
slightly
 slower at evaluation, than a plain for-loop. The effect is 
similar to that
 of a for-loop with a check-pointing instruction after 
each chunk of 
iterations with size base.

Set base to -1 to unroll all the way; no gains in memory efficiency 
here.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L698

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L477-L483

";

%feature("docstring") casadi::Function::map "

[INTERNAL]

>  Function casadi::Function::map(casadi_int n, const std::string &parallelization, casadi_int max_num_threads) const
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::Function::generate_in "

[INTERNAL] 
Export an input file that can be passed to generate C code with 
a 
main.

See: 
 generate_out

See: 
 convert_in to convert between dict/map and vector

Extra doc: https://github.com/casadi/casadi/wiki/L_1wx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L855

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1176-L1186

>  std::vector< DM > casadi::Function::generate_in(const std::string &fname)
------------------------------------------------------------------------
[INTERNAL] 
Export an input file that can be passed to generate C code with a 
main.

See: 
 generate_out

See: 
 convert_in to convert between dict/map and vector

Extra doc: https://github.com/casadi/casadi/wiki/L_1wx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L855

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1176-L1186

";

";

%feature("docstring") casadi::Function::generate_out "

[INTERNAL] 
Export an output file that can be checked with generated C code
 
output.

See: 
 generate_in

See: 
 convert_out to convert between dict/map and vector

Extra doc: https://github.com/casadi/casadi/wiki/L_1wy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L866

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1188-L1198

>  std::vector< DM > casadi::Function::generate_out(const std::string &fname)
------------------------------------------------------------------------
[INTERNAL] 
Export an output file that can be checked with generated C code 
output.

See: 
 generate_in

See: 
 convert_out to convert between dict/map and vector

Extra doc: https://github.com/casadi/casadi/wiki/L_1wy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L866

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1188-L1198

";

";

%feature("docstring") casadi::Function::export_code "

[INTERNAL] 
Export function in specific language.

Only allowed for (a subset of) SX/MX Functions

Extra doc: https://github.com/casadi/casadi/wiki/L_1wz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L904

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1200-L1203

>  void casadi::Function::export_code(const std::string &lang, std::ostream &stream, const Dict &options=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Export function in specific language.

Only allowed for (a subset of) SX/MX Functions

Extra doc: https://github.com/casadi/casadi/wiki/L_1wz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L904

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1200-L1203

";

";

%feature("docstring") casadi::Function::serialize "

[INTERNAL] 
Serialize.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L893

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1217-L1221

>  std::string casadi::Function::serialize(const Dict &opts=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Serialize.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L893

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1217-L1221

";

";

%feature("docstring") casadi::Function::save "

[INTERNAL] 
Save  Function to a file.

See: 
 load

Extra doc: https://github.com/casadi/casadi/wiki/L_240

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L900

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1212-L1215

";

%feature("docstring") casadi::Function::sx_in "

[INTERNAL] 
Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L944

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1511-L1517

>  const std::vector< SX > casadi::Function::sx_in() const
------------------------------------------------------------------------
[INTERNAL] 
Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L944

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1511-L1517

";

";

%feature("docstring") casadi::Function::mx_in "

[INTERNAL] 
Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L949

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1535-L1537

>  const std::vector< MX > casadi::Function::mx_in() const
------------------------------------------------------------------------
[INTERNAL] 
Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L949

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1535-L1537

";

";

%feature("docstring") casadi::Function::sx_out "

[INTERNAL] 
Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L962

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1519-L1525

>  const std::vector< SX > casadi::Function::sx_out() const
------------------------------------------------------------------------
[INTERNAL] 
Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L962

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1519-L1525

";

";

%feature("docstring") casadi::Function::mx_out "

[INTERNAL] 
Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L967

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1539-L1541

>  const std::vector< MX > casadi::Function::mx_out() const
------------------------------------------------------------------------
[INTERNAL] 
Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L967

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1539-L1541

";

";

%feature("docstring") casadi::Function::nz_from_in "

[INTERNAL] 
Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L974

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1543-L1545

";

%feature("docstring") casadi::Function::nz_from_out "

[INTERNAL] 
Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L975

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1547-L1549

";

%feature("docstring") casadi::Function::nz_to_in "

[INTERNAL] 
Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L976

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1551-L1553

";

%feature("docstring") casadi::Function::nz_to_out "

[INTERNAL] 
Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L977

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1555-L1557

";

%feature("docstring") casadi::Function::convert_in "

[INTERNAL] 
Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L996

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1595-L1597

>  std::vector< MX > casadi::Function::convert_in(const MXDict &arg) const
------------------------------------------------------------------------
[INTERNAL] 
Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L996

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1595-L1597

";

";

%feature("docstring") casadi::Function::convert_out "

[INTERNAL] 
Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L998

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1603-L1605

>  std::vector< MX > casadi::Function::convert_out(const MXDict &arg) const
------------------------------------------------------------------------
[INTERNAL] 
Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L998

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1603-L1605

";

";

%feature("docstring") casadi::Function::has_spfwd "

[INTERNAL] 
Is the class able to propagate seeds through the algorithm?

Extra doc: https://github.com/casadi/casadi/wiki/L_1xl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1078

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1627-L1629

";

%feature("docstring") casadi::Function::has_sprev "

[INTERNAL] 
Is the class able to propagate seeds through the algorithm?

Extra doc: https://github.com/casadi/casadi/wiki/L_1xl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1079

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1631-L1633

";

%feature("docstring") casadi::Function::~Function "

[INTERNAL] 
Destructor.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L59-L60

";

%feature("docstring") casadi::Function::n_in "

[INTERNAL] 
Get the number of function inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L228

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L787-L789

";

%feature("docstring") casadi::Function::n_out "

[INTERNAL] 
Get the number of function outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L791-L793

";

%feature("docstring") casadi::Function::name_in "

[INTERNAL] 
Get input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L320

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L943-L949

>  const std::string & casadi::Function::name_in(casadi_int ind) const
------------------------------------------------------------------------
[INTERNAL] 
Get input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L320

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L943-L949

";

";

%feature("docstring") casadi::Function::name_out "

[INTERNAL] 
Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L325

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L951-L957

>  const std::string & casadi::Function::name_out(casadi_int ind) const
------------------------------------------------------------------------
[INTERNAL] 
Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L325

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L951-L957

";

";

%feature("docstring") casadi::Function::index_in "

[INTERNAL] 
Find the index for a string describing a particular entry of an 
input 
scheme.

example: schemeEntry(\"x_opt\") -> returns NLPSOL_X if 
FunctionInternal 
adheres to SCHEME_NLPINput

Extra doc: https://github.com/casadi/casadi/wiki/L_1vk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L333

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L927-L933

";

%feature("docstring") casadi::Function::index_out "

[INTERNAL] 
Find the index for a string describing a particular entry of an 
output
 scheme.

example: schemeEntry(\"x_opt\") -> returns NLPSOL_X if 
FunctionInternal 
adheres to SCHEME_NLPINput

Extra doc: https://github.com/casadi/casadi/wiki/L_1vl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L341

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L935-L941

";

%feature("docstring") casadi::Function::default_in "

[INTERNAL] 
Get default input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L346

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1423-L1425

";

%feature("docstring") casadi::Function::max_in "

[INTERNAL] 
Get largest input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L351

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1427-L1429

";

%feature("docstring") casadi::Function::min_in "

[INTERNAL] 
Get smallest input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L356

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1431-L1433

";

%feature("docstring") casadi::Function::nominal_in "

[INTERNAL] 
Get nominal input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L361

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1435-L1437

";

%feature("docstring") casadi::Function::nominal_out "

[INTERNAL] 
Get nominal output value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L366

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1439-L1441

";

%feature("docstring") casadi::Function::factory "

[INTERNAL] ";

%feature("docstring") casadi::Function::oracle "

[INTERNAL] 
Get oracle.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L407

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1846-L1852

";

%feature("docstring") casadi::Function::wrap "

[INTERNAL] 
Wrap in an  Function instance consisting of only one  MX call.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L412

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1854-L1856

";

%feature("docstring") casadi::Function::wrap_as_needed "

[INTERNAL] 
Wrap in a  Function with options.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L417

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1858-L1860

";

%feature("docstring") casadi::Function::which_depends "

[INTERNAL] 
Which variables enter with some order.

Parameters:
-----------

order: 
Only 1 (linear) and 2 (nonlinear) allowed

tr: 
Flip the relationship. Return which expressions contain the variables

Extra doc: https://github.com/casadi/casadi/wiki/L_1vx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L425

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1764-L1771

";

%feature("docstring") casadi::Function::print_dimensions "

[INTERNAL] 
Print dimensions of inputs and outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L432

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1095-L1097

";

%feature("docstring") casadi::Function::print_options "

[INTERNAL] 
Print options to a stream.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L437

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1099-L1101

";

%feature("docstring") casadi::Function::print_option "

[INTERNAL] 
Print all information there is to know about a certain option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L442

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1103-L1105

";

%feature("docstring") casadi::Function::has_option "

[INTERNAL] 
Does a particular option exist.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L447

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1107-L1114

";

%feature("docstring") casadi::Function::change_option "

[INTERNAL] 
Change option after object creation for debugging.

This is only possible for a selected number of options that do not 
change 
the numerical results of the computation, e.g. to enable a more
 verbose 
output or saving to file.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L455

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1116-L1126

";

%feature("docstring") casadi::Function::uses_output "

[INTERNAL] 
Do the derivative functions need nondifferentiated outputs?

Extra doc: https://github.com/casadi/casadi/wiki/L_1w3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L460

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L851-L853

";

%feature("docstring") casadi::Function::jacobian_old "

[DEPRECATED] Replaced by  Function::factory.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L466

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L856-L862

";

%feature("docstring") casadi::Function::hessian_old "

[DEPRECATED] Replaced by  Function::factory.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L471

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L864-L872

";

%feature("docstring") casadi::Function::jacobian "

[INTERNAL] 
Calculate all Jacobian blocks.

Generates a function that takes all non-differentiated inputs and 
outputs 
and calculates all Jacobian blocks. Inputs that are not needed
 by the 
routine are all-zero sparse matrices with the correct 
dimensions.  Output 
blocks that are not calculated, e.g. if the corresponding input or 
output 
is marked non-differentiated are also all-zero sparse. The 
Jacobian blocks 
are sorted starting by all the blocks for the first 
output, then all the 
blocks for the second output and so on. E.g. f : 
(x, y) -> (r, s) results 
in the function jac_f : (x, y, out_r, out_s) 
-> (jac_r_x, jac_r_y, jac_s_x,
 jac_s_y)

This function is cached.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L508

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L884-L890

";

%feature("docstring") casadi::Function::rev "

[INTERNAL] 
Propagate sparsity backward with temporary memory allocation.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L633

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L460-L475

>  int casadi::Function::rev(std::vector< bvec_t * > arg, std::vector< bvec_t * > res) const
------------------------------------------------------------------------
[INTERNAL] 
Propagate sparsity backward with temporary memory allocation.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L633

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L460-L475

";

";

%feature("docstring") casadi::Function::mapsum "

[INTERNAL] 
Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization: 
Type of parallelization used: unroll|serial|openmp

Extra doc: https://github.com/casadi/casadi/wiki/L_1wh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L642

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L724-L731

";

%feature("docstring") casadi::Function::slice "

[INTERNAL] 
returns a new function with a selection of inputs/outputs of the
 
original

Extra doc: https://github.com/casadi/casadi/wiki/L_1wl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L754

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L715-L722

";

%feature("docstring") casadi::Function::forward "

[INTERNAL] 
Get a function that calculates  nfwd forward derivatives.



::

     Returns a function with <tt>n_in + n_out + n_in</tt> inputs
     and <tt>nfwd</tt> outputs.
     The first <tt>n_in</tt> inputs correspond to nondifferentiated inputs.
     The next <tt>n_out</tt> inputs correspond to nondifferentiated outputs.
     and the last <tt>n_in</tt> inputs correspond to forward seeds,
     stacked horizontally
     The  <tt>n_out</tt> outputs correspond to forward sensitivities,
     stacked horizontally.     *
     <tt>(n_in = n_in(), n_out = n_out())</tt>
  
    The functions returned are cached, meaning that if called multiple timed
    with the same value, then multiple references to the same function will be returned.
  

Extra doc: https://github.com/casadi/casadi/wiki/L_1wq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L800

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1079-L1085

";

%feature("docstring") casadi::Function::reverse "

[INTERNAL] 
Get a function that calculates  nadj adjoint derivatives.



::

     Returns a function with <tt>n_in + n_out + n_out</tt> inputs
     and <tt>n_in</tt> outputs.
     The first <tt>n_in</tt> inputs correspond to nondifferentiated inputs.
     The next <tt>n_out</tt> inputs correspond to nondifferentiated outputs.
     and the last <tt>n_out</tt> inputs correspond to adjoint seeds,
     stacked horizontally
     The  <tt>n_in</tt> outputs correspond to adjoint sensitivities,
     stacked horizontally.     *
     <tt>(n_in = n_in(), n_out = n_out())</tt>
  
     <tt>(n_in = n_in(), n_out = n_out())</tt>
  
    The functions returned are cached, meaning that if called multiple timed
    with the same value, then multiple references to the same function will be returned.
  

Extra doc: https://github.com/casadi/casadi/wiki/L_1wr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L820

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1087-L1093

";

%feature("docstring") casadi::Function::jac_sparsity "

[INTERNAL] 
Get, if necessary generate, the sparsity of a single Jacobian 
block.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L830

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L911-L917

>  Sparsity casadi::Function::jac_sparsity(casadi_int oind, casadi_int iind, bool compact=false) const
------------------------------------------------------------------------
[INTERNAL] 
Get, if necessary generate, the sparsity of a single Jacobian block.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L830

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L911-L917

";

";

%feature("docstring") casadi::Function::generate "

[INTERNAL] 
Export / Generate C code for the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L840

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1132-L1134

>  std::string casadi::Function::generate(const Dict &opts=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Export / Generate C code for the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L840

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1132-L1134

";

";

%feature("docstring") casadi::Function::generate_dependencies "

[INTERNAL] 
Export / Generate C code for the dependency function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ww

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L845

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1142-L1144

";

%feature("docstring") casadi::Function::stats "

[INTERNAL] 
Get all statistics obtained at the end of the last evaluate 
call.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L932

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L896-L898

";

%feature("docstring") casadi::Function::has_free "

[INTERNAL] 
Does the function have free variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1004

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1635-L1637

";

%feature("docstring") casadi::Function::get_free "

[INTERNAL] 
Get free variables as a string.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1009

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1128-L1130

";

%feature("docstring") casadi::Function::free_sx "

[INTERNAL] 
Get all the free variables of the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1014

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1611-L1617

";

%feature("docstring") casadi::Function::free_mx "

[INTERNAL] 
Get all the free variables of the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1019

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1619-L1625

";

%feature("docstring") casadi::Function::generate_lifted "

[INTERNAL] 
Extract the functions needed for the Lifted  Newton method.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1024

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1639-L1645

";

%feature("docstring") casadi::Function::n_nodes "

[INTERNAL] 
Number of nodes in the algorithm.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1030

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1703-L1709

";

%feature("docstring") casadi::Function::n_instructions "

[INTERNAL] 
Number of instruction in the algorithm (SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xe

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1035

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1647-L1653

";

%feature("docstring") casadi::Function::instruction_id "

[INTERNAL] 
Identifier index of the instruction (SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1040

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1671-L1677

";

%feature("docstring") casadi::Function::instruction_input "

[INTERNAL] 
Locations in the work vector for the inputs of the instruction.

(SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1047

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1679-L1685

";

%feature("docstring") casadi::Function::instruction_constant "

[INTERNAL] 
Get the floating point output argument of an instruction 
(SXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1052

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1687-L1693

";

%feature("docstring") casadi::Function::instruction_output "

[INTERNAL] 
Location in the work vector for the output of the instruction.

(SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1059

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1695-L1701

";

%feature("docstring") casadi::Function::instruction_MX "

[INTERNAL] 
Get the  MX node corresponding to an instruction (MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1064

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1655-L1661

";

%feature("docstring") casadi::Function::instructions_sx "

[INTERNAL] 
Get the SX node corresponding to all instructions (SXFunction)

Note: input and output instructions have no SX representation. This 
method 
returns nan for those instructions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1072

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1663-L1669

";

%feature("docstring") casadi::Function::sz_arg "

[INTERNAL] 
Get required length of arg field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1085

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1027-L1027

";

%feature("docstring") casadi::Function::sz_res "

[INTERNAL] 
Get required length of res field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1090

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1029-L1029

";

%feature("docstring") casadi::Function::sz_iw "

[INTERNAL] 
Get required length of iw field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1095

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1031-L1031

";

%feature("docstring") casadi::Function::sz_w "

[INTERNAL] 
Get required length of w field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1033-L1033

";

%feature("docstring") casadi::Function::sz_work "

[INTERNAL] 
Get number of temporary variables needed.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1106

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1023-L1025

";

%feature("docstring") casadi::Function::set_work "

[INTERNAL] 
Set the (persistent) work vectors.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1111

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1052-L1059

";

%feature("docstring") casadi::Function::set_temp "

[INTERNAL] 
Set the (temporary) work vectors.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1117

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1061-L1068

";

%feature("docstring") casadi::Function::setup "

[INTERNAL] 
Set the (persistent and temporary) work vectors.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1123

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1070-L1077

";

%feature("docstring") casadi::Function::name "

[INTERNAL] 
Name of the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1137

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1250-L1257

";

%feature("docstring") casadi::Function::is_a "

[INTERNAL] 
Check if the function is of a particular type.

Optionally check if name matches one of the base classes (default 
true)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1144

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1607-L1609

";

%feature("docstring") casadi::Function::assert_size_in "

[INTERNAL] 
Assert that an input dimension is equal so some given value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1723-L1729

";

%feature("docstring") casadi::Function::assert_size_out "

[INTERNAL] 
Assert that an output dimension is equal so some given value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1192

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1731-L1736

";

%feature("docstring") casadi::Function::assert_sparsity_out "

[INTERNAL] 
Assert that an output sparsity is a multiple of some given 
sparsity.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1738-L1747

";

%feature("docstring") casadi::Function::checkout "

[INTERNAL] 
Checkout a memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1199

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1711-L1713

";

%feature("docstring") casadi::Function::release "

[INTERNAL] 
Release a memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1202

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1715-L1717

";

%feature("docstring") casadi::Function::memory "

[INTERNAL] 
Get memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1206

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1719-L1721

";

%feature("docstring") casadi::Function::cache "

[INTERNAL] 
Get all functions in the cache.

Extra doc: https://github.com/casadi/casadi/wiki/L_26i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1212

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1773-L1780

";

%feature("docstring") casadi::Function::get_function "

[INTERNAL] 
Get a dependency function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1791-L1797

>  Function casadi::Function::get_function(const std::string &name) const
------------------------------------------------------------------------
[INTERNAL] 
Get a dependency function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1791-L1797

";

";

%feature("docstring") casadi::Function::has_function "

[INTERNAL] 
Check if a particular dependency exists.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1227

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1799-L1806

";

%feature("docstring") casadi::Function::find_functions "

[INTERNAL] 
Get all functions embedded in the expression graphs.

Parameters:
-----------

max_depth: 
Maximum depth - a negative number indicates no maximum

Extra doc: https://github.com/casadi/casadi/wiki/L_1y6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1234

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1808-L1824

";

%feature("docstring") casadi::Function::find_function "

[INTERNAL] 
Get a specific function embedded in the expression graphs.

Parameters:
-----------

name: 
Name of function needed

max_depth: 
Maximum depth - a negative number indicates no maximum

Extra doc: https://github.com/casadi/casadi/wiki/L_1y7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1826-L1843

";

%feature("docstring") casadi::Function::info "

[INTERNAL] 
Obtain information about function

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1245

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1872-L1874

";

%feature("docstring") casadi::Function::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L99-L101

";

%feature("docstring") casadi::Function::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L103-L109

";

%feature("docstring") casadi::Function::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::Function::is_null "

[INTERNAL] 
Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L73-L75

";

%feature("docstring") casadi::Function::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L129-L131

";


// File: classcasadi_1_1FunctionBuffer.xml
%feature("docstring") casadi::FunctionBuffer "

[INTERNAL] 
Class to achieve minimal overhead function evaluations.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y9

C++ includes: function.hpp
";

%feature("docstring") casadi::FunctionBuffer::FunctionBuffer "

[INTERNAL]

>  casadi::FunctionBuffer::FunctionBuffer(const FunctionBuffer &f)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::FunctionBuffer::~FunctionBuffer "

[INTERNAL] ";

%feature("docstring") casadi::FunctionBuffer::set_arg "

[INTERNAL] 
Set input buffer for input i.

mem.set_arg(0, memoryview(a))

Note that CasADi uses 'fortran' order: column-by-column

Extra doc: https://github.com/casadi/casadi/wiki/L_1yb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1323

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1912-L1917

";

%feature("docstring") casadi::FunctionBuffer::set_res "

[INTERNAL] 
Set output buffer for ouput i.

mem.set_res(0, memoryview(a))

Note that CasADi uses 'fortran' order: column-by-column

Extra doc: https://github.com/casadi/casadi/wiki/L_1yc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1332

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1918-L1923

";

%feature("docstring") casadi::FunctionBuffer::ret "

[INTERNAL] 
Get last return value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1334

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1931-L1933

";

%feature("docstring") casadi::FunctionBuffer::_eval "

[INTERNAL] ";

%feature("docstring") casadi::FunctionBuffer::_self "

[INTERNAL] ";


// File: classcasadi_1_1GenericExpression.xml
%feature("docstring") casadi::GenericExpressionCommon "

[INTERNAL] 
Expression interface.

This is a common base class for SX,  MX and Matrix<>, introducing a uniform 
syntax and implementing common 
functionality using the curiously recurring 
template pattern (CRTP) 
idiom.
Joel Andersson

C++ includes: generic_expression.hpp
";

%feature("docstring") casadi::GenericExpressionCommon::minus "

[INTERNAL] 
Subtraction: (x,y) -> x - y.

Extra doc: https://github.com/casadi/casadi/wiki/L_oo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L86

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L86-L88

>  ExType casadi::GenericExpression::minus(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Subtraction: (x,y) -> x - y.

Extra doc: https://github.com/casadi/casadi/wiki/L_oo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L86

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L86-L88

";

";

%feature("docstring") casadi::GenericExpressionCommon::times "

[INTERNAL] 
Elementwise multiplication: (x,y) -> x .* y.

Extra doc: https://github.com/casadi/casadi/wiki/L_op

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L102

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L102-L104

>  ExType casadi::GenericExpression::times(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Elementwise multiplication: (x,y) -> x .* y.

Extra doc: https://github.com/casadi/casadi/wiki/L_op

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L102

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L102-L104

";

";

%feature("docstring") casadi::GenericExpressionCommon::rdivide "

[INTERNAL] 
Elementwise division: (x,y) -> x ./ y.

Extra doc: https://github.com/casadi/casadi/wiki/L_oq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L118

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L118-L120

>  ExType casadi::GenericExpression::rdivide(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Elementwise division: (x,y) -> x ./ y.

Extra doc: https://github.com/casadi/casadi/wiki/L_oq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L118

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L118-L120

";

";

%feature("docstring") casadi::GenericExpressionCommon::lt "

[INTERNAL] 
Logical less than: (x,y) -> x < y.

Extra doc: https://github.com/casadi/casadi/wiki/L_or

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L134

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L134-L136

>  ExType casadi::GenericExpression::lt(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Logical less than: (x,y) -> x < y.

Extra doc: https://github.com/casadi/casadi/wiki/L_or

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L134

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L134-L136

";

";

%feature("docstring") casadi::GenericExpressionCommon::le "

[INTERNAL] 
Logical less or equal to: (x,y) -> x <= y.

Extra doc: https://github.com/casadi/casadi/wiki/L_os

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L149

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L149-L151

>  ExType casadi::GenericExpression::le(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Logical less or equal to: (x,y) -> x <= y.

Extra doc: https://github.com/casadi/casadi/wiki/L_os

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L149

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L149-L151

";

";

%feature("docstring") casadi::GenericExpressionCommon::gt "

[INTERNAL] 
Logical greater than: (x,y) -> x > y.

Extra doc: https://github.com/casadi/casadi/wiki/L_ot

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L164

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L164-L166

>  ExType casadi::GenericExpression::gt(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Logical greater than: (x,y) -> x > y.

Extra doc: https://github.com/casadi/casadi/wiki/L_ot

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L164

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L164-L166

";

";

%feature("docstring") casadi::GenericExpressionCommon::ge "

[INTERNAL] 
Logical greater or equal to: (x,y) -> x <= y.

Extra doc: https://github.com/casadi/casadi/wiki/L_ou

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L179

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L179-L181

>  ExType casadi::GenericExpression::ge(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Logical greater or equal to: (x,y) -> x <= y.

Extra doc: https://github.com/casadi/casadi/wiki/L_ou

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L179

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L179-L181

";

";

%feature("docstring") casadi::GenericExpressionCommon::eq "

[INTERNAL] 
Logical equal to: (x,y) -> x == y.

Extra doc: https://github.com/casadi/casadi/wiki/L_ov

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L194

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L194-L196

>  ExType casadi::GenericExpression::eq(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Logical equal to: (x,y) -> x == y.

Extra doc: https://github.com/casadi/casadi/wiki/L_ov

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L194

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L194-L196

";

";

%feature("docstring") casadi::GenericExpressionCommon::ne "

[INTERNAL] 
Logical not equal to: (x,y) -> x != y.

Extra doc: https://github.com/casadi/casadi/wiki/L_ow

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L209

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L209-L211

>  ExType casadi::GenericExpression::ne(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Logical not equal to: (x,y) -> x != y.

Extra doc: https://github.com/casadi/casadi/wiki/L_ow

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L209

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L209-L211

";

";

%feature("docstring") casadi::GenericExpressionCommon::logic_and "

[INTERNAL] 
Logical  and

Returns (an expression evaluating to) 1 if both expressions are 
nonzero and
 0 otherwise

Extra doc: https://github.com/casadi/casadi/wiki/L_ox

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L227

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L227-L229

>  ExType casadi::GenericExpression::logic_and(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Logical  and

Returns (an expression evaluating to) 1 if both expressions are 
nonzero and
 0 otherwise

Extra doc: https://github.com/casadi/casadi/wiki/L_ox

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L227

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L227-L229

";

";

%feature("docstring") casadi::GenericExpressionCommon::logic_or "

[INTERNAL] 
Logical  or

returns (an expression evaluating to) 1 if at least one expression is 

nonzero and 0 otherwise

Extra doc: https://github.com/casadi/casadi/wiki/L_oy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L245

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L245-L247

>  ExType casadi::GenericExpression::logic_or(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Logical  or

returns (an expression evaluating to) 1 if at least one expression is 

nonzero and 0 otherwise

Extra doc: https://github.com/casadi/casadi/wiki/L_oy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L245

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L245-L247

";

";

%feature("docstring") casadi::GenericExpressionCommon::logic_not "

[INTERNAL] 
Logical  not x -> !x.

Returns (an expression evaluating to) 1 if expression is zero and 0 

otherwise

Extra doc: https://github.com/casadi/casadi/wiki/L_oz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L263

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L263-L265

>  ExType casadi::GenericExpression::logic_not(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Logical  not x -> !x.

Returns (an expression evaluating to) 1 if expression is zero and 0 

otherwise

Extra doc: https://github.com/casadi/casadi/wiki/L_oz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L263

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L263-L265

";

";

%feature("docstring") casadi::GenericExpressionCommon::abs "

[INTERNAL] 
Absolute value: x -> abs(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L278

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L278-L280

>  ExType casadi::GenericExpression::abs(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Absolute value: x -> abs(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L278

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L278-L280

";

";

%feature("docstring") casadi::GenericExpressionCommon::fabs "

[INTERNAL] 
Absolute value: x -> abs(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L281

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L281-L283

";

%feature("docstring") casadi::GenericExpressionCommon::sqrt "

[INTERNAL] 
Square root: x -> sqrt(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L293-L295

>  ExType casadi::GenericExpression::sqrt(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Square root: x -> sqrt(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L293-L295

";

";

%feature("docstring") casadi::GenericExpressionCommon::sq "

[INTERNAL] 
Square: x -> x^2.

Extra doc: https://github.com/casadi/casadi/wiki/L_p2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L305

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L305-L307

>  ExType casadi::GenericExpression::sq(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Square: x -> x^2.

Extra doc: https://github.com/casadi/casadi/wiki/L_p2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L305

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L305-L307

";

";

%feature("docstring") casadi::GenericExpressionCommon::sin "

[INTERNAL] 
Sine: x -> sin(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L317

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L317-L319

>  ExType casadi::GenericExpression::sin(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Sine: x -> sin(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L317

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L317-L319

";

";

%feature("docstring") casadi::GenericExpressionCommon::cos "

[INTERNAL] 
Cosine: x -> cos(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L329

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L329-L331

>  ExType casadi::GenericExpression::cos(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Cosine: x -> cos(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L329

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L329-L331

";

";

%feature("docstring") casadi::GenericExpressionCommon::tan "

[INTERNAL] 
Tangent: x -> tan(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L341

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L341-L343

>  ExType casadi::GenericExpression::tan(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Tangent: x -> tan(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L341

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L341-L343

";

";

%feature("docstring") casadi::GenericExpressionCommon::atan "

[INTERNAL] 
Arc tangent: x -> atan(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L353

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L353-L355

>  ExType casadi::GenericExpression::atan(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Arc tangent: x -> atan(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L353

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L353-L355

";

";

%feature("docstring") casadi::GenericExpressionCommon::asin "

[INTERNAL] 
Arc sine: x -> asin(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L365

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L365-L367

>  ExType casadi::GenericExpression::asin(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Arc sine: x -> asin(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L365

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L365-L367

";

";

%feature("docstring") casadi::GenericExpressionCommon::acos "

[INTERNAL] 
Arc cosine: x -> acos(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L377

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L377-L379

>  ExType casadi::GenericExpression::acos(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Arc cosine: x -> acos(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L377

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L377-L379

";

";

%feature("docstring") casadi::GenericExpressionCommon::tanh "

[INTERNAL] 
Hyperbolic tangent: x -> tanh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L389

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L389-L391

>  ExType casadi::GenericExpression::tanh(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Hyperbolic tangent: x -> tanh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L389

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L389-L391

";

";

%feature("docstring") casadi::GenericExpressionCommon::sinh "

[INTERNAL] 
Hyperbolic sin: x -> sinh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L401

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L401-L403

>  ExType casadi::GenericExpression::sinh(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Hyperbolic sin: x -> sinh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L401

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L401-L403

";

";

%feature("docstring") casadi::GenericExpressionCommon::cosh "

[INTERNAL] 
Hyperbolic cosine: x -> cosh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L413

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L413-L415

>  ExType casadi::GenericExpression::cosh(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Hyperbolic cosine: x -> cosh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L413

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L413-L415

";

";

%feature("docstring") casadi::GenericExpressionCommon::atanh "

[INTERNAL] 
Inverse hyperbolic tangent: x -> atanh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L425

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L425-L427

>  ExType casadi::GenericExpression::atanh(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Inverse hyperbolic tangent: x -> atanh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L425

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L425-L427

";

";

%feature("docstring") casadi::GenericExpressionCommon::asinh "

[INTERNAL] 
Inverse hyperbolic sin: x -> asinh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L437

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L437-L439

>  ExType casadi::GenericExpression::asinh(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Inverse hyperbolic sin: x -> asinh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L437

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L437-L439

";

";

%feature("docstring") casadi::GenericExpressionCommon::acosh "

[INTERNAL] 
Inverse hyperbolic cosine: x -> acosh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pe

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L449

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L449-L451

>  ExType casadi::GenericExpression::acosh(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Inverse hyperbolic cosine: x -> acosh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pe

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L449

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L449-L451

";

";

%feature("docstring") casadi::GenericExpressionCommon::exp "

[INTERNAL] 
Elementwise exponential: x -> exp(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L461

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L461-L463

>  ExType casadi::GenericExpression::exp(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Elementwise exponential: x -> exp(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L461

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L461-L463

";

";

%feature("docstring") casadi::GenericExpressionCommon::log "

[INTERNAL] 
Natural logarithm: x -> log(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L473

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L473-L475

>  ExType casadi::GenericExpression::log(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Natural logarithm: x -> log(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L473

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L473-L475

";

";

%feature("docstring") casadi::GenericExpressionCommon::log10 "

[INTERNAL] 
Base-10 logarithm: x -> log10(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_ph

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L485

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L485-L487

>  ExType casadi::GenericExpression::log10(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Base-10 logarithm: x -> log10(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_ph

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L485

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L485-L487

";

";

%feature("docstring") casadi::GenericExpressionCommon::log1p "

[INTERNAL] 
Precision variant for natural logarithm: x -> log(x+1)

Extra doc: https://github.com/casadi/casadi/wiki/L_pi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L497

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L497-L499

>  ExType casadi::GenericExpression::log1p(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Precision variant for natural logarithm: x -> log(x+1)

Extra doc: https://github.com/casadi/casadi/wiki/L_pi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L497

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L497-L499

";

";

%feature("docstring") casadi::GenericExpressionCommon::expm1 "

[INTERNAL] 
Precision variant for elementwise exponential: x -> exp(x)-1.

Extra doc: https://github.com/casadi/casadi/wiki/L_pj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L509

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L509-L511

>  ExType casadi::GenericExpression::expm1(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Precision variant for elementwise exponential: x -> exp(x)-1.

Extra doc: https://github.com/casadi/casadi/wiki/L_pj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L509

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L509-L511

";

";

%feature("docstring") casadi::GenericExpressionCommon::floor "

[INTERNAL] 
Round down to nearest integer: x -> floor(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L521

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L521-L523

>  ExType casadi::GenericExpression::floor(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Round down to nearest integer: x -> floor(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L521

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L521-L523

";

";

%feature("docstring") casadi::GenericExpressionCommon::ceil "

[INTERNAL] 
Round up to nearest integer: x -> ceil(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L533

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L533-L535

>  ExType casadi::GenericExpression::ceil(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Round up to nearest integer: x -> ceil(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L533

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L533-L535

";

";

%feature("docstring") casadi::GenericExpressionCommon::erf "

[INTERNAL] 
Error function: x -> erf(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L545

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L545-L547

>  ExType casadi::GenericExpression::erf(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Error function: x -> erf(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L545

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L545-L547

";

";

%feature("docstring") casadi::GenericExpressionCommon::erfinv "

[INTERNAL] 
Inverse error function: x -> erfinv(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L557

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L557-L559

>  ExType casadi::GenericExpression::erfinv(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Inverse error function: x -> erfinv(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L557

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L557-L559

";

";

%feature("docstring") casadi::GenericExpressionCommon::sign "

[INTERNAL] 
Sign function:

sign(x) := -1 for x<0 sign(x) := 1 for x>0, sign(0) := 0 sign(NaN) := 
NaN

Extra doc: https://github.com/casadi/casadi/wiki/L_po

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L574

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L574-L576

>  ExType casadi::GenericExpression::sign(const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Sign function:

sign(x) := -1 for x<0 sign(x) := 1 for x>0, sign(0) := 0 sign(NaN) := 
NaN

Extra doc: https://github.com/casadi/casadi/wiki/L_po

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L574

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L574-L576

";

";

%feature("docstring") casadi::GenericExpressionCommon::pow "

[INTERNAL] 
Elementwise power: (x,y) -> x.^y.

Extra doc: https://github.com/casadi/casadi/wiki/L_pp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L586

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L586-L588

>  ExType casadi::GenericExpression::pow(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Elementwise power: (x,y) -> x.^y.

Extra doc: https://github.com/casadi/casadi/wiki/L_pp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L586

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L586-L588

";

";

%feature("docstring") casadi::GenericExpressionCommon::mod "

[INTERNAL] 
Remainder after division: (x,y) -> fmod(x,y)

This  Function follows the convention of 
https://en.cppreference.com/w/c/numeric/math/fmod

Notably:
fmod(5,3) -> 2

fmod(5,-3) -> 2

fmod(-5,3) -> -2

fmod(-5,-3) -> -2

This is equivalent to Python's numpy.fmod and Matlab's rem.

\\\\seealso remainder

::

  Extra doc: https://github.com/casadi/casadi/wiki/L_pq 
  



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L610

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L610-L612

>  ExType casadi::GenericExpression::mod(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Remainder after division: (x,y) -> fmod(x,y)

This  Function follows the convention of 
https://en.cppreference.com/w/c/numeric/math/fmod

Notably:
fmod(5,3) -> 2

fmod(5,-3) -> 2

fmod(-5,3) -> -2

fmod(-5,-3) -> -2

This is equivalent to Python's numpy.fmod and Matlab's rem.

\\\\seealso remainder

::

  Extra doc: https://github.com/casadi/casadi/wiki/L_pq 
  



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L610

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L610-L612

";

";

%feature("docstring") casadi::GenericExpressionCommon::fmod "

[INTERNAL] 
Remainder after division: (x,y) -> fmod(x,y)

This  Function follows the convention of 
https://en.cppreference.com/w/c/numeric/math/fmod

Notably:
fmod(5,3) -> 2

fmod(5,-3) -> 2

fmod(-5,3) -> -2

fmod(-5,-3) -> -2

This is equivalent to Python's numpy.fmod and Matlab's rem.

\\\\seealso remainder

::

  Extra doc: https://github.com/casadi/casadi/wiki/L_pq 
  



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L613

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L613-L615

";

%feature("docstring") casadi::GenericExpressionCommon::remainder "

[INTERNAL] 
Remainder after division: (x,y) -> remainder(x,y)

This  Function follows the convention of 
https://en.cppreference.com/w/c/numeric/math/remainder

Notably:
remainder(5,3) -> -1

remainder(5,-3) -> -1

remainder(-5,3) -> 1

remainder(-5,-3) -> 1

This is equivalent to Python's math.remainder. There is no equivalence
 in 
Matlab.

\\\\seealso fmod

::

  Extra doc: https://github.com/casadi/casadi/wiki/L_24x 
  



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L637

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L637-L639

>  ExType casadi::GenericExpression::remainder(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Remainder after division: (x,y) -> remainder(x,y)

This  Function follows the convention of 
https://en.cppreference.com/w/c/numeric/math/remainder

Notably:
remainder(5,3) -> -1

remainder(5,-3) -> -1

remainder(-5,3) -> 1

remainder(-5,-3) -> 1

This is equivalent to Python's math.remainder. There is no equivalence
 in 
Matlab.

\\\\seealso fmod

::

  Extra doc: https://github.com/casadi/casadi/wiki/L_24x 
  



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L637

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L637-L639

";

";

%feature("docstring") casadi::GenericExpressionCommon::atan2 "

[INTERNAL] 
Two argument arc tangent: (y,x) -> atan2(y,x)

theta = atan2(y,x) corresponds to x = r cos(theta), y = r sin(theta)

Extra doc: https://github.com/casadi/casadi/wiki/L_pr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L651

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L651-L653

>  ExType casadi::GenericExpression::atan2(const ExType &y, const ExType &x)
------------------------------------------------------------------------
[INTERNAL] 
Two argument arc tangent: (y,x) -> atan2(y,x)

theta = atan2(y,x) corresponds to x = r cos(theta), y = r sin(theta)

Extra doc: https://github.com/casadi/casadi/wiki/L_pr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L651

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L651-L653

";

";

%feature("docstring") casadi::GenericExpressionCommon::if_else_zero "

[INTERNAL] 
Conditional assignment: (x,y) -> x ? y : 0.

Extra doc: https://github.com/casadi/casadi/wiki/L_ps

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L663

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L663-L665

>  ExType casadi::GenericExpression::if_else_zero(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Conditional assignment: (x,y) -> x ? y : 0.

Extra doc: https://github.com/casadi/casadi/wiki/L_ps

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L663

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L663-L665

";

";

%feature("docstring") casadi::GenericExpressionCommon::fmin "

[INTERNAL] 
Smallest of two values: (x,y) -> min(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L675

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L675-L677

>  ExType casadi::GenericExpression::fmin(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Smallest of two values: (x,y) -> min(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L675

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L675-L677

";

";

%feature("docstring") casadi::GenericExpressionCommon::fmax "

[INTERNAL] 
Largest of two values: (x,y) -> max(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L687

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L687-L689

>  ExType casadi::GenericExpression::fmax(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Largest of two values: (x,y) -> max(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L687

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L687-L689

";

";

%feature("docstring") casadi::GenericExpressionCommon::copysign "

[INTERNAL] 
Copy sign

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L713

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L713-L715

>  ExType casadi::GenericExpression::copysign(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Copy sign

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L713

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L713-L715

";

";

%feature("docstring") casadi::GenericExpressionCommon::constpow "

[INTERNAL] 
Elementwise power with const power

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L723

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L723-L725

>  ExType casadi::GenericExpression::constpow(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Elementwise power with const power

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L723

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L723-L725

";

";

%feature("docstring") casadi::GenericExpressionCommon::printme "

[INTERNAL] 
Debug printing

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L733

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L733-L735

>  ExType casadi::GenericExpression::printme(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Debug printing

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L733

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L733-L735

";

";

%feature("docstring") casadi::GenericExpressionCommon::hypot "

[INTERNAL] 
Precision variant for 2 norm: (x,y) -> sqrt(x^2+y^2)

Extra doc: https://github.com/casadi/casadi/wiki/L_pw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L745

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L745-L747

>  ExType casadi::GenericExpression::hypot(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Precision variant for 2 norm: (x,y) -> sqrt(x^2+y^2)

Extra doc: https://github.com/casadi/casadi/wiki/L_pw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L745

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L745-L747

";

";

%feature("docstring") casadi::GenericExpressionCommon::is_equal "

[INTERNAL] 
Check if two nodes are equivalent up to a given depth.

Depth=0 checks if the expressions are identical, i.e. points to the 
same 
node.

a = x*x b = x*x

is_equal(a,b,0) will return false, but a.is_equal(a,b,1) will return 
true

Extra doc: https://github.com/casadi/casadi/wiki/L_pv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L703

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L703-L705

";

%feature("docstring") casadi::GenericExpressionCommon::plus "

[INTERNAL] ";


// File: classcasadi_1_1GenericMatrix.xml


/*
 Construct symbolic primitives 
*/

/*
The \"sym\" function is intended to work in a similar way as \"sym\" 

used in the Symbolic Toolbox for Matlab but instead creating a CasADi 

symbolic primitive.

*/
%feature("docstring") casadi::GenericMatrixCommon "

[INTERNAL] 
 Matrix base class.

This is a common base class for  MX and Matrix<>, introducing a uniform 
syntax and implementing common 
functionality using the curiously recurring 
template pattern (CRTP) 
idiom.
 The class is designed with the idea that 
\"everything is a matrix\",
 that is, also scalars and vectors.
This 
philosophy makes it easy to use and to interface in particularly
 with 
Python and Matlab/Octave.
 The syntax tries to stay as close as possible to 
the ublas syntax 
when it comes to vector/matrix operations.
 Index starts 
with 0.
Index vec happens as follows: (rr, cc) -> k = rr+cc*size1()
Vectors 
are column vectors.
 The storage format is Compressed Column Storage (CCS), 
similar to 
that used for sparse matrices in Matlab, 
but unlike this 
format, we do allow for elements to be structurally 
non-zero but 
numerically zero.
 The sparsity pattern, which is reference counted and 
cached, can be 
accessed with Sparsity&  sparsity()
Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_1am

C++ includes: generic_matrix.hpp
";

%feature("docstring") casadi::GenericMatrixCommon::get_row "

[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L194

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L194-L194

";

%feature("docstring") casadi::GenericMatrixCommon::get_colind "

[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L195-L195

";

%feature("docstring") casadi::GenericMatrixCommon::row "

[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200-L200

>  casadi_int casadi::GenericMatrix< MatType >::row(casadi_int el) const
------------------------------------------------------------------------
[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200-L200

";

";

%feature("docstring") casadi::GenericMatrixCommon::colind "

[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201-L201

>  casadi_int casadi::GenericMatrix< MatType >::colind(casadi_int col) const
------------------------------------------------------------------------
[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201-L201

";

";

%feature("docstring") casadi::GenericMatrixCommon::interp1d "

[INTERNAL] 
Performs 1d linear interpolation.

The data-points to be interpolated are given as (x[i], v[i]). xq[j] is
 used
 as interplating value

Extra doc: https://github.com/casadi/casadi/wiki/L_1bh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L311

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L311-L314

>  MatType casadi::GenericMatrix::interp1d(const std::vector< double > &x, const MatType &v, const std::vector< double > &xq, const std::string &mode, bool equidistant=false)
------------------------------------------------------------------------
[INTERNAL] 
Performs 1d linear interpolation.

The data-points to be interpolated are given as (x[i], v[i]). xq[j] is
 used
 as interplating value

Extra doc: https://github.com/casadi/casadi/wiki/L_1bh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L311

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L311-L314

";

";

%feature("docstring") casadi::GenericMatrixCommon::sprank "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L215

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L215-L215

";

%feature("docstring") casadi::GenericMatrixCommon::norm_0_mul "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L216

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L216-L218

";

%feature("docstring") casadi::GenericMatrixCommon::tril "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L219

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L219-L221

";

%feature("docstring") casadi::GenericMatrixCommon::triu "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L222-L224

";

%feature("docstring") casadi::GenericMatrixCommon::sumsqr "

[INTERNAL] 
Calculate sum of squares: sum_ij X_ij^2.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L429

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L429-L431

>  MatType casadi::GenericMatrix::sumsqr(const MatType &x)
------------------------------------------------------------------------
[INTERNAL] 
Calculate sum of squares: sum_ij X_ij^2.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L429

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L429-L431

";

";

%feature("docstring") casadi::GenericMatrixCommon::linspace "

[INTERNAL] 
Matlab's  linspace command.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L460

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L460-L462

>  MatType casadi::GenericMatrix::linspace(const MatType &a, const MatType &b, casadi_int nsteps)
------------------------------------------------------------------------
[INTERNAL] 
Matlab's  linspace command.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L460

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L460-L462

";

";

%feature("docstring") casadi::GenericMatrixCommon::cross "

[INTERNAL] 
Matlab's  cross command.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L467

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L467-L469

>  MatType casadi::GenericMatrix::cross(const MatType &a, const MatType &b, casadi_int dim=-1)
------------------------------------------------------------------------
[INTERNAL] 
Matlab's  cross command.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L467

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L467-L469

";

";

%feature("docstring") casadi::GenericMatrixCommon::skew "

[INTERNAL] 
Generate a skew symmetric matrix from a 3-vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L474

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L474-L476

>  MatType casadi::GenericMatrix::skew(const MatType &a)
------------------------------------------------------------------------
[INTERNAL] 
Generate a skew symmetric matrix from a 3-vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L474

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L474-L476

";

";

%feature("docstring") casadi::GenericMatrixCommon::inv_skew "

[INTERNAL] 
Generate the 3-vector progenitor of a skew symmetric matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L481

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L481-L483

>  MatType casadi::GenericMatrix::inv_skew(const MatType &a)
------------------------------------------------------------------------
[INTERNAL] 
Generate the 3-vector progenitor of a skew symmetric matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L481

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L481-L483

";

";

%feature("docstring") casadi::GenericMatrixCommon::tril2symm "

[INTERNAL] 
Convert a lower triangular matrix to a symmetric one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L519

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L519-L519

>  MatType casadi::GenericMatrix::tril2symm(const MatType &a)
------------------------------------------------------------------------
[INTERNAL] 
Convert a lower triangular matrix to a symmetric one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L519

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L519-L519

";

";

%feature("docstring") casadi::GenericMatrixCommon::triu2symm "

[INTERNAL] 
Convert a upper triangular matrix to a symmetric one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L524

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L524-L524

>  MatType casadi::GenericMatrix::triu2symm(const MatType &a)
------------------------------------------------------------------------
[INTERNAL] 
Convert a upper triangular matrix to a symmetric one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L524

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L524-L524

";

";

%feature("docstring") casadi::GenericMatrixCommon::repsum "

[INTERNAL] 
Given a repeated matrix, computes the sum of repeated parts.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L980

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L980-L982

>  MatType casadi::GenericMatrix::repsum(const MatType &A, casadi_int n, casadi_int m=1)
------------------------------------------------------------------------
[INTERNAL] 
Given a repeated matrix, computes the sum of repeated parts.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L980

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L980-L982

";

";

%feature("docstring") casadi::GenericMatrixCommon::diff "

[INTERNAL] 
Returns difference (n-th order) along given axis (MATLAB 
convention)

Extra doc: https://github.com/casadi/casadi/wiki/L_1c8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L549

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L549-L551

>  MatType casadi::GenericMatrix::diff(const MatType &x, casadi_int n=1, casadi_int axis=-1)
------------------------------------------------------------------------
[INTERNAL] 
Returns difference (n-th order) along given axis (MATLAB convention)

Extra doc: https://github.com/casadi/casadi/wiki/L_1c8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L549

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L549-L551

";

";

%feature("docstring") casadi::GenericMatrixCommon::is_linear "

[INTERNAL] 
Is expr linear in var?

False negatives are possible (an expression may not be recognised as 
linear
 while it really is), false positives not.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L895

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L895-L897

>  bool casadi::GenericMatrix::is_linear(const MatType &expr, const MatType &var)
------------------------------------------------------------------------
[INTERNAL] 
Is expr linear in var?

False negatives are possible (an expression may not be recognised as 
linear
 while it really is), false positives not.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L895

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L895-L897

";

";

%feature("docstring") casadi::GenericMatrixCommon::is_quadratic "

[INTERNAL] 
Is expr quadratic in var?

False negatives are possible (an expression may not be recognised as 

quadratic while it really is), false positives not.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L906

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L906-L908

>  bool casadi::GenericMatrix::is_quadratic(const MatType &expr, const MatType &var)
------------------------------------------------------------------------
[INTERNAL] 
Is expr quadratic in var?

False negatives are possible (an expression may not be recognised as 

quadratic while it really is), false positives not.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L906

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L906-L908

";

";

%feature("docstring") casadi::GenericMatrixCommon::quadratic_coeff "

[INTERNAL] 
Recognizes quadratic form in scalar expression.

1/2*x' A x + b' x + c

e = 0.5*bilin(A,x,x)+dot(b,x)+c

Parameters:
-----------

check[in]: 
When true (default), A is checked to be independent of x. 
Provided to 
deal with false positive dependency checks.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L920

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L920-L923

>  void casadi::GenericMatrix::quadratic_coeff(const MatType &expr, const MatType &var, MatType &A, MatType &b, MatType &c, bool check=true)
------------------------------------------------------------------------
[INTERNAL] 
Recognizes quadratic form in scalar expression.

1/2*x' A x + b' x + c

e = 0.5*bilin(A,x,x)+dot(b,x)+c

Parameters:
-----------

check[in]: 
When true (default), A is checked to be independent of x. 
Provided to 
deal with false positive dependency checks.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L920

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L920-L923

";

";

%feature("docstring") casadi::GenericMatrixCommon::linear_coeff "

[INTERNAL] 
Recognizes linear form in vector expression.

A x + b

Parameters:
-----------

check[in]: 
When true (default)m, A is checked to be independent of x. 
Provided to
 deal with false positive dependency checks.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L933

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L933-L936

>  void casadi::GenericMatrix::linear_coeff(const MatType &expr, const MatType &var, MatType &A, MatType &b, bool check=true)
------------------------------------------------------------------------
[INTERNAL] 
Recognizes linear form in vector expression.

A x + b

Parameters:
-----------

check[in]: 
When true (default)m, A is checked to be independent of x. 
Provided to
 deal with false positive dependency checks.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L933

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L933-L936

";

";

%feature("docstring") casadi::GenericMatrixCommon::bilin "

[INTERNAL] 
Calculate bilinear/quadratic form x^T A y.

Parameters:
-----------

y: 
can be omitted, in which case x^T A x is calculated

Extra doc: https://github.com/casadi/casadi/wiki/L_1bo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L406

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L406-L408

>  MatType casadi::GenericMatrix::bilin(const MatType &A, const MatType &x)
------------------------------------------------------------------------
[INTERNAL] 
Calculate bilinear/quadratic form x^T A y.

Parameters:
-----------

y: 
can be omitted, in which case x^T A x is calculated

Extra doc: https://github.com/casadi/casadi/wiki/L_1bo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L406

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L406-L408

";

";

%feature("docstring") casadi::GenericMatrixCommon::rank1 "

[INTERNAL] 
Make a rank-1 update to a matrix A.

Calculates A + 1/2 * alpha * x*y'

Extra doc: https://github.com/casadi/casadi/wiki/L_1bp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L418

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L418-L421

>  MatType casadi::GenericMatrix< MatType >::rank1(const MatType &A, const MatType &alpha, const MatType &x, const MatType &y)
------------------------------------------------------------------------
[INTERNAL] 
Make a rank-1 update to a matrix A.

Calculates A + 1/2 * alpha * x*y'

Extra doc: https://github.com/casadi/casadi/wiki/L_1bp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L418

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L418-L421

";

";

%feature("docstring") casadi::GenericMatrixCommon::jtimes "

[INTERNAL] 
Calculate the Jacobian and multiply by a vector from the right.

This is equivalent to  mul(jacobian(ex, arg), v) or  mul(jacobian(ex, 
arg).T, v) for tr set to false and true respectively. If contrast to these 

expressions, it will use directional derivatives which is typically 
(but 
not necessarily) more efficient if the complete Jacobian is not 
needed and 
v has few rows.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L832

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L832-L835

>  MatType casadi::GenericMatrix< MatType >::jtimes(const MatType &ex, const MatType &arg, const MatType &v, bool tr=false, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Calculate the Jacobian and multiply by a vector from the right.

This is equivalent to  mul(jacobian(ex, arg), v) or  mul(jacobian(ex, 
arg).T, v) for tr set to false and true respectively. If contrast to these 

expressions, it will use directional derivatives which is typically 
(but 
not necessarily) more efficient if the complete Jacobian is not 
needed and 
v has few rows.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L832

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L832-L835

";

";

%feature("docstring") casadi::GenericMatrixCommon::gradient "

[INTERNAL] 
Calculate the gradient of an expression.

Parameters:
-----------

ex[in]: 
Scalar expression to take the gradient of

arg[in]: 
Vector expression of symbols

opts[in]: 
Options

Dense column vector

Extra doc: https://github.com/casadi/casadi/wiki/L_23x

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L812

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L812-L814

>  MatType casadi::GenericMatrix< MatType >::gradient(const MatType &ex, const MatType &arg, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Calculate the gradient of an expression.

Parameters:
-----------

ex[in]: 
Scalar expression to take the gradient of

arg[in]: 
Vector expression of symbols

opts[in]: 
Options

Dense column vector

Extra doc: https://github.com/casadi/casadi/wiki/L_23x

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L812

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L812-L814

";

";

%feature("docstring") casadi::GenericMatrixCommon::tangent "

[INTERNAL] 
Calculate the tangent of an expression.

Extra doc: https://github.com/casadi/casadi/wiki/L_23y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L819

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L819-L821

>  MatType casadi::GenericMatrix< MatType >::tangent(const MatType &ex, const MatType &arg, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Calculate the tangent of an expression.

Extra doc: https://github.com/casadi/casadi/wiki/L_23y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L819

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L819-L821

";

";

%feature("docstring") casadi::GenericMatrixCommon::linearize "

[INTERNAL] 
Linearize an expression.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L740

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L740-L743

>  MatType casadi::GenericMatrix< MatType >::linearize(const MatType &f, const MatType &x, const MatType &x0, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Linearize an expression.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L740

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L740-L743

";

";

%feature("docstring") casadi::GenericMatrixCommon::mpower "

[INTERNAL] 
 Matrix power x^n.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L319

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L319-L321

>  MatType casadi::GenericMatrix< MatType >::mpower(const MatType &x, const MatType &n)
------------------------------------------------------------------------
[INTERNAL] 
 Matrix power x^n.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L319

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L319-L321

";

";

%feature("docstring") casadi::GenericMatrixCommon::soc "

[INTERNAL] 
Construct second-order-convex.

Parameters:
-----------

x: 
vector expression of size n

y: 
scalar expression

soc(x,y) computes [y*eye(n) x; x' y]

soc(x,y) positive semi definite <=> || x ||_2 <= y

Extra doc: https://github.com/casadi/casadi/wiki/L_1bj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L334

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L334-L336

>  MatType casadi::GenericMatrix< MatType >::soc(const MatType &x, const MatType &y)
------------------------------------------------------------------------
[INTERNAL] 
Construct second-order-convex.

Parameters:
-----------

x: 
vector expression of size n

y: 
scalar expression

soc(x,y) computes [y*eye(n) x; x' y]

soc(x,y) positive semi definite <=> || x ||_2 <= y

Extra doc: https://github.com/casadi/casadi/wiki/L_1bj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L334

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L334-L336

";

";

%feature("docstring") casadi::GenericMatrixCommon::sym "

[INTERNAL] 
Create a vector of length r of vectors of length p.

with nrow-by-ncol symbolic primitives

Extra doc: https://github.com/casadi/casadi/wiki/L_1dg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1074

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1074-L1076

>  static std::vector<std::vector<MatType> > casadi::GenericMatrix< MatType >::sym(const std::string &name, casadi_int nrow, casadi_int ncol, casadi_int p, casadi_int r)
------------------------------------------------------------------------
[INTERNAL] 
Create a vector of length r of vectors of length p.

with nrow-by-ncol symbolic primitives

Extra doc: https://github.com/casadi/casadi/wiki/L_1dg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1074

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1074-L1076

";

";

%feature("docstring") casadi::GenericMatrixCommon::zeros "

[INTERNAL] 
Create a dense matrix or a matrix with specified sparsity with 
all 
entries zero.

Extra doc: https://github.com/casadi/casadi/wiki/L_1dh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1087

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1087-L1089

>  static MatType casadi::GenericMatrix< MatType >::zeros(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
Create a dense matrix or a matrix with specified sparsity with all 
entries zero.

Extra doc: https://github.com/casadi/casadi/wiki/L_1dh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1087

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1087-L1089

";

";

%feature("docstring") casadi::GenericMatrixCommon::ones "

[INTERNAL] 
Create a dense matrix or a matrix with specified sparsity with 
all 
entries one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1di

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100-L1102

>  static MatType casadi::GenericMatrix< MatType >::ones(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
Create a dense matrix or a matrix with specified sparsity with all 
entries one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1di

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100-L1102

";

";

%feature("docstring") casadi::GenericMatrixCommon::einstein "

[INTERNAL] 
Compute any contraction of two dense tensors, using 
index/einstein 
notation.

einstein(A, B, a, b, c) -> C

Given two tensors, A and B, computes a third tensor C such that:

C_c = A_a * B_b

With a, b, c representing einstein indices. Instead of the classical 
index 
labels i,j,k,... we employ -1,-2,-3,...

A, B, C are represented as CasADi vectors, with dim_a, dim_b, dim_c 

indictating theire tensorial dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L364

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L364-L370

>  MatType casadi::GenericMatrix::einstein(const MatType &A, const MatType &B, const std::vector< casadi_int > &dim_a, const std::vector< casadi_int > &dim_b, const std::vector< casadi_int > &dim_c, const std::vector< casadi_int > &a, const std::vector< casadi_int > &b, const std::vector< casadi_int > &c)
------------------------------------------------------------------------
[INTERNAL] 
Compute any contraction of two dense tensors, using index/einstein 
notation.

einstein(A, B, a, b, c) -> C

Given two tensors, A and B, computes a third tensor C such that:

C_c = A_a * B_b

With a, b, c representing einstein indices. Instead of the classical 
index 
labels i,j,k,... we employ -1,-2,-3,...

A, B, C are represented as CasADi vectors, with dim_a, dim_b, dim_c 

indictating theire tensorial dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L364

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L364-L370

";

";

%feature("docstring") casadi::GenericMatrixCommon::hessian "

[INTERNAL] 
Hessian and (optionally) gradient.

Extra doc: https://github.com/casadi/casadi/wiki/L_23z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L865

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L865-L868

>  MatType casadi::GenericMatrix::hessian(const MatType &ex, const MatType &arg, MatType &output_g, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Hessian and (optionally) gradient.

Extra doc: https://github.com/casadi/casadi/wiki/L_23z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L865

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L865-L868

";

";

%feature("docstring") casadi::GenericMatrixCommon::mmin "

[INTERNAL] 
Smallest element in a matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L988

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L988-L990

";

%feature("docstring") casadi::GenericMatrixCommon::mmax "

[INTERNAL] 
Largest element in a matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L997

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L997-L999

";

%feature("docstring") casadi::GenericMatrixCommon::nnz "

[INTERNAL] 
Get the number of (structural) non-zero elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1an

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L84

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1119-L1121

";

%feature("docstring") casadi::GenericMatrixCommon::nnz_lower "

[INTERNAL] 
Get the number of non-zeros in the lower triangular half.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ao

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L89

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1124-L1126

";

%feature("docstring") casadi::GenericMatrixCommon::nnz_upper "

[INTERNAL] 
Get the number of non-zeros in the upper triangular half.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ap

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L94

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1129-L1131

";

%feature("docstring") casadi::GenericMatrixCommon::nnz_diag "

[INTERNAL] 
Get get the number of non-zeros on the diagonal.

Extra doc: https://github.com/casadi/casadi/wiki/L_1aq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L99

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1134-L1136

";

%feature("docstring") casadi::GenericMatrixCommon::numel "

[INTERNAL] 
Get the number of elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ar

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L104

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1139-L1141

";

%feature("docstring") casadi::GenericMatrixCommon::size1 "

[INTERNAL] 
Get the first dimension (i.e. number of rows)

Extra doc: https://github.com/casadi/casadi/wiki/L_1as

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L109

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1144-L1146

";

%feature("docstring") casadi::GenericMatrixCommon::rows "

[INTERNAL] 
Get the number of rows, Octave-style syntax.

Extra doc: https://github.com/casadi/casadi/wiki/L_1at

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L114

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L114-L114

";

%feature("docstring") casadi::GenericMatrixCommon::size2 "

[INTERNAL] 
Get the second dimension (i.e. number of columns)

Extra doc: https://github.com/casadi/casadi/wiki/L_1au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L119

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1149-L1151

";

%feature("docstring") casadi::GenericMatrixCommon::columns "

[INTERNAL] 
Get the number of columns, Octave-style syntax.

Extra doc: https://github.com/casadi/casadi/wiki/L_1av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L124

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L124-L124

";

%feature("docstring") casadi::GenericMatrixCommon::dim "

[INTERNAL] 
Get string representation of dimensions.

The representation is e.g. \"4x5\" or \"4x5,10nz\"

Extra doc: https://github.com/casadi/casadi/wiki/L_1aw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L131

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1164-L1166

";

%feature("docstring") casadi::GenericMatrixCommon::size "

[INTERNAL] 
Get the size along a particular dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ay

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L141

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1159-L1161

>  casadi_int casadi::GenericMatrix< MatType >::size(casadi_int axis) const
------------------------------------------------------------------------
[INTERNAL] 
Get the size along a particular dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ay

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L141

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1159-L1161

";

";

%feature("docstring") casadi::GenericMatrixCommon::is_empty "

[INTERNAL] 
Check if the sparsity is empty, i.e. if one of the dimensions is
 zero.

(or optionally both dimensions)

Extra doc: https://github.com/casadi/casadi/wiki/L_1az

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L148

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L148-L148

";

%feature("docstring") casadi::GenericMatrixCommon::is_dense "

[INTERNAL] 
Check if the matrix expression is dense.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L153-L153

";

%feature("docstring") casadi::GenericMatrixCommon::is_scalar "

[INTERNAL] 
Check if the matrix expression is scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L158

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1169-L1171

";

%feature("docstring") casadi::GenericMatrixCommon::is_square "

[INTERNAL] 
Check if the matrix expression is square.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L163

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L163-L163

";

%feature("docstring") casadi::GenericMatrixCommon::is_vector "

[INTERNAL] 
Check if the matrix is a row or column vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L168

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L168-L168

";

%feature("docstring") casadi::GenericMatrixCommon::is_row "

[INTERNAL] 
Check if the matrix is a row vector (i.e.  size1()==1)

Extra doc: https://github.com/casadi/casadi/wiki/L_1b4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L173

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L173-L173

";

%feature("docstring") casadi::GenericMatrixCommon::is_column "

[INTERNAL] 
Check if the matrix is a column vector (i.e.  size2()==1)

Extra doc: https://github.com/casadi/casadi/wiki/L_1b5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L178

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L178-L178

";

%feature("docstring") casadi::GenericMatrixCommon::is_triu "

[INTERNAL] 
Check if the matrix is upper triangular.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L183

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L183-L183

";

%feature("docstring") casadi::GenericMatrixCommon::is_tril "

[INTERNAL] 
Check if the matrix is lower triangular.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L188

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L188-L188

";

%feature("docstring") casadi::GenericMatrixCommon::sparsity "

[INTERNAL] 
Get the sparsity pattern.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L207

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1114-L1116

";

%feature("docstring") casadi::GenericMatrixCommon::nz "

[INTERNAL] 
Access vector nonzero or slice of nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258-L260

>  NonZeros<MatType, K> casadi::GenericMatrix< MatType >::nz(const K &k)
------------------------------------------------------------------------
[INTERNAL] 
Access vector nonzero or slice of nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258-L260

";

";

%feature("docstring") casadi::GenericMatrixCommon::mrdivide "

[INTERNAL] 
 Matrix divide (cf. slash '/' in MATLAB)

Extra doc: https://github.com/casadi/casadi/wiki/L_1bl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L376

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L376-L378

";

%feature("docstring") casadi::GenericMatrixCommon::mldivide "

[INTERNAL] 
 Matrix divide (cf. backslash '\\\\' in MATLAB)

Extra doc: https://github.com/casadi/casadi/wiki/L_1bm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L383

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L383-L385

";

%feature("docstring") casadi::GenericMatrixCommon::symvar "

[INTERNAL] 
Get all symbols contained in the supplied expression.

Get all symbols on which the supplied expression depends 
See: 

SXFunction::getFree(), MXFunction::getFree()

Extra doc: https://github.com/casadi/casadi/wiki/L_1bn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L393

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L393-L395

";

%feature("docstring") casadi::GenericMatrixCommon::logsumexp "

[INTERNAL] 
Scaled version of logsumexp.

Scaled such that max(x) <= logsumexp(x, margin) <= max(x)+margin

Extra doc: https://github.com/casadi/casadi/wiki/L_1bs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L451

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L451-L454

>  MatType casadi::GenericMatrix::logsumexp(const MatType &x, const MatType &margin)
------------------------------------------------------------------------
[INTERNAL] 
Scaled version of logsumexp.

Scaled such that max(x) <= logsumexp(x, margin) <= max(x)+margin

Extra doc: https://github.com/casadi/casadi/wiki/L_1bs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L451

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L451-L454

";

";

%feature("docstring") casadi::GenericMatrixCommon::det "

[INTERNAL] 
 Matrix determinant (experimental)

Extra doc: https://github.com/casadi/casadi/wiki/L_1bx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L488

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L488-L488

";

%feature("docstring") casadi::GenericMatrixCommon::inv_minor "

[INTERNAL] 
 Matrix inverse (experimental)

Extra doc: https://github.com/casadi/casadi/wiki/L_1by

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L493

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L493-L493

";

%feature("docstring") casadi::GenericMatrixCommon::inv "

[INTERNAL] 
 Matrix inverse.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L505

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L505-L509

>  MatType casadi::GenericMatrix::inv(const MatType &A, const std::string &lsolver, const Dict &options=Dict())
------------------------------------------------------------------------
[INTERNAL] 
 Matrix inverse.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L505

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L505-L509

";

";

%feature("docstring") casadi::GenericMatrixCommon::trace "

[INTERNAL] 
 Matrix trace.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L514

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L514-L514

";

%feature("docstring") casadi::GenericMatrixCommon::norm_fro "

[INTERNAL] 
Frobenius norm.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L529

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L529-L529

";

%feature("docstring") casadi::GenericMatrixCommon::norm_2 "

[INTERNAL] 
2-norm

Extra doc: https://github.com/casadi/casadi/wiki/L_1c5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L534

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L534-L534

";

%feature("docstring") casadi::GenericMatrixCommon::norm_1 "

[INTERNAL] 
1-norm

Extra doc: https://github.com/casadi/casadi/wiki/L_1c6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L539

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L539-L539

";

%feature("docstring") casadi::GenericMatrixCommon::norm_inf "

[INTERNAL] 
Infinity-norm.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L544

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L544-L544

";

%feature("docstring") casadi::GenericMatrixCommon::cumsum "

[INTERNAL] 
Returns cumulative sum along given axis (MATLAB convention)

Extra doc: https://github.com/casadi/casadi/wiki/L_1c9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L556

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L556-L558

";

%feature("docstring") casadi::GenericMatrixCommon::dot "

[INTERNAL] 
Inner product of two matrices.

with x and y matrices of the same dimension

Extra doc: https://github.com/casadi/casadi/wiki/L_1ca

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L565

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L565-L567

";

%feature("docstring") casadi::GenericMatrixCommon::nullspace "

[INTERNAL] 
Computes the nullspace of a matrix A.

Finds Z m-by-(m-n) such that AZ = 0 with A n-by-m with m > n

Assumes A is full rank

Inspired by Numerical Methods in Scientific Computing by Ake Bjorck

Extra doc: https://github.com/casadi/casadi/wiki/L_1cb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L579

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L579-L581

";

%feature("docstring") casadi::GenericMatrixCommon::polyval "

[INTERNAL] 
Evaluate a polynomial with coefficients p in x.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L586

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L586-L588

";

%feature("docstring") casadi::GenericMatrixCommon::diag "

[INTERNAL] 
Get the diagonal of a matrix or construct a diagonal.

When the input is square, the diagonal elements are returned. If the 
input 
is vector-like, a diagonal matrix is constructed with it.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L596

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L596-L598

";

%feature("docstring") casadi::GenericMatrixCommon::unite "

[INTERNAL] 
Unite two matrices no overlapping sparsity.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ce

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L603

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L603-L605

";

%feature("docstring") casadi::GenericMatrixCommon::densify "

[INTERNAL] 
Make the matrix dense and assign nonzeros to a value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L617

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L617-L619

>  MatType casadi::GenericMatrix::densify(const MatType &x, const MatType &val)
------------------------------------------------------------------------
[INTERNAL] 
Make the matrix dense and assign nonzeros to a value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L617

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L617-L619

";

";

%feature("docstring") casadi::GenericMatrixCommon::project "

[INTERNAL] 
Create a new matrix with a given sparsity pattern but with the.

nonzeros taken from an existing matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_1ch

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L626

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L626-L629

";

%feature("docstring") casadi::GenericMatrixCommon::if_else "

[INTERNAL] 
Branching on  MX nodes.

Ternary operator, \"cond ? if_true : if_false\"

Extra doc: https://github.com/casadi/casadi/wiki/L_1ci

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L636

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L636-L639

";

%feature("docstring") casadi::GenericMatrixCommon::conditional "

[INTERNAL] 
Create a switch.

If the condition

Parameters:
-----------

ind: 
evaluates to the integer k, where 0<=k<f.size(), then x[k] will be 

returned, otherwise

x_default: 
will be returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L647

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L647-L650

";

%feature("docstring") casadi::GenericMatrixCommon::depends_on "

[INTERNAL] 
Check if expression depends on the argument.

The argument must be symbolic

Extra doc: https://github.com/casadi/casadi/wiki/L_1ck

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L657

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L657-L659

";

%feature("docstring") casadi::GenericMatrixCommon::substitute "

[INTERNAL] 
Substitute variable var with expression expr in multiple 
expressions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L673

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L673-L676

>  std::vector<MatType> casadi::GenericMatrix::substitute(const std::vector< MatType > &ex, const std::vector< MatType > &v, const std::vector< MatType > &vdef)
------------------------------------------------------------------------
[INTERNAL] 
Substitute variable var with expression expr in multiple expressions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L673

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L673-L676

";

";

%feature("docstring") casadi::GenericMatrixCommon::substitute_inplace "

[INTERNAL] 
Inplace substitution with piggyback expressions.

Substitute variables v out of the expressions vdef sequentially, as 
well as
 out of a number of other expressions piggyback

Extra doc: https://github.com/casadi/casadi/wiki/L_1cn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L685

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L685-L689

";

%feature("docstring") casadi::GenericMatrixCommon::cse "

[INTERNAL] 
Common subexpression elimination.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L702

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L702-L704

>  std::vector<MatType> casadi::GenericMatrix::cse(const std::vector< MatType > &e)
------------------------------------------------------------------------
[INTERNAL] 
Common subexpression elimination.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L702

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L702-L704

";

";

%feature("docstring") casadi::GenericMatrixCommon::solve "

[INTERNAL] 
 Solve a system of equations: A*x = b.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L731

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L731-L735

>  MatType casadi::GenericMatrix::solve(const MatType &A, const MatType &b, const std::string &lsolver, const Dict &dict=Dict())
------------------------------------------------------------------------
[INTERNAL] 
 Solve a system of equations: A*x = b.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L731

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L731-L735

";

";

%feature("docstring") casadi::GenericMatrixCommon::pinv "

[INTERNAL] 
Computes the Moore-Penrose pseudo-inverse.

If the matrix A is fat (size1>size2), mul(A, pinv(A)) is unity. If the
 
matrix A is slender (size2<size1), mul(pinv(A), A) is unity.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L767

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L767-L770

>  MatType casadi::GenericMatrix::pinv(const MatType &A, const std::string &lsolver, const Dict &dict=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Computes the Moore-Penrose pseudo-inverse.

If the matrix A is fat (size1>size2), mul(A, pinv(A)) is unity. If the
 
matrix A is slender (size2<size1), mul(pinv(A), A) is unity.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L767

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L767-L770

";

";

%feature("docstring") casadi::GenericMatrixCommon::expm_const "

[INTERNAL] 
Calculate  Matrix exponential.

Computes expm(A*t) with A constant

Parameters:
-----------

A[in]: 
Square matrix

t[in]: 
Scalar

Extra doc: https://github.com/casadi/casadi/wiki/L_23v

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L782

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L782-L784

";

%feature("docstring") casadi::GenericMatrixCommon::expm "

[INTERNAL] 
Calculate  Matrix exponential.

Extra doc: https://github.com/casadi/casadi/wiki/L_23w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L790

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L790-L792

";

%feature("docstring") casadi::GenericMatrixCommon::jacobian "

[INTERNAL] 
Calculate Jacobian.

Sparse matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_1cv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L799

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L799-L802

";

%feature("docstring") casadi::GenericMatrixCommon::forward "

[INTERNAL] 
Forward directional derivative.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L841

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L841-L845

";

%feature("docstring") casadi::GenericMatrixCommon::reverse "

[INTERNAL] 
Reverse directional derivative.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L851

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L851-L855

";

%feature("docstring") casadi::GenericMatrixCommon::which_depends "

[INTERNAL] 
Find out which variables enter with some order.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L874

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L874-L877

";

%feature("docstring") casadi::GenericMatrixCommon::jacobian_sparsity "

[INTERNAL] 
Get the sparsity pattern of a jacobian.

Equivalent to, but cheaper to compute than, jacobian(f,x). sparsity()

Extra doc: https://github.com/casadi/casadi/wiki/L_259

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L884

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L884-L886

";

%feature("docstring") casadi::GenericMatrixCommon::n_nodes "

[INTERNAL] 
Count number of nodes

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L939

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L939-L941

";

%feature("docstring") casadi::GenericMatrixCommon::simplify "

[INTERNAL] 
Simplify an expression.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L944

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L944-L946

";

%feature("docstring") casadi::GenericMatrixCommon::print_operator "

[INTERNAL] 
Get a string representation for a binary MatType, using custom 

arguments.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L952

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L952-L954

";

%feature("docstring") casadi::GenericMatrixCommon::extract "

[INTERNAL] 
Introduce intermediate variables for selected nodes in a graph.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L959

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L959-L964

";

%feature("docstring") casadi::GenericMatrixCommon::shared "

[INTERNAL] 
Extract shared subexpressions from an set of expressions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L969

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L969-L975

";


// File: classcasadi_1_1GenericType.xml
%feature("docstring") casadi::GenericType "

[INTERNAL] 
Generic data type, can hold different types such as bool, 
casadi_int, 
std::string etc.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_17n

C++ includes: generic_type.hpp
";

%feature("docstring") casadi::GenericType::is_bool "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L161

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L134-L136

";

%feature("docstring") casadi::GenericType::is_int "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L162

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L138-L140

";

%feature("docstring") casadi::GenericType::is_double "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L163

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L142-L144

";

%feature("docstring") casadi::GenericType::is_string "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L164

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L146-L148

";

%feature("docstring") casadi::GenericType::is_empty_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L165

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L150-L157

";

%feature("docstring") casadi::GenericType::is_int_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L166

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L159-L161

";

%feature("docstring") casadi::GenericType::is_int_vector_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L167

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L167-L169

";

%feature("docstring") casadi::GenericType::is_double_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L168

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L171-L173

";

%feature("docstring") casadi::GenericType::is_double_vector_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L169

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L175-L177

";

%feature("docstring") casadi::GenericType::is_bool_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L170

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L163-L165

";

%feature("docstring") casadi::GenericType::is_string_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L171

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L179-L181

";

%feature("docstring") casadi::GenericType::is_dict "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L172

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L195-L197

";

%feature("docstring") casadi::GenericType::is_function "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L173

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L183-L185

";

%feature("docstring") casadi::GenericType::is_function_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L174

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L187-L189

";

%feature("docstring") casadi::GenericType::is_void_pointer "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L175

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L191-L193

";

%feature("docstring") casadi::GenericType::as_bool "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L182

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L262-L265

";

%feature("docstring") casadi::GenericType::as_int "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L183

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L267-L270

";

%feature("docstring") casadi::GenericType::as_double "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L184

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L272-L275

";

%feature("docstring") casadi::GenericType::as_string "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L185

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L277-L280

";

%feature("docstring") casadi::GenericType::as_int_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L186

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L282-L285

";

%feature("docstring") casadi::GenericType::as_bool_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L187

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L287-L290

";

%feature("docstring") casadi::GenericType::as_int_vector_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L188

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L292-L295

";

%feature("docstring") casadi::GenericType::as_double_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L297-L300

";

%feature("docstring") casadi::GenericType::as_double_vector_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L190

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L302-L305

";

%feature("docstring") casadi::GenericType::as_string_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L191

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L307-L310

";

%feature("docstring") casadi::GenericType::as_dict "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L192

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L312-L315

";

%feature("docstring") casadi::GenericType::as_function "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L193

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L317-L320

";

%feature("docstring") casadi::GenericType::as_function_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L194

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L322-L325

";

%feature("docstring") casadi::GenericType::as_void_pointer "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L327-L330

";

%feature("docstring") casadi::GenericType::to_bool "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L332-L341

";

%feature("docstring") casadi::GenericType::to_int "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L343-L352

";

%feature("docstring") casadi::GenericType::to_double "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L202

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L354-L361

";

%feature("docstring") casadi::GenericType::to_string "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L203

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L363-L366

";

%feature("docstring") casadi::GenericType::to_int_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L204

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L373-L376

";

%feature("docstring") casadi::GenericType::to_bool_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L205

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L384-L393

";

%feature("docstring") casadi::GenericType::to_int_vector_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L206

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L395-L398

";

%feature("docstring") casadi::GenericType::to_double_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L207

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L400-L408

";

%feature("docstring") casadi::GenericType::to_double_vector_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L208

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L410-L421

";

%feature("docstring") casadi::GenericType::to_string_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L209

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L423-L439

";

%feature("docstring") casadi::GenericType::to_dict "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L210

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L441-L444

";

%feature("docstring") casadi::GenericType::to_function "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L211

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L446-L449

";

%feature("docstring") casadi::GenericType::to_function_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L212

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L451-L454

";

%feature("docstring") casadi::GenericType::to_void_pointer "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L213

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L456-L464

";

%feature("docstring") casadi::GenericType::to_int_type_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L214

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L368-L371

";

%feature("docstring") casadi::GenericType::GenericType "

[INTERNAL]

>  casadi::GenericType::GenericType(casadi_int i)

>  casadi::GenericType::GenericType(int i)

>  casadi::GenericType::GenericType(double d)

>  casadi::GenericType::GenericType(const std::string &s)

>  casadi::GenericType::GenericType(const std::vector< bool > &iv)

>  casadi::GenericType::GenericType(const std::vector< casadi_int > &iv)

>  casadi::GenericType::GenericType(const std::vector< int > &iv)

>  casadi::GenericType::GenericType(const std::vector< std::vector< casadi_int > > &ivv)

>  casadi::GenericType::GenericType(const std::vector< double > &dv)

>  casadi::GenericType::GenericType(const std::vector< std::vector< double > > &dv)

>  casadi::GenericType::GenericType(const std::vector< std::string > &sv)

>  casadi::GenericType::GenericType(const char s[])

>  casadi::GenericType::GenericType(const Function &f)

>  casadi::GenericType::GenericType(const std::vector< Function > &f)

>  casadi::GenericType::GenericType(const Dict &dict)

>  casadi::GenericType::GenericType(void *ptr)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::GenericType::get_description "

[INTERNAL] 
Get a description of the object's type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L120

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L120-L120

";

%feature("docstring") casadi::GenericType::can_cast_to "

[INTERNAL] ";

%feature("docstring") casadi::GenericType::getType "

[INTERNAL] ";

%feature("docstring") casadi::GenericType::serialize "

[INTERNAL] 
Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_17r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L225

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L547-L550

";

%feature("docstring") casadi::GenericType::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L99-L101

";

%feature("docstring") casadi::GenericType::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L103-L109

";

%feature("docstring") casadi::GenericType::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::GenericType::is_null "

[INTERNAL] 
Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L73-L75

";

%feature("docstring") casadi::GenericType::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L129-L131

";


// File: classcasadi_1_1GlobalOptions.xml
%feature("docstring") casadi::GlobalOptions "

[INTERNAL] 
Collects global CasADi options.

Note to developers: 
use sparingly. Global options are - in general - a 
rather bad idea

this class must never be instantiated. Access its static members 
directly 

Joris Gillis

Extra doc: https://github.com/casadi/casadi/wiki/L_23m

C++ includes: global_options.hpp
";


// File: classcasadi_1_1ImplicitFixedStepIntegrator.xml
%feature("docstring") casadi::ImplicitFixedStepIntegrator "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1ImplicitToNlp.xml
%feature("docstring") casadi::ImplicitToNlp "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Importer.xml
%feature("docstring") casadi::Importer "

[INTERNAL] 
 Importer.

Just-in-time compilation of code
General informationList of plugins
- clang

- shell

Note: some of the plugins in this list might not be available on your 

system.  Also, there might be extra plugins available to you that are 
not 
listed here. You can obtain their documentation with   
Importer.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

clang
-----



Interface to the JIT compiler CLANG

Extra doc: https://github.com/casadi/casadi/wiki/L_21v

>List of available options

+--------------+-----------------+-----------------------------------------+
|      Id      |      Type       |               Description               |
+==============+=================+=========================================+
| flags        | OT_STRINGVECTOR | Compile flags for the JIT compiler.     |
|              |                 | Default: None                           |
+--------------+-----------------+-----------------------------------------+
| include_path | OT_STRING       | Include paths for the JIT compiler. The |
|              |                 | include directory shipped with CasADi   |
|              |                 | will be automatically appended.         |
+--------------+-----------------+-----------------------------------------+



--------------------------------------------------------------------------------

shell
-----



Interface to the JIT compiler SHELL

Extra doc: https://github.com/casadi/casadi/wiki/L_22w

>List of available options

+----------------------+-----------------+---------------------------------+
|          Id          |      Type       |           Description           |
+======================+=================+=================================+
| cleanup              | OT_BOOL         | Cleanup temporary files when    |
|                      |                 | unloading. Default: true        |
+----------------------+-----------------+---------------------------------+
| compiler             | OT_STRING       | Compiler command                |
+----------------------+-----------------+---------------------------------+
| compiler_flags       | OT_STRINGVECTOR | Alias for 'compiler_flags'      |
+----------------------+-----------------+---------------------------------+
| compiler_output_flag | OT_STRING       | Compiler flag to denote object  |
|                      |                 | output. Default: '-o '          |
+----------------------+-----------------+---------------------------------+
| compiler_setup       | OT_STRING       | Compiler setup command.         |
|                      |                 | Intended to be fixed. The       |
|                      |                 | 'flag' option is the prefered   |
|                      |                 | way to set custom flags.        |
+----------------------+-----------------+---------------------------------+
| directory            | OT_STRING       | Directory to put temporary      |
|                      |                 | objects in. Must end with a     |
|                      |                 | file separator.                 |
+----------------------+-----------------+---------------------------------+
| extra_suffixes       | OT_STRINGVECTOR | List of suffixes for extra      |
|                      |                 | files that the compiler may     |
|                      |                 | generate. Default: None         |
+----------------------+-----------------+---------------------------------+
| flags                | OT_STRINGVECTOR | Compile flags for the JIT       |
|                      |                 | compiler. Default: None         |
+----------------------+-----------------+---------------------------------+
| linker               | OT_STRING       | Linker command                  |
+----------------------+-----------------+---------------------------------+
| linker_flags         | OT_STRINGVECTOR | Linker flags for the JIT        |
|                      |                 | compiler. Default: None         |
+----------------------+-----------------+---------------------------------+
| linker_output_flag   | OT_STRING       | Linker flag to denote shared    |
|                      |                 | library output. Default: '-o '  |
+----------------------+-----------------+---------------------------------+
| linker_setup         | OT_STRING       | Linker setup command. Intended  |
|                      |                 | to be fixed. The 'flag' option  |
|                      |                 | is the prefered way to set      |
|                      |                 | custom flags.                   |
+----------------------+-----------------+---------------------------------+
| name                 | OT_STRING       | The file name used to write out |
|                      |                 | compiled objects/libraries. The |
|                      |                 | actual file names used depend   |
|                      |                 | on 'temp_suffix' and include    |
|                      |                 | extensions. Default:            |
|                      |                 | 'tmp_casadi_compiler_shell'     |
+----------------------+-----------------+---------------------------------+
| temp_suffix          | OT_BOOL         | Use a temporary (seemingly      |
|                      |                 | random) filename suffix for     |
|                      |                 | file names. This is desired for |
|                      |                 | thread-safety. This behaviour   |
|                      |                 | may defeat caching compiler     |
|                      |                 | wrappers. Default: true         |
+----------------------+-----------------+---------------------------------+

Joris Gillis

Extra doc: https://github.com/casadi/casadi/wiki/L_161

C++ includes: importer.hpp
";

%feature("docstring") casadi::Importer::Importer "

[INTERNAL] 
 Importer factory.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L97

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L34-L45

>  casadi::Importer::Importer(const std::string &name, const std::string &compiler, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
 Importer factory.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L97

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L34-L45

";

";

%feature("docstring") casadi::Importer::plugin_name "

[INTERNAL] 
Query plugin name.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L118

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L71-L73

";

%feature("docstring") casadi::Importer::has_function "

[INTERNAL] ";

%feature("docstring") casadi::Importer::get_function "

[INTERNAL] 
Get a function pointer for numerical evaluation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L125

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L79-L81

";

%feature("docstring") casadi::Importer::has_meta "

[INTERNAL] 
Does a meta entry exist?

Extra doc: https://github.com/casadi/casadi/wiki/L_165

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L145

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L83-L85

";

%feature("docstring") casadi::Importer::get_meta "

[INTERNAL] 
Get entry as a text.

Extra doc: https://github.com/casadi/casadi/wiki/L_166

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L87-L89

";

%feature("docstring") casadi::Importer::inlined "

[INTERNAL] 
Check if a function is inlined.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L91-L93

";

%feature("docstring") casadi::Importer::body "

[INTERNAL] 
Get the function body, if inlined.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L156

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L95-L97

";

%feature("docstring") casadi::Importer::library "

[INTERNAL] 
Get library name.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L159

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L99-L101

";

%feature("docstring") casadi::Importer::to "

[INTERNAL] 
Convert to a type.

Extra doc: https://github.com/casadi/casadi/wiki/L_167

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L173

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L173-L175

";

%feature("docstring") casadi::Importer::meta_string "

[INTERNAL] 
Get entry as a string.

Extra doc: https://github.com/casadi/casadi/wiki/L_168

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L180

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L180-L182

";

%feature("docstring") casadi::Importer::meta_vector "

[INTERNAL] 
Get entry as a vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_169

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L188

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L188-L190

";

%feature("docstring") casadi::Importer::meta_set "

[INTERNAL] 
Get entry as a set.

Extra doc: https://github.com/casadi/casadi/wiki/L_16a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L196

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L196-L198

";

%feature("docstring") casadi::Importer::meta_int "

[INTERNAL] 
Get entry as an integer.

Extra doc: https://github.com/casadi/casadi/wiki/L_16b

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L203

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L203-L205

";

%feature("docstring") casadi::Importer::serialize "

[INTERNAL] 
Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_16c

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L211

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L103-L105

";

%feature("docstring") casadi::Importer::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L99-L101

";

%feature("docstring") casadi::Importer::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L103-L109

";

%feature("docstring") casadi::Importer::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::Importer::is_null "

[INTERNAL] 
Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L73-L75

";

%feature("docstring") casadi::Importer::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L129-L131

";


// File: classcasadi_1_1IncrementalSerializer.xml
%feature("docstring") casadi::IncrementalSerializer "

[INTERNAL] ";


// File: classcasadi_1_1Input.xml


// File: classcasadi_1_1Integrator.xml
%feature("docstring") casadi::Integrator "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Interpolant.xml
%feature("docstring") casadi::Interpolant "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1IOInstruction.xml


// File: classcasadi_1_1Ipqp.xml
%feature("docstring") casadi::Ipqp "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1KeyboardInterruptException.xml
%feature("docstring") casadi::KeyboardInterruptException "

[INTERNAL] C++ includes: exception.hpp
";

%feature("docstring") 
casadi::KeyboardInterruptException::KeyboardInterruptException "

[INTERNAL] 
Default constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L101

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L101-L101

";

%feature("docstring") 
casadi::KeyboardInterruptException::~KeyboardInterruptException "

[INTERNAL]  throw ()
Destructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L103

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L103-L103

";

%feature("docstring") casadi::KeyboardInterruptException::what "

[INTERNAL]  throw ()
Display error.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L90

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L90-L92

";


// File: classcasadi_1_1LapackLu.xml
%feature("docstring") casadi::LapackLu "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1LapackQr.xml
%feature("docstring") casadi::LapackQr "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1LinearInterpolant.xml
%feature("docstring") casadi::LinearInterpolant "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Linsol.xml
%feature("docstring") casadi::Linsol "

[INTERNAL] 
Linear solver.

Create a solver for linear systems of equations Solves the linear 
system 
A*X = B or A^T*X = B for X with A square and non-singular

If A is structurally singular, an error will be thrown during init. If
 A is
 numerically singular, the prepare step will fail.
General informationList 
of plugins
- csparsecholesky

- csparse

- ma27

- lapacklu

- lapackqr

- mumps

- ldl

- qr

- tridiag

- symbolicqr

Note: some of the plugins in this list might not be available on your 

system.  Also, there might be extra plugins available to you that are 
not 
listed here. You can obtain their documentation with   
Linsol.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

csparsecholesky
---------------



Linsol with CSparseCholesky Interface

Extra doc: https://github.com/casadi/casadi/wiki/L_21u



--------------------------------------------------------------------------------

csparse
-------



Linsol with CSparse Interface

Extra doc: https://github.com/casadi/casadi/wiki/L_21t



--------------------------------------------------------------------------------

ma27
----



Interface to the sparse direct linear solver MA27 Works for symmetric
 
indefinite systems Partly adopted from qpOASES 3.2 
Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_229



--------------------------------------------------------------------------------

lapacklu
--------



This class solves the linear system  A.x=b by making an LU factorization of 
A:  A = L.U, with L lower and U upper triangular

Extra doc: https://github.com/casadi/casadi/wiki/L_22h

>List of available options

+-----------------------------+---------+----------------------------------+
|             Id              |  Type   |           Description            |
+=============================+=========+==================================+
| allow_equilibration_failure | OT_BOOL | Non-fatal error when             |
|                             |         | equilibration fails              |
+-----------------------------+---------+----------------------------------+
| equilibration               | OT_BOOL | Equilibrate the matrix           |
+-----------------------------+---------+----------------------------------+



--------------------------------------------------------------------------------

lapackqr
--------



This class solves the linear system  A.x=b by making an QR factorization of 
A:  A = Q.R, with Q orthogonal and R upper triangular

Extra doc: https://github.com/casadi/casadi/wiki/L_22g

>List of available options

+----------+--------+------------------------------------------------------+
|    Id    |  Type  |                     Description                      |
+==========+========+======================================================+
| max_nrhs | OT_INT | Maximum number of right-hand-sides that get          |
|          |        | processed in a single pass [default:10].             |
+----------+--------+------------------------------------------------------+



--------------------------------------------------------------------------------

mumps
-----



Interface to the sparse direct linear solver MUMPS Works for 
symmetric 
indefinite systems 
Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_22t

>List of available options

+-----------+---------+-------------------+
|    Id     |  Type   |    Description    |
+===========+=========+===================+
| posdef    | OT_BOOL | Positive definite |
+-----------+---------+-------------------+
| symmetric | OT_BOOL | Symmetric matrix  |
+-----------+---------+-------------------+



--------------------------------------------------------------------------------

ldl
---



Linear solver using sparse direct LDL factorization

Extra doc: https://github.com/casadi/casadi/wiki/L_233



--------------------------------------------------------------------------------

qr
--



Linear solver using sparse direct QR factorization

Extra doc: https://github.com/casadi/casadi/wiki/L_22z



--------------------------------------------------------------------------------

tridiag
-------



Linear solver for tridiagonal matrices

Extra doc: https://github.com/casadi/casadi/wiki/L_22v



--------------------------------------------------------------------------------

symbolicqr
----------



Linear solver for sparse least-squares problems Inspired from 
https://github.com/scipy/scipy/blob/v0.14.0/scipy/sparse/linalg/isolve/lsqr.py#L96

Extra doc: https://github.com/casadi/casadi/wiki/L_230

Linsol based on QR factorization with sparsity pattern based reordering  
without partial pivoting

Extra doc: https://github.com/casadi/casadi/wiki/L_231

>List of available options

+-------+---------+----------------------------------------------------+
|  Id   |  Type   |                    Description                     |
+=======+=========+====================================================+
| fopts | OT_DICT | Options to be passed to generated function objects |
+-------+---------+----------------------------------------------------+

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_1kh

C++ includes: linsol.hpp
";

%feature("docstring") casadi::Linsol::solve "

[INTERNAL] 
Low-level API

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L136

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L187-L194

>  int casadi::Linsol::solve(const double *A, double *x, casadi_int nrhs=1, bool tr=false, int mem=0) const
------------------------------------------------------------------------
[INTERNAL] 
Low-level API

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L136

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L187-L194

";

";

%feature("docstring") casadi::Linsol::sfact "

[INTERNAL] 
Symbolic factorization of the linear system, e.g. selecting 
pivots.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L103

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L104-L107

>  void casadi::Linsol::sfact(const DM &A) const
------------------------------------------------------------------------
[INTERNAL] 
Symbolic factorization of the linear system, e.g. selecting pivots.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L103

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L104-L107

";

";

%feature("docstring") casadi::Linsol::nfact "

[INTERNAL] 
Numeric factorization of the linear system.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L106

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L126-L129

>  void casadi::Linsol::nfact(const DM &A) const
------------------------------------------------------------------------
[INTERNAL] 
Numeric factorization of the linear system.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L106

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L126-L129

";

";

%feature("docstring") casadi::Linsol::neig "

[INTERNAL] 
Number of negative eigenvalues.

Not available for all solvers

Extra doc: https://github.com/casadi/casadi/wiki/L_1kk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L119

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L165-L170

>  casadi_int casadi::Linsol::neig(const DM &A) const
------------------------------------------------------------------------
[INTERNAL] 
Number of negative eigenvalues.

Not available for all solvers

Extra doc: https://github.com/casadi/casadi/wiki/L_1kk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L119

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L165-L170

";

";

%feature("docstring") casadi::Linsol::rank "

[INTERNAL] 
 Matrix rank.

Not available for all solvers

Extra doc: https://github.com/casadi/casadi/wiki/L_1kl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L126

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L176-L181

>  casadi_int casadi::Linsol::rank(const DM &A) const
------------------------------------------------------------------------
[INTERNAL] 
 Matrix rank.

Not available for all solvers

Extra doc: https://github.com/casadi/casadi/wiki/L_1kl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L126

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L176-L181

";

";

%feature("docstring") casadi::Linsol::Linsol "

[INTERNAL] 
Constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L66

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L34-L38

>  casadi::Linsol::Linsol(const std::string &name, const std::string &solver, const Sparsity &sp, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L66

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L34-L38

";

";

%feature("docstring") casadi::Linsol::plugin_name "

[INTERNAL] 
Query plugin name.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L97

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L64-L66

";

%feature("docstring") casadi::Linsol::sparsity "

[INTERNAL] 
Get linear system sparsity.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L68-L70

";

%feature("docstring") casadi::Linsol::stats "

[INTERNAL] 
Get all statistics obtained at the end of the last evaluate 
call.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L129

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L217-L219

";

%feature("docstring") casadi::Linsol::checkout "

[INTERNAL] 
Checkout a memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L142

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L196-L198

";

%feature("docstring") casadi::Linsol::release "

[INTERNAL] 
Release a memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L145

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L200-L202

";

%feature("docstring") casadi::Linsol::serialize "

[INTERNAL] 
Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_1km

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L221-L224

";

%feature("docstring") casadi::Linsol::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L99-L101

";

%feature("docstring") casadi::Linsol::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L103-L109

";

%feature("docstring") casadi::Linsol::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::Linsol::is_null "

[INTERNAL] 
Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L73-L75

";

%feature("docstring") casadi::Linsol::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L129-L131

";


// File: classcasadi_1_1LinsolCall.xml


// File: classcasadi_1_1LinsolLdl.xml
%feature("docstring") casadi::LinsolLdl "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1LinsolQr.xml
%feature("docstring") casadi::LinsolQr "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Logger.xml
%feature("docstring") casadi::Logger "

[INTERNAL] 
Keeps track of logging output to screen and/or files.

All printout from CasADi routines should go through this files.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_23u

C++ includes: casadi_logger.hpp
";


// File: classcasadi_1_1Matrix.xml


/*
 Construct symbolic primitives 
*/

/*
The \"sym\" function is intended to work in a similar way as \"sym\" 

used in the Symbolic Toolbox for Matlab but instead creating a CasADi 

symbolic primitive.

*/
%feature("docstring") casadi::MatrixCommon "

[INTERNAL] 
Sparse matrix class. SX and DM are specializations.

General sparse matrix class that is designed with the idea that 

\"everything is a matrix\", that is, also scalars and vectors.
This 
philosophy makes it easy to use and to interface in particularly
 with 
Python and Matlab/Octave.
 Index starts with 0.
Index vec happens as 
follows: (rr, cc) -> k = rr+cc*size1()
Vectors are column vectors.
 The 
storage format is Compressed Column Storage (CCS), similar to 
that used for
 sparse matrices in Matlab, 
but unlike this format, we do allow for 
elements to be structurally 
non-zero but numerically zero.
 Matrix<Scalar> 
is polymorphic with a std::vector<Scalar> that 
contain all non-identical-
zero elements.
The sparsity can be accessed with Sparsity&  sparsity()
Joel 
Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_18e

C++ includes: matrix_decl.hpp
";

%feature("docstring") casadi::MatrixCommon::get "

[INTERNAL] 
Get a submatrix, two arguments

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L239

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L111-L131

>  void casadi::Matrix< Scalar >::get(Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &rr, const Matrix< casadi_int > &cc) const
------------------------------------------------------------------------
[INTERNAL] 
Get a submatrix, two arguments

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L239

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L111-L131

";

";

%feature("docstring") casadi::MatrixCommon::set "

[INTERNAL] 
Set a submatrix, two arguments

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L255

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L218-L282

>  void casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &rr, const Matrix< casadi_int > &cc)
------------------------------------------------------------------------
[INTERNAL] 
Set a submatrix, two arguments

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L255

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L218-L282

";

";

%feature("docstring") casadi::MatrixCommon::get_nz "

[INTERNAL] 
Get a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L262

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L406-L432

>  void casadi::Matrix< Scalar >::get_nz(Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &k) const
------------------------------------------------------------------------
[INTERNAL] 
Get a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L262

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L406-L432

";

";

%feature("docstring") casadi::MatrixCommon::set_nz "

[INTERNAL] 
Set a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L268

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L447-L489

>  void casadi::Matrix< Scalar >::set_nz(const Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &k)
------------------------------------------------------------------------
[INTERNAL] 
Set a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L268

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L447-L489

";

";

%feature("docstring") casadi::MatrixCommon::nonzeros "

[INTERNAL] 
Access the non-zero elements

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L960

>  const std::vector<Scalar>& casadi::Matrix< Scalar >::nonzeros() const
------------------------------------------------------------------------
[INTERNAL] 
Access the non-zero elements

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L960

";

";

%feature("docstring") casadi::MatrixCommon::ptr "

[INTERNAL] 
Get a pointer to the data

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L966

>  const Scalar* casadi::Matrix< Scalar >::ptr() const
------------------------------------------------------------------------
[INTERNAL] 
Get a pointer to the data

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L966

";

";

%feature("docstring") casadi::MatrixCommon::get_ptr "

[INTERNAL] 
Get a pointer to the data

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L968

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L968-L968

>  const Scalar* casadi::Matrix::get_ptr(const Matrix< Scalar > &v)
------------------------------------------------------------------------
[INTERNAL] 
Get a pointer to the data

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L968

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L968-L968

";

";

%feature("docstring") casadi::MatrixCommon::get_row "

[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L194

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L194-L194

";

%feature("docstring") casadi::MatrixCommon::get_colind "

[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L195-L195

";

%feature("docstring") casadi::MatrixCommon::row "

[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200-L200

>  casadi_int casadi::GenericMatrix< Matrix< Scalar >  >::row(casadi_int el) const
------------------------------------------------------------------------
[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200-L200

";

";

%feature("docstring") casadi::MatrixCommon::colind "

[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201-L201

>  casadi_int casadi::GenericMatrix< Matrix< Scalar >  >::colind(casadi_int col) const
------------------------------------------------------------------------
[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201-L201

";

";

%feature("docstring") casadi::MatrixCommon::is_equal "

[INTERNAL] 
Functions called by friend functions defined for  
GenericExpression

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L294

";

%feature("docstring") casadi::MatrixCommon::mmin "

[INTERNAL] 
Functions called by friend functions defined for  
GenericExpression

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L295

";

%feature("docstring") casadi::MatrixCommon::mmax "

[INTERNAL] 
Functions called by friend functions defined for  
GenericExpression

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L296

";

%feature("docstring") casadi::MatrixCommon::simplify "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L301

";

%feature("docstring") casadi::MatrixCommon::jacobian "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L302

";

%feature("docstring") casadi::MatrixCommon::jacobian_sparsity "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L304

";

%feature("docstring") casadi::MatrixCommon::hessian "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L307

>  static Matrix<Scalar> casadi::Matrix< Scalar >::hessian(const Matrix< Scalar > &f, const Matrix< Scalar > &x, Matrix< Scalar > &g, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L307

";

";

%feature("docstring") casadi::MatrixCommon::substitute "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L314

>  static std::vector<Matrix<Scalar> > casadi::Matrix< Scalar >::substitute(const std::vector< Matrix< Scalar > > &ex, const std::vector< Matrix< Scalar > > &v, const std::vector< Matrix< Scalar > > &vdef)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L314

";

";

%feature("docstring") casadi::MatrixCommon::substitute_inplace "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L317

";

%feature("docstring") casadi::MatrixCommon::pinv "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L322

>  static Matrix<Scalar> casadi::Matrix< Scalar >::pinv(const Matrix< Scalar > &A, const std::string &lsolver, const Dict &opts)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L322

";

";

%feature("docstring") casadi::MatrixCommon::expm_const "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L324

";

%feature("docstring") casadi::MatrixCommon::expm "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L325

";

%feature("docstring") casadi::MatrixCommon::solve "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L327

>  static Matrix<Scalar> casadi::Matrix< Scalar >::solve(const Matrix< Scalar > &A, const Matrix< Scalar > &b, const std::string &lsolver, const Dict &opts)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L327

";

";

%feature("docstring") casadi::MatrixCommon::inv "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L330

>  static Matrix<Scalar> casadi::Matrix< Scalar >::inv(const Matrix< Scalar > &A, const std::string &lsolver, const Dict &opts)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L330

";

";

%feature("docstring") casadi::MatrixCommon::n_nodes "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L333

";

%feature("docstring") casadi::MatrixCommon::print_operator "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L334

";

%feature("docstring") casadi::MatrixCommon::extract "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L336

";

%feature("docstring") casadi::MatrixCommon::shared "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L338

";

%feature("docstring") casadi::MatrixCommon::_bilin "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L343

";

%feature("docstring") casadi::MatrixCommon::_rank1 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L346

";

%feature("docstring") casadi::MatrixCommon::if_else "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L350

";

%feature("docstring") casadi::MatrixCommon::conditional "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L354

";

%feature("docstring") casadi::MatrixCommon::depends_on "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L358

";

%feature("docstring") casadi::MatrixCommon::mrdivide "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L359

";

%feature("docstring") casadi::MatrixCommon::mldivide "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L360

";

%feature("docstring") casadi::MatrixCommon::symvar "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L361

";

%feature("docstring") casadi::MatrixCommon::det "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L362

";

%feature("docstring") casadi::MatrixCommon::inv_minor "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L363

";

%feature("docstring") casadi::MatrixCommon::trace "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L364

";

%feature("docstring") casadi::MatrixCommon::norm_1 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L365

";

%feature("docstring") casadi::MatrixCommon::norm_2 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L366

";

%feature("docstring") casadi::MatrixCommon::norm_fro "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L367

";

%feature("docstring") casadi::MatrixCommon::norm_inf "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L368

";

%feature("docstring") casadi::MatrixCommon::sum2 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L369

";

%feature("docstring") casadi::MatrixCommon::sum1 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L370

";

%feature("docstring") casadi::MatrixCommon::dot "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L371

";

%feature("docstring") casadi::MatrixCommon::nullspace "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L372

";

%feature("docstring") casadi::MatrixCommon::diag "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L373

";

%feature("docstring") casadi::MatrixCommon::unite "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L374

";

%feature("docstring") casadi::MatrixCommon::project "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L375

";

%feature("docstring") casadi::MatrixCommon::polyval "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L377

";

%feature("docstring") casadi::MatrixCommon::densify "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L379

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L492-L494

>  Matrix< Scalar > casadi::Matrix< Scalar >::densify(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L379

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L492-L494

";

";

%feature("docstring") casadi::MatrixCommon::einstein "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L387

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L566-L573

>  Matrix< Scalar > casadi::Matrix< Scalar >::einstein(const Matrix< Scalar > &A, const Matrix< Scalar > &B, const std::vector< casadi_int > &dim_a, const std::vector< casadi_int > &dim_b, const std::vector< casadi_int > &dim_c, const std::vector< casadi_int > &a, const std::vector< casadi_int > &b, const std::vector< casadi_int > &c)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L387

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L566-L573

";

";

%feature("docstring") casadi::MatrixCommon::cumsum "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L392

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L527-L538

";

%feature("docstring") casadi::MatrixCommon::_logsumexp "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L393

";

%feature("docstring") casadi::MatrixCommon::cse "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L394

";

%feature("docstring") casadi::MatrixCommon::blockcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L399

";

%feature("docstring") casadi::MatrixCommon::horzcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L400

";

%feature("docstring") casadi::MatrixCommon::horzsplit "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L402

";

%feature("docstring") casadi::MatrixCommon::vertcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L404

";

%feature("docstring") casadi::MatrixCommon::vertsplit "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L406

";

%feature("docstring") casadi::MatrixCommon::diagsplit "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L409

";

%feature("docstring") casadi::MatrixCommon::reshape "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L413

>  static Matrix<Scalar> casadi::Matrix< Scalar >::reshape(const Matrix< Scalar > &x, const Sparsity &sp)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L413

";

";

%feature("docstring") casadi::MatrixCommon::sparsity_cast "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L414

";

%feature("docstring") casadi::MatrixCommon::kron "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L415

";

%feature("docstring") casadi::MatrixCommon::mtimes "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L416

";

%feature("docstring") casadi::MatrixCommon::mac "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L417

";

%feature("docstring") casadi::MatrixCommon::sparsify "

[INTERNAL] 
Make a matrix sparse by removing numerical zeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_191

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L613

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L613-L615

>  Matrix<Scalar> casadi::Matrix::sparsify(const Matrix< Scalar > &A, double tol=0)
------------------------------------------------------------------------
[INTERNAL] 
Make a matrix sparse by removing numerical zeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_191

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L613

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L613-L615

";

";

%feature("docstring") casadi::MatrixCommon::expand "

[INTERNAL] 
Expand the expression as a weighted sum (with constant weights)

Extra doc: https://github.com/casadi/casadi/wiki/L_192

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L620

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L620-L623

>  void casadi::Matrix::expand(const Matrix< Scalar > &ex, Matrix< Scalar > &weights, Matrix< Scalar > &terms)
------------------------------------------------------------------------
[INTERNAL] 
Expand the expression as a weighted sum (with constant weights)

Extra doc: https://github.com/casadi/casadi/wiki/L_192

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L620

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L620-L623

";

";

%feature("docstring") casadi::MatrixCommon::pw_const "

[INTERNAL] 
Create a piecewise constant function.

Create a piecewise constant function with n=val.size() intervals

Inputs:

Parameters:
-----------

t: 
a scalar variable (e.g. time)

tval: 
vector with the discrete values of t at the interval transitions 

(length n-1)

val: 
vector with the value of the function for each interval (length n)

Extra doc: https://github.com/casadi/casadi/wiki/L_193

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L635

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L635-L639

>  Matrix<Scalar> casadi::Matrix::pw_const(const Matrix< Scalar > &t, const Matrix< Scalar > &tval, const Matrix< Scalar > &val)
------------------------------------------------------------------------
[INTERNAL] 
Create a piecewise constant function.

Create a piecewise constant function with n=val.size() intervals

Inputs:

Parameters:
-----------

t: 
a scalar variable (e.g. time)

tval: 
vector with the discrete values of t at the interval transitions 

(length n-1)

val: 
vector with the value of the function for each interval (length n)

Extra doc: https://github.com/casadi/casadi/wiki/L_193

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L635

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L635-L639

";

";

%feature("docstring") casadi::MatrixCommon::pw_lin "

[INTERNAL] 
t a scalar variable (e.g. time)

Create a piecewise linear function

Create a piecewise linear function:

Inputs:

tval vector with the the discrete values of t (monotonically 
increasing)

val vector with the corresponding function values (same length as 
tval)

::

                                                                                                      Extra doc: https://github.com/casadi/casadi/wiki/L_194 
  



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L652

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L652-L655

>  Matrix<Scalar> casadi::Matrix::pw_lin(const Matrix< Scalar > &t, const Matrix< Scalar > &tval, const Matrix< Scalar > &val)
------------------------------------------------------------------------
[INTERNAL] 
t a scalar variable (e.g. time)

Create a piecewise linear function

Create a piecewise linear function:

Inputs:

tval vector with the the discrete values of t (monotonically 
increasing)

val vector with the corresponding function values (same length as 
tval)

::

                                                                                                      Extra doc: https://github.com/casadi/casadi/wiki/L_194 
  



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L652

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L652-L655

";

";

%feature("docstring") casadi::MatrixCommon::heaviside "

[INTERNAL] 
Heaviside function.

\\\\[ \\\\begin {cases} H(x) = 0 & x<0 \\\\\\\\ H(x) = 1/2 & x=0 
\\\\\\\\ 
H(x) = 1 & x>0 \\\\\\\\ \\\\end {cases} \\\\]

Extra doc: https://github.com/casadi/casadi/wiki/L_195

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L668

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L668-L670

>  Matrix<Scalar> casadi::Matrix::heaviside(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
Heaviside function.

\\\\[ \\\\begin {cases} H(x) = 0 & x<0 \\\\\\\\ H(x) = 1/2 & x=0 
\\\\\\\\ 
H(x) = 1 & x>0 \\\\\\\\ \\\\end {cases} \\\\]

Extra doc: https://github.com/casadi/casadi/wiki/L_195

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L668

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L668-L670

";

";

%feature("docstring") casadi::MatrixCommon::rectangle "

[INTERNAL] 
rectangle function

\\\\[ \\\\begin {cases} \\\\Pi(x) = 1 & |x| < 1/2 \\\\\\\\ \\\\Pi(x) =
 1/2 
& |x| = 1/2 \\\\\\\\ \\\\Pi(x) = 0 & |x| > 1/2 \\\\\\\\ \\\\end 
{cases} 
\\\\]

Also called: gate function, block function, band function, pulse 
function, 
window function

Extra doc: https://github.com/casadi/casadi/wiki/L_23n

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L685

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L685-L687

>  Matrix<Scalar> casadi::Matrix::rectangle(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
rectangle function

\\\\[ \\\\begin {cases} \\\\Pi(x) = 1 & |x| < 1/2 \\\\\\\\ \\\\Pi(x) =
 1/2 
& |x| = 1/2 \\\\\\\\ \\\\Pi(x) = 0 & |x| > 1/2 \\\\\\\\ \\\\end 
{cases} 
\\\\]

Also called: gate function, block function, band function, pulse 
function, 
window function

Extra doc: https://github.com/casadi/casadi/wiki/L_23n

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L685

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L685-L687

";

";

%feature("docstring") casadi::MatrixCommon::triangle "

[INTERNAL] 
triangle function

\\\\[ \\\\begin {cases} \\\\Lambda(x) = 0 & |x| >= 1 \\\\\\\\ 
\\\\Lambda(x)
 = 1-|x| & |x| < 1 \\\\end {cases} \\\\]

Extra doc: https://github.com/casadi/casadi/wiki/L_23o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L700

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L700-L702

>  Matrix<Scalar> casadi::Matrix::triangle(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
triangle function

\\\\[ \\\\begin {cases} \\\\Lambda(x) = 0 & |x| >= 1 \\\\\\\\ 
\\\\Lambda(x)
 = 1-|x| & |x| < 1 \\\\end {cases} \\\\]

Extra doc: https://github.com/casadi/casadi/wiki/L_23o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L700

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L700-L702

";

";

%feature("docstring") casadi::MatrixCommon::ramp "

[INTERNAL] 
ramp function

\\\\[ \\\\begin {cases} R(x) = 0 & x <= 1 \\\\\\\\ R(x) = x & x > 1 

\\\\\\\\ \\\\end {cases} \\\\]

Also called: slope function

Extra doc: https://github.com/casadi/casadi/wiki/L_23p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L717

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L717-L719

>  Matrix<Scalar> casadi::Matrix::ramp(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
ramp function

\\\\[ \\\\begin {cases} R(x) = 0 & x <= 1 \\\\\\\\ R(x) = x & x > 1 

\\\\\\\\ \\\\end {cases} \\\\]

Also called: slope function

Extra doc: https://github.com/casadi/casadi/wiki/L_23p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L717

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L717-L719

";

";

%feature("docstring") casadi::MatrixCommon::gauss_quadrature "

[INTERNAL] 
Integrate f from a to b using Gaussian quadrature with n points.

Extra doc: https://github.com/casadi/casadi/wiki/L_196

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L732

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L732-L736

>  Matrix<Scalar> casadi::Matrix::gauss_quadrature(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Matrix< Scalar > &a, const Matrix< Scalar > &b, casadi_int order, const Matrix< Scalar > &w)
------------------------------------------------------------------------
[INTERNAL] 
Integrate f from a to b using Gaussian quadrature with n points.

Extra doc: https://github.com/casadi/casadi/wiki/L_196

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L732

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L732-L736

";

";

%feature("docstring") casadi::MatrixCommon::forward "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L444

";

%feature("docstring") casadi::MatrixCommon::reverse "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L449

";

%feature("docstring") casadi::MatrixCommon::which_depends "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L453

";

%feature("docstring") casadi::MatrixCommon::taylor "

[INTERNAL] 
univariate Taylor series expansion

Calculate the Taylor expansion of expression 'ex' up to order 'order' 
with 
respect to variable 'x' around the point 'a'

$(x)=f(a)+f'(a)(x-a)+f''(a)\\\\frac 

{(x-a)^2}{2!}+f'''(a)\\\\frac{(x-a)^3}{3!}+\\\\ldots$

Example usage:

::

>>> taylor(sin(x), x)

::

>>   x



Extra doc: https://github.com/casadi/casadi/wiki/L_23q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L758

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L758-L760

>  Matrix<Scalar> casadi::Matrix::taylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
univariate Taylor series expansion

Calculate the Taylor expansion of expression 'ex' up to order 'order' 
with 
respect to variable 'x' around the point 'a'

$(x)=f(a)+f'(a)(x-a)+f''(a)\\\\frac 

{(x-a)^2}{2!}+f'''(a)\\\\frac{(x-a)^3}{3!}+\\\\ldots$

Example usage:

::

>>> taylor(sin(x), x)

::

>>   x



Extra doc: https://github.com/casadi/casadi/wiki/L_23q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L758

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L758-L760

";

";

%feature("docstring") casadi::MatrixCommon::mtaylor "

[INTERNAL] 
multivariate Taylor series expansion

Do Taylor expansions until the aggregated order of a term is equal to 

'order'. The aggregated order of  $x^n y^m$ equals  $n+m$.

The argument order_contributions can denote how match each variable 

contributes to the aggregated order. If x=[x, y] and 

order_contributions=[1, 2], then the aggregated order of  $x^n y^m$ equals  
$1n+2m$.

Example usage



::

>>> taylor(sin(x+y),[x, y],[a, b], 1)
 $ \\\\sin(b+a)+\\\\cos(b+a)(x-a)+\\\\cos(b+a)(y-b) $

::

>>> taylor(sin(x+y),[x, y],[0, 0], 4)
 $ y+x-(x^3+3y x^2+3 y^2 x+y^3)/6 $

::

>>> taylor(sin(x+y),[x, y],[0, 0], 4,[1, 2])
 $ (-3 x^2 y-x^3)/6+y+x $

Extra doc: https://github.com/casadi/casadi/wiki/L_23s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L801

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L801-L805

>  Matrix<Scalar> casadi::Matrix::mtaylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x, const Matrix< Scalar > &a, casadi_int order, const std::vector< casadi_int > &order_contributions)
------------------------------------------------------------------------
[INTERNAL] 
multivariate Taylor series expansion

Do Taylor expansions until the aggregated order of a term is equal to 

'order'. The aggregated order of  $x^n y^m$ equals  $n+m$.

The argument order_contributions can denote how match each variable 

contributes to the aggregated order. If x=[x, y] and 

order_contributions=[1, 2], then the aggregated order of  $x^n y^m$ equals  
$1n+2m$.

Example usage



::

>>> taylor(sin(x+y),[x, y],[a, b], 1)
 $ \\\\sin(b+a)+\\\\cos(b+a)(x-a)+\\\\cos(b+a)(y-b) $

::

>>> taylor(sin(x+y),[x, y],[0, 0], 4)
 $ y+x-(x^3+3y x^2+3 y^2 x+y^3)/6 $

::

>>> taylor(sin(x+y),[x, y],[0, 0], 4,[1, 2])
 $ (-3 x^2 y-x^3)/6+y+x $

Extra doc: https://github.com/casadi/casadi/wiki/L_23s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L801

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L801-L805

";

";

%feature("docstring") casadi::MatrixCommon::poly_coeff "

[INTERNAL] 
extracts polynomial coefficients from an expression

Parameters:
-----------

ex: 
Scalar expression that represents a polynomial

x: 
Scalar symbol that the polynomial is build up with

Extra doc: https://github.com/casadi/casadi/wiki/L_197

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L813

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L813-L816

>  Matrix<Scalar> casadi::Matrix::poly_coeff(const Matrix< Scalar > &f, const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
extracts polynomial coefficients from an expression

Parameters:
-----------

ex: 
Scalar expression that represents a polynomial

x: 
Scalar symbol that the polynomial is build up with

Extra doc: https://github.com/casadi/casadi/wiki/L_197

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L813

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L813-L816

";

";

%feature("docstring") casadi::MatrixCommon::poly_roots "

[INTERNAL] 
Attempts to find the roots of a polynomial.

This will only work for polynomials up to order 3 It is assumed that 
the 
roots are real.

Extra doc: https://github.com/casadi/casadi/wiki/L_198

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L824

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L824-L826

>  Matrix<Scalar> casadi::Matrix::poly_roots(const Matrix< Scalar > &p)
------------------------------------------------------------------------
[INTERNAL] 
Attempts to find the roots of a polynomial.

This will only work for polynomials up to order 3 It is assumed that 
the 
roots are real.

Extra doc: https://github.com/casadi/casadi/wiki/L_198

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L824

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L824-L826

";

";

%feature("docstring") casadi::MatrixCommon::eig_symbolic "

[INTERNAL] 
Attempts to find the eigenvalues of a symbolic matrix.

This will only work for up to 3x3 matrices

Extra doc: https://github.com/casadi/casadi/wiki/L_199

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L833

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L833-L835

>  Matrix<Scalar> casadi::Matrix::eig_symbolic(const Matrix< Scalar > &m)
------------------------------------------------------------------------
[INTERNAL] 
Attempts to find the eigenvalues of a symbolic matrix.

This will only work for up to 3x3 matrices

Extra doc: https://github.com/casadi/casadi/wiki/L_199

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L833

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L833-L835

";

";

%feature("docstring") casadi::MatrixCommon::evalf "

[INTERNAL] 
Evaluates the expression numerically.

An error is raised when the expression contains symbols

Extra doc: https://github.com/casadi/casadi/wiki/L_19a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L843

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L843-L845

>  Matrix<double> casadi::Matrix::evalf(const Matrix< Scalar > &expr)
------------------------------------------------------------------------
[INTERNAL] 
Evaluates the expression numerically.

An error is raised when the expression contains symbols

Extra doc: https://github.com/casadi/casadi/wiki/L_19a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L843

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L843-L845

";

";

%feature("docstring") casadi::MatrixCommon::qr_sparse "

[INTERNAL] 
Sparse direct QR factorization.

See T. Davis: Direct Methods for Sparse Linear Systems

Extra doc: https://github.com/casadi/casadi/wiki/L_18t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L539

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L539-L543

>  void casadi::Matrix::qr_sparse(const Matrix< Scalar > &A, Matrix< Scalar > &V, Matrix< Scalar > &R, Matrix< Scalar > &beta, std::vector< casadi_int > &prinv, std::vector< casadi_int > &pc, bool amd=true)
------------------------------------------------------------------------
[INTERNAL] 
Sparse direct QR factorization.

See T. Davis: Direct Methods for Sparse Linear Systems

Extra doc: https://github.com/casadi/casadi/wiki/L_18t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L539

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L539-L543

";

";

%feature("docstring") casadi::MatrixCommon::qr_solve "

[INTERNAL] 
 Solve using a sparse QR factorization.

Extra doc: https://github.com/casadi/casadi/wiki/L_18u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L549

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L549-L554

>  Matrix<Scalar> casadi::Matrix::qr_solve(const Matrix< Scalar > &b, const Matrix< Scalar > &v, const Matrix< Scalar > &r, const Matrix< Scalar > &beta, const std::vector< casadi_int > &prinv, const std::vector< casadi_int > &pc, bool tr=false)
------------------------------------------------------------------------
[INTERNAL] 
 Solve using a sparse QR factorization.

Extra doc: https://github.com/casadi/casadi/wiki/L_18u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L549

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L549-L554

";

";

%feature("docstring") casadi::MatrixCommon::qr "

[INTERNAL] 
QR factorization using the modified Gram-Schmidt algorithm.

More stable than the classical Gram-Schmidt, but may break down if the
 rows
 of A are nearly linearly dependent See J. Demmel: Applied 
Numerical Linear
 Algebra (algorithm 3.1.). Note that in SWIG, Q and R 
are returned by 
value.

Extra doc: https://github.com/casadi/casadi/wiki/L_18s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L530

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L530-L532

>  void casadi::Matrix::qr(const Matrix< Scalar > &A, Matrix< Scalar > &Q, Matrix< Scalar > &R)
------------------------------------------------------------------------
[INTERNAL] 
QR factorization using the modified Gram-Schmidt algorithm.

More stable than the classical Gram-Schmidt, but may break down if the
 rows
 of A are nearly linearly dependent See J. Demmel: Applied 
Numerical Linear
 Algebra (algorithm 3.1.). Note that in SWIG, Q and R 
are returned by 
value.

Extra doc: https://github.com/casadi/casadi/wiki/L_18s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L530

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L530-L532

";

";

%feature("docstring") casadi::MatrixCommon::ldl "

[INTERNAL] 
Sparse LDL^T factorization.

Returns D and the strictly upper triangular entries of L^T I.e. ones 
on the
 diagonal are ignored. Only guarenteed to work for positive 
definite 
matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_18w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L573

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L573-L576

>  void casadi::Matrix::ldl(const Matrix< Scalar > &A, Matrix< Scalar > &D, Matrix< Scalar > &LT, std::vector< casadi_int > &p, bool amd=true)
------------------------------------------------------------------------
[INTERNAL] 
Sparse LDL^T factorization.

Returns D and the strictly upper triangular entries of L^T I.e. ones 
on the
 diagonal are ignored. Only guarenteed to work for positive 
definite 
matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_18w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L573

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L573-L576

";

";

%feature("docstring") casadi::MatrixCommon::ldl_solve "

[INTERNAL] 
 Solve using a sparse LDL^T factorization.

Extra doc: https://github.com/casadi/casadi/wiki/L_18x

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L582

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L582-L585

>  Matrix<Scalar> casadi::Matrix::ldl_solve(const Matrix< Scalar > &b, const Matrix< Scalar > &D, const Matrix< Scalar > &LT, const std::vector< casadi_int > &p)
------------------------------------------------------------------------
[INTERNAL] 
 Solve using a sparse LDL^T factorization.

Extra doc: https://github.com/casadi/casadi/wiki/L_18x

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L582

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L582-L585

";

";

%feature("docstring") casadi::MatrixCommon::all "

[INTERNAL] 
Returns true only if every element in the matrix is true.

Extra doc: https://github.com/casadi/casadi/wiki/L_18z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L597

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L597-L599

>  Matrix<Scalar> casadi::Matrix::all(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
Returns true only if every element in the matrix is true.

Extra doc: https://github.com/casadi/casadi/wiki/L_18z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L597

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L597-L599

";

";

%feature("docstring") casadi::MatrixCommon::any "

[INTERNAL] 
Returns true only if any element in the matrix is true.

Extra doc: https://github.com/casadi/casadi/wiki/L_18y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L590

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L590-L592

>  Matrix<Scalar> casadi::Matrix::any(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
Returns true only if any element in the matrix is true.

Extra doc: https://github.com/casadi/casadi/wiki/L_18y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L590

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L590-L592

";

";

%feature("docstring") casadi::MatrixCommon::adj "

[INTERNAL] 
 Matrix adjoint.

Extra doc: https://github.com/casadi/casadi/wiki/L_18p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L504

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L504-L506

>  Matrix<Scalar> casadi::Matrix::adj(const Matrix< Scalar > &A)
------------------------------------------------------------------------
[INTERNAL] 
 Matrix adjoint.

Extra doc: https://github.com/casadi/casadi/wiki/L_18p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L504

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L504-L506

";

";

%feature("docstring") casadi::MatrixCommon::minor "

[INTERNAL] 
Get the (i,j) minor matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_18q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L511

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L511-L513

>  Matrix<Scalar> casadi::Matrix::minor(const Matrix< Scalar > &x, casadi_int i, casadi_int j)
------------------------------------------------------------------------
[INTERNAL] 
Get the (i,j) minor matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_18q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L511

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L511-L513

";

";

%feature("docstring") casadi::MatrixCommon::cofactor "

[INTERNAL] 
Get the (i,j) cofactor matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_18r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L518

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L518-L520

>  Matrix<Scalar> casadi::Matrix::cofactor(const Matrix< Scalar > &x, casadi_int i, casadi_int j)
------------------------------------------------------------------------
[INTERNAL] 
Get the (i,j) cofactor matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_18r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L518

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L518-L520

";

";

%feature("docstring") casadi::MatrixCommon::chol "

[INTERNAL] 
Obtain a Cholesky factorisation of a matrix.

Performs and LDL transformation [L,D] = ldl(A) and returns 
diag(sqrt(D))*L'

Extra doc: https://github.com/casadi/casadi/wiki/L_18v

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L562

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L562-L564

>  Matrix<Scalar> casadi::Matrix::chol(const Matrix< Scalar > &A)
------------------------------------------------------------------------
[INTERNAL] 
Obtain a Cholesky factorisation of a matrix.

Performs and LDL transformation [L,D] = ldl(A) and returns 
diag(sqrt(D))*L'

Extra doc: https://github.com/casadi/casadi/wiki/L_18v

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L562

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L562-L564

";

";

%feature("docstring") casadi::MatrixCommon::norm_inf_mul "

[INTERNAL] 
Inf-norm of a Matrix-Matrix product.

Extra doc: https://github.com/casadi/casadi/wiki/L_190

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L605

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L605-L607

>  Matrix<Scalar> casadi::Matrix::norm_inf_mul(const Matrix< Scalar > &x, const Matrix< Scalar > &y)
------------------------------------------------------------------------
[INTERNAL] 
Inf-norm of a Matrix-Matrix product.

Extra doc: https://github.com/casadi/casadi/wiki/L_190

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L605

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L605-L607

";

";

%feature("docstring") casadi::MatrixCommon::diagcat "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L486

";

%feature("docstring") casadi::MatrixCommon::triplet "

[INTERNAL] 
Construct a sparse matrix from triplet form.

Default matrix size is max(col) x max(row)

Extra doc: https://github.com/casadi/casadi/wiki/L_23t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L993

>  static Matrix<Scalar> casadi::Matrix< Scalar >::triplet(const std::vector< casadi_int > &row, const std::vector< casadi_int > &col, const Matrix< Scalar > &d, const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
Construct a sparse matrix from triplet form.

Default matrix size is max(col) x max(row)

Extra doc: https://github.com/casadi/casadi/wiki/L_23t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L993

";

";

%feature("docstring") casadi::MatrixCommon::inf "

[INTERNAL] 
create a matrix with all inf

Extra doc: https://github.com/casadi/casadi/wiki/L_19k

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1005

>  static Matrix<Scalar> casadi::Matrix< Scalar >::inf(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
create a matrix with all inf

Extra doc: https://github.com/casadi/casadi/wiki/L_19k

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1005

";

";

%feature("docstring") casadi::MatrixCommon::nan "

[INTERNAL] 
create a matrix with all nan

Extra doc: https://github.com/casadi/casadi/wiki/L_19l

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1014

>  static Matrix<Scalar> casadi::Matrix< Scalar >::nan(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
create a matrix with all nan

Extra doc: https://github.com/casadi/casadi/wiki/L_19l

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1014

";

";

%feature("docstring") casadi::MatrixCommon::rand "

[INTERNAL] 
Create a matrix with uniformly distributed random numbers.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ab

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1201

>  static Matrix<Scalar> casadi::Matrix< Scalar >::rand(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
Create a matrix with uniformly distributed random numbers.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ab

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1201

";

";

%feature("docstring") casadi::MatrixCommon::interp1d "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L213

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1252-L1308

";

%feature("docstring") casadi::MatrixCommon::sprank "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L215

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L215-L215

";

%feature("docstring") casadi::MatrixCommon::norm_0_mul "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L216

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L216-L218

";

%feature("docstring") casadi::MatrixCommon::tril "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L219

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L219-L221

";

%feature("docstring") casadi::MatrixCommon::triu "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L222-L224

";

%feature("docstring") casadi::MatrixCommon::sumsqr "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L225

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L225-L225

";

%feature("docstring") casadi::MatrixCommon::linspace "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L226

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1202-L1212

";

%feature("docstring") casadi::MatrixCommon::cross "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L227

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1215-L1246

";

%feature("docstring") casadi::MatrixCommon::skew "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L228

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1311-L1319

";

%feature("docstring") casadi::MatrixCommon::inv_skew "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L229

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1322-L1327

";

%feature("docstring") casadi::MatrixCommon::tril2symm "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L230

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1331-L1337

";

%feature("docstring") casadi::MatrixCommon::triu2symm "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L231

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1355-L1361

";

%feature("docstring") casadi::MatrixCommon::repsum "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L232

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1340-L1352

";

%feature("docstring") casadi::MatrixCommon::diff "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1567-L1592

";

%feature("docstring") casadi::MatrixCommon::is_linear "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L235

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1536-L1538

";

%feature("docstring") casadi::MatrixCommon::is_quadratic "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L236

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1541-L1543

";

%feature("docstring") casadi::MatrixCommon::quadratic_coeff "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L237

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1546-L1554

";

%feature("docstring") casadi::MatrixCommon::linear_coeff "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L239

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1557-L1564

";

%feature("docstring") casadi::MatrixCommon::mpower "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1010

";

%feature("docstring") casadi::MatrixCommon::soc "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1011

";

%feature("docstring") casadi::MatrixCommon::linearize "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1008

";

%feature("docstring") casadi::MatrixCommon::gradient "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1006

";

%feature("docstring") casadi::MatrixCommon::tangent "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1007

";

%feature("docstring") casadi::MatrixCommon::jtimes "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1004

";

%feature("docstring") casadi::MatrixCommon::bilin "

[INTERNAL] 
Calculate bilinear/quadratic form x^T A y.

Parameters:
-----------

y: 
can be omitted, in which case x^T A x is calculated

Extra doc: https://github.com/casadi/casadi/wiki/L_1bo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L409

";

%feature("docstring") casadi::MatrixCommon::rank1 "

[INTERNAL] 
Make a rank-1 update to a matrix A.

Calculates A + 1/2 * alpha * x*y'

Extra doc: https://github.com/casadi/casadi/wiki/L_1bp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L422

";

%feature("docstring") casadi::MatrixCommon::sym "

[INTERNAL] 
Create a vector of length r of vectors of length p.

with nrow-by-ncol symbolic primitives

Extra doc: https://github.com/casadi/casadi/wiki/L_1dg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1074

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1074-L1076

>  static std::vector<std::vector<Matrix< Scalar > > > casadi::GenericMatrix< Matrix< Scalar >  >::sym(const std::string &name, casadi_int nrow, casadi_int ncol, casadi_int p, casadi_int r)
------------------------------------------------------------------------
[INTERNAL] 
Create a vector of length r of vectors of length p.

with nrow-by-ncol symbolic primitives

Extra doc: https://github.com/casadi/casadi/wiki/L_1dg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1074

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1074-L1076

";

";

%feature("docstring") casadi::MatrixCommon::zeros "

[INTERNAL] 
Create a dense matrix or a matrix with specified sparsity with 
all 
entries zero.

Extra doc: https://github.com/casadi/casadi/wiki/L_1dh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1087

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1087-L1089

>  static Matrix< Scalar >  casadi::GenericMatrix< Matrix< Scalar >  >::zeros(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
Create a dense matrix or a matrix with specified sparsity with all 
entries zero.

Extra doc: https://github.com/casadi/casadi/wiki/L_1dh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1087

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1087-L1089

";

";

%feature("docstring") casadi::MatrixCommon::ones "

[INTERNAL] 
Create a dense matrix or a matrix with specified sparsity with 
all 
entries one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1di

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100-L1102

>  static Matrix< Scalar >  casadi::GenericMatrix< Matrix< Scalar >  >::ones(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
Create a dense matrix or a matrix with specified sparsity with all 
entries one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1di

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100-L1102

";

";

%feature("docstring") casadi::MatrixCommon::printme "

[INTERNAL]

>  Matrix<Scalar> casadi::Matrix< Scalar >::printme(const Matrix< Scalar > &y) const
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::MatrixCommon::MatrixCommon "

[INTERNAL] 
Sparse matrix with a given sparsity and non-zero elements.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1283

>  casadi::Matrix< Scalar >::Matrix(const Sparsity &sp, const std::vector< Scalar > &d, bool dummy)
------------------------------------------------------------------------
[INTERNAL] 
Sparse matrix with a given sparsity and non-zero elements.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1283

";

";

%feature("docstring") casadi::MatrixCommon::scalar "

[INTERNAL] 
Convert to scalar type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L180

";

%feature("docstring") casadi::MatrixCommon::has_nz "

[INTERNAL] 
Returns true if the matrix has a non-zero at location rr, cc.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L219

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L65-L67

";

%feature("docstring") casadi::MatrixCommon::__nonzero__ "

[INTERNAL] 
Returns the truth value of a  Matrix.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L70-L76

";

%feature("docstring") casadi::MatrixCommon::T "

[INTERNAL] 
Transpose the matrix.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L494

";

%feature("docstring") casadi::MatrixCommon::print_split "

[INTERNAL] 
Get strings corresponding to the nonzeros and the 
interdependencies.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L873

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L696-L700

";

%feature("docstring") casadi::MatrixCommon::disp "

[INTERNAL] 
Print a representation of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L877

";

%feature("docstring") casadi::MatrixCommon::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L880

";

%feature("docstring") casadi::MatrixCommon::print_scalar "

[INTERNAL] 
Print scalar.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L883

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L609-L633

";

%feature("docstring") casadi::MatrixCommon::print_vector "

[INTERNAL] 
Print vector-style.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L886

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L636-L638

";

%feature("docstring") casadi::MatrixCommon::print_dense "

[INTERNAL] 
Print dense matrix-stype.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L889

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L686-L688

";

%feature("docstring") casadi::MatrixCommon::print_sparse "

[INTERNAL] 
Print sparse matrix style.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L892

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L691-L693

";

%feature("docstring") casadi::MatrixCommon::clear "

[INTERNAL] ";

%feature("docstring") casadi::MatrixCommon::resize "

[INTERNAL] ";

%feature("docstring") casadi::MatrixCommon::reserve "

[INTERNAL] ";

%feature("docstring") casadi::MatrixCommon::erase "

[INTERNAL] 
Erase a submatrix (leaving structural zeros in its place)

Erase elements of a matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_19g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L937

>  void casadi::Matrix< Scalar >::erase(const std::vector< casadi_int > &rr, bool ind1=false)
------------------------------------------------------------------------
[INTERNAL] 
Erase a submatrix (leaving structural zeros in its place)

Erase elements of a matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_19g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L937

";

";

%feature("docstring") casadi::MatrixCommon::remove "

[INTERNAL] 
Remove columns and rows.

Remove/delete rows and/or columns of a matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_19h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L944

";

%feature("docstring") casadi::MatrixCommon::enlarge "

[INTERNAL] 
Enlarge matrix.

Make the matrix larger by inserting empty rows and columns, keeping 
the 
existing non-zeros

Extra doc: https://github.com/casadi/casadi/wiki/L_19i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L952

";

%feature("docstring") casadi::MatrixCommon::sparsity "

[INTERNAL] 
Const access the sparsity - reference to data member.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L972

";

%feature("docstring") casadi::MatrixCommon::get_sparsity "

[INTERNAL] 
Get an owning reference to the sparsity pattern.

Extra doc: https://github.com/casadi/casadi/wiki/L_19j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L979

";

%feature("docstring") casadi::MatrixCommon::element_hash "

[INTERNAL] 
Returns a number that is unique for a given symbolic scalar.

Only defined if symbolic scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_19n

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1027

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L125-L127

";

%feature("docstring") casadi::MatrixCommon::is_regular "

[INTERNAL] 
Checks if expression does not contain NaN or Inf.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1030

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L99-L112

";

%feature("docstring") casadi::MatrixCommon::is_smooth "

[INTERNAL] 
Check if smooth.

Extra doc: https://github.com/casadi/casadi/wiki/L_19o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1035

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L115-L122

";

%feature("docstring") casadi::MatrixCommon::is_leaf "

[INTERNAL] 
Check if SX is a leaf of the SX graph.

Only defined if symbolic scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_19p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1042

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L130-L132

";

%feature("docstring") casadi::MatrixCommon::is_commutative "

[INTERNAL] 
Check whether a binary SX is commutative.

Only defined if symbolic scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_19q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1049

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L135-L137

";

%feature("docstring") casadi::MatrixCommon::is_symbolic "

[INTERNAL] 
Check if symbolic (Dense)

Sparse matrices invariable return false

Extra doc: https://github.com/casadi/casadi/wiki/L_19r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1056

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L149-L155

";

%feature("docstring") casadi::MatrixCommon::is_valid_input "

[INTERNAL] 
Check if matrix can be used to define function inputs.

Sparse matrices can return true if all non-zero elements are symbolic

Extra doc: https://github.com/casadi/casadi/wiki/L_19s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1063

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L140-L146

";

%feature("docstring") casadi::MatrixCommon::is_constant "

[INTERNAL] 
Check if the matrix is constant (note that false negative 
answers are 
possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_19v

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1085

";

%feature("docstring") casadi::MatrixCommon::is_integer "

[INTERNAL] 
Check if the matrix is integer-valued.

(note that false negative answers are possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_19w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1092

";

%feature("docstring") casadi::MatrixCommon::is_zero "

[INTERNAL] 
check if the matrix is 0 (note that false negative answers are 

possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_19x

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1097

";

%feature("docstring") casadi::MatrixCommon::is_one "

[INTERNAL] 
check if the matrix is 1 (note that false negative answers are 

possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_19y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1102

";

%feature("docstring") casadi::MatrixCommon::is_minus_one "

[INTERNAL] 
check if the matrix is -1 (note that false negative answers are
 
possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_19z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1107

";

%feature("docstring") casadi::MatrixCommon::is_eye "

[INTERNAL] 
check if the matrix is an identity matrix (note that false 
negative 
answers

are possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_1a0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1114

";

%feature("docstring") casadi::MatrixCommon::op "

[INTERNAL] 
Get operation type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1117

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L158-L160

";

%feature("docstring") casadi::MatrixCommon::is_op "

[INTERNAL] 
Is it a certain operation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1120

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L163-L165

";

%feature("docstring") casadi::MatrixCommon::has_zeros "

[INTERNAL] 
Check if the matrix has any zero entries which are not 
structural 
zeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1a1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1125

";

%feature("docstring") casadi::MatrixCommon::get_nonzeros "

[INTERNAL] 
Get all nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1a4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1142

>  std::vector<A> casadi::Matrix< Scalar >::get_nonzeros() const
------------------------------------------------------------------------
[INTERNAL] 
Get all nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1a4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1142

";

";

%feature("docstring") casadi::MatrixCommon::get_elements "

[INTERNAL] 
Get all elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1a3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1135

";

%feature("docstring") casadi::MatrixCommon::name "

[INTERNAL] 
Get name (only if symbolic scalar)

Extra doc: https://github.com/casadi/casadi/wiki/L_1a8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1166

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L187-L189

";

%feature("docstring") casadi::MatrixCommon::dep "

[INTERNAL] 
Get expressions of the children of the expression.

Only defined if symbolic scalar. Wraps  SXElem SXElem::dep(casadi_int ch=0) 
const.

Extra doc: https://github.com/casadi/casadi/wiki/L_1a9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1174

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L192-L194

";

%feature("docstring") casadi::MatrixCommon::n_dep "

[INTERNAL] 
Get the number of dependencies of a binary  SXElem.

Only defined if symbolic scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_1aa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1181

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L197-L199

";

%feature("docstring") casadi::MatrixCommon::export_code "

[INTERNAL] 
Export matrix in specific language.

lang: only 'matlab' supported for now

::

  * options:
  *   inline: Indicates if you want everything on a single line (default: False)
  *   name: Name of exported variable (default: 'm')
  *   indent_level: Level of indentation (default: 0)
  *   spoof_zero: Replace numerical zero by a 1e-200 (default: false)
  *               might be needed for matlab sparse construct,
  *               which doesn't allow numerical zero
  * 



Extra doc: https://github.com/casadi/casadi/wiki/L_1ac

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1219

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dm_instantiator.cpp#L94-L192

";

%feature("docstring") casadi::MatrixCommon::info "

[INTERNAL] 
Obtain information about sparsity

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1223

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dm_instantiator.cpp#L195-L197

";

%feature("docstring") casadi::MatrixCommon::serialize "

[INTERNAL] 
Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ah

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1249

>  void casadi::Matrix< Scalar >::serialize(SerializingStream &s) const
------------------------------------------------------------------------
[INTERNAL] 
Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ah

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1249

";

";

%feature("docstring") casadi::MatrixCommon::to_file "

[INTERNAL] 
Export numerical matrix to file

Supported formats:



::

  *   - .mtx   Matrix Market (sparse)
  *   - .txt   Ascii full precision representation (sparse)
  *            Whitespace separated, aligned.
  *            Comments with # % or /
  *            Uses C locale
  *            Structural zeros represented by 00
  *            Does not scale well for large sparse matrices
  * 



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1269

";

%feature("docstring") casadi::MatrixCommon::nnz "

[INTERNAL] 
Get the number of (structural) non-zero elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1an

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L84

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1119-L1121

";

%feature("docstring") casadi::MatrixCommon::nnz_lower "

[INTERNAL] 
Get the number of non-zeros in the lower triangular half.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ao

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L89

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1124-L1126

";

%feature("docstring") casadi::MatrixCommon::nnz_upper "

[INTERNAL] 
Get the number of non-zeros in the upper triangular half.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ap

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L94

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1129-L1131

";

%feature("docstring") casadi::MatrixCommon::nnz_diag "

[INTERNAL] 
Get get the number of non-zeros on the diagonal.

Extra doc: https://github.com/casadi/casadi/wiki/L_1aq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L99

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1134-L1136

";

%feature("docstring") casadi::MatrixCommon::numel "

[INTERNAL] 
Get the number of elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ar

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L104

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1139-L1141

";

%feature("docstring") casadi::MatrixCommon::size1 "

[INTERNAL] 
Get the first dimension (i.e. number of rows)

Extra doc: https://github.com/casadi/casadi/wiki/L_1as

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L109

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1144-L1146

";

%feature("docstring") casadi::MatrixCommon::rows "

[INTERNAL] 
Get the number of rows, Octave-style syntax.

Extra doc: https://github.com/casadi/casadi/wiki/L_1at

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L114

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L114-L114

";

%feature("docstring") casadi::MatrixCommon::size2 "

[INTERNAL] 
Get the second dimension (i.e. number of columns)

Extra doc: https://github.com/casadi/casadi/wiki/L_1au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L119

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1149-L1151

";

%feature("docstring") casadi::MatrixCommon::columns "

[INTERNAL] 
Get the number of columns, Octave-style syntax.

Extra doc: https://github.com/casadi/casadi/wiki/L_1av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L124

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L124-L124

";

%feature("docstring") casadi::MatrixCommon::dim "

[INTERNAL] 
Get string representation of dimensions.

The representation is e.g. \"4x5\" or \"4x5,10nz\"

Extra doc: https://github.com/casadi/casadi/wiki/L_1aw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L131

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1164-L1166

";

%feature("docstring") casadi::MatrixCommon::size "

[INTERNAL] 
Get the size along a particular dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ay

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L141

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1159-L1161

>  casadi_int casadi::GenericMatrix< Matrix< Scalar >  >::size(casadi_int axis) const
------------------------------------------------------------------------
[INTERNAL] 
Get the size along a particular dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ay

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L141

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1159-L1161

";

";

%feature("docstring") casadi::MatrixCommon::is_empty "

[INTERNAL] 
Check if the sparsity is empty, i.e. if one of the dimensions is
 zero.

(or optionally both dimensions)

Extra doc: https://github.com/casadi/casadi/wiki/L_1az

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L148

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L148-L148

";

%feature("docstring") casadi::MatrixCommon::is_dense "

[INTERNAL] 
Check if the matrix expression is dense.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L153-L153

";

%feature("docstring") casadi::MatrixCommon::is_scalar "

[INTERNAL] 
Check if the matrix expression is scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L158

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1169-L1171

";

%feature("docstring") casadi::MatrixCommon::is_square "

[INTERNAL] 
Check if the matrix expression is square.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L163

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L163-L163

";

%feature("docstring") casadi::MatrixCommon::is_vector "

[INTERNAL] 
Check if the matrix is a row or column vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L168

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L168-L168

";

%feature("docstring") casadi::MatrixCommon::is_row "

[INTERNAL] 
Check if the matrix is a row vector (i.e.  size1()==1)

Extra doc: https://github.com/casadi/casadi/wiki/L_1b4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L173

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L173-L173

";

%feature("docstring") casadi::MatrixCommon::is_column "

[INTERNAL] 
Check if the matrix is a column vector (i.e.  size2()==1)

Extra doc: https://github.com/casadi/casadi/wiki/L_1b5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L178

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L178-L178

";

%feature("docstring") casadi::MatrixCommon::is_triu "

[INTERNAL] 
Check if the matrix is upper triangular.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L183

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L183-L183

";

%feature("docstring") casadi::MatrixCommon::is_tril "

[INTERNAL] 
Check if the matrix is lower triangular.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L188

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L188-L188

";

%feature("docstring") casadi::MatrixCommon::nz "

[INTERNAL] 
Access vector nonzero or slice of nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258-L260

>  NonZeros<Matrix< Scalar > , K> casadi::GenericMatrix< Matrix< Scalar >  >::nz(const K &k)
------------------------------------------------------------------------
[INTERNAL] 
Access vector nonzero or slice of nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258-L260

";

";


// File: classcasadi_1_1MX.xml


/*
 Construct symbolic primitives 
*/

/*
The \"sym\" function is intended to work in a similar way as \"sym\" 

used in the Symbolic Toolbox for Matlab but instead creating a CasADi 

symbolic primitive.

*/
%feature("docstring") casadi::MX "

[INTERNAL] 
 MX -  Matrix expression.

The  MX class is used to build up trees made up from MXNodes. It is a more 

general graph representation than the scalar expression, SX, and much 
less 
efficient for small objects. On the other hand, the class allows 
much more 
general operations than does SX, in particular matrix valued
 operations and
 calls to arbitrary differentiable functions.

The  MX class is designed to have identical syntax with the Matrix<> 
template
 class, and uses DM (i.e. Matrix<double>) as its internal 

representation of the values at a node. By keeping the syntaxes 
identical, 
it is possible to switch from one class to the other, as 
well as inlining  
MX functions to  SXElem functions.

Note that an operation is always \"lazy\", making a matrix 
multiplication 
will create a matrix multiplication node, not perform 
the actual 
multiplication.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_py

C++ includes: mx.hpp
";

%feature("docstring") casadi::MX::binary "

[INTERNAL] 
Create nodes by their ID.

Extra doc: https://github.com/casadi/casadi/wiki/L_r1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L398

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L507-L532

";

%feature("docstring") casadi::MX::unary "

[INTERNAL] 
Create nodes by their ID.

Extra doc: https://github.com/casadi/casadi/wiki/L_r1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L399

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L534-L536

";

%feature("docstring") casadi::MX::inf "

[INTERNAL] 
create a matrix with all inf

Extra doc: https://github.com/casadi/casadi/wiki/L_r2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L408

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L554-L556

>  MX casadi::MX::inf(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
create a matrix with all inf

Extra doc: https://github.com/casadi/casadi/wiki/L_r2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L408

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L554-L556

";

";

%feature("docstring") casadi::MX::nan "

[INTERNAL] 
create a matrix with all nan

Extra doc: https://github.com/casadi/casadi/wiki/L_r3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L417

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L566-L568

>  MX casadi::MX::nan(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
create a matrix with all nan

Extra doc: https://github.com/casadi/casadi/wiki/L_r3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L417

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L566-L568

";

";

%feature("docstring") casadi::MX::einstein "

[INTERNAL] 
Computes an einstein dense tensor contraction.

Computes the product: C_c = A_a + B_b where a b c are index/einstein 

notation in an encoded form

For example, an matrix-matrix product may be written as: C_ij = A_ik 
B_kj

The encoded form uses strictly negative numbers to indicate labels. 
For the
 above example, we would have: a {-1, -3} b {-3, -2} c {-1 -2}

Extra doc: https://github.com/casadi/casadi/wiki/L_r5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L520

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L646-L652

>  MX casadi::MX::einstein(const MX &A, const MX &B, const std::vector< casadi_int > &dim_a, const std::vector< casadi_int > &dim_b, const std::vector< casadi_int > &dim_c, const std::vector< casadi_int > &a, const std::vector< casadi_int > &b, const std::vector< casadi_int > &c)
------------------------------------------------------------------------
[INTERNAL] 
Computes an einstein dense tensor contraction.

Computes the product: C_c = A_a + B_b where a b c are index/einstein 

notation in an encoded form

For example, an matrix-matrix product may be written as: C_ij = A_ik 
B_kj

The encoded form uses strictly negative numbers to indicate labels. 
For the
 above example, we would have: a {-1, -3} b {-3, -2} c {-1 -2}

Extra doc: https://github.com/casadi/casadi/wiki/L_r5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L520

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L646-L652

";

";

%feature("docstring") casadi::MX::is_equal "

[INTERNAL] 
Functions called by friend functions defined for  
GenericExpression

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L531

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L810-L812

";

%feature("docstring") casadi::MX::mmin "

[INTERNAL] 
Functions called by friend functions defined for  
GenericExpression

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L532

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L814-L816

";

%feature("docstring") casadi::MX::mmax "

[INTERNAL] 
Functions called by friend functions defined for  
GenericExpression

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L533

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L818-L820

";

%feature("docstring") casadi::MX::horzcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L538

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L993-L1027

";

%feature("docstring") casadi::MX::diagcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L539

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1029-L1052

";

%feature("docstring") casadi::MX::vertcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L540

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1054-L1093

";

%feature("docstring") casadi::MX::horzsplit "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L541

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1095-L1110

";

%feature("docstring") casadi::MX::diagsplit "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L542

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1112-L1127

";

%feature("docstring") casadi::MX::vertsplit "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L544

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1129-L1150

";

%feature("docstring") casadi::MX::blockcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L545

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1152-L1171

";

%feature("docstring") casadi::MX::mtimes "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L546

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L628-L636

";

%feature("docstring") casadi::MX::mac "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L547

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L668-L689

";

%feature("docstring") casadi::MX::reshape "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L549

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1205-L1213

>  MX casadi::MX::reshape(const MX &x, const Sparsity &sp)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L549

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1205-L1213

";

";

%feature("docstring") casadi::MX::sparsity_cast "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L550

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1215-L1225

";

%feature("docstring") casadi::MX::kron "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L551

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1852-L1865

";

%feature("docstring") casadi::MX::repmat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L552

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1867-L1879

";

%feature("docstring") casadi::MX::jacobian "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L557

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1700-L1710

";

%feature("docstring") casadi::MX::hessian "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L559

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1717-L1726

>  MX casadi::MX::hessian(const MX &f, const MX &x, MX &g, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L559

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1717-L1726

";

";

%feature("docstring") casadi::MX::forward "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L561

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1729-L1757

";

%feature("docstring") casadi::MX::reverse "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L566

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1760-L1790

";

%feature("docstring") casadi::MX::which_depends "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L570

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1792-L1794

";

%feature("docstring") casadi::MX::jacobian_sparsity "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L572

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1796-L1798

";

%feature("docstring") casadi::MX::substitute "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L574

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1383-L1403

>  std::vector< MX > casadi::MX::substitute(const std::vector< MX > &ex, const std::vector< MX > &v, const std::vector< MX > &vdef)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L574

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1383-L1403

";

";

%feature("docstring") casadi::MX::substitute_inplace "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L577

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1357-L1377

";

%feature("docstring") casadi::MX::solve "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L581

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1903-L1907

>  MX casadi::MX::solve(const MX &a, const MX &b, const std::string &lsolver, const Dict &dict=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L581

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1903-L1907

";

";

%feature("docstring") casadi::MX::inv_minor "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L583

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1808-L1810

";

%feature("docstring") casadi::MX::inv_node "

[INTERNAL] 
Inverse node.

Extra doc: https://github.com/casadi/casadi/wiki/L_re

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L793

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L793-L795

>  MX casadi::MX::inv_node(const MX &x)
------------------------------------------------------------------------
[INTERNAL] 
Inverse node.

Extra doc: https://github.com/casadi/casadi/wiki/L_re

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L793

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L793-L795

";

";

%feature("docstring") casadi::MX::inv "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L585

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1812-L1814

";

%feature("docstring") casadi::MX::pinv "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L586

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1909-L1915

";

%feature("docstring") casadi::MX::expm_const "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L588

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1917-L1922

";

%feature("docstring") casadi::MX::expm "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L589

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1924-L1927

";

%feature("docstring") casadi::MX::n_nodes "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L590

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1329-L1333

";

%feature("docstring") casadi::MX::print_operator "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L591

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1353-L1355

";

%feature("docstring") casadi::MX::extract "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L592

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1520-L1691

";

%feature("docstring") casadi::MX::shared "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L594

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1693-L1698

";

%feature("docstring") casadi::MX::if_else "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L596

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1227-L1247

";

%feature("docstring") casadi::MX::conditional "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L598

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1249-L1281

";

%feature("docstring") casadi::MX::depends_on "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L600

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1935-L1951

";

%feature("docstring") casadi::MX::simplify "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L601

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1193-L1195

";

%feature("docstring") casadi::MX::dot "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L602

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L691-L693

";

%feature("docstring") casadi::MX::mrdivide "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L603

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L720-L723

";

%feature("docstring") casadi::MX::mldivide "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L604

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L725-L728

";

%feature("docstring") casadi::MX::norm_2 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L605

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1173-L1179

";

%feature("docstring") casadi::MX::norm_fro "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L606

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1181-L1183

";

%feature("docstring") casadi::MX::norm_1 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L607

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1185-L1187

";

%feature("docstring") casadi::MX::norm_inf "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L608

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1189-L1191

";

%feature("docstring") casadi::MX::unite "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L609

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1283-L1307

";

%feature("docstring") casadi::MX::trace "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L610

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1309-L1316

";

%feature("docstring") casadi::MX::diag "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L611

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1318-L1327

";

%feature("docstring") casadi::MX::sum2 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L612

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1335-L1337

";

%feature("docstring") casadi::MX::sum1 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L613

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1339-L1341

";

%feature("docstring") casadi::MX::polyval "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L614

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1343-L1351

";

%feature("docstring") casadi::MX::det "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L615

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1800-L1802

";

%feature("docstring") casadi::MX::symvar "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L616

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1816-L1819

";

%feature("docstring") casadi::MX::nullspace "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L617

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1929-L1933

";

%feature("docstring") casadi::MX::repsum "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L232

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1340-L1352

>  MX  casadi::GenericMatrix< MX  >::repsum(const MX &x, casadi_int n, casadi_int m=1)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L232

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1340-L1352

";

";

%feature("docstring") casadi::MX::densify "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L619

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L866-L877

";

%feature("docstring") casadi::MX::_bilin "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L620

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2275-L2277

";

%feature("docstring") casadi::MX::_rank1 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L621

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2279-L2281

";

%feature("docstring") casadi::MX::project "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L622

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L849-L864

";

%feature("docstring") casadi::MX::cumsum "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L623

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L654-L666

";

%feature("docstring") casadi::MX::_logsumexp "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L624

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2283-L2285

";

%feature("docstring") casadi::MX::cse "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L625

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2024-L2091

";

%feature("docstring") casadi::MX::find "

[INTERNAL] 
Find first nonzero, returned as row index.

If failed, returns the number of rows

Extra doc: https://github.com/casadi/casadi/wiki/L_r7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L690

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L690-L692

>  MX casadi::MX::find(const MX &x)
------------------------------------------------------------------------
[INTERNAL] 
Find first nonzero, returned as row index.

If failed, returns the number of rows

Extra doc: https://github.com/casadi/casadi/wiki/L_r7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L690

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L690-L692

";

";

%feature("docstring") casadi::MX::low "

[INTERNAL] 
Find first nonzero.

If failed, returns the number of rows

Extra doc: https://github.com/casadi/casadi/wiki/L_r8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L699

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L699-L701

>  MX casadi::MX::low(const MX &v, const MX &p, const Dict &options=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Find first nonzero.

If failed, returns the number of rows

Extra doc: https://github.com/casadi/casadi/wiki/L_r8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L699

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L699-L701

";

";

%feature("docstring") casadi::MX::graph_substitute "

[INTERNAL] 
Substitute multiple expressions in graph.

Substitute variable var with expression expr in multiple expressions, 

preserving nodes

Extra doc: https://github.com/casadi/casadi/wiki/L_ra

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L720

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L720-L724

>  std::vector<MX> casadi::MX::graph_substitute(const std::vector< MX > &ex, const std::vector< MX > &v, const std::vector< MX > &vdef)
------------------------------------------------------------------------
[INTERNAL] 
Substitute multiple expressions in graph.

Substitute variable var with expression expr in multiple expressions, 

preserving nodes

Extra doc: https://github.com/casadi/casadi/wiki/L_ra

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L720

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L720-L724

";

";

%feature("docstring") casadi::MX::matrix_expand "

[INTERNAL] 
Expand  MX graph to SXFunction call.

Expand the given expression e, optionally supplying expressions 
contained 
in it at which expansion should stop.

Extra doc: https://github.com/casadi/casadi/wiki/L_rc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L745

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L745-L749

>  std::vector<MX> casadi::MX::matrix_expand(const std::vector< MX > &e, const std::vector< MX > &boundary=std::vector< MX >(), const Dict &options=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Expand  MX graph to SXFunction call.

Expand the given expression e, optionally supplying expressions 
contained 
in it at which expansion should stop.

Extra doc: https://github.com/casadi/casadi/wiki/L_rc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L745

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L745-L749

";

";

%feature("docstring") casadi::MX::lift "

[INTERNAL] 
Lift the expression.

Experimental feature

Extra doc: https://github.com/casadi/casadi/wiki/L_rd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L786

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L786-L788

>  MX casadi::MX::lift(const MX &x, const MX &x_guess)
------------------------------------------------------------------------
[INTERNAL] 
Lift the expression.

Experimental feature

Extra doc: https://github.com/casadi/casadi/wiki/L_rd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L786

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L786-L788

";

";

%feature("docstring") casadi::MX::evalf "

[INTERNAL] 
Evaluates the expression numerically.

An error is raised when the expression contains symbols

Extra doc: https://github.com/casadi/casadi/wiki/L_rf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L802

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L802-L804

>  DM casadi::MX::evalf(const MX &expr)
------------------------------------------------------------------------
[INTERNAL] 
Evaluates the expression numerically.

An error is raised when the expression contains symbols

Extra doc: https://github.com/casadi/casadi/wiki/L_rf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L802

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L802-L804

";

";

%feature("docstring") casadi::MX::bspline "

[INTERNAL]

>  MX casadi::MX::bspline(const MX &x, const DM &coeffs, const std::vector< std::vector< double > > &knots, const std::vector< casadi_int > &degree, casadi_int m, const Dict &opts=Dict())

>  MX casadi::MX::bspline(const MX &x, const MX &coeffs, const std::vector< std::vector< double > > &knots, const std::vector< casadi_int > &degree, casadi_int m, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::MX::convexify "

[INTERNAL]

>  MX casadi::MX::convexify(const MX &H, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::MX::stop_diff "

[INTERNAL] 
Stop derivatives of an expression wrt to a select set of 
symbolic 
variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_25o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L835

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L835-L837

>  MX casadi::MX::stop_diff(const MX &expr, const MX &var, casadi_int order)
------------------------------------------------------------------------
[INTERNAL] 
Stop derivatives of an expression wrt to a select set of symbolic 
variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_25o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L835

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L835-L837

";

";

%feature("docstring") casadi::MX::interp1d "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L213

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1252-L1308

";

%feature("docstring") casadi::MX::sprank "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L215

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L215-L215

";

%feature("docstring") casadi::MX::norm_0_mul "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L216

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L216-L218

";

%feature("docstring") casadi::MX::tril "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L219

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L219-L221

";

%feature("docstring") casadi::MX::triu "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L222-L224

";

%feature("docstring") casadi::MX::sumsqr "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L225

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L225-L225

";

%feature("docstring") casadi::MX::linspace "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L226

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1202-L1212

";

%feature("docstring") casadi::MX::cross "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L227

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1215-L1246

";

%feature("docstring") casadi::MX::skew "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L228

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1311-L1319

";

%feature("docstring") casadi::MX::inv_skew "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L229

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1322-L1327

";

%feature("docstring") casadi::MX::tril2symm "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L230

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1331-L1337

";

%feature("docstring") casadi::MX::triu2symm "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L231

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1355-L1361

";

%feature("docstring") casadi::MX::diff "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1567-L1592

";

%feature("docstring") casadi::MX::is_linear "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L235

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1536-L1538

";

%feature("docstring") casadi::MX::is_quadratic "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L236

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1541-L1543

";

%feature("docstring") casadi::MX::quadratic_coeff "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L237

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1546-L1554

";

%feature("docstring") casadi::MX::linear_coeff "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L239

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1557-L1564

";

%feature("docstring") casadi::MX::mpower "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1010

";

%feature("docstring") casadi::MX::soc "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1011

";

%feature("docstring") casadi::MX::linearize "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1008

";

%feature("docstring") casadi::MX::gradient "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1006

";

%feature("docstring") casadi::MX::tangent "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1007

";

%feature("docstring") casadi::MX::jtimes "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1004

";

%feature("docstring") casadi::MX::bilin "

[INTERNAL] 
Calculate bilinear/quadratic form x^T A y.

Parameters:
-----------

y: 
can be omitted, in which case x^T A x is calculated

Extra doc: https://github.com/casadi/casadi/wiki/L_1bo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L409

";

%feature("docstring") casadi::MX::rank1 "

[INTERNAL] 
Make a rank-1 update to a matrix A.

Calculates A + 1/2 * alpha * x*y'

Extra doc: https://github.com/casadi/casadi/wiki/L_1bp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L422

";

%feature("docstring") casadi::MX::sym "

[INTERNAL] 
Create a vector of length r of vectors of length p.

with nrow-by-ncol symbolic primitives

Extra doc: https://github.com/casadi/casadi/wiki/L_1dg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1074

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1074-L1076

>  static std::vector<std::vector<MX > > casadi::GenericMatrix< MX  >::sym(const std::string &name, casadi_int nrow, casadi_int ncol, casadi_int p, casadi_int r)
------------------------------------------------------------------------
[INTERNAL] 
Create a vector of length r of vectors of length p.

with nrow-by-ncol symbolic primitives

Extra doc: https://github.com/casadi/casadi/wiki/L_1dg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1074

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1074-L1076

";

";

%feature("docstring") casadi::MX::zeros "

[INTERNAL] 
Create a dense matrix or a matrix with specified sparsity with 
all 
entries zero.

Extra doc: https://github.com/casadi/casadi/wiki/L_1dh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1087

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1087-L1089

>  static MX  casadi::GenericMatrix< MX  >::zeros(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
Create a dense matrix or a matrix with specified sparsity with all 
entries zero.

Extra doc: https://github.com/casadi/casadi/wiki/L_1dh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1087

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1087-L1089

";

";

%feature("docstring") casadi::MX::ones "

[INTERNAL] 
Create a dense matrix or a matrix with specified sparsity with 
all 
entries one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1di

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100-L1102

>  static MX  casadi::GenericMatrix< MX  >::ones(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
Create a dense matrix or a matrix with specified sparsity with all 
entries one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1di

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100-L1102

";

";

%feature("docstring") casadi::MX::get "

[INTERNAL] 
Get a const pointer to the node.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L427

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L538-L540

>  MXNode * casadi::MX::get() const
------------------------------------------------------------------------
[INTERNAL] 
Get a const pointer to the node.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L427

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L538-L540

";

";

%feature("docstring") casadi::MX::set "

[INTERNAL] 
Set a submatrix, two arguments

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L475

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L238-L296

>  void casadi::MX::set(const MX &m, bool ind1, const Matrix< casadi_int > &rr, const Matrix< casadi_int > &cc)
------------------------------------------------------------------------
[INTERNAL] 
Set a submatrix, two arguments

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L475

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L238-L296

";

";

%feature("docstring") casadi::MX::get_nz "

[INTERNAL] 
Get a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L488

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L427-L430

>  void casadi::MX::get_nz(MX &m, bool ind1, const MX &inner, const MX &outer) const
------------------------------------------------------------------------
[INTERNAL] 
Get a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L488

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L427-L430

";

";

%feature("docstring") casadi::MX::set_nz "

[INTERNAL] 
Set a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L496

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L496-L496

>  void casadi::MX::set_nz(const MX &m, bool ind1, casadi_int kk)
------------------------------------------------------------------------
[INTERNAL] 
Set a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L496

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L496-L496

";

";

%feature("docstring") casadi::MX::ad_forward "

[INTERNAL] 
Called from MXFunction.

Extra doc: https://github.com/casadi/casadi/wiki/L_ro

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L904

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2296-L2303

";

%feature("docstring") casadi::MX::ad_reverse "

[INTERNAL] 
Called from MXFunction.

Extra doc: https://github.com/casadi/casadi/wiki/L_ro

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L906

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2305-L2312

";

%feature("docstring") casadi::MX::get_row "

[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L194

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L194-L194

";

%feature("docstring") casadi::MX::get_colind "

[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L195-L195

";

%feature("docstring") casadi::MX::row "

[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200-L200

>  casadi_int casadi::GenericMatrix< MX  >::row(casadi_int el) const
------------------------------------------------------------------------
[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200-L200

";

";

%feature("docstring") casadi::MX::colind "

[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201-L201

>  casadi_int casadi::GenericMatrix< MX  >::colind(casadi_int col) const
------------------------------------------------------------------------
[INTERNAL] 
Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201-L201

";

";

%feature("docstring") casadi::MX::printme "

[INTERNAL]

>  MX casadi::MX::printme(const MX &b) const
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::MX::MX "

[INTERNAL] 
Construct constant matrix with a given sparsity and values.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L911

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L114-L116

>  casadi::MX::MX(const Sparsity &sp, double val, bool dummy)
------------------------------------------------------------------------
[INTERNAL] 
Construct constant matrix with a given sparsity and values.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L911

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L114-L116

";

";

%feature("docstring") casadi::MX::sparsity "

[INTERNAL] 
Get the sparsity pattern.

Extra doc: https://github.com/casadi/casadi/wiki/L_qc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L165

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L586-L588

";

%feature("docstring") casadi::MX::__nonzero__ "

[INTERNAL] 
Returns the truth value of an  MX expression.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L184

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L137-L139

";

%feature("docstring") casadi::MX::get_sparsity "

[INTERNAL] 
Get an owning reference to the sparsity pattern.

Extra doc: https://github.com/casadi/casadi/wiki/L_qd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L189-L189

";

%feature("docstring") casadi::MX::erase "

[INTERNAL] 
Erase a submatrix (leaving structural zeros in its place)

Erase elements of a matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_qf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L204

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L604-L616

>  void casadi::MX::erase(const std::vector< casadi_int > &rr, bool ind1=false)
------------------------------------------------------------------------
[INTERNAL] 
Erase a submatrix (leaving structural zeros in its place)

Erase elements of a matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_qf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L204

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L604-L616

";

";

%feature("docstring") casadi::MX::enlarge "

[INTERNAL] 
Enlarge matrix.

Make the matrix larger by inserting empty rows and columns, keeping 
the 
existing non-zeros

Extra doc: https://github.com/casadi/casadi/wiki/L_qg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L211

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L618-L626

";

%feature("docstring") casadi::MX::dep "

[INTERNAL] 
Get the nth dependency as  MX.

Extra doc: https://github.com/casadi/casadi/wiki/L_qj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L236

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L730-L732

";

%feature("docstring") casadi::MX::n_out "

[INTERNAL] 
Number of outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_qk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L241

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L841-L843

";

%feature("docstring") casadi::MX::get_output "

[INTERNAL] 
Get an output.

Extra doc: https://github.com/casadi/casadi/wiki/L_ql

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L246

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L845-L847

";

%feature("docstring") casadi::MX::n_dep "

[INTERNAL] 
Get the number of dependencies of a binary  SXElem.

Extra doc: https://github.com/casadi/casadi/wiki/L_qm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L251

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L734-L736

";

%feature("docstring") casadi::MX::name "

[INTERNAL] 
Get the name.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L254

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L738-L740

";

%feature("docstring") casadi::MX::is_symbolic "

[INTERNAL] 
Check if symbolic.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L263

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L742-L744

";

%feature("docstring") casadi::MX::is_constant "

[INTERNAL] 
Check if constant.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L266

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L746-L748

";

%feature("docstring") casadi::MX::is_call "

[INTERNAL] 
Check if evaluation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L269

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L750-L752

";

%feature("docstring") casadi::MX::which_function "

[INTERNAL] 
Get function - only valid when  is_call() is true.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L272

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L754-L756

";

%feature("docstring") casadi::MX::is_output "

[INTERNAL] 
Check if evaluation output.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L275

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L758-L760

";

%feature("docstring") casadi::MX::which_output "

[INTERNAL] 
Get the index of evaluation output - only valid when  
is_output() is true.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L278

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L762-L764

";

%feature("docstring") casadi::MX::is_op "

[INTERNAL] 
Is it a certain operation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L281

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L766-L768

";

%feature("docstring") casadi::MX::is_multiplication "

[INTERNAL] 
Check if multiplication.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L284

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L770-L772

";

%feature("docstring") casadi::MX::is_commutative "

[INTERNAL] 
Check if commutative operation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L287

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L822-L827

";

%feature("docstring") casadi::MX::is_norm "

[INTERNAL] 
Check if norm.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L290

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L774-L776

";

%feature("docstring") casadi::MX::is_valid_input "

[INTERNAL] 
Check if matrix can be used to define function inputs.

Valid inputs for MXFunctions are combinations of Reshape, 
concatenations 
and SymbolicMX

Extra doc: https://github.com/casadi/casadi/wiki/L_qn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L297

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L897-L899

";

%feature("docstring") casadi::MX::n_primitives "

[INTERNAL] 
Get the number of primitives for MXFunction inputs/outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_qo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L302

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L901-L903

";

%feature("docstring") casadi::MX::primitives "

[INTERNAL] 
Get primitives.

Extra doc: https://github.com/casadi/casadi/wiki/L_qp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L307

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L905-L911

";

%feature("docstring") casadi::MX::split_primitives "

[INTERNAL] 
Split up an expression along symbolic primitives.

Extra doc: https://github.com/casadi/casadi/wiki/L_qq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L312

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L913-L919

";

%feature("docstring") casadi::MX::join_primitives "

[INTERNAL] 
Join an expression along symbolic primitives.

Extra doc: https://github.com/casadi/casadi/wiki/L_qr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L317

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L921-L927

";

%feature("docstring") casadi::MX::is_eye "

[INTERNAL] 
check if identity

Extra doc: https://github.com/casadi/casadi/wiki/L_qu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L339

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L937-L939

";

%feature("docstring") casadi::MX::is_zero "

[INTERNAL] 
check if zero (note that false negative answers are possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_qv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L344

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L941-L947

";

%feature("docstring") casadi::MX::is_one "

[INTERNAL] 
check if zero (note that false negative answers are possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_qw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L349

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L949-L951

";

%feature("docstring") casadi::MX::is_minus_one "

[INTERNAL] 
check if zero (note that false negative answers are possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_qx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L354

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L953-L955

";

%feature("docstring") casadi::MX::is_transpose "

[INTERNAL] 
Is the expression a transpose?

Extra doc: https://github.com/casadi/casadi/wiki/L_qy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L359

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L957-L959

";

%feature("docstring") casadi::MX::is_regular "

[INTERNAL] 
Checks if expression does not contain NaN or Inf.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L362

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L961-L967

";

%feature("docstring") casadi::MX::is_binary "

[INTERNAL] 
Is binary operation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L365

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L786-L788

";

%feature("docstring") casadi::MX::is_unary "

[INTERNAL] 
Is unary operation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L368

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L790-L792

";

%feature("docstring") casadi::MX::op "

[INTERNAL] 
Get operation type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L371

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L794-L796

";

%feature("docstring") casadi::MX::info "

[INTERNAL] 
Obtain information about node

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L374

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L798-L800

";

%feature("docstring") casadi::MX::serialize "

[INTERNAL] 
Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_qz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L379

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L802-L804

";

%feature("docstring") casadi::MX::attachAssert "

[INTERNAL] 
returns itself, but with an assertion attached

If y does not evaluate to 1, a runtime error is raised

Extra doc: https://github.com/casadi/casadi/wiki/L_rg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L851

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L699-L704

";

%feature("docstring") casadi::MX::monitor "

[INTERNAL] 
Monitor an expression.

Returns itself, but with the side effect of printing the nonzeros 
along 
with a comment

Extra doc: https://github.com/casadi/casadi/wiki/L_rh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L858

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L706-L708

";

%feature("docstring") casadi::MX::T "

[INTERNAL] 
Transpose the matrix.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L861

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L969-L971

";

%feature("docstring") casadi::MX::mapping "

[INTERNAL] 
Get an IM representation of a GetNonzeros or SetNonzeros node.

Extra doc: https://github.com/casadi/casadi/wiki/L_ri

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L866

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L829-L831

";

%feature("docstring") casadi::MX::eval_mx "

[INTERNAL] 
Evaluate the  MX node with new symbolic dependencies.

Extra doc: https://github.com/casadi/casadi/wiki/L_rn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L897

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2288-L2294

";

%feature("docstring") casadi::MX::nnz "

[INTERNAL] 
Get the number of (structural) non-zero elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1an

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L84

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1119-L1121

";

%feature("docstring") casadi::MX::nnz_lower "

[INTERNAL] 
Get the number of non-zeros in the lower triangular half.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ao

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L89

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1124-L1126

";

%feature("docstring") casadi::MX::nnz_upper "

[INTERNAL] 
Get the number of non-zeros in the upper triangular half.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ap

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L94

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1129-L1131

";

%feature("docstring") casadi::MX::nnz_diag "

[INTERNAL] 
Get get the number of non-zeros on the diagonal.

Extra doc: https://github.com/casadi/casadi/wiki/L_1aq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L99

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1134-L1136

";

%feature("docstring") casadi::MX::numel "

[INTERNAL] 
Get the number of elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ar

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L104

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1139-L1141

";

%feature("docstring") casadi::MX::size1 "

[INTERNAL] 
Get the first dimension (i.e. number of rows)

Extra doc: https://github.com/casadi/casadi/wiki/L_1as

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L109

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1144-L1146

";

%feature("docstring") casadi::MX::rows "

[INTERNAL] 
Get the number of rows, Octave-style syntax.

Extra doc: https://github.com/casadi/casadi/wiki/L_1at

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L114

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L114-L114

";

%feature("docstring") casadi::MX::size2 "

[INTERNAL] 
Get the second dimension (i.e. number of columns)

Extra doc: https://github.com/casadi/casadi/wiki/L_1au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L119

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1149-L1151

";

%feature("docstring") casadi::MX::columns "

[INTERNAL] 
Get the number of columns, Octave-style syntax.

Extra doc: https://github.com/casadi/casadi/wiki/L_1av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L124

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L124-L124

";

%feature("docstring") casadi::MX::dim "

[INTERNAL] 
Get string representation of dimensions.

The representation is e.g. \"4x5\" or \"4x5,10nz\"

Extra doc: https://github.com/casadi/casadi/wiki/L_1aw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L131

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1164-L1166

";

%feature("docstring") casadi::MX::size "

[INTERNAL] 
Get the size along a particular dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ay

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L141

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1159-L1161

>  casadi_int casadi::GenericMatrix< MX  >::size(casadi_int axis) const
------------------------------------------------------------------------
[INTERNAL] 
Get the size along a particular dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ay

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L141

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1159-L1161

";

";

%feature("docstring") casadi::MX::is_empty "

[INTERNAL] 
Check if the sparsity is empty, i.e. if one of the dimensions is
 zero.

(or optionally both dimensions)

Extra doc: https://github.com/casadi/casadi/wiki/L_1az

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L148

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L148-L148

";

%feature("docstring") casadi::MX::is_dense "

[INTERNAL] 
Check if the matrix expression is dense.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L153-L153

";

%feature("docstring") casadi::MX::is_scalar "

[INTERNAL] 
Check if the matrix expression is scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L158

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1169-L1171

";

%feature("docstring") casadi::MX::is_square "

[INTERNAL] 
Check if the matrix expression is square.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L163

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L163-L163

";

%feature("docstring") casadi::MX::is_vector "

[INTERNAL] 
Check if the matrix is a row or column vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L168

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L168-L168

";

%feature("docstring") casadi::MX::is_row "

[INTERNAL] 
Check if the matrix is a row vector (i.e.  size1()==1)

Extra doc: https://github.com/casadi/casadi/wiki/L_1b4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L173

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L173-L173

";

%feature("docstring") casadi::MX::is_column "

[INTERNAL] 
Check if the matrix is a column vector (i.e.  size2()==1)

Extra doc: https://github.com/casadi/casadi/wiki/L_1b5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L178

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L178-L178

";

%feature("docstring") casadi::MX::is_triu "

[INTERNAL] 
Check if the matrix is upper triangular.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L183

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L183-L183

";

%feature("docstring") casadi::MX::is_tril "

[INTERNAL] 
Check if the matrix is lower triangular.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L188

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L188-L188

";

%feature("docstring") casadi::MX::nz "

[INTERNAL] 
Access vector nonzero or slice of nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258-L260

>  NonZeros<MX , K> casadi::GenericMatrix< MX  >::nz(const K &k)
------------------------------------------------------------------------
[INTERNAL] 
Access vector nonzero or slice of nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258-L260

";

";

%feature("docstring") casadi::MX::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L99-L101

";

%feature("docstring") casadi::MX::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L103-L109

";

%feature("docstring") casadi::MX::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::MX::is_null "

[INTERNAL] 
Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L73-L75

";

%feature("docstring") casadi::MX::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L129-L131

";

%feature("docstring") casadi::MX::bspline_dual "

[INTERNAL] ";

%feature("docstring") casadi::MX::no_grad "

[INTERNAL] 
Stop first derivatives of an expression wrt to all its symbolic
 
variables.

\\\\seealso stop_diff

Extra doc: https://github.com/casadi/casadi/wiki/L_25m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L818

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L818-L820

";

%feature("docstring") casadi::MX::no_hess "

[INTERNAL] 
Stop second derivatives of an expression wrt to all its symbolic
 
variables.

\\\\seealso stop_diff

Extra doc: https://github.com/casadi/casadi/wiki/L_25n

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L827

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L827-L829

";

%feature("docstring") casadi::MX::difference "

[INTERNAL] 
\\\\bried Return all elements of a that do not occur in b, 
preserving 
order

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L844

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2165-L2179

";


// File: classcasadi_1_1Newton.xml
%feature("docstring") casadi::Newton "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1NlImporter.xml
%feature("docstring") casadi::NlImporter "

[INTERNAL] 
\\\\Helper class for .nl import The .nl format is described in 

\"Writing .nl Files\" paper by David M. Gay (2005) 
Joel Andersson

C++ includes: nlp_builder.hpp
";

%feature("docstring") casadi::NlImporter::NlImporter "

[INTERNAL] ";

%feature("docstring") casadi::NlImporter::~NlImporter "

[INTERNAL] ";


// File: classcasadi_1_1NlpBuilder.xml


/*
 Symbolic representation of the NLP 
*/

/*
Data members

*/
%feature("docstring") casadi::NlpBuilder "

[INTERNAL] 
A symbolic NLP representation.

Joel Andersson

::

  Extra doc: https://github.com/casadi/casadi/wiki/L_1e2 
  



C++ includes: nlp_builder.hpp
";

%feature("docstring") casadi::NlpBuilder::import_nl "

[INTERNAL] 
Import an .nl file.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.hpp#L74

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.cpp#L32-L35

";

%feature("docstring") casadi::NlpBuilder::type_name "

[INTERNAL] 
Readable name of the class.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.hpp#L77

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.hpp#L77-L77

";

%feature("docstring") casadi::NlpBuilder::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.hpp#L80

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.cpp#L37-L45

";

%feature("docstring") casadi::NlpBuilder::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.hpp#L83

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.hpp#L83-L87

";


// File: classcasadi_1_1Nlpsol.xml
%feature("docstring") casadi::Nlpsol "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1NonZeros.xml
%feature("docstring") casadi::NonZeros "

[INTERNAL] 
Access to a set of nonzeros.

NonZeros class for  Matrix NonZeros is the return type for operator[] of the
  Matrix class, it allows access to the value as well as changing the parent
 
object 
Joel Andersson

C++ includes: nonzeros.hpp
";

%feature("docstring") casadi::NonZeros::NonZeros "

[INTERNAL] 
Default copy constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nonzeros.hpp#L46

>  casadi::NonZeros< M, K >::NonZeros(const NonZeros< M, K > &y)=default
------------------------------------------------------------------------
[INTERNAL] 
Default copy constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nonzeros.hpp#L46

";

";


// File: classcasadi_1_1null__ptr__on__copy.xml
%feature("docstring") casadi::null_ptr_on_copy "

[INTERNAL] 
Pointer that gets set to null when copied.

C++ includes: optistack_internal.hpp
";

%feature("docstring") casadi::null_ptr_on_copy::null_ptr_on_copy "

[INTERNAL] ";


// File: classstd_1_1numeric__limits_3_01casadi_1_1SXElem_01_4.xml
%feature("docstring") std::numeric_limits< casadi::SXElem > "

[INTERNAL] C++ includes: sx_elem.hpp
";


// File: classcasadi_1_1Opti.xml
%feature("docstring") casadi::Opti "

[INTERNAL] 
A simplified interface for NLP modeling/solving.

This class offers a view with model description facilities The API is 

guaranteed to be stable.

Example NLP:

::

    opti = casadi.Opti();
  
    x = opti.variable();
    y = opti.variable();
  
    opti.minimize(  (y-x^2)^2   );
    opti.subject_to( x^2+y^2==1 );
    opti.subject_to(     x+y>=1 );
  
    opti.solver('ipopt');
    sol = opti.solve();
  
    sol.value(x)
    sol.value(y)



Example parametric NLP:

::

    opti = casadi.Opti();
  
    x = opti.variable(2,1);
    p = opti.parameter();
  
    opti.minimize(  (p*x(2)-x(1)^2)^2   );
    opti.subject_to( 1<=sum(x)<=2 );
  
    opti.solver('ipopt');
  
    opti.set_value(p, 3);
    sol = opti.solve();
    sol.value(x)
  
    opti.set_value(p, 5);
    sol = opti.solve();
    sol.value(x)



Joris Gillis, Erik Lambrechts, Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_16

C++ includes: optistack.hpp
";

%feature("docstring") casadi::Opti::subject_to "

[INTERNAL] 
Clear constraints.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L166

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L102-L108

>  void casadi::Opti::subject_to()
------------------------------------------------------------------------
[INTERNAL] 
Clear constraints.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L166

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L102-L108

";

";

%feature("docstring") casadi::Opti::set_initial "

[INTERNAL] 
Set initial guess for decision variables

::

  * opti.set_initial(x, 2)
  * opti.set_initial(10*x(1), 2)
  * 



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L190

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L128-L134

>  void casadi::Opti::set_initial(const std::vector< MX > &assignments)
------------------------------------------------------------------------
[INTERNAL] 
Set initial guess for decision variables

::

  * opti.set_initial(x, 2)
  * opti.set_initial(10*x(1), 2)
  * 



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L190

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L128-L134

";

";

%feature("docstring") casadi::Opti::set_value "

[INTERNAL] 
Set value of parameter.

Each parameter must be given a value before 'solve' can be called

Extra doc: https://github.com/casadi/casadi/wiki/L_1d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L145-L151

>  void casadi::Opti::set_value(const std::vector< MX > &assignments)
------------------------------------------------------------------------
[INTERNAL] 
Set value of parameter.

Each parameter must be given a value before 'solve' can be called

Extra doc: https://github.com/casadi/casadi/wiki/L_1d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L145-L151

";

";

%feature("docstring") casadi::Opti::value "

[INTERNAL] 
Obtain value of expression at the current value

In regular mode, teh current value is the converged solution In debug 
mode,
 the value can be non-converged

Parameters:
-----------

values: 
Optional assignment expressions (e.g. x==3) to overrule the current
 
value

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L225

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L186-L192

>  DM casadi::Opti::value(const SX &x, const std::vector< MX > &values=std::vector< MX >()) const
------------------------------------------------------------------------
[INTERNAL] 
Obtain value of expression at the current value

In regular mode, teh current value is the converged solution In debug 
mode,
 the value can be non-converged

Parameters:
-----------

values: 
Optional assignment expressions (e.g. x==3) to overrule the current
 
value

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L225

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L186-L192

";

";

%feature("docstring") casadi::Opti::to_function "

[INTERNAL] 
Create a CasADi  Function from the  Opti solver.

Parameters:
-----------

name: 
Name of the resulting CasADi  Function

args: 
List of parameters and decision/dual variables (which can be given an
 
initial guess) with the resulting  Function

res: 
List of expressions that will get evaluated at the optimal solution

opts: 
Standard CasADi Funcion options

Extra doc: https://github.com/casadi/casadi/wiki/L_1j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L336

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L341-L361

>  Function casadi::Opti::to_function(const std::string &name, const std::map< std::string, MX > &dict, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Create a CasADi  Function from the  Opti solver.

Parameters:
-----------

name: 
Name of the resulting CasADi  Function

args: 
List of parameters and decision/dual variables (which can be given an
 
initial guess) with the resulting  Function

res: 
List of expressions that will get evaluated at the optimal solution

opts: 
Standard CasADi Funcion options

Extra doc: https://github.com/casadi/casadi/wiki/L_1j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L336

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L341-L361

";

";

%feature("docstring") casadi::Opti::callback_class "

[INTERNAL] 
Helper methods for callback()

Do not use directly.

Extra doc: https://github.com/casadi/casadi/wiki/L_1p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L409

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L371-L377

>  void casadi::Opti::callback_class()
------------------------------------------------------------------------
[INTERNAL] 
Helper methods for callback()

Do not use directly.

Extra doc: https://github.com/casadi/casadi/wiki/L_1p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L409

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L371-L377

";

";

%feature("docstring") casadi::Opti::Opti "

[INTERNAL]

>  casadi::Opti::Opti(const Opti &x)

>  casadi::Opti::Opti(OptiNode *node)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::Opti::variable "

[INTERNAL] 
Create a decision variable (symbol)

The order of creation matters. The order will be reflected in the 

optimization problem. It is not required for decision variables to 
actualy 
appear in the optimization problem.

Parameters:
-----------

n: 
number of rows (default 1)

m: 
number of columnss (default 1)

attribute: 
'full' (default) or 'symmetric'

Extra doc: https://github.com/casadi/casadi/wiki/L_18

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L112

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L51-L57

";

%feature("docstring") casadi::Opti::parameter "

[INTERNAL] 
Create a parameter (symbol); fixed during optimization.

The order of creation does not matter. It is not required for 
parameter to 
actualy appear in the optimization problem. Parameters 
that do appear, must
 be given a value before the problem can be 
solved.

Parameters:
-----------

n: 
number of rows (default 1)

m: 
number of columnss (default 1)

attribute: 
'full' (default) or 'symmetric'

Extra doc: https://github.com/casadi/casadi/wiki/L_19

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L125

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L74-L80

";

%feature("docstring") casadi::Opti::minimize "

[INTERNAL] 
Set objective.

Objective must be a scalar. Default objective: 0 When method is called
 
multiple times, the last call takes effect

Extra doc: https://github.com/casadi/casadi/wiki/L_1a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L133

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L82-L88

";

%feature("docstring") casadi::Opti::solver "

[INTERNAL] 
Set a solver.

Parameters:
-----------

solver: 
any of the nlpsol plugins can be used here In practice, not all 
nlpsol
 plugins may be supported yet

options: 
passed on to nlpsol plugin No stability can be guaranteed about 
this 
part of the API

options: 
to be passed to nlpsol solver No stability can be guaranteed about
 
this part of the API

Extra doc: https://github.com/casadi/casadi/wiki/L_1c

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L178

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L111-L119

";

%feature("docstring") casadi::Opti::solve "

[INTERNAL] 
Crunch the numbers; solve the problem.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L204

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L153-L159

";

%feature("docstring") casadi::Opti::solve_limited "

[INTERNAL] 
Crunch the numbers; solve the problem.

Allows the solver to return without error when an iteration or time 
limit 
is reached

Extra doc: https://github.com/casadi/casadi/wiki/L_1e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L212

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L161-L167

";

%feature("docstring") casadi::Opti::stats "

[INTERNAL] 
Get statistics.

nlpsol stats are passed as-is. No stability can be guaranteed about 
this 
part of the API

Extra doc: https://github.com/casadi/casadi/wiki/L_1f

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L234

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L194-L200

";

%feature("docstring") casadi::Opti::return_status "

[INTERNAL] 
Get return status of solver.



::

     passed as-is from nlpsol
  

No stability can be guaranteed about this part of the API

Extra doc: https://github.com/casadi/casadi/wiki/L_1g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L202-L208

";

%feature("docstring") casadi::Opti::initial "

[INTERNAL] 
get assignment expressions for initial values

Extra doc: https://github.com/casadi/casadi/wiki/L_266

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L247

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L210-L216

";

%feature("docstring") casadi::Opti::value_variables "

[INTERNAL] 
get assignment expressions for latest values

Extra doc: https://github.com/casadi/casadi/wiki/L_267

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L252

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L218-L224

";

%feature("docstring") casadi::Opti::value_parameters "

[INTERNAL] ";

%feature("docstring") casadi::Opti::dual "

[INTERNAL] 
get the dual variable

m must be a constraint expression. The returned value is still a 
symbolic 
expression. Use  value on it to obtain the numerical value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L262

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L234-L240

";

%feature("docstring") casadi::Opti::nx "

[INTERNAL] 
Number of (scalarised) decision variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_268

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L267

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L242-L248

";

%feature("docstring") casadi::Opti::np "

[INTERNAL] 
Number of (scalarised) parameters.

Extra doc: https://github.com/casadi/casadi/wiki/L_269

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L272

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L250-L256

";

%feature("docstring") casadi::Opti::ng "

[INTERNAL] 
Number of (scalarised) constraints.

Extra doc: https://github.com/casadi/casadi/wiki/L_26a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L277

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L258-L264

";

%feature("docstring") casadi::Opti::x "

[INTERNAL] 
Get all (scalarised) decision variables as a symbolic column 
vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_26b

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L282

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L266-L272

";

%feature("docstring") casadi::Opti::p "

[INTERNAL] 
Get all (scalarised) parameters as a symbolic column vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_26c

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L287

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L274-L280

";

%feature("docstring") casadi::Opti::g "

[INTERNAL] 
Get all (scalarised) constraint expressions as a column vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_26d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L292

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L282-L288

";

%feature("docstring") casadi::Opti::f "

[INTERNAL] 
Get objective expression.

Extra doc: https://github.com/casadi/casadi/wiki/L_26e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L297

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L290-L296

";

%feature("docstring") casadi::Opti::lbg "

[INTERNAL] 
Get all (scalarised) bounds on constraints as a column vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_26f

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L302

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L298-L304

";

%feature("docstring") casadi::Opti::ubg "

[INTERNAL] ";

%feature("docstring") casadi::Opti::lam_g "

[INTERNAL] 
Get all (scalarised) dual variables as a symbolic column vector.

Useful for obtaining the Lagrange Hessian:

::

  * sol.value(hessian(opti.f+opti.lam_g'*opti.g,opti.x)) % MATLAB
  * sol.value(hessian(opti.f+dot(opti.lam_g,opti.g),opti.x)[0]) # Python
  * 



Extra doc: https://github.com/casadi/casadi/wiki/L_1i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L314

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L315-L321

";

%feature("docstring") casadi::Opti::debug "

[INTERNAL] 
Get a copy with advanced functionality.

You get access to more methods, but you have no guarantees about API 

stability

The copy is effectively a deep copy: Updating the state of the copy 
does 
not update the original.

Extra doc: https://github.com/casadi/casadi/wiki/L_1l

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L362

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L733-L735

";

%feature("docstring") casadi::Opti::advanced "

[INTERNAL] 
Get a copy with advanced functionality.

You get access to more methods, but you have no guarantees about API 

stability

The copy is effectively a deep copy: Updating the state of the copy 
does 
not update the original.

Extra doc: https://github.com/casadi/casadi/wiki/L_1m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L372

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L736-L738

";

%feature("docstring") casadi::Opti::copy "

[INTERNAL] 
Get a copy of the.

The copy is effectively a deep copy: Updating the state of the copy 
does 
not update the original.

Extra doc: https://github.com/casadi/casadi/wiki/L_1n

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L380

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L739-L741

";

%feature("docstring") casadi::Opti::update_user_dict "

[INTERNAL]

>  void casadi::Opti::update_user_dict(const std::vector< MX > &m, const Dict &meta)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::Opti::user_dict "

[INTERNAL] 
Get user data.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L391

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L392-L398

";

%feature("docstring") casadi::Opti::type_name "

[INTERNAL] 
Readable name of the class.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L394

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L394-L394

";

%feature("docstring") casadi::Opti::disp "

[INTERNAL] 
Print representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L397

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L636-L656

";

%feature("docstring") casadi::Opti::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L400

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L658-L662

";

%feature("docstring") casadi::Opti::~Opti "

[INTERNAL] 
Destructor.

Extra doc: https://github.com/casadi/casadi/wiki/L_1q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L418

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L418-L418

";

%feature("docstring") casadi::Opti::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L99-L101

";

%feature("docstring") casadi::Opti::is_null "

[INTERNAL] 
Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L73-L75

";

%feature("docstring") casadi::Opti::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L129-L131

";


// File: classcasadi_1_1OptiAdvanced.xml
%feature("docstring") casadi::OptiAdvanced "

[INTERNAL] C++ includes: optistack.hpp
";

%feature("docstring") casadi::OptiAdvanced::symvar "

[INTERNAL] 
Get symbols present in expression.

Returned vector is ordered according to the order of  variable()/parameter()
 calls used to create the variables

Extra doc: https://github.com/casadi/casadi/wiki/L_1u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L525

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L432-L438

>  std::vector< MX > casadi::OptiAdvanced::symvar(const MX &expr, VariableType type) const
------------------------------------------------------------------------
[INTERNAL] 
Get symbols present in expression.

Returned vector is ordered according to the order of  variable()/parameter()
 calls used to create the variables

Extra doc: https://github.com/casadi/casadi/wiki/L_1u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L525

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L432-L438

";

";

%feature("docstring") casadi::OptiAdvanced::subject_to "

[INTERNAL] 
Clear constraints.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L166

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L102-L108

>  void casadi::Opti::subject_to()
------------------------------------------------------------------------
[INTERNAL] 
Clear constraints.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L166

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L102-L108

";

";

%feature("docstring") casadi::OptiAdvanced::set_initial "

[INTERNAL] 
Set initial guess for decision variables

::

  * opti.set_initial(x, 2)
  * opti.set_initial(10*x(1), 2)
  * 



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L190

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L128-L134

>  void casadi::Opti::set_initial(const std::vector< MX > &assignments)
------------------------------------------------------------------------
[INTERNAL] 
Set initial guess for decision variables

::

  * opti.set_initial(x, 2)
  * opti.set_initial(10*x(1), 2)
  * 



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L190

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L128-L134

";

";

%feature("docstring") casadi::OptiAdvanced::set_value "

[INTERNAL] 
Set value of parameter.

Each parameter must be given a value before 'solve' can be called

Extra doc: https://github.com/casadi/casadi/wiki/L_1d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L145-L151

>  void casadi::Opti::set_value(const std::vector< MX > &assignments)
------------------------------------------------------------------------
[INTERNAL] 
Set value of parameter.

Each parameter must be given a value before 'solve' can be called

Extra doc: https://github.com/casadi/casadi/wiki/L_1d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L145-L151

";

";

%feature("docstring") casadi::OptiAdvanced::value "

[INTERNAL] 
Obtain value of expression at the current value

In regular mode, teh current value is the converged solution In debug 
mode,
 the value can be non-converged

Parameters:
-----------

values: 
Optional assignment expressions (e.g. x==3) to overrule the current
 
value

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L225

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L186-L192

>  DM casadi::Opti::value(const SX &x, const std::vector< MX > &values=std::vector< MX >()) const
------------------------------------------------------------------------
[INTERNAL] 
Obtain value of expression at the current value

In regular mode, teh current value is the converged solution In debug 
mode,
 the value can be non-converged

Parameters:
-----------

values: 
Optional assignment expressions (e.g. x==3) to overrule the current
 
value

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L225

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L186-L192

";

";

%feature("docstring") casadi::OptiAdvanced::to_function "

[INTERNAL] 
Create a CasADi  Function from the  Opti solver.

Parameters:
-----------

name: 
Name of the resulting CasADi  Function

args: 
List of parameters and decision/dual variables (which can be given an
 
initial guess) with the resulting  Function

res: 
List of expressions that will get evaluated at the optimal solution

opts: 
Standard CasADi Funcion options

Extra doc: https://github.com/casadi/casadi/wiki/L_1j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L336

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L341-L361

>  Function casadi::Opti::to_function(const std::string &name, const std::map< std::string, MX > &dict, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Create a CasADi  Function from the  Opti solver.

Parameters:
-----------

name: 
Name of the resulting CasADi  Function

args: 
List of parameters and decision/dual variables (which can be given an
 
initial guess) with the resulting  Function

res: 
List of expressions that will get evaluated at the optimal solution

opts: 
Standard CasADi Funcion options

Extra doc: https://github.com/casadi/casadi/wiki/L_1j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L336

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L341-L361

";

";

%feature("docstring") casadi::OptiAdvanced::callback_class "

[INTERNAL] 
Helper methods for callback()

Do not use directly.

Extra doc: https://github.com/casadi/casadi/wiki/L_1p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L409

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L371-L377

>  void casadi::Opti::callback_class()
------------------------------------------------------------------------
[INTERNAL] 
Helper methods for callback()

Do not use directly.

Extra doc: https://github.com/casadi/casadi/wiki/L_1p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L409

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L371-L377

";

";

%feature("docstring") casadi::OptiAdvanced::OptiAdvanced "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::~OptiAdvanced "

[INTERNAL] 
Destructor.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L507

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L507-L507

";

%feature("docstring") casadi::OptiAdvanced::casadi_solver "

[INTERNAL] 
Get the underlying CasADi solver of the  Opti stack.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L511

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L400-L406

";

%feature("docstring") casadi::OptiAdvanced::is_parametric "

[INTERNAL] 
return true if expression is only dependant on  Opti parameters,
 not variables

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L514

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L408-L414

";

%feature("docstring") casadi::OptiAdvanced::canon_expr "

[INTERNAL] 
Interpret an expression (for internal use only)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L529

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L440-L446

";

%feature("docstring") casadi::OptiAdvanced::get_meta "

[INTERNAL] 
Get meta-data of symbol (for internal use only)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L532

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L448-L454

";

%feature("docstring") casadi::OptiAdvanced::get_meta_con "

[INTERNAL] 
Get meta-data of symbol (for internal use only)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L535

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L456-L462

";

%feature("docstring") casadi::OptiAdvanced::set_meta "

[INTERNAL] 
Set meta-data of an expression.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L538

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L464-L470

";

%feature("docstring") casadi::OptiAdvanced::set_meta_con "

[INTERNAL] 
Set meta-data of an expression.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L541

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L472-L478

";

%feature("docstring") casadi::OptiAdvanced::assert_active_symbol "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::active_symvar "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::active_values "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::x_lookup "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::g_lookup "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::x_describe "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::g_describe "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::describe "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::show_infeasibilities "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::solve_prepare "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::solve_actual "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::arg "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::res "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::constraints "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::objective "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::baked_copy "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::assert_empty "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::bake "

[INTERNAL] 
Fix the structure of the optimization problem.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L572

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L664-L670

";

%feature("docstring") casadi::OptiAdvanced::mark_problem_dirty "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::problem_dirty "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::mark_solver_dirty "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::solver_dirty "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::mark_solved "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::solved "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::assert_solved "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::assert_baked "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::instance_number "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::variable "

[INTERNAL] 
Create a decision variable (symbol)

The order of creation matters. The order will be reflected in the 

optimization problem. It is not required for decision variables to 
actualy 
appear in the optimization problem.

Parameters:
-----------

n: 
number of rows (default 1)

m: 
number of columnss (default 1)

attribute: 
'full' (default) or 'symmetric'

Extra doc: https://github.com/casadi/casadi/wiki/L_18

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L112

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L51-L57

";

%feature("docstring") casadi::OptiAdvanced::parameter "

[INTERNAL] 
Create a parameter (symbol); fixed during optimization.

The order of creation does not matter. It is not required for 
parameter to 
actualy appear in the optimization problem. Parameters 
that do appear, must
 be given a value before the problem can be 
solved.

Parameters:
-----------

n: 
number of rows (default 1)

m: 
number of columnss (default 1)

attribute: 
'full' (default) or 'symmetric'

Extra doc: https://github.com/casadi/casadi/wiki/L_19

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L125

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L74-L80

";

%feature("docstring") casadi::OptiAdvanced::minimize "

[INTERNAL] 
Set objective.

Objective must be a scalar. Default objective: 0 When method is called
 
multiple times, the last call takes effect

Extra doc: https://github.com/casadi/casadi/wiki/L_1a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L133

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L82-L88

";

%feature("docstring") casadi::OptiAdvanced::solver "

[INTERNAL] 
Set a solver.

Parameters:
-----------

solver: 
any of the nlpsol plugins can be used here In practice, not all 
nlpsol
 plugins may be supported yet

options: 
passed on to nlpsol plugin No stability can be guaranteed about 
this 
part of the API

options: 
to be passed to nlpsol solver No stability can be guaranteed about
 
this part of the API

Extra doc: https://github.com/casadi/casadi/wiki/L_1c

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L178

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L111-L119

";

%feature("docstring") casadi::OptiAdvanced::solve "

[INTERNAL] 
Crunch the numbers; solve the problem.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L204

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L153-L159

";

%feature("docstring") casadi::OptiAdvanced::solve_limited "

[INTERNAL] 
Crunch the numbers; solve the problem.

Allows the solver to return without error when an iteration or time 
limit 
is reached

Extra doc: https://github.com/casadi/casadi/wiki/L_1e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L212

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L161-L167

";

%feature("docstring") casadi::OptiAdvanced::stats "

[INTERNAL] 
Get statistics.

nlpsol stats are passed as-is. No stability can be guaranteed about 
this 
part of the API

Extra doc: https://github.com/casadi/casadi/wiki/L_1f

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L234

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L194-L200

";

%feature("docstring") casadi::OptiAdvanced::return_status "

[INTERNAL] 
Get return status of solver.



::

     passed as-is from nlpsol
  

No stability can be guaranteed about this part of the API

Extra doc: https://github.com/casadi/casadi/wiki/L_1g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L202-L208

";

%feature("docstring") casadi::OptiAdvanced::initial "

[INTERNAL] 
get assignment expressions for initial values

Extra doc: https://github.com/casadi/casadi/wiki/L_266

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L247

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L210-L216

";

%feature("docstring") casadi::OptiAdvanced::value_variables "

[INTERNAL] 
get assignment expressions for latest values

Extra doc: https://github.com/casadi/casadi/wiki/L_267

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L252

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L218-L224

";

%feature("docstring") casadi::OptiAdvanced::value_parameters "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::dual "

[INTERNAL] 
get the dual variable

m must be a constraint expression. The returned value is still a 
symbolic 
expression. Use  value on it to obtain the numerical value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L262

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L234-L240

";

%feature("docstring") casadi::OptiAdvanced::nx "

[INTERNAL] 
Number of (scalarised) decision variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_268

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L267

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L242-L248

";

%feature("docstring") casadi::OptiAdvanced::np "

[INTERNAL] 
Number of (scalarised) parameters.

Extra doc: https://github.com/casadi/casadi/wiki/L_269

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L272

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L250-L256

";

%feature("docstring") casadi::OptiAdvanced::ng "

[INTERNAL] 
Number of (scalarised) constraints.

Extra doc: https://github.com/casadi/casadi/wiki/L_26a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L277

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L258-L264

";

%feature("docstring") casadi::OptiAdvanced::x "

[INTERNAL] 
Get all (scalarised) decision variables as a symbolic column 
vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_26b

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L282

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L266-L272

";

%feature("docstring") casadi::OptiAdvanced::p "

[INTERNAL] 
Get all (scalarised) parameters as a symbolic column vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_26c

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L287

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L274-L280

";

%feature("docstring") casadi::OptiAdvanced::g "

[INTERNAL] 
Get all (scalarised) constraint expressions as a column vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_26d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L292

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L282-L288

";

%feature("docstring") casadi::OptiAdvanced::f "

[INTERNAL] 
Get objective expression.

Extra doc: https://github.com/casadi/casadi/wiki/L_26e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L297

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L290-L296

";

%feature("docstring") casadi::OptiAdvanced::lbg "

[INTERNAL] 
Get all (scalarised) bounds on constraints as a column vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_26f

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L302

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L298-L304

";

%feature("docstring") casadi::OptiAdvanced::ubg "

[INTERNAL] ";

%feature("docstring") casadi::OptiAdvanced::lam_g "

[INTERNAL] 
Get all (scalarised) dual variables as a symbolic column vector.

Useful for obtaining the Lagrange Hessian:

::

  * sol.value(hessian(opti.f+opti.lam_g'*opti.g,opti.x)) % MATLAB
  * sol.value(hessian(opti.f+dot(opti.lam_g,opti.g),opti.x)[0]) # Python
  * 



Extra doc: https://github.com/casadi/casadi/wiki/L_1i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L314

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L315-L321

";

%feature("docstring") casadi::OptiAdvanced::debug "

[INTERNAL] 
Get a copy with advanced functionality.

You get access to more methods, but you have no guarantees about API 

stability

The copy is effectively a deep copy: Updating the state of the copy 
does 
not update the original.

Extra doc: https://github.com/casadi/casadi/wiki/L_1l

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L362

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L733-L735

";

%feature("docstring") casadi::OptiAdvanced::advanced "

[INTERNAL] 
Get a copy with advanced functionality.

You get access to more methods, but you have no guarantees about API 

stability

The copy is effectively a deep copy: Updating the state of the copy 
does 
not update the original.

Extra doc: https://github.com/casadi/casadi/wiki/L_1m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L372

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L736-L738

";

%feature("docstring") casadi::OptiAdvanced::copy "

[INTERNAL] 
Get a copy of the.

The copy is effectively a deep copy: Updating the state of the copy 
does 
not update the original.

Extra doc: https://github.com/casadi/casadi/wiki/L_1n

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L380

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L739-L741

";

%feature("docstring") casadi::OptiAdvanced::update_user_dict "

[INTERNAL]

>  void casadi::Opti::update_user_dict(const std::vector< MX > &m, const Dict &meta)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::OptiAdvanced::user_dict "

[INTERNAL] 
Get user data.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L391

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L392-L398

";

%feature("docstring") casadi::OptiAdvanced::type_name "

[INTERNAL] 
Readable name of the class.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L394

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L394-L394

";

%feature("docstring") casadi::OptiAdvanced::disp "

[INTERNAL] 
Print representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L397

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L636-L656

";

%feature("docstring") casadi::OptiAdvanced::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L400

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L658-L662

";

%feature("docstring") casadi::OptiAdvanced::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L99-L101

";

%feature("docstring") casadi::OptiAdvanced::is_null "

[INTERNAL] 
Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L73-L75

";

%feature("docstring") casadi::OptiAdvanced::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L129-L131

";


// File: classcasadi_1_1OptiCallback.xml
%feature("docstring") casadi::OptiCallback "

[INTERNAL] C++ includes: optistack.hpp
";

%feature("docstring") casadi::OptiCallback::OptiCallback "

[INTERNAL] ";

%feature("docstring") casadi::OptiCallback::call "

[INTERNAL] ";

%feature("docstring") casadi::OptiCallback::~OptiCallback "

[INTERNAL] ";


// File: classcasadi_1_1OptiSol.xml
%feature("docstring") casadi::OptiSol "

[INTERNAL] 
A simplified interface for NLP modeling/solving.

This class offers a view with solution retrieval facilities The API is
 
guaranteed to be stable.

Joris Gillis, Erik Lambrechts

Extra doc: https://github.com/casadi/casadi/wiki/L_1v

C++ includes: optistack.hpp
";

%feature("docstring") casadi::OptiSol::value "

[INTERNAL] 
Obtain value of expression at the current value

In regular mode, teh current value is the converged solution In debug 
mode,
 the value can be non-converged

Parameters:
-----------

values: 
Optional assignment expressions (e.g. x==3) to overrule the current
 
value

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L622

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L760-L762

>  DM casadi::OptiSol::value(const SX &x, const std::vector< MX > &values=std::vector< MX >()) const
------------------------------------------------------------------------
[INTERNAL] 
Obtain value of expression at the current value

In regular mode, teh current value is the converged solution In debug 
mode,
 the value can be non-converged

Parameters:
-----------

values: 
Optional assignment expressions (e.g. x==3) to overrule the current
 
value

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L622

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L760-L762

";

";

%feature("docstring") casadi::OptiSol::type_name "

[INTERNAL] ";

%feature("docstring") casadi::OptiSol::disp "

[INTERNAL] ";

%feature("docstring") casadi::OptiSol::get_str "

[INTERNAL] ";

%feature("docstring") casadi::OptiSol::value_variables "

[INTERNAL] 
get assignment expressions for the optimal solution

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L626

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L764-L766

";

%feature("docstring") casadi::OptiSol::value_parameters "

[INTERNAL] ";

%feature("docstring") casadi::OptiSol::stats "

[INTERNAL] 
Get statistics.

nlpsol stats are passed as-is. No stability can be guaranteed about 
this 
part of the API

Extra doc: https://github.com/casadi/casadi/wiki/L_1w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L635

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L772-L774

";

%feature("docstring") casadi::OptiSol::opti "

[INTERNAL] ";


// File: classcasadi_1_1Output.xml


// File: classcasadi_1_1Polynomial.xml
%feature("docstring") casadi::Polynomial "

[INTERNAL] 
Helper class for differentiating and integrating polynomials.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_8y

C++ includes: polynomial.hpp
";

%feature("docstring") casadi::Polynomial::Polynomial "

[INTERNAL] 
Construct from a vector of polynomial coefficients.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L56

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L56-L56

>  casadi::Polynomial::Polynomial(const std::vector< T > &coeff)
------------------------------------------------------------------------
[INTERNAL] 
Construct from a vector of polynomial coefficients.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L56

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L56-L56

";

";

%feature("docstring") casadi::Polynomial::coeff "

[INTERNAL] 
Coefficients of the polynomial.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L71

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L71-L71

";

%feature("docstring") casadi::Polynomial::degree "

[INTERNAL] 
Degree of the polynomial.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L74

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.cpp#L72-L74

";

%feature("docstring") casadi::Polynomial::scalar "

[INTERNAL] 
Get scalar value (error if  degree()!=0)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L77

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.cpp#L76-L79

";

%feature("docstring") casadi::Polynomial::derivative "

[INTERNAL] 
Create a new polynomial for the derivative.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L80

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.cpp#L141-L147

";

%feature("docstring") casadi::Polynomial::anti_derivative "

[INTERNAL] 
Create a new polynomial for the anti-derivative (primitive 
function)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L83

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.cpp#L149-L156

";

%feature("docstring") casadi::Polynomial::trim "

[INTERNAL] 
Remove excess zeros.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L86

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.cpp#L133-L139

";

%feature("docstring") casadi::Polynomial::type_name "

[INTERNAL] 
Readable name of the class.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L89

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L89-L89

";

%feature("docstring") casadi::Polynomial::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L92

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.cpp#L56-L70

";


// File: classcasadi_1_1Printable.xml
%feature("docstring") casadi::Printable "

[INTERNAL] 
Base class for objects that have a natural string 
representation.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_7i

C++ includes: printable.hpp
";

%feature("docstring") casadi::Printable::repr "

[INTERNAL] 
Get string representation with type information.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/printable.hpp#L64

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/printable.hpp#L64-L66

";


// File: classcasadi_1_1QpToNlp.xml
%feature("docstring") casadi::QpToNlp "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Qrqp.xml
%feature("docstring") casadi::Qrqp "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Qrsqp.xml
%feature("docstring") casadi::Qrsqp "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Rootfinder.xml
%feature("docstring") casadi::Rootfinder "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1scoped__checkout.xml
%feature("docstring") casadi::scoped_checkout "

[INTERNAL] C++ includes: casadi_misc.hpp
";

%feature("docstring") casadi::scoped_checkout::scoped_checkout "

[INTERNAL] ";

%feature("docstring") casadi::scoped_checkout::~scoped_checkout "

[INTERNAL] ";


// File: classcasadi_1_1Scpgen.xml
%feature("docstring") casadi::Scpgen "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1SerializerBase.xml
%feature("docstring") casadi::SerializerBase "

[INTERNAL] C++ includes: serializer.hpp
";

%feature("docstring") casadi::SerializerBase::SerializerBase "

[INTERNAL] ";

%feature("docstring") casadi::SerializerBase::~SerializerBase "

[INTERNAL] ";

%feature("docstring") casadi::SerializerBase::pack "

[INTERNAL] ";

%feature("docstring") casadi::SerializerBase::connect "

[INTERNAL] ";

%feature("docstring") casadi::SerializerBase::reset "

[INTERNAL] ";


// File: classcasadi_1_1SerializingStream.xml
%feature("docstring") casadi::SerializingStream "

[INTERNAL] 
Helper class for Serialization.

Joris Gillis

Extra doc: https://github.com/casadi/casadi/wiki/L_ao

C++ includes: serializing_stream.hpp
";

%feature("docstring") casadi::SerializingStream::SerializingStream "

[INTERNAL]

>  casadi::SerializingStream::SerializingStream(std::ostream &out, const Dict &opts)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::SerializingStream::pack "

[INTERNAL]

>  void casadi::SerializingStream::pack(const MX &e)

>  void casadi::SerializingStream::pack(const SXElem &e)

>  void casadi::SerializingStream::pack(const Linsol &e)

>  void casadi::SerializingStream::pack(const Matrix< T > &e)

>  void casadi::SerializingStream::pack(const Function &e)

>  void casadi::SerializingStream::pack(const Importer &e)

>  void casadi::SerializingStream::pack(const Slice &e)

>  void casadi::SerializingStream::pack(const GenericType &e)

>  void casadi::SerializingStream::pack(std::istream &s)

>  void casadi::SerializingStream::pack(int e)

>  void casadi::SerializingStream::pack(bool e)

>  void casadi::SerializingStream::pack(casadi_int e)

>  void casadi::SerializingStream::pack(size_t e)

>  void casadi::SerializingStream::pack(double e)

>  void casadi::SerializingStream::pack(const std::string &e)

>  void casadi::SerializingStream::pack(char e)

>  void casadi::SerializingStream::pack(const std::vector< T > &e)

>  void casadi::SerializingStream::pack(const std::map< K, V > &e)

>  void casadi::SerializingStream::pack(const std::pair< A, B > &e)

>  void casadi::SerializingStream::pack(const std::string &descr, const T &e)

>  void casadi::SerializingStream::pack(const std::string &descr, T &e)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::SerializingStream::version "

[INTERNAL] ";

%feature("docstring") casadi::SerializingStream::connect "

[INTERNAL] ";

%feature("docstring") casadi::SerializingStream::reset "

[INTERNAL] ";


// File: classcasadi_1_1SharedObject.xml
%feature("docstring") casadi::SharedObject "

[INTERNAL] 
 SharedObject implements a reference counting framework similar 
for efficient and.

easily-maintained memory management.

To use the class, both the  SharedObject class (the public class), and the 
SharedObjectInternal class (the 
internal class) must be inherited from. It 
can be done in two 
different files and together with memory management, 
this approach 
provides a clear distinction of which methods of the class 
are to be 
considered \"public\", i.e. methods for public use that can be 

considered to remain over time with small changes, and the internal 
memory.

When interfacing a software, which typically includes including some 
header
 file, this is best done only in the file where the internal 
class is 
defined, to avoid polluting the global namespace and other 
side effects.

The default constructor always means creating a null pointer to an 
internal
 class only. To allocate an internal class (this works only 
when the 
internal class isn't abstract), use the constructor with 
arguments.

The copy constructor and the assignment operator perform shallow 
copies 
only, to make a deep copy you must use the clone method 
explicitly. This 
will give a shared pointer instance.

In an inheritance hierarchy, you can cast down automatically, e.g. 

(SXFunction is a child class of  Function): SXFunction derived(...);  
Function base = derived;

To cast up, use the shared_cast template function, which works 
analogously 
to dynamic_cast, static_cast, const_cast etc, e.g.: 
SXFunction 
derived(...);  Function base = derived; SXFunction derived_from_base = 

shared_cast<SXFunction>(base);

A failed shared_cast will result in a null pointer (cf. dynamic_cast)

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_as

C++ includes: shared_object.hpp
";

%feature("docstring") casadi::SharedObject::SharedObject "

[INTERNAL] 
Copy constructor (shallow copy)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L96

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L37-L40

>  casadi::SharedObject::SharedObject(const SharedObject &ref)
------------------------------------------------------------------------
[INTERNAL] 
Copy constructor (shallow copy)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L96

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L37-L40

";

";

%feature("docstring") casadi::SharedObject::~SharedObject "

[INTERNAL] 
Destructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L99

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L42-L44

";

%feature("docstring") casadi::SharedObject::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L99-L101

";

%feature("docstring") casadi::SharedObject::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L103-L109

";

%feature("docstring") casadi::SharedObject::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::SharedObject::is_null "

[INTERNAL] 
Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L73-L75

";

%feature("docstring") casadi::SharedObject::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L129-L131

";


// File: classcasadi_1_1ShellCompiler.xml
%feature("docstring") casadi::ShellCompiler "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Slice.xml
%feature("docstring") casadi::Slice "

[INTERNAL] 
Class representing a  Slice.

Note that Python or Octave do not need to use this class. They can 
just use
 slicing utility from the host language ( M[0:6] in Python, 
M(1:7) )

Extra doc: https://github.com/casadi/casadi/wiki/L_13

C++ includes: slice.hpp
";

%feature("docstring") casadi::Slice::Slice "

[INTERNAL]

>  casadi::Slice::Slice(int start, int stop, int step=1)

>  casadi::Slice::Slice(int start, casadi_int stop, int step=1)

>  casadi::Slice::Slice(casadi_int start, int stop, int step=1)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::Slice::all "

[INTERNAL] 
Get a vector of indices (nested slice)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L75

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L129-L137

>  std::vector< casadi_int > casadi::Slice::all(const Slice &outer, casadi_int len) const
------------------------------------------------------------------------
[INTERNAL] 
Get a vector of indices (nested slice)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L75

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L129-L137

";

";

%feature("docstring") casadi::Slice::size "

[INTERNAL] 
Get number of elements.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L78

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L104-L109

";

%feature("docstring") casadi::Slice::is_empty "

[INTERNAL] 
Check if slice is empty.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L81

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L111-L113

";

%feature("docstring") casadi::Slice::is_scalar "

[INTERNAL] 
Is the slice a scalar.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L84

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L139-L144

";

%feature("docstring") casadi::Slice::scalar "

[INTERNAL] 
Get scalar (if is_scalar)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L87

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L146-L150

";

%feature("docstring") casadi::Slice::apply "

[INTERNAL] 
Apply concrete length.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L98

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L66-L88

";

%feature("docstring") casadi::Slice::type_name "

[INTERNAL] 
Get name of the class.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L107

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L107-L107

";

%feature("docstring") casadi::Slice::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L110

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L115-L127

";

%feature("docstring") casadi::Slice::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L113

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L113-L117

";

%feature("docstring") casadi::Slice::info "

[INTERNAL] 
Obtain information

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L120

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L120-L122

";

%feature("docstring") casadi::Slice::serialize "

[INTERNAL] 
Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_14

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L127

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L292-L296

";


// File: classcasadi_1_1SlicotDple.xml
%feature("docstring") casadi::SlicotDple "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Solve.xml


// File: classcasadi_1_1SolveUnity.xml


// File: classcasadi_1_1SparsityInterface.xml
%feature("docstring") casadi::SparsityInterfaceCommon "

[INTERNAL] 
Sparsity interface class.

This is a common base class for  GenericMatrix (i.e.  MX and Matrix<>) and 
Sparsity, introducing a uniform syntax and 
implementing common 
functionality using the curiously recurring 
template pattern (CRTP) idiom.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_3d

C++ includes: sparsity_interface.hpp
";

%feature("docstring") casadi::SparsityInterfaceCommon::horzcat "

[INTERNAL] 
Concatenate horizontally, six matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_4e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L502

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L502-L505

>  MatType casadi::SparsityInterface::horzcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w, const MatType &v, const MatType &u)
------------------------------------------------------------------------
[INTERNAL] 
Concatenate horizontally, six matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_4e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L502

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L502-L505

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::vertcat "

[INTERNAL] 
Concatenate vertically, six matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_4j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L540

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L540-L543

>  MatType casadi::SparsityInterface::vertcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w, const MatType &v, const MatType &u)
------------------------------------------------------------------------
[INTERNAL] 
Concatenate vertically, six matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_4j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L540

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L540-L543

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::horzsplit "

[INTERNAL] 
split horizontally, retaining fixed-sized groups of columns

Parameters:
-----------

incr: 
Size (width) of each group of columns

horzcat(horzsplit(x, ...)) = x

\\\\seealso horzsplit_n

Extra doc: https://github.com/casadi/casadi/wiki/L_3h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L134

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L134-L136

>  std::vector<MatType > casadi::SparsityInterface::horzsplit(const MatType &x, casadi_int incr=1)
------------------------------------------------------------------------
[INTERNAL] 
split horizontally, retaining fixed-sized groups of columns

Parameters:
-----------

incr: 
Size (width) of each group of columns

horzcat(horzsplit(x, ...)) = x

\\\\seealso horzsplit_n

Extra doc: https://github.com/casadi/casadi/wiki/L_3h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L134

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L134-L136

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::horzsplit_n "

[INTERNAL] 
split horizontally, retaining fixed-sized groups of columns

Parameters:
-----------

n: 
Number of groups of columns

Will error when the number of columns is not a multiple of n

horzcat(horzsplit(x, ...)) = x

\\\\seealso horzsplit

Extra doc: https://github.com/casadi/casadi/wiki/L_277

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L149

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L149-L151

";

%feature("docstring") casadi::SparsityInterfaceCommon::vertsplit "

[INTERNAL] 
split vertically, retaining fixed-sized groups of rows

Parameters:
-----------

incr: 
Size of each group of rows

vertcat(vertsplit(x, ...)) = x



::

  >>> print vertsplit(SX.sym(\"a\",4))
  [SX(a_0), SX(a_1), SX(a_2), SX(a_3)]
  





::

  >>> print vertsplit(SX.sym(\"a\",4),2)
  [SX([a_0, a_1]), SX([a_2, a_3])]
  



If the number of rows is not a multiple of  incr, the last entry returned 
will have a size smaller than  incr.



::

  >>> print vertsplit(DM([0,1,2,3,4]),2)
  [DM([0, 1]), DM([2, 3]), DM(4)]
  



\\\\seealso vertsplit_n

Extra doc: https://github.com/casadi/casadi/wiki/L_3k

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L204

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L204-L206

>  std::vector<MatType > casadi::SparsityInterface::vertsplit(const MatType &x, casadi_int incr=1)
------------------------------------------------------------------------
[INTERNAL] 
split vertically, retaining fixed-sized groups of rows

Parameters:
-----------

incr: 
Size of each group of rows

vertcat(vertsplit(x, ...)) = x



::

  >>> print vertsplit(SX.sym(\"a\",4))
  [SX(a_0), SX(a_1), SX(a_2), SX(a_3)]
  





::

  >>> print vertsplit(SX.sym(\"a\",4),2)
  [SX([a_0, a_1]), SX([a_2, a_3])]
  



If the number of rows is not a multiple of  incr, the last entry returned 
will have a size smaller than  incr.



::

  >>> print vertsplit(DM([0,1,2,3,4]),2)
  [DM([0, 1]), DM([2, 3]), DM(4)]
  



\\\\seealso vertsplit_n

Extra doc: https://github.com/casadi/casadi/wiki/L_3k

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L204

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L204-L206

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::offset "

[INTERNAL] 
Helper function, get offsets corresponding to a vector of 
matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_3j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L169

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L169-L171

";

%feature("docstring") casadi::SparsityInterfaceCommon::vertsplit_n "

[INTERNAL] 
split vertically, retaining fixed-sized groups of rows

Parameters:
-----------

n: 
Number of groups of rows

Will error when the number of rows is not a multiple of n

vertcat(vertsplit(x, ...)) = x

\\\\seealso vertsplit

Extra doc: https://github.com/casadi/casadi/wiki/L_278

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L219

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L219-L221

";

%feature("docstring") casadi::SparsityInterfaceCommon::blockcat "

[INTERNAL] 
Construct a matrix from 4 blocks.

Extra doc: https://github.com/casadi/casadi/wiki/L_3m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L234

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L234-L236

>  MatType casadi::SparsityInterface::blockcat(const MatType &A, const MatType &B, const MatType &C, const MatType &D)
------------------------------------------------------------------------
[INTERNAL] 
Construct a matrix from 4 blocks.

Extra doc: https://github.com/casadi/casadi/wiki/L_3m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L234

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L234-L236

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::blocksplit "

[INTERNAL] 
chop up into blocks

Parameters:
-----------

vert_incr: 
Defines the increment for block boundaries in row dimension

horz_incr: 
Defines the increment for block boundaries in column dimension

blockcat(blocksplit(x,..., ...)) = x

Extra doc: https://github.com/casadi/casadi/wiki/L_3o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L262

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L262-L264

>  std::vector< std::vector< MatType > > casadi::SparsityInterface::blocksplit(const MatType &x, casadi_int vert_incr=1, casadi_int horz_incr=1)
------------------------------------------------------------------------
[INTERNAL] 
chop up into blocks

Parameters:
-----------

vert_incr: 
Defines the increment for block boundaries in row dimension

horz_incr: 
Defines the increment for block boundaries in column dimension

blockcat(blocksplit(x,..., ...)) = x

Extra doc: https://github.com/casadi/casadi/wiki/L_3o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L262

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L262-L264

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::diagcat "

[INTERNAL] 
Concatenate along diagonal, six matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_4o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L578

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L578-L581

>  MatType casadi::SparsityInterface::diagcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w, const MatType &v, const MatType &u)
------------------------------------------------------------------------
[INTERNAL] 
Concatenate along diagonal, six matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_4o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L578

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L578-L581

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::diagsplit "

[INTERNAL] 
split diagonally, retaining fixed-sized matrices

Parameters:
-----------

incr1: 
Row dimension of each matrix

incr2: 
Column dimension of each matrix

diagsplit(diagsplit(x, ...)) = x

Extra doc: https://github.com/casadi/casadi/wiki/L_3t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L324

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L324-L326

>  std::vector< MatType > casadi::SparsityInterface::diagsplit(const MatType &x, casadi_int incr1, casadi_int incr2)
------------------------------------------------------------------------
[INTERNAL] 
split diagonally, retaining fixed-sized matrices

Parameters:
-----------

incr1: 
Row dimension of each matrix

incr2: 
Column dimension of each matrix

diagsplit(diagsplit(x, ...)) = x

Extra doc: https://github.com/casadi/casadi/wiki/L_3t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L324

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L324-L326

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::veccat "

[INTERNAL] 
concatenate vertically while vectorizing all arguments with vec

Extra doc: https://github.com/casadi/casadi/wiki/L_3u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L331

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L331-L333

";

%feature("docstring") casadi::SparsityInterfaceCommon::mtimes "

[INTERNAL] 
 Matrix product of n matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_3w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L345

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L345-L347

>  MatType casadi::SparsityInterface::mtimes(const std::vector< MatType > &args)
------------------------------------------------------------------------
[INTERNAL] 
 Matrix product of n matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_3w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L345

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L345-L347

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::mac "

[INTERNAL] 
Multiply-accumulate operation.

Matrix product of two matrices (x and y), adding the result to a third 

matrix z. The result has the same sparsity pattern as C meaning that 
other 
entries of (x*y) are ignored. The operation is equivalent to: 

z+mtimes(x,y).project(z.sparsity()).

Extra doc: https://github.com/casadi/casadi/wiki/L_3x

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L358

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L358-L360

";

%feature("docstring") casadi::SparsityInterfaceCommon::transpose "

[INTERNAL] 
Transpose.

Extra doc: https://github.com/casadi/casadi/wiki/L_3y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L365

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L365-L367

";

%feature("docstring") casadi::SparsityInterfaceCommon::vec "

[INTERNAL] 
make a vector

Reshapes/vectorizes the matrix such that the shape becomes 
(expr.numel(), 
1). Columns are stacked on top of each other. Same as 
reshape(expr, 
expr.numel(), 1)

a c 
b d 
 turns into

a 
b 
c 
d 
 Extra doc: https://github.com/casadi/casadi/wiki/L_3z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L386

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L386-L388

";

%feature("docstring") casadi::SparsityInterfaceCommon::reshape "

[INTERNAL] 
Reshape the matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_42

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L407

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L407-L409

>  MatType casadi::SparsityInterface::reshape(const MatType &x, const Sparsity &sp)
------------------------------------------------------------------------
[INTERNAL] 
Reshape the matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_42

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L407

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L407-L409

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::sparsity_cast "

[INTERNAL] 
Cast matrix nonzeros to different Sparsity.

Extra doc: https://github.com/casadi/casadi/wiki/L_24z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L414

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L414-L416

";

%feature("docstring") casadi::SparsityInterfaceCommon::sprank "

[INTERNAL] 
Obtain the structural rank of a sparsity-pattern.

Extra doc: https://github.com/casadi/casadi/wiki/L_43

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L421

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L421-L423

";

%feature("docstring") casadi::SparsityInterfaceCommon::norm_0_mul "

[INTERNAL] 
0-norm (nonzero count) of a Matrix-matrix product

Extra doc: https://github.com/casadi/casadi/wiki/L_44

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L428

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L428-L430

";

%feature("docstring") casadi::SparsityInterfaceCommon::triu "

[INTERNAL] 
Get the upper triangular part of a matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_45

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L435

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L435-L437

";

%feature("docstring") casadi::SparsityInterfaceCommon::tril "

[INTERNAL] 
Get the lower triangular part of a matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_46

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L442

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L442-L444

";

%feature("docstring") casadi::SparsityInterfaceCommon::kron "

[INTERNAL] 
Kronecker tensor product.

Creates a block matrix in which each element (i, j) is a_ij*b

Extra doc: https://github.com/casadi/casadi/wiki/L_47

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L451

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L451-L453

";

%feature("docstring") casadi::SparsityInterfaceCommon::repmat "

[INTERNAL] 
Repeat matrix A n times vertically and m times horizontally.

Extra doc: https://github.com/casadi/casadi/wiki/L_49

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L465

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L465-L467

>  MatType casadi::SparsityInterface::repmat(const MatType &A, const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
Repeat matrix A n times vertically and m times horizontally.

Extra doc: https://github.com/casadi/casadi/wiki/L_49

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L465

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L465-L467

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::sum1 "

[INTERNAL] 
Return a row-wise summation of elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_4p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L586

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L586-L586

";

%feature("docstring") casadi::SparsityInterfaceCommon::sum2 "

[INTERNAL] 
Return a column-wise summation of elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_4q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L591

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L591-L591

";


// File: classcasadi_1_1Sqpmethod.xml
%feature("docstring") casadi::Sqpmethod "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Logger_1_1Stream.xml
%feature("docstring") casadi::Logger::Stream "

[INTERNAL] C++ includes: casadi_logger.hpp
";

%feature("docstring") casadi::Logger::Stream::Stream "

[INTERNAL] ";


// File: classcasadi_1_1Logger_1_1Streambuf.xml
%feature("docstring") casadi::Logger::Streambuf "

[INTERNAL] C++ includes: casadi_logger.hpp
";

%feature("docstring") casadi::Logger::Streambuf::Streambuf "

[INTERNAL] ";


// File: classcasadi_1_1StringDeserializer.xml
%feature("docstring") casadi::StringDeserializer "

[INTERNAL] C++ includes: serializer.hpp
";

%feature("docstring") casadi::StringDeserializer::StringDeserializer "

[INTERNAL] 
Advanced deserialization of CasADi objects.

See: 
 StringDeserializer

Extra doc: https://github.com/casadi/casadi/wiki/L_7r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.hpp#L233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.cpp#L114-L117

";

%feature("docstring") casadi::StringDeserializer::~StringDeserializer "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::decode "

[INTERNAL] 
Sets the string to deserialize objects from.

Extra doc: https://github.com/casadi/casadi/wiki/L_7s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.hpp#L240

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.cpp#L91-L96

";

%feature("docstring") casadi::StringDeserializer::pop_type "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_sparsity "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_mx "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_dm "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_sx "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_linsol "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_function "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_generictype "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_int "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_double "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_string "

[INTERNAL] ";

%feature("docstring") 
casadi::StringDeserializer::blind_unpack_sparsity_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_mx_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_dm_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_sx_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_linsol_vector
 "

[INTERNAL] ";

%feature("docstring") 
casadi::StringDeserializer::blind_unpack_function_vector "

[INTERNAL] ";

%feature("docstring") 
casadi::StringDeserializer::blind_unpack_generictype_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_int_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_double_vector
 "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::blind_unpack_string_vector
 "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_sparsity "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_mx "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_dm "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_sx "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_linsol "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_function "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_generictype "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_int "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_double "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_string "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_sparsity_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_mx_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_dm_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_sx_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_linsol_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_function_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_generictype_vector 
"

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_int_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_double_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::unpack_string_vector "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::connect "

[INTERNAL] ";

%feature("docstring") casadi::StringDeserializer::reset "

[INTERNAL] ";


// File: classcasadi_1_1StringSerializer.xml
%feature("docstring") casadi::StringSerializer "

[INTERNAL] C++ includes: serializer.hpp
";

%feature("docstring") casadi::StringSerializer::StringSerializer "

[INTERNAL] 
Advanced serialization of CasADi objects.

This class is intended for advanced users that want to circumvent the 

restrictions of standard pickling/matlab save load, ie no raw SX/MX 
symbols
 allowed.



::

  x = SX.sym('x');
  s = StringSerializer();
  s.pack(x);
  s.pack(sin(x));
   
  data = s.encode();
  
  s = StringDeserializer(data);
  a = s.unpack();
  b = s.unpack();
  



Note: Saving SX/MX objects individually has a substantial overhead 
(both 
time and length of encoded string). You are encouraged to use 
the 
vector/list variants of 'save' for SX/MX to reduce the overhead.

See: 
 Function::save,  Function::serialize,  StringDeserializer,  
FileSerializer

Extra doc: https://github.com/casadi/casadi/wiki/L_7o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.hpp#L203

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.cpp#L37-L39

";

%feature("docstring") casadi::StringSerializer::~StringSerializer "

[INTERNAL] ";

%feature("docstring") casadi::StringSerializer::encode "

[INTERNAL] 
Returns a string that holds the serialized objects.

As a side effect, this method clears the internal buffer

Extra doc: https://github.com/casadi/casadi/wiki/L_7p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.hpp#L211

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.cpp#L85-L90

";

%feature("docstring") casadi::StringSerializer::pack "

[INTERNAL] ";

%feature("docstring") casadi::StringSerializer::connect "

[INTERNAL] ";

%feature("docstring") casadi::StringSerializer::reset "

[INTERNAL] ";


// File: classcasadi_1_1SubIndex.xml
%feature("docstring") casadi::SubIndex "

[INTERNAL] 
 SubIndex class for  Matrix Same as the above class but for 
single argument return for operator()
 
Joel Andersson

C++ includes: submatrix.hpp
";

%feature("docstring") casadi::SubIndex::SubIndex "

[INTERNAL] 
Default copy constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/submatrix.hpp#L113

>  casadi::SubIndex< M, I >::SubIndex(const SubIndex< M, I > &y)=default
------------------------------------------------------------------------
[INTERNAL] 
Default copy constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/submatrix.hpp#L113

";

";


// File: classcasadi_1_1SubMatrix.xml
%feature("docstring") casadi::SubMatrix "

[INTERNAL] 
 SubMatrix class for  Matrix SubMatrix is the return type for 
operator() of the  Matrix class, it allows access to the value as well as 
changing the parent 
object 
Joel Andersson

C++ includes: submatrix.hpp
";

%feature("docstring") casadi::SubMatrix::SubMatrix "

[INTERNAL] 
Default copy constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/submatrix.hpp#L53

>  casadi::SubMatrix< M, I, J >::SubMatrix(const SubMatrix< M, I, J > &y)=default
------------------------------------------------------------------------
[INTERNAL] 
Default copy constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/submatrix.hpp#L53

";

";


// File: classcasadi_1_1SXElem.xml
%feature("docstring") casadi::SXElem "

[INTERNAL] 
The basic scalar symbolic class of CasADi.

SXElem is exposed only as an empty struct to SWIG

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_10i

C++ includes: sx_elem.hpp
";

%feature("docstring") casadi::SXElem::minus "

[INTERNAL] 
Subtraction: (x,y) -> x - y.

Extra doc: https://github.com/casadi/casadi/wiki/L_oo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L83

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L83-L85

";

%feature("docstring") casadi::SXElem::times "

[INTERNAL] 
Elementwise multiplication: (x,y) -> x .* y.

Extra doc: https://github.com/casadi/casadi/wiki/L_op

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L99

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L99-L101

";

%feature("docstring") casadi::SXElem::rdivide "

[INTERNAL] 
Elementwise division: (x,y) -> x ./ y.

Extra doc: https://github.com/casadi/casadi/wiki/L_oq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L115

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L115-L117

";

%feature("docstring") casadi::SXElem::logic_not "

[INTERNAL] 
Logical  not x -> !x.

Returns (an expression evaluating to) 1 if expression is zero and 0 

otherwise

Extra doc: https://github.com/casadi/casadi/wiki/L_oz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L260

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L260-L262

";

%feature("docstring") casadi::SXElem::lt "

[INTERNAL] 
Logical less than: (x,y) -> x < y.

Extra doc: https://github.com/casadi/casadi/wiki/L_or

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L131

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L131-L133

";

%feature("docstring") casadi::SXElem::le "

[INTERNAL] 
Logical less or equal to: (x,y) -> x <= y.

Extra doc: https://github.com/casadi/casadi/wiki/L_os

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L146

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L146-L148

";

%feature("docstring") casadi::SXElem::gt "

[INTERNAL] 
Logical greater than: (x,y) -> x > y.

Extra doc: https://github.com/casadi/casadi/wiki/L_ot

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L161

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L161-L163

";

%feature("docstring") casadi::SXElem::ge "

[INTERNAL] 
Logical greater or equal to: (x,y) -> x <= y.

Extra doc: https://github.com/casadi/casadi/wiki/L_ou

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L176

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L176-L178

";

%feature("docstring") casadi::SXElem::eq "

[INTERNAL] 
Logical equal to: (x,y) -> x == y.

Extra doc: https://github.com/casadi/casadi/wiki/L_ov

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L191

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L191-L193

";

%feature("docstring") casadi::SXElem::ne "

[INTERNAL] 
Logical not equal to: (x,y) -> x != y.

Extra doc: https://github.com/casadi/casadi/wiki/L_ow

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L206

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L206-L208

";

%feature("docstring") casadi::SXElem::logic_and "

[INTERNAL] 
Logical  and

Returns (an expression evaluating to) 1 if both expressions are 
nonzero and
 0 otherwise

Extra doc: https://github.com/casadi/casadi/wiki/L_ox

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L224

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L224-L226

";

%feature("docstring") casadi::SXElem::logic_or "

[INTERNAL] 
Logical  or

returns (an expression evaluating to) 1 if at least one expression is 

nonzero and 0 otherwise

Extra doc: https://github.com/casadi/casadi/wiki/L_oy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L242-L244

";

%feature("docstring") casadi::SXElem::abs "

[INTERNAL] 
Absolute value: x -> abs(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L275

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L275-L277

";

%feature("docstring") casadi::SXElem::sqrt "

[INTERNAL] 
Square root: x -> sqrt(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L290

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L290-L292

";

%feature("docstring") casadi::SXElem::sq "

[INTERNAL] 
Square: x -> x^2.

Extra doc: https://github.com/casadi/casadi/wiki/L_p2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L302

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L302-L304

";

%feature("docstring") casadi::SXElem::sin "

[INTERNAL] 
Sine: x -> sin(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L314

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L314-L316

";

%feature("docstring") casadi::SXElem::cos "

[INTERNAL] 
Cosine: x -> cos(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L326

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L326-L328

";

%feature("docstring") casadi::SXElem::tan "

[INTERNAL] 
Tangent: x -> tan(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L338

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L338-L340

";

%feature("docstring") casadi::SXElem::atan "

[INTERNAL] 
Arc tangent: x -> atan(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L350

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L350-L352

";

%feature("docstring") casadi::SXElem::asin "

[INTERNAL] 
Arc sine: x -> asin(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L362

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L362-L364

";

%feature("docstring") casadi::SXElem::acos "

[INTERNAL] 
Arc cosine: x -> acos(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L374

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L374-L376

";

%feature("docstring") casadi::SXElem::tanh "

[INTERNAL] 
Hyperbolic tangent: x -> tanh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_p9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L386

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L386-L388

";

%feature("docstring") casadi::SXElem::sinh "

[INTERNAL] 
Hyperbolic sin: x -> sinh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L398

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L398-L400

";

%feature("docstring") casadi::SXElem::cosh "

[INTERNAL] 
Hyperbolic cosine: x -> cosh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L410

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L410-L412

";

%feature("docstring") casadi::SXElem::atanh "

[INTERNAL] 
Inverse hyperbolic tangent: x -> atanh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L422

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L422-L424

";

%feature("docstring") casadi::SXElem::asinh "

[INTERNAL] 
Inverse hyperbolic sin: x -> asinh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L434

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L434-L436

";

%feature("docstring") casadi::SXElem::acosh "

[INTERNAL] 
Inverse hyperbolic cosine: x -> acosh(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pe

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L446

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L446-L448

";

%feature("docstring") casadi::SXElem::exp "

[INTERNAL] 
Elementwise exponential: x -> exp(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L458

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L458-L460

";

%feature("docstring") casadi::SXElem::log "

[INTERNAL] 
Natural logarithm: x -> log(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L470

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L470-L472

";

%feature("docstring") casadi::SXElem::log10 "

[INTERNAL] 
Base-10 logarithm: x -> log10(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_ph

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L482

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L482-L484

";

%feature("docstring") casadi::SXElem::log1p "

[INTERNAL] 
Precision variant for natural logarithm: x -> log(x+1)

Extra doc: https://github.com/casadi/casadi/wiki/L_pi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L494

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L494-L496

";

%feature("docstring") casadi::SXElem::expm1 "

[INTERNAL] 
Precision variant for elementwise exponential: x -> exp(x)-1.

Extra doc: https://github.com/casadi/casadi/wiki/L_pj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L506

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L506-L508

";

%feature("docstring") casadi::SXElem::floor "

[INTERNAL] 
Round down to nearest integer: x -> floor(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L518

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L518-L520

";

%feature("docstring") casadi::SXElem::ceil "

[INTERNAL] 
Round up to nearest integer: x -> ceil(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L530

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L530-L532

";

%feature("docstring") casadi::SXElem::erf "

[INTERNAL] 
Error function: x -> erf(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L542

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L542-L544

";

%feature("docstring") casadi::SXElem::erfinv "

[INTERNAL] 
Inverse error function: x -> erfinv(x)

Extra doc: https://github.com/casadi/casadi/wiki/L_pn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L554

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L554-L556

";

%feature("docstring") casadi::SXElem::sign "

[INTERNAL] 
Sign function:

sign(x) := -1 for x<0 sign(x) := 1 for x>0, sign(0) := 0 sign(NaN) := 
NaN

Extra doc: https://github.com/casadi/casadi/wiki/L_po

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L571

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L571-L573

";

%feature("docstring") casadi::SXElem::pow "

[INTERNAL] 
Elementwise power: (x,y) -> x.^y.

Extra doc: https://github.com/casadi/casadi/wiki/L_pp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L583

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L583-L585

";

%feature("docstring") casadi::SXElem::mod "

[INTERNAL] 
Remainder after division: (x,y) -> fmod(x,y)

This Function follows the convention of 
https://en.cppreference.com/w/c/numeric/math/fmod

Notably:
fmod(5,3) -> 2

fmod(5,-3) -> 2

fmod(-5,3) -> -2

fmod(-5,-3) -> -2

This is equivalent to Python's numpy.fmod and Matlab's rem.

\\\\seealso remainder

::

  Extra doc: https://github.com/casadi/casadi/wiki/L_pq 
  



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L607

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L607-L609

";

%feature("docstring") casadi::SXElem::remainder "

[INTERNAL] 
Remainder after division: (x,y) -> remainder(x,y)

This Function follows the convention of 
https://en.cppreference.com/w/c/numeric/math/remainder

Notably:
remainder(5,3) -> -1

remainder(5,-3) -> -1

remainder(-5,3) -> 1

remainder(-5,-3) -> 1

This is equivalent to Python's math.remainder. There is no equivalence
 in 
Matlab.

\\\\seealso fmod

::

  Extra doc: https://github.com/casadi/casadi/wiki/L_24x 
  



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L634

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L634-L636

";

%feature("docstring") casadi::SXElem::atan2 "

[INTERNAL] 
Two argument arc tangent: (y,x) -> atan2(y,x)

theta = atan2(y,x) corresponds to x = r cos(theta), y = r sin(theta)

Extra doc: https://github.com/casadi/casadi/wiki/L_pr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L648

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L648-L650

";

%feature("docstring") casadi::SXElem::if_else_zero "

[INTERNAL] 
Conditional assignment: (x,y) -> x ? y : 0.

Extra doc: https://github.com/casadi/casadi/wiki/L_ps

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L660

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L660-L662

";

%feature("docstring") casadi::SXElem::fmin "

[INTERNAL] 
Smallest of two values: (x,y) -> min(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L672

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L672-L674

";

%feature("docstring") casadi::SXElem::fmax "

[INTERNAL] 
Largest of two values: (x,y) -> max(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L684

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L684-L686

";

%feature("docstring") casadi::SXElem::copysign "

[INTERNAL] 
Copy sign

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L710

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L710-L712

";

%feature("docstring") casadi::SXElem::constpow "

[INTERNAL] 
Elementwise power with const power

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L720

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L720-L722

";

%feature("docstring") casadi::SXElem::printme "

[INTERNAL] 
Debug printing

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L730

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L730-L732

";

%feature("docstring") casadi::SXElem::hypot "

[INTERNAL] 
Precision variant for 2 norm: (x,y) -> sqrt(x^2+y^2)

Extra doc: https://github.com/casadi/casadi/wiki/L_pw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L742

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L742-L744

";

%feature("docstring") casadi::SXElem::if_else "

[INTERNAL] 
Ternary if_else: x ? y : z.

Extra doc: https://github.com/casadi/casadi/wiki/L_113

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L265

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L265-L267

";

%feature("docstring") casadi::SXElem::SXElem "

[INTERNAL] 
Copy constructor.

Extra doc: https://github.com/casadi/casadi/wiki/L_10m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L107

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L59-L62

>  casadi::SXElem::SXElem(const SXElem &scalar)
------------------------------------------------------------------------
[INTERNAL] 
Copy constructor.

Extra doc: https://github.com/casadi/casadi/wiki/L_10m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L107

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L59-L62

";

";

%feature("docstring") casadi::SXElem::~SXElem "

[INTERNAL] 
Destructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L110

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L87-L89

";

%feature("docstring") casadi::SXElem::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L129

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L133-L135

";

%feature("docstring") casadi::SXElem::__nonzero__ "

[INTERNAL] 
Check the truth value of this node.

Introduced to catch bool(x) situations in python

Extra doc: https://github.com/casadi/casadi/wiki/L_10q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L156

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L150-L153

";

%feature("docstring") casadi::SXElem::is_leaf "

[INTERNAL] 
check if this  SXElem is a leaf of the SX graph

An  SXElem qualifies as leaf when it has no dependencies.

Extra doc: https://github.com/casadi/casadi/wiki/L_10r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L163

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L435-L438

";

%feature("docstring") casadi::SXElem::is_constant "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::is_integer "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::is_symbolic "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::is_commutative "

[INTERNAL] 
Check whether a binary  SXElem is commutative.

Extra doc: https://github.com/casadi/casadi/wiki/L_10s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L170

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L440-L443

";

%feature("docstring") casadi::SXElem::is_zero "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::is_almost_zero "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::is_one "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::is_minus_one "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::is_nan "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::is_inf "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::is_minus_inf "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::name "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::op "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::is_op "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::is_regular "

[INTERNAL] 
Checks if expression does not contain NaN or Inf.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L183

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L611-L617

";

%feature("docstring") casadi::SXElem::is_nonnegative "

[INTERNAL] 
Check if a value is always nonnegative (false negatives are 
allowed)

Extra doc: https://github.com/casadi/casadi/wiki/L_10t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L188

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L508-L516

";

%feature("docstring") casadi::SXElem::dep "

[INTERNAL] ";

%feature("docstring") casadi::SXElem::is_doubled "

[INTERNAL] 
Check if the node is the sum of two equal expressions.

Extra doc: https://github.com/casadi/casadi/wiki/L_10u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L155-L157

";

%feature("docstring") casadi::SXElem::n_dep "

[INTERNAL] 
Get the number of dependencies of a binary  SXElem.

Extra doc: https://github.com/casadi/casadi/wiki/L_10v

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L205

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L531-L533

";

%feature("docstring") casadi::SXElem::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given SXNode.

If the  SXElem does not point to any node, 0 is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_10w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L212

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L535-L537

";

%feature("docstring") casadi::SXElem::inv "

[INTERNAL] 
Element-wise inverse.

Extra doc: https://github.com/casadi/casadi/wiki/L_10y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L159-L165

";

%feature("docstring") casadi::SXElem::is_null "

[INTERNAL] 
 SXElem nodes are not allowed to be null.

Extra doc: https://github.com/casadi/casadi/wiki/L_112

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L260

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L260-L260

";

%feature("docstring") casadi::SXElem::serialize "

[INTERNAL] 
Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_114

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L272

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L619-L621

";


// File: classcasadi_1_1SymbolicQr.xml
%feature("docstring") casadi::SymbolicQr "

[INTERNAL] 
Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1TrilSolve.xml


// File: classcasadi_1_1TrilSolveUnity.xml


// File: classcasadi_1_1TriuSolve.xml


// File: classcasadi_1_1TriuSolveUnity.xml


// File: classcasadi_1_1WeakRef.xml
%feature("docstring") casadi::WeakRef "

[INTERNAL] 
Weak reference type.

A weak reference to a  SharedObject
Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_ax

C++ includes: shared_object.hpp
";

%feature("docstring") casadi::WeakRef::WeakRef "

[INTERNAL] 
Construct from a shared object (also implicit type conversion)

Extra doc: https://github.com/casadi/casadi/wiki/L_az

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L193

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L157-L159

>  casadi::WeakRef::WeakRef(SharedObject shared)
------------------------------------------------------------------------
[INTERNAL] 
Construct from a shared object (also implicit type conversion)

Extra doc: https://github.com/casadi/casadi/wiki/L_az

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L193

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L157-L159

";

";

%feature("docstring") casadi::WeakRef::shared "

[INTERNAL] 
Get a shared (owning) reference.

Extra doc: https://github.com/casadi/casadi/wiki/L_b0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L198

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L141-L147

";

%feature("docstring") casadi::WeakRef::alive "

[INTERNAL] 
Check if alive.

Extra doc: https://github.com/casadi/casadi/wiki/L_b1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L203

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L137-L139

";

%feature("docstring") casadi::WeakRef::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L99-L101

";

%feature("docstring") casadi::WeakRef::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L103-L109

";

%feature("docstring") casadi::WeakRef::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::WeakRef::is_null "

[INTERNAL] 
Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L73-L75

";

%feature("docstring") casadi::WeakRef::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L129-L131

";


// File: classcasadi_1_1XmlFile.xml
%feature("docstring") casadi::XmlFile "

[INTERNAL] 
XML parser.

Can be used for parsing XML files into CasADi data structures.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_7k

C++ includes: xml_file.hpp
";

%feature("docstring") casadi::XmlFile::XmlFile "

[INTERNAL] ";

%feature("docstring") casadi::XmlFile::~XmlFile "

[INTERNAL] ";

%feature("docstring") casadi::XmlFile::parse "

[INTERNAL] ";

%feature("docstring") casadi::XmlFile::dump "

[INTERNAL] ";

%feature("docstring") casadi::XmlFile::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L99-L101

";

%feature("docstring") casadi::XmlFile::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L103-L109

";

%feature("docstring") casadi::XmlFile::get_str "

[INTERNAL] 
Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::XmlFile::is_null "

[INTERNAL] 
Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L73-L75

";

%feature("docstring") casadi::XmlFile::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L129-L131

";


// File: namespace_0d364.xml


// File: namespacecasadi.xml
%feature("docstring") casadi::IndexRecution::collocation_points "
[INTERNAL] 
Obtain collocation points of specific order and scheme.

Parameters:
-----------

order: 
Which order (1 to 9 supported)

scheme: 
'radau' or 'legendre'

Extra doc: https://github.com/casadi/casadi/wiki/L_1so

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L120

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L120-L122

";

%feature("docstring") casadi::IndexRecution::collocation_pointsL "

[INTERNAL] 
Obtain collocation points of specific order and scheme.

Parameters:
-----------

order: 
Which order (1 to 9 supported)

scheme: 
'radau' or 'legendre'

Extra doc: https://github.com/casadi/casadi/wiki/L_1so

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L124

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L124-L126

";

%feature("docstring") casadi::IndexRecution::dae_reduce_index "

[INTERNAL] 
Reduce index.

Index reduction leads to a new set of variables and equations.

In the process, a set of constraints (algebraic equations or 
derivatives) 
a.k.a invariants is constructed that are invariant to the
 problem: whenever
 an initial point satisfies these constraints, the 
boundary-value-problem 
outcome will keep satisfying those constraints 
automatically, even though 
they are  not part of the reduced DAE.

For any practical numerical integration method, there will be 
numerical 
drift away from satisfaction of those constraints. In other 
words, you will
 see the value of invariants slowly moving away from 
original zero.

A classic mitigation technique is Baumgarte stabilization: you add 
these 
invariants to the reduced DAE as a correction term that acts in 
a way to 
make small (numerical) perturbations to the invariants decay 
to the origin 
as a dampened linear system.

in which a certain set of constraints (algebraic equations or 
derivatives) 
has been dropped in favour of

Parameters:
-----------

dae: 
Expression dictionary describing the DAE

Each value must be a dense column vector.

keys:
x_impl: symbol for implicit differential states

dx_impl: symbol for implicit differential state derivatives

z: symbol for algebraic variables

alg: expression for algebraic equations

t: symbol for time

p: symbol for parameters

Parameters:
-----------

opts: 
Option dictionary

'baumgarte_pole': double Poles (inverse time constants) of the 
Baumgarte 
invariant correction term. Must be <0 to dampen out 
perturbations 0 
(default) amounts to no correction. Corresponds to 
-gamma of equation (1.5)
 in Ascher, Uri M., Hongsheng Chin, and 
Sebastian Reich. \"Stabilization of
 DAEs and invariant manifolds.\" 
Numerische Mathematik 67.2 (1994): 
131-149.

Parameters:
-----------

stats: 
Statistics

Expression dictionary describing the reduced DAE

In addition the fields allowed in the input DAE, the following keys 
occur:

x: symbol for explicit differential states

ode: expression for right-hand-side of explicit differential states

I: expression for invariants

Extra doc: https://github.com/casadi/casadi/wiki/L_23h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L1060

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L1060-L1062

>  SXDict casadi::dae_reduce_index(const SXDict &dae, Dict &stats, const Dict &opts)
------------------------------------------------------------------------
[INTERNAL] 
Reduce index.

Index reduction leads to a new set of variables and equations.

In the process, a set of constraints (algebraic equations or 
derivatives) 
a.k.a invariants is constructed that are invariant to the
 problem: whenever
 an initial point satisfies these constraints, the 
boundary-value-problem 
outcome will keep satisfying those constraints 
automatically, even though 
they are  not part of the reduced DAE.

For any practical numerical integration method, there will be 
numerical 
drift away from satisfaction of those constraints. In other 
words, you will
 see the value of invariants slowly moving away from 
original zero.

A classic mitigation technique is Baumgarte stabilization: you add 
these 
invariants to the reduced DAE as a correction term that acts in 
a way to 
make small (numerical) perturbations to the invariants decay 
to the origin 
as a dampened linear system.

in which a certain set of constraints (algebraic equations or 
derivatives) 
has been dropped in favour of

Parameters:
-----------

dae: 
Expression dictionary describing the DAE

Each value must be a dense column vector.

keys:
x_impl: symbol for implicit differential states

dx_impl: symbol for implicit differential state derivatives

z: symbol for algebraic variables

alg: expression for algebraic equations

t: symbol for time

p: symbol for parameters

Parameters:
-----------

opts: 
Option dictionary

'baumgarte_pole': double Poles (inverse time constants) of the 
Baumgarte 
invariant correction term. Must be <0 to dampen out 
perturbations 0 
(default) amounts to no correction. Corresponds to 
-gamma of equation (1.5)
 in Ascher, Uri M., Hongsheng Chin, and 
Sebastian Reich. \"Stabilization of
 DAEs and invariant manifolds.\" 
Numerische Mathematik 67.2 (1994): 
131-149.

Parameters:
-----------

stats: 
Statistics

Expression dictionary describing the reduced DAE

In addition the fields allowed in the input DAE, the following keys 
occur:

x: symbol for explicit differential states

ode: expression for right-hand-side of explicit differential states

I: expression for invariants

Extra doc: https://github.com/casadi/casadi/wiki/L_23h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L1060

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L1060-L1062

";

";

%feature("docstring") casadi::IndexRecution::dae_map_semi_expl "

[INTERNAL] 
Turn a reduced DAE into a semi explicit form suitable for CasADi
 
integrator.

Parameters:
-----------

dae: 
Original (unreduced) DAE structure

dae_red: 
Reduced DAE (see dae_reduce_index)

state_to_orig: 
A mapping of integrator (semi explicit) states to states of 
the 
original DAE

phi: 
A function to compute the invariants of the reduced DAE Inputs:
x and 
z: (semi explicit) integrator states; typically integrator 
outputs xf and 
zf

p: parameters

t: time

Semi explicit DAE dictionary, suitable to pass to a CasADi integrator

See: 
 dae_reduce_index

Extra doc: https://github.com/casadi/casadi/wiki/L_1su

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L1205

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L1205-L1208

>  SXDict casadi::dae_map_semi_expl(const SXDict &dae, const SXDict &dae_red, Function &state_to_orig, Function &phi)
------------------------------------------------------------------------
[INTERNAL] 
Turn a reduced DAE into a semi explicit form suitable for CasADi 
integrator.

Parameters:
-----------

dae: 
Original (unreduced) DAE structure

dae_red: 
Reduced DAE (see dae_reduce_index)

state_to_orig: 
A mapping of integrator (semi explicit) states to states of 
the 
original DAE

phi: 
A function to compute the invariants of the reduced DAE Inputs:
x and 
z: (semi explicit) integrator states; typically integrator 
outputs xf and 
zf

p: parameters

t: time

Semi explicit DAE dictionary, suitable to pass to a CasADi integrator

See: 
 dae_reduce_index

Extra doc: https://github.com/casadi/casadi/wiki/L_1su

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L1205

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L1205-L1208

";

";

%feature("docstring") casadi::IndexRecution::dae_init_gen "

[INTERNAL] 
Obtain a generator  Function for producing consistent initial 
guesses of a reduced DAE.

Parameters:
-----------

dae: 
Original (unreduced) DAE structure

dae_red: 
Reduced DAE (see dae_reduce_index)

init_solver: 
NLP solver plugin name for nlpsol used to construct an initial
 guess

init_strength: 
Influence the nature of the NLP Structure with keys x_impl, 
dx_impl, z
 corresponding to inputs of init_gen Each key maps to a DM that 
should
 match the variable size corresponding to that key. For each variable
 
the meaning of the corresponding DM value is as follows: When >=0, 

indicates that the provided initial guess is used in a quadratic 
penalty 
(value used as weight) When -1, indicates that the provided 
initial guess 
must be observed (simple bound on variable)

init_solver_options: 
NLP solver options to be passed to nlpsol

init_gen A function to generate a consistent initial guess that can be
 used
 to pass to an integrator constructed from a semi explict reduced
 DAE 
Inputs:
x_impl, dx_impl, z: initial guesses in the original DAE space

p: parameters

t: time Outputs:

x0, z0: (semi explicit) integrator states and algebraic variables; 

typically used as input for integrators

See: 
 dae_reduce_index

Extra doc: https://github.com/casadi/casadi/wiki/L_1sv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L1215

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L1215-L1218

>  Function casadi::dae_init_gen(const SXDict &dae, const SXDict &dae_red, const std::string &init_solver, const DMDict &init_strength, const Dict &init_solver_options)
------------------------------------------------------------------------
[INTERNAL] 
Obtain a generator  Function for producing consistent initial guesses of a reduced DAE.

Parameters:
-----------

dae: 
Original (unreduced) DAE structure

dae_red: 
Reduced DAE (see dae_reduce_index)

init_solver: 
NLP solver plugin name for nlpsol used to construct an initial
 guess

init_strength: 
Influence the nature of the NLP Structure with keys x_impl, 
dx_impl, z
 corresponding to inputs of init_gen Each key maps to a DM that 
should
 match the variable size corresponding to that key. For each variable
 
the meaning of the corresponding DM value is as follows: When >=0, 

indicates that the provided initial guess is used in a quadratic 
penalty 
(value used as weight) When -1, indicates that the provided 
initial guess 
must be observed (simple bound on variable)

init_solver_options: 
NLP solver options to be passed to nlpsol

init_gen A function to generate a consistent initial guess that can be
 used
 to pass to an integrator constructed from a semi explict reduced
 DAE 
Inputs:
x_impl, dx_impl, z: initial guesses in the original DAE space

p: parameters

t: time Outputs:

x0, z0: (semi explicit) integrator states and algebraic variables; 

typically used as input for integrators

See: 
 dae_reduce_index

Extra doc: https://github.com/casadi/casadi/wiki/L_1sv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L1215

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L1215-L1218

";

";

%feature("docstring") casadi::IndexRecution::to_string "

[INTERNAL] 
Convert to string

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L43

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L43-L51

>  std::string casadi::to_string(DynOut v)
------------------------------------------------------------------------
[INTERNAL] 
Convert to string

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L43

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L43-L51

";

";

%feature("docstring") casadi::IndexRecution::matrixName "

[INTERNAL] 
Get typename.

Extra doc: https://github.com/casadi/casadi/wiki/L_18d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L56

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L56-L57

";

%feature("docstring") casadi::IndexRecution::matrixName< double > "

[INTERNAL] 
Get typename.

Extra doc: https://github.com/casadi/casadi/wiki/L_18d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L58

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L58-L58

";

%feature("docstring") casadi::IndexRecution::matrixName< casadi_int > "

[INTERNAL] 
Get typename.

Extra doc: https://github.com/casadi/casadi/wiki/L_18d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L59

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L59-L59

";

%feature("docstring") casadi::IndexRecution::nlpsol_default_in "

[INTERNAL] 
Default input for an NLP solver.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L222-L226

>  std::vector< double > casadi::nlpsol_default_in()
------------------------------------------------------------------------
[INTERNAL] 
Default input for an NLP solver.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L222-L226

";

";

%feature("docstring") casadi::IndexRecution::str "

[INTERNAL] 
String representation of a dictionary.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_common.hpp#L300

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_common.hpp#L300-L310

>  std::string casadi::str(const std::map< std::string, T2 > &p, bool more=false)
------------------------------------------------------------------------
[INTERNAL] 
String representation of a dictionary.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_common.hpp#L300

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_common.hpp#L300-L310

";

";

%feature("docstring") casadi::IndexRecution::strvec "

[INTERNAL] 
Create a list of strings from  VA_ARGS, six arguments.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_common.hpp#L225

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_common.hpp#L225-L228

>  std::vector<std::string> casadi::strvec(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6)
------------------------------------------------------------------------
[INTERNAL] 
Create a list of strings from  VA_ARGS, six arguments.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_common.hpp#L225

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_common.hpp#L225-L228

";

";

%feature("docstring") casadi::IndexRecution::fmtstr "

[INTERNAL] 
Create a string from a formatted string.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_common.hpp#L231

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_common.hpp#L231-L239

";

%feature("docstring") casadi::IndexRecution::is_zero "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::uout "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::uerr "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::to_int "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::all "

[INTERNAL] 
Check if all arguments are true.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L76

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L76-L81

";

%feature("docstring") casadi::IndexRecution::any "

[INTERNAL] 
Check if any arguments are true.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L83

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L83-L88

";

%feature("docstring") casadi::IndexRecution::is_range "

[INTERNAL] 
Check if a vector matches a range.

Extra doc: https://github.com/casadi/casadi/wiki/L_1l7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L90

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L90-L100

";

%feature("docstring") casadi::IndexRecution::range "

[INTERNAL] 
Range function.

Parameters:
-----------

stop:

list [0, 1, 2...stop-1]

Extra doc: https://github.com/casadi/casadi/wiki/L_1l8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L132-L134

>  std::vector< casadi_int > casadi::range(casadi_int stop)
------------------------------------------------------------------------
[INTERNAL] 
Range function.

Parameters:
-----------

stop:

list [0, 1, 2...stop-1]

Extra doc: https://github.com/casadi/casadi/wiki/L_1l8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L132-L134

";

";

%feature("docstring") casadi::IndexRecution::is_equally_spaced "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::complement "

[INTERNAL] 
Returns the list of all i in [0, size[ not found in supplied 
list.

The supplied vector may contain duplicates and may be non-monotonous 
The 
supplied vector will be checked for bounds The result vector is 
guaranteed 
to be monotonously increasing

Extra doc: https://github.com/casadi/casadi/wiki/L_1lf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L136

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L136-L152

";

%feature("docstring") casadi::IndexRecution::lookupvector "

[INTERNAL]

>  std::vector< casadi_int > casadi::lookupvector(const std::vector< casadi_int > &v)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::IndexRecution::is_permutation "

[INTERNAL] 
Does the list represent a permutation?

Extra doc: https://github.com/casadi/casadi/wiki/L_250

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L170

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L170-L175

";

%feature("docstring") casadi::IndexRecution::invert_permutation "

[INTERNAL] 
inverse a permutation vector

Extra doc: https://github.com/casadi/casadi/wiki/L_251

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L177

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L177-L184

";

%feature("docstring") casadi::IndexRecution::tensor_permute_mapping "

[INTERNAL] 
Computes a mapping for a (dense) tensor permutation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L187

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L187-L252

";

%feature("docstring") casadi::IndexRecution::join "

[INTERNAL] 
Join three lists.

Extra doc: https://github.com/casadi/casadi/wiki/L_1lc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L544

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L544-L549

>  std::vector< T > casadi::join(const std::vector< T > &a, const std::vector< T > &b, const std::vector< T > &c)
------------------------------------------------------------------------
[INTERNAL] 
Join three lists.

Extra doc: https://github.com/casadi/casadi/wiki/L_1lc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L544

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L544-L549

";

";

%feature("docstring") casadi::IndexRecution::startswith "

[INTERNAL] 
Checks if s starts with p.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L280

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L280-L286

";

%feature("docstring") casadi::IndexRecution::replace "

[INTERNAL] 
Replace all occurences of p with r in s.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L288

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L288-L297

";

%feature("docstring") casadi::IndexRecution::temporary_file "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::boolvec_not "

[INTERNAL] 
Invert all entries.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L371

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L371-L376

";

%feature("docstring") casadi::IndexRecution::boolvec_and "

[INTERNAL] 
And operation on boolean vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L378

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L378-L384

";

%feature("docstring") casadi::IndexRecution::boolvec_or "

[INTERNAL] 
Or operation on boolean vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L386

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L386-L392

";

%feature("docstring") casadi::IndexRecution::boolvec_to_index "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::normalized_setup "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::str_bvec "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::vector_static_cast "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::vector_slice "

[INTERNAL] 
Slicing vector.

Parameters:
-----------

v: 
Vector to slice

i: 
List of indices

Extra doc: https://github.com/casadi/casadi/wiki/L_1l9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L513

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L513-L527

";

%feature("docstring") casadi::IndexRecution::reverse "

[INTERNAL] 
Reverse a list.

Extra doc: https://github.com/casadi/casadi/wiki/L_1la

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L530

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L530-L534

";

%feature("docstring") casadi::IndexRecution::permute "

[INTERNAL] 
permute a list

Extra doc: https://github.com/casadi/casadi/wiki/L_1ld

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L552

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L552-L556

";

%feature("docstring") casadi::IndexRecution::find "

[INTERNAL] 
find nonzeros

Extra doc: https://github.com/casadi/casadi/wiki/L_1le

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L559

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L559-L565

";

%feature("docstring") casadi::IndexRecution::in_range "

[INTERNAL] 
Check if for each element of v holds: lower <= v_i < upper.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L613

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L613-L619

>  bool casadi::in_range(const std::vector< T > &v, casadi_int lower, casadi_int upper)
------------------------------------------------------------------------
[INTERNAL] 
Check if for each element of v holds: lower <= v_i < upper.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L613

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L613-L619

";

";

%feature("docstring") casadi::IndexRecution::flatten_nested_vector "

[INTERNAL] 
Flatten a nested std::vector tot a single flattened vector.

Contents of nested[i] ends up in 
flat[indices[i]]..flat[indices[i+1]-1]

Extra doc: https://github.com/casadi/casadi/wiki/L_1li

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L639

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L639-L654

>  void casadi::flatten_nested_vector(const std::vector< std::vector< T > > &nested, std::vector< S > &flat, std::vector< I > &indices)
------------------------------------------------------------------------
[INTERNAL] 
Flatten a nested std::vector tot a single flattened vector.

Contents of nested[i] ends up in 
flat[indices[i]]..flat[indices[i+1]-1]

Extra doc: https://github.com/casadi/casadi/wiki/L_1li

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L639

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L639-L654

";

";

%feature("docstring") casadi::IndexRecution::is_increasing "

[INTERNAL] 
Check if the vector is strictly increasing.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L663

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L663-L671

";

%feature("docstring") casadi::IndexRecution::is_decreasing "

[INTERNAL] 
Check if the vector is strictly decreasing.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L674

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L674-L682

";

%feature("docstring") casadi::IndexRecution::is_nonincreasing "

[INTERNAL] 
Check if the vector is non-increasing.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L685

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L685-L693

";

%feature("docstring") casadi::IndexRecution::is_nondecreasing "

[INTERNAL] 
Check if the vector is non-decreasing.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L696

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L696-L704

";

%feature("docstring") casadi::IndexRecution::is_monotone "

[INTERNAL] 
Check if the vector is monotone.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L707

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L707-L709

";

%feature("docstring") casadi::IndexRecution::is_strictly_monotone "

[INTERNAL] 
Check if the vector is strictly monotone.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L712

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L712-L714

";

%feature("docstring") casadi::IndexRecution::has_negative "

[INTERNAL] 
Check if the vector has negative entries.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L717

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L717-L722

";

%feature("docstring") casadi::IndexRecution::write_matlab "

[INTERNAL] 
Print matrix, matlab style.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L730

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L730-L735

>  void casadi::write_matlab(std::ostream &stream, const std::vector< std::vector< T > > &v)
------------------------------------------------------------------------
[INTERNAL] 
Print matrix, matlab style.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L730

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L730-L735

";

";

%feature("docstring") casadi::IndexRecution::read_matlab "

[INTERNAL] 
Read matrix, matlab style.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L758

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L758-L781

>  void casadi::read_matlab(std::ifstream &file, std::vector< std::vector< T > > &v)
------------------------------------------------------------------------
[INTERNAL] 
Read matrix, matlab style.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L758

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L758-L781

";

";

%feature("docstring") casadi::IndexRecution::linspace "

[INTERNAL] 
Matlab's linspace.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L784

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L784-L795

";

%feature("docstring") casadi::IndexRecution::sort "

[INTERNAL] 
Sort the data in a vector.

Parameters:
-----------

values: 
the vector that needs sorting

sorted_values: 
the sorted vector

indices: 
The indices such that 'sorted_values= values[indices]'

invert_indices: 
 Output indices such that 'sorted_values[indices=values'

Extra doc: https://github.com/casadi/casadi/wiki/L_1lj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L822

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L822-L847

";

%feature("docstring") casadi::IndexRecution::product "

[INTERNAL] 
product

Extra doc: https://github.com/casadi/casadi/wiki/L_1lk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L850

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L850-L854

";

%feature("docstring") casadi::IndexRecution::sum "

[INTERNAL] 
sum

Extra doc: https://github.com/casadi/casadi/wiki/L_1ll

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L857

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L857-L861

";

%feature("docstring") casadi::IndexRecution::cumsum "

[INTERNAL] 
cumulative sum

Extra doc: https://github.com/casadi/casadi/wiki/L_1lm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L864

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L864-L872

";

%feature("docstring") casadi::IndexRecution::diff "

[INTERNAL] 
diff

Extra doc: https://github.com/casadi/casadi/wiki/L_1ln

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L886

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L886-L893

";

%feature("docstring") casadi::IndexRecution::cumsum0 "

[INTERNAL] 
cumulative sum, starting with zero

Extra doc: https://github.com/casadi/casadi/wiki/L_1lo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L875

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L875-L883

";

%feature("docstring") casadi::IndexRecution::is_regular "

[INTERNAL] 
Checks if array does not contain NaN or Inf.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L393

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L393-L399

";

%feature("docstring") casadi::IndexRecution::normalized_out "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::normalized_in "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::applymap "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::copy_vector "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::assign_vector "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::init_vector "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::isUnique "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::get_ptr "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::dot "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::norm_inf "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::norm_1 "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::norm_2 "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::get_bvec_t "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::has_conic "

[INTERNAL] 
Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L31

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L31-L33

";

%feature("docstring") casadi::IndexRecution::load_conic "

[INTERNAL] 
Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L35

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L35-L37

";

%feature("docstring") casadi::IndexRecution::doc_conic "

[INTERNAL] 
Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L39

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L39-L41

";

%feature("docstring") casadi::IndexRecution::conic "

[INTERNAL]

Create a QP solver Solves the following strictly convex problem:



::

  min          1/2 x' H x + g' x
  x
  
  subject to
  LBA <= A x <= UBA
  LBX <= x   <= UBX
  
  resize(Q x, np, np) + P >= 0 (psd)
  
  with :
  H sparse (n x n) positive definite
  g dense  (n x 1)
  A sparse (nc x n)
  Q sparse symmetric (np^2 x n)
  P sparse symmetric (np x nq)
  
  n: number of decision variables (x)
  nc: number of constraints (A)
  nq: shape of psd constraint matrix



If H is not positive-definite, the solver should throw an error.

Second-order cone constraints can be added as psd constraints through 
a 
helper function 'soc':

x in R^n y in R

|| x ||_2 <= y

<=>

soc(x, y) psd

This can be proven with soc(x, y)=[y*I x; x' y] using the Shur 
complement.

General information

>List of available options

+------------------+-----------------+------------------+------------------+
|        Id        |      Type       |   Description    |     Used in      |
+==================+=================+==================+==================+
| ad_weight        | OT_DOUBLE       | Weighting factor | casadi::Function |
|                  |                 | for derivative   | Internal         |
|                  |                 | calculation.When |                  |
|                  |                 | there is an      |                  |
|                  |                 | option of either |                  |
|                  |                 | using forward or |                  |
|                  |                 | reverse mode     |                  |
|                  |                 | directional      |                  |
|                  |                 | derivatives, the |                  |
|                  |                 | condition ad_wei |                  |
|                  |                 | ght*nf<=(1-      |                  |
|                  |                 | ad_weight)*na is |                  |
|                  |                 | used where nf    |                  |
|                  |                 | and na are       |                  |
|                  |                 | estimates of the |                  |
|                  |                 | number of        |                  |
|                  |                 | forward/reverse  |                  |
|                  |                 | mode directional |                  |
|                  |                 | derivatives      |                  |
|                  |                 | needed. By       |                  |
|                  |                 | default,         |                  |
|                  |                 | ad_weight is     |                  |
|                  |                 | calculated       |                  |
|                  |                 | automatically,   |                  |
|                  |                 | but this can be  |                  |
|                  |                 | overridden by    |                  |
|                  |                 | setting this     |                  |
|                  |                 | option. In       |                  |
|                  |                 | particular, 0    |                  |
|                  |                 | means forcing    |                  |
|                  |                 | forward mode and |                  |
|                  |                 | 1 forcing        |                  |
|                  |                 | reverse mode.    |                  |
|                  |                 | Leave unset for  |                  |
|                  |                 | (class specific) |                  |
|                  |                 | heuristics.      |                  |
+------------------+-----------------+------------------+------------------+
| ad_weight_sp     | OT_DOUBLE       | Weighting factor | casadi::Function |
|                  |                 | for sparsity     | Internal         |
|                  |                 | pattern          |                  |
|                  |                 | calculation calc |                  |
|                  |                 | ulation.Override |                  |
|                  |                 | s default        |                  |
|                  |                 | behavior. Set to |                  |
|                  |                 | 0 and 1 to force |                  |
|                  |                 | forward and      |                  |
|                  |                 | reverse mode     |                  |
|                  |                 | respectively.    |                  |
|                  |                 | Cf. option       |                  |
|                  |                 | \"ad_weight\".     |                  |
|                  |                 | When set to -1,  |                  |
|                  |                 | sparsity is      |                  |
|                  |                 | completely       |                  |
|                  |                 | ignored and      |                  |
|                  |                 | dense matrices   |                  |
|                  |                 | are used.        |                  |
+------------------+-----------------+------------------+------------------+
| always_inline    | OT_BOOL         | Force inlining.  | casadi::Function |
|                  |                 |                  | Internal         |
+------------------+-----------------+------------------+------------------+
| cache            | OT_DICT         | Prepopulate the  | casadi::Function |
|                  |                 | function cache.  | Internal         |
|                  |                 | Default: empty   |                  |
+------------------+-----------------+------------------+------------------+
| compiler         | OT_STRING       | Just-in-time     | casadi::Function |
|                  |                 | compiler plugin  | Internal         |
|                  |                 | to be used.      |                  |
+------------------+-----------------+------------------+------------------+
| custom_jacobian  | OT_FUNCTION     | Override         | casadi::Function |
|                  |                 | CasADi's AD. Use | Internal         |
|                  |                 | together with    |                  |
|                  |                 | 'jac_penalty':   |                  |
|                  |                 | 0. Note: Highly  |                  |
|                  |                 | experimental.    |                  |
|                  |                 | Syntax may break |                  |
|                  |                 | often.           |                  |
+------------------+-----------------+------------------+------------------+
| der_options      | OT_DICT         | Default options  | casadi::Function |
|                  |                 | to be used to    | Internal         |
|                  |                 | populate         |                  |
|                  |                 | forward_options, |                  |
|                  |                 | reverse_options, |                  |
|                  |                 | and              |                  |
|                  |                 | jacobian_options |                  |
|                  |                 | before those     |                  |
|                  |                 | options are      |                  |
|                  |                 | merged in.       |                  |
+------------------+-----------------+------------------+------------------+
| derivative_of    | OT_FUNCTION     | The function is  | casadi::Function |
|                  |                 | a derivative of  | Internal         |
|                  |                 | another          |                  |
|                  |                 | function. The    |                  |
|                  |                 | type of          |                  |
|                  |                 | derivative       |                  |
|                  |                 | (directional     |                  |
|                  |                 | derivative,      |                  |
|                  |                 | Jacobian) is     |                  |
|                  |                 | inferred from    |                  |
|                  |                 | the function     |                  |
|                  |                 | name.            |                  |
+------------------+-----------------+------------------+------------------+
| discrete         | OT_BOOLVECTOR   | Indicates which  | casadi::Conic    |
|                  |                 | of the variables |                  |
|                  |                 | are discrete,    |                  |
|                  |                 | i.e. integer-    |                  |
|                  |                 | valued           |                  |
+------------------+-----------------+------------------+------------------+
| dump             | OT_BOOL         | Dump function to | casadi::Function |
|                  |                 | file upon first  | Internal         |
|                  |                 | evaluation.      |                  |
|                  |                 | [false]          |                  |
+------------------+-----------------+------------------+------------------+
| dump_dir         | OT_STRING       | Directory to     | casadi::Function |
|                  |                 | dump             | Internal         |
|                  |                 | inputs/outputs   |                  |
|                  |                 | to. Make sure    |                  |
|                  |                 | the directory    |                  |
|                  |                 | exists [.]       |                  |
+------------------+-----------------+------------------+------------------+
| dump_format      | OT_STRING       | Choose file      | casadi::Function |
|                  |                 | format to dump   | Internal         |
|                  |                 | matrices. See    |                  |
|                  |                 | DM.from_file     |                  |
|                  |                 | [mtx]            |                  |
+------------------+-----------------+------------------+------------------+
| dump_in          | OT_BOOL         | Dump numerical   | casadi::Function |
|                  |                 | values of inputs | Internal         |
|                  |                 | to file          |                  |
|                  |                 | (readable with   |                  |
|                  |                 | DM.from_file )   |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| dump_out         | OT_BOOL         | Dump numerical   | casadi::Function |
|                  |                 | values of        | Internal         |
|                  |                 | outputs to file  |                  |
|                  |                 | (readable with   |                  |
|                  |                 | DM.from_file )   |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| enable_fd        | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation by   |                  |
|                  |                 | finite           |                  |
|                  |                 | differencing.    |                  |
|                  |                 | [default:        |                  |
|                  |                 | false]]          |                  |
+------------------+-----------------+------------------+------------------+
| enable_forward   | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | Jacobian-times-  |                  |
|                  |                 | vector products  |                  |
|                  |                 | - typically      |                  |
|                  |                 | using forward    |                  |
|                  |                 | mode AD - if     |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| enable_jacobian  | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | Jacobians of all |                  |
|                  |                 | differentiable   |                  |
|                  |                 | outputs with     |                  |
|                  |                 | respect to all   |                  |
|                  |                 | differentiable   |                  |
|                  |                 | inputs - if      |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| enable_reverse   | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | transposed       |                  |
|                  |                 | Jacobian-times-  |                  |
|                  |                 | vector products  |                  |
|                  |                 | - typically      |                  |
|                  |                 | using reverse    |                  |
|                  |                 | mode AD - if     |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| error_on_fail    | OT_BOOL         | Throw exceptions | casadi::ProtoFun |
|                  |                 | when function    | ction            |
|                  |                 | evaluation fails |                  |
|                  |                 | (default true).  |                  |
+------------------+-----------------+------------------+------------------+
| fd_method        | OT_STRING       | Method for       | casadi::Function |
|                  |                 | finite           | Internal         |
|                  |                 | differencing     |                  |
|                  |                 | [default         |                  |
|                  |                 | 'central']       |                  |
+------------------+-----------------+------------------+------------------+
| fd_options       | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to the    | Internal         |
|                  |                 | finite           |                  |
|                  |                 | difference       |                  |
|                  |                 | instance         |                  |
+------------------+-----------------+------------------+------------------+
| forward_options  | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | forward mode     |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| gather_stats     | OT_BOOL         | Deprecated       | casadi::Function |
|                  |                 | option           | Internal         |
|                  |                 | (ignored):       |                  |
|                  |                 | Statistics are   |                  |
|                  |                 | now always       |                  |
|                  |                 | collected.       |                  |
+------------------+-----------------+------------------+------------------+
| input_scheme     | OT_STRINGVECTOR | Deprecated       | casadi::Function |
|                  |                 | option (ignored) | Internal         |
+------------------+-----------------+------------------+------------------+
| inputs_check     | OT_BOOL         | Throw exceptions | casadi::Function |
|                  |                 | when the         | Internal         |
|                  |                 | numerical values |                  |
|                  |                 | of the inputs    |                  |
|                  |                 | don't make sense |                  |
+------------------+-----------------+------------------+------------------+
| is_diff_in       | OT_BOOLVECTOR   | Indicate for     | casadi::Function |
|                  |                 | each input if it | Internal         |
|                  |                 | should be        |                  |
|                  |                 | differentiable.  |                  |
+------------------+-----------------+------------------+------------------+
| is_diff_out      | OT_BOOLVECTOR   | Indicate for     | casadi::Function |
|                  |                 | each output if   | Internal         |
|                  |                 | it should be     |                  |
|                  |                 | differentiable.  |                  |
+------------------+-----------------+------------------+------------------+
| jac_penalty      | OT_DOUBLE       | When requested   | casadi::Function |
|                  |                 | for a number of  | Internal         |
|                  |                 | forward/reverse  |                  |
|                  |                 | directions, it   |                  |
|                  |                 | may be cheaper   |                  |
|                  |                 | to compute first |                  |
|                  |                 | the full         |                  |
|                  |                 | jacobian and     |                  |
|                  |                 | then multiply    |                  |
|                  |                 | with seeds,      |                  |
|                  |                 | rather than      |                  |
|                  |                 | obtain the       |                  |
|                  |                 | requested        |                  |
|                  |                 | directions in a  |                  |
|                  |                 | straightforward  |                  |
|                  |                 | manner. Casadi   |                  |
|                  |                 | uses a heuristic |                  |
|                  |                 | to decide which  |                  |
|                  |                 | is cheaper. A    |                  |
|                  |                 | high value of    |                  |
|                  |                 | 'jac_penalty'    |                  |
|                  |                 | makes it less    |                  |
|                  |                 | likely for the   |                  |
|                  |                 | heurstic to      |                  |
|                  |                 | chose the full   |                  |
|                  |                 | Jacobian         |                  |
|                  |                 | strategy. The    |                  |
|                  |                 | special value -1 |                  |
|                  |                 | indicates never  |                  |
|                  |                 | to use the full  |                  |
|                  |                 | Jacobian         |                  |
|                  |                 | strategy         |                  |
+------------------+-----------------+------------------+------------------+
| jacobian_options | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | Jacobian         |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| jit              | OT_BOOL         | Use just-in-time | casadi::Function |
|                  |                 | compiler to      | Internal         |
|                  |                 | speed up the     |                  |
|                  |                 | evaluation       |                  |
+------------------+-----------------+------------------+------------------+
| jit_cleanup      | OT_BOOL         | Cleanup up the   | casadi::Function |
|                  |                 | temporary source | Internal         |
|                  |                 | file that jit    |                  |
|                  |                 | creates.         |                  |
|                  |                 | Default: true    |                  |
+------------------+-----------------+------------------+------------------+
| jit_name         | OT_STRING       | The file name    | casadi::Function |
|                  |                 | used to write    | Internal         |
|                  |                 | out code. The    |                  |
|                  |                 | actual file      |                  |
|                  |                 | names used       |                  |
|                  |                 | depend on 'jit_t |                  |
|                  |                 | emp_suffix' and  |                  |
|                  |                 | include          |                  |
|                  |                 | extensions.      |                  |
|                  |                 | Default:         |                  |
|                  |                 | 'jit_tmp'        |                  |
+------------------+-----------------+------------------+------------------+
| jit_options      | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to the    | Internal         |
|                  |                 | jit compiler.    |                  |
+------------------+-----------------+------------------+------------------+
| jit_serialize    | OT_STRING       | Specify          | casadi::Function |
|                  |                 | behaviour when   | Internal         |
|                  |                 | serializing a    |                  |
|                  |                 | jitted function: |                  |
|                  |                 | SOURCE|link|embe |                  |
|                  |                 | d.               |                  |
+------------------+-----------------+------------------+------------------+
| jit_temp_suffix  | OT_BOOL         | Use a temporary  | casadi::Function |
|                  |                 | (seemingly       | Internal         |
|                  |                 | random) filename |                  |
|                  |                 | suffix for       |                  |
|                  |                 | generated code   |                  |
|                  |                 | and libraries.   |                  |
|                  |                 | This is desired  |                  |
|                  |                 | for thread-      |                  |
|                  |                 | safety. This     |                  |
|                  |                 | behaviour may    |                  |
|                  |                 | defeat caching   |                  |
|                  |                 | compiler         |                  |
|                  |                 | wrappers.        |                  |
|                  |                 | Default: true    |                  |
+------------------+-----------------+------------------+------------------+
| max_io           | OT_INT          | Acceptable       | casadi::Function |
|                  |                 | number of inputs | Internal         |
|                  |                 | and outputs.     |                  |
|                  |                 | Warn if          |                  |
|                  |                 | exceeded.        |                  |
+------------------+-----------------+------------------+------------------+
| max_num_dir      | OT_INT          | Specify the      | casadi::Function |
|                  |                 | maximum number   | Internal         |
|                  |                 | of directions    |                  |
|                  |                 | for derivative   |                  |
|                  |                 | functions.       |                  |
|                  |                 | Overrules the    |                  |
|                  |                 | builtin optimize |                  |
|                  |                 | d_num_dir.       |                  |
+------------------+-----------------+------------------+------------------+
| never_inline     | OT_BOOL         | Forbid inlining. | casadi::Function |
|                  |                 |                  | Internal         |
+------------------+-----------------+------------------+------------------+
| output_scheme    | OT_STRINGVECTOR | Deprecated       | casadi::Function |
|                  |                 | option (ignored) | Internal         |
+------------------+-----------------+------------------+------------------+
| post_expand      | OT_BOOL         | After            | casadi::Function |
|                  |                 | construction,    | Internal         |
|                  |                 | expand this      |                  |
|                  |                 | Function .       |                  |
|                  |                 | Default: False   |                  |
+------------------+-----------------+------------------+------------------+
| post_expand_opti | OT_DICT         | Options to be    | casadi::Function |
| ons              |                 | passed to post-  | Internal         |
|                  |                 | construction     |                  |
|                  |                 | expansion.       |                  |
|                  |                 | Default: empty   |                  |
+------------------+-----------------+------------------+------------------+
| print_in         | OT_BOOL         | Print numerical  | casadi::Function |
|                  |                 | values of inputs | Internal         |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| print_out        | OT_BOOL         | Print numerical  | casadi::Function |
|                  |                 | values of        | Internal         |
|                  |                 | outputs          |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| print_problem    | OT_BOOL         | Print a numeric  | casadi::Conic    |
|                  |                 | description of   |                  |
|                  |                 | the problem      |                  |
+------------------+-----------------+------------------+------------------+
| print_time       | OT_BOOL         | print            | casadi::ProtoFun |
|                  |                 | information      | ction            |
|                  |                 | about execution  |                  |
|                  |                 | time. Implies    |                  |
|                  |                 | record_time.     |                  |
+------------------+-----------------+------------------+------------------+
| record_time      | OT_BOOL         | record           | casadi::ProtoFun |
|                  |                 | information      | ction            |
|                  |                 | about execution  |                  |
|                  |                 | time, for        |                  |
|                  |                 | retrieval with   |                  |
|                  |                 | stats().         |                  |
+------------------+-----------------+------------------+------------------+
| regularity_check | OT_BOOL         | Throw exceptions | casadi::ProtoFun |
|                  |                 | when NaN or Inf  | ction            |
|                  |                 | appears during   |                  |
|                  |                 | evaluation       |                  |
+------------------+-----------------+------------------+------------------+
| reverse_options  | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | reverse mode     |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| user_data        | OT_VOIDPTR      | A user-defined   | casadi::Function |
|                  |                 | field that can   | Internal         |
|                  |                 | be used to       |                  |
|                  |                 | identify the     |                  |
|                  |                 | function or pass |                  |
|                  |                 | additional       |                  |
|                  |                 | information      |                  |
+------------------+-----------------+------------------+------------------+
| verbose          | OT_BOOL         | Verbose          | casadi::ProtoFun |
|                  |                 | evaluation  for  | ction            |
|                  |                 | debugging        |                  |
+------------------+-----------------+------------------+------------------+

>Input scheme: casadi::ConicInput (CONIC_NUM_IN = 12)

+--------------+--------+--------------------------------------------------+
|  Full name   | Short  |                   Description                    |
+==============+========+==================================================+
| CONIC_H      | h      | The square matrix H: sparse, (n x n). Only the   |
|              |        | lower triangular part is actually used. The      |
|              |        | matrix is assumed to be symmetrical.             |
+--------------+--------+--------------------------------------------------+
| CONIC_G      | g      | The vector g: dense, (n x 1)                     |
+--------------+--------+--------------------------------------------------+
| CONIC_A      | a      | The matrix A: sparse, (nc x n) - product with x  |
|              |        | must be dense.                                   |
+--------------+--------+--------------------------------------------------+
| CONIC_LBA    | lba    | dense, (nc x 1)                                  |
+--------------+--------+--------------------------------------------------+
| CONIC_UBA    | uba    | dense, (nc x 1)                                  |
+--------------+--------+--------------------------------------------------+
| CONIC_LBX    | lbx    | dense, (n x 1)                                   |
+--------------+--------+--------------------------------------------------+
| CONIC_UBX    | ubx    | dense, (n x 1)                                   |
+--------------+--------+--------------------------------------------------+
| CONIC_X0     | x0     | dense, (n x 1)                                   |
+--------------+--------+--------------------------------------------------+
| CONIC_LAM_X0 | lam_x0 | dense                                            |
+--------------+--------+--------------------------------------------------+
| CONIC_LAM_A0 | lam_a0 | dense                                            |
+--------------+--------+--------------------------------------------------+
| CONIC_Q      | q      | The matrix Q: sparse symmetric, (np^2 x n)       |
+--------------+--------+--------------------------------------------------+
| CONIC_P      | p      | The matrix P: sparse symmetric, (np x np)        |
+--------------+--------+--------------------------------------------------+

>Output scheme: casadi::ConicOutput (CONIC_NUM_OUT = 4)

+-------------+-------+---------------------------------------------------+
|  Full name  | Short |                    Description                    |
+=============+=======+===================================================+
| CONIC_X     | x     | The primal solution.                              |
+-------------+-------+---------------------------------------------------+
| CONIC_COST  | cost  | The optimal cost.                                 |
+-------------+-------+---------------------------------------------------+
| CONIC_LAM_A | lam_a | The dual solution corresponding to linear bounds. |
+-------------+-------+---------------------------------------------------+
| CONIC_LAM_X | lam_x | The dual solution corresponding to simple bounds. |
+-------------+-------+---------------------------------------------------+

List of plugins
- cbc

- clp

- cplex

- gurobi

- highs

- hpipm

- hpmpc

- ooqp

- osqp

- proxqp

- qpoases

- sqic

- superscs

- ipqp

- nlpsol

- qrqp

Note: some of the plugins in this list might not be available on your 

system.  Also, there might be extra plugins available to you that are 
not 
listed here. You can obtain their documentation with  
Conic.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

cbc
---



Interface to Cbc solver for sparse Quadratic Programs

Extra doc: https://github.com/casadi/casadi/wiki/L_221

>List of available options

+-------------+-----------------------+------------------------------------+
|     Id      |         Type          |            Description             |
+=============+=======================+====================================+
| cbc         | OT_DICT               | Options to be passed to CBC.Three  |
|             |                       | sets of options are supported. The |
|             |                       | first can be found in              |
|             |                       | OsiSolverParameters.hpp. The       |
|             |                       | second can be found in             |
|             |                       | CbcModel.hpp. The third are        |
|             |                       | options that can be passed to      |
|             |                       | CbcMain1.                          |
+-------------+-----------------------+------------------------------------+
| hot_start   | OT_BOOL               | Hot start with x0 [Default false]. |
+-------------+-----------------------+------------------------------------+
| sos_groups  | OT_INTVECTORVECTOR    | Definition of SOS groups by        |
|             |                       | indices.                           |
+-------------+-----------------------+------------------------------------+
| sos_types   | OT_INTVECTOR          | Specify 1 or 2 for each SOS group. |
+-------------+-----------------------+------------------------------------+
| sos_weights | OT_DOUBLEVECTORVECTOR | Weights corresponding to SOS       |
|             |                       | entries.                           |
+-------------+-----------------------+------------------------------------+



--------------------------------------------------------------------------------

clp
---



Interface to Clp solver for sparse Quadratic Programs

Extra doc: https://github.com/casadi/casadi/wiki/L_22d

>List of available options

+-----+---------+----------------------------------------------------------+
| Id  |  Type   |                       Description                        |
+=====+=========+==========================================================+
| clp | OT_DICT | Options to be passed to CLP. A first set of options can  |
|     |         | be found in ClpParameters.hpp. eg. 'PrimalTolerance'.    |
|     |         | There are other options in additions. 'AutomaticScaling' |
|     |         | (bool) is recognised. 'initial_solve' (default off)      |
|     |         | activates the use of Clp's initialSolve.                 |
|     |         | 'initial_solve_options' takes a dictionary with          |
|     |         | following keys (see ClpSolve.hpp): SolveType (string),   |
|     |         | PresolveType (string), NumberPasses, SpecialOptions      |
|     |         | (intvectorvector), IndependentOptions (intvectorvector). |
+-----+---------+----------------------------------------------------------+



--------------------------------------------------------------------------------

cplex
-----



Interface to Cplex solver for sparse Quadratic Programs

Extra doc: https://github.com/casadi/casadi/wiki/L_22a

>List of available options

+----------------+-----------------------+---------------------------------+
|       Id       |         Type          |           Description           |
+================+=======================+=================================+
| cplex          | OT_DICT               | Options to be passed to CPLEX   |
+----------------+-----------------------+---------------------------------+
| dep_check      | OT_INT                | Detect redundant constraints.   |
+----------------+-----------------------+---------------------------------+
| dump_filename  | OT_STRING             | The filename to dump to.        |
+----------------+-----------------------+---------------------------------+
| dump_to_file   | OT_BOOL               | Dumps QP to file in CPLEX       |
|                |                       | format.                         |
+----------------+-----------------------+---------------------------------+
| mip_start      | OT_BOOL               | Hot start integers with x0      |
|                |                       | [Default false].                |
+----------------+-----------------------+---------------------------------+
| qp_method      | OT_INT                | Determines which CPLEX          |
|                |                       | algorithm to use.               |
+----------------+-----------------------+---------------------------------+
| sos_groups     | OT_INTVECTORVECTOR    | Definition of SOS groups by     |
|                |                       | indices.                        |
+----------------+-----------------------+---------------------------------+
| sos_types      | OT_INTVECTOR          | Specify 1 or 2 for each SOS     |
|                |                       | group.                          |
+----------------+-----------------------+---------------------------------+
| sos_weights    | OT_DOUBLEVECTORVECTOR | Weights corresponding to SOS    |
|                |                       | entries.                        |
+----------------+-----------------------+---------------------------------+
| tol            | OT_DOUBLE             | Tolerance of solver             |
+----------------+-----------------------+---------------------------------+
| version_suffix | OT_STRING             | Specify version of cplex to     |
|                |                       | load. We will attempt to load l |
|                |                       | ibcplex<version_suffix>.[so|dll |
|                |                       | |dylib]. Default value is taken |
|                |                       | from CPLEX_VERSION env          |
|                |                       | variable.                       |
+----------------+-----------------------+---------------------------------+
| warm_start     | OT_BOOL               | Use warm start with simplex     |
|                |                       | methods (affects only the       |
|                |                       | simplex methods).               |
+----------------+-----------------------+---------------------------------+



--------------------------------------------------------------------------------

gurobi
------



Interface to the GUROBI Solver for quadratic programming

Extra doc: https://github.com/casadi/casadi/wiki/L_22q

>List of available options

+-------------+-----------------------+------------------------------------+
|     Id      |         Type          |            Description             |
+=============+=======================+====================================+
| gurobi      | OT_DICT               | Options to be passed to gurobi.    |
+-------------+-----------------------+------------------------------------+
| sos_groups  | OT_INTVECTORVECTOR    | Definition of SOS groups by        |
|             |                       | indices.                           |
+-------------+-----------------------+------------------------------------+
| sos_types   | OT_INTVECTOR          | Specify 1 or 2 for each SOS group. |
+-------------+-----------------------+------------------------------------+
| sos_weights | OT_DOUBLEVECTORVECTOR | Weights corresponding to SOS       |
|             |                       | entries.                           |
+-------------+-----------------------+------------------------------------+
| vtype       | OT_STRINGVECTOR       | Type of variables: [CONTINUOUS|bin |
|             |                       | ary|integer|semicont|semiint]      |
+-------------+-----------------------+------------------------------------+



--------------------------------------------------------------------------------

highs
-----



Interface to HiGHS solver for sparse Quadratic Programs, see 
highs.dev for 
more information and https://www.maths.ed.ac.uk/hall/HiGHS/HighsOptions.html
  for a list of options.

Extra doc: https://github.com/casadi/casadi/wiki/L_22f

>List of available options

+-------+---------+--------------------------------+
|  Id   |  Type   |          Description           |
+=======+=========+================================+
| highs | OT_DICT | Options to be passed to HiGHS. |
+-------+---------+--------------------------------+



--------------------------------------------------------------------------------

hpipm
-----



Interface to HMPC Solver

In order to use this interface, you must:

Decision variables must only by state and control, and the variable 

ordering must be [x0 u0 x1 u1 ...]

The constraints must be in order: [ gap0 lincon0 gap1 lincon1 ]

gap: Ak+1 = Ak xk + Bk uk lincon: yk= Ck xk + Dk uk



::

         A0 B0 -I
         C0 D0
                A1 B1 -I
                C1 D1



where I must be a diagonal sparse matrix
Either supply all of N, nx, ng, nu 
options or rely on automatic 
detection

Extra doc: https://github.com/casadi/casadi/wiki/L_242

>List of available options

+-------+--------------+--------------------------------------------------+
|  Id   |     Type     |                   Description                    |
+=======+==============+==================================================+
| N     | OT_INT       | OCP horizon                                      |
+-------+--------------+--------------------------------------------------+
| hpipm | OT_DICT      | Options to be passed to hpipm                    |
+-------+--------------+--------------------------------------------------+
| inf   | OT_DOUBLE    | Replace infinities by this amount [default: 1e8] |
+-------+--------------+--------------------------------------------------+
| ng    | OT_INTVECTOR | Number of non-dynamic constraints, length N+1    |
+-------+--------------+--------------------------------------------------+
| nu    | OT_INTVECTOR | Number of controls, length N                     |
+-------+--------------+--------------------------------------------------+
| nx    | OT_INTVECTOR | Number of states, length N+1                     |
+-------+--------------+--------------------------------------------------+



--------------------------------------------------------------------------------

hpmpc
-----



Interface to HMPC Solver

In order to use this interface, you must:

Decision variables must only by state and control, and the variable 

ordering must be [x0 u0 x1 u1 ...]

The constraints must be in order: [ gap0 lincon0 gap1 lincon1 ]

gap: Ak+1 = Ak xk + Bk uk lincon: yk= Ck xk + Dk uk



::

         A0 B0 -I
         C0 D0
                A1 B1 -I
                C1 D1



where I must be a diagonal sparse matrix
Either supply all of N, nx, ng, nu 
options or rely on automatic 
detection

Extra doc: https://github.com/casadi/casadi/wiki/L_22p

>List of available options

+----------------+--------------+------------------------------------------+
|       Id       |     Type     |               Description                |
+================+==============+==========================================+
| N              | OT_INT       | OCP horizon                              |
+----------------+--------------+------------------------------------------+
| blasfeo_target | OT_STRING    | hpmpc target                             |
+----------------+--------------+------------------------------------------+
| inf            | OT_DOUBLE    | HPMPC cannot handle infinities.          |
|                |              | Infinities will be replaced by this      |
|                |              | option's value.                          |
+----------------+--------------+------------------------------------------+
| max_iter       | OT_INT       | Max number of iterations                 |
+----------------+--------------+------------------------------------------+
| mu0            | OT_DOUBLE    | Max element in cost function as estimate |
|                |              | of max multiplier                        |
+----------------+--------------+------------------------------------------+
| ng             | OT_INTVECTOR | Number of non-dynamic constraints,       |
|                |              | length N+1                               |
+----------------+--------------+------------------------------------------+
| nu             | OT_INTVECTOR | Number of controls, length N             |
+----------------+--------------+------------------------------------------+
| nx             | OT_INTVECTOR | Number of states, length N+1             |
+----------------+--------------+------------------------------------------+
| print_level    | OT_INT       | Amount of diagnostic printing [Default:  |
|                |              | 1].                                      |
+----------------+--------------+------------------------------------------+
| target         | OT_STRING    | hpmpc target                             |
+----------------+--------------+------------------------------------------+
| tol            | OT_DOUBLE    | Tolerance in the duality measure         |
+----------------+--------------+------------------------------------------+
| warm_start     | OT_BOOL      | Use warm-starting                        |
+----------------+--------------+------------------------------------------+



--------------------------------------------------------------------------------

ooqp
----



Interface to the OOQP Solver for quadratic programming The current 

implementation assumes that OOQP is configured with the MA27 sparse 
linear 
solver.

NOTE: when doing multiple calls to evaluate(), check if you need to 

reInit();

Extra doc: https://github.com/casadi/casadi/wiki/L_222

>List of available options

+-------------+-----------+------------------------------------------------+
|     Id      |   Type    |                  Description                   |
+=============+===========+================================================+
| artol       | OT_DOUBLE | tolerance as provided with setArTol to OOQP    |
+-------------+-----------+------------------------------------------------+
| mutol       | OT_DOUBLE | tolerance as provided with setMuTol to OOQP    |
+-------------+-----------+------------------------------------------------+
| print_level | OT_INT    | Print level. OOQP listens to print_level 0, 10 |
|             |           | and 100                                        |
+-------------+-----------+------------------------------------------------+



--------------------------------------------------------------------------------

osqp
----



Interface to the OSQP Solver for quadratic programming

Extra doc: https://github.com/casadi/casadi/wiki/L_220

Interface to the PROXQP Solver for quadratic programming

Extra doc: https://github.com/casadi/casadi/wiki/L_243

>List of available options

+-------------------+---------+--------------------------------------------+
|        Id         |  Type   |                Description                 |
+===================+=========+============================================+
| osqp              | OT_DICT | const Options to be passed to osqp.        |
+-------------------+---------+--------------------------------------------+
| warm_start_dual   | OT_BOOL | Use lam_a0 and lam_x0 input to warmstart   |
|                   |         | [Default: truw].                           |
+-------------------+---------+--------------------------------------------+
| warm_start_primal | OT_BOOL | Use x0 input to warmstart [Default: true]. |
+-------------------+---------+--------------------------------------------+



--------------------------------------------------------------------------------

proxqp
------



>List of available options

+-------------------+---------+--------------------------------------------+
|        Id         |  Type   |                Description                 |
+===================+=========+============================================+
| proxqp            | OT_DICT | const proxqp options.                      |
+-------------------+---------+--------------------------------------------+
| warm_start_dual   | OT_BOOL | Use y and z input to warmstart [Default:   |
|                   |         | true].                                     |
+-------------------+---------+--------------------------------------------+
| warm_start_primal | OT_BOOL | Use x input to warmstart [Default: true].  |
+-------------------+---------+--------------------------------------------+



--------------------------------------------------------------------------------

qpoases
-------



Interface to QPOases Solver for quadratic programming

::

  Extra doc: https://github.com/casadi/casadi/wiki/L_22o 
  



>List of available options

+-------------------------------+-----------+------------------------------+
|              Id               |   Type    |         Description          |
+===============================+===========+==============================+
| CPUtime                       | OT_DOUBLE | The maximum allowed CPU time |
|                               |           | in seconds for the whole     |
|                               |           | initialisation (and the      |
|                               |           | actually required one on     |
|                               |           | output). Disabled if unset.  |
+-------------------------------+-----------+------------------------------+
| boundRelaxation               | OT_DOUBLE | Initial relaxation of bounds |
|                               |           | to start homotopy and        |
|                               |           | initial value for far        |
|                               |           | bounds.                      |
+-------------------------------+-----------+------------------------------+
| boundTolerance                | OT_DOUBLE | If upper and lower bounds    |
|                               |           | differ less than this        |
|                               |           | tolerance, they are regarded |
|                               |           | equal, i.e. as equality      |
|                               |           | constraint.                  |
+-------------------------------+-----------+------------------------------+
| enableCholeskyRefactorisation | OT_INT    | Specifies the frequency of a |
|                               |           | full re-factorisation of     |
|                               |           | projected Hessian matrix: 0: |
|                               |           | turns them off, 1: uses them |
|                               |           | at each iteration etc.       |
+-------------------------------+-----------+------------------------------+
| enableDriftCorrection         | OT_INT    | Specifies the frequency of   |
|                               |           | drift corrections: 0: turns  |
|                               |           | them off.                    |
+-------------------------------+-----------+------------------------------+
| enableEqualities              | OT_BOOL   | Specifies whether equalities |
|                               |           | should be treated as always  |
|                               |           | active (True) or not (False) |
+-------------------------------+-----------+------------------------------+
| enableFarBounds               | OT_BOOL   | Enables the use of far       |
|                               |           | bounds.                      |
+-------------------------------+-----------+------------------------------+
| enableFlippingBounds          | OT_BOOL   | Enables the use of flipping  |
|                               |           | bounds.                      |
+-------------------------------+-----------+------------------------------+
| enableFullLITests             | OT_BOOL   | Enables condition-hardened   |
|                               |           | (but more expensive) LI      |
|                               |           | test.                        |
+-------------------------------+-----------+------------------------------+
| enableInertiaCorrection       | OT_BOOL   | Should working set be        |
|                               |           | repaired when negative       |
|                               |           | curvature is discovered      |
|                               |           | during hotstart.             |
+-------------------------------+-----------+------------------------------+
| enableNZCTests                | OT_BOOL   | Enables nonzero curvature    |
|                               |           | tests.                       |
+-------------------------------+-----------+------------------------------+
| enableRamping                 | OT_BOOL   | Enables ramping.             |
+-------------------------------+-----------+------------------------------+
| enableRegularisation          | OT_BOOL   | Enables automatic Hessian    |
|                               |           | regularisation.              |
+-------------------------------+-----------+------------------------------+
| epsDen                        | OT_DOUBLE | Denominator tolerance for    |
|                               |           | ratio tests.                 |
+-------------------------------+-----------+------------------------------+
| epsFlipping                   | OT_DOUBLE | Tolerance of squared         |
|                               |           | Cholesky diagonal factor     |
|                               |           | which triggers flipping      |
|                               |           | bound.                       |
+-------------------------------+-----------+------------------------------+
| epsIterRef                    | OT_DOUBLE | Early termination tolerance  |
|                               |           | for iterative refinement.    |
+-------------------------------+-----------+------------------------------+
| epsLITests                    | OT_DOUBLE | Tolerance for linear         |
|                               |           | independence tests.          |
+-------------------------------+-----------+------------------------------+
| epsNZCTests                   | OT_DOUBLE | Tolerance for nonzero        |
|                               |           | curvature tests.             |
+-------------------------------+-----------+------------------------------+
| epsNum                        | OT_DOUBLE | Numerator tolerance for      |
|                               |           | ratio tests.                 |
+-------------------------------+-----------+------------------------------+
| epsRegularisation             | OT_DOUBLE | Scaling factor of identity   |
|                               |           | matrix used for Hessian      |
|                               |           | regularisation.              |
+-------------------------------+-----------+------------------------------+
| finalRamping                  | OT_DOUBLE | Final value for ramping      |
|                               |           | strategy.                    |
+-------------------------------+-----------+------------------------------+
| growFarBounds                 | OT_DOUBLE | Factor to grow far bounds.   |
+-------------------------------+-----------+------------------------------+
| hessian_type                  | OT_STRING | Type of Hessian - see        |
|                               |           | qpOASES documentation [UNKNO |
|                               |           | WN|posdef|semidef|indef|zero |
|                               |           | |identity]]                  |
+-------------------------------+-----------+------------------------------+
| initialFarBounds              | OT_DOUBLE | Initial size for far bounds. |
+-------------------------------+-----------+------------------------------+
| initialRamping                | OT_DOUBLE | Start value for ramping      |
|                               |           | strategy.                    |
+-------------------------------+-----------+------------------------------+
| initialStatusBounds           | OT_STRING | Initial status of bounds at  |
|                               |           | first iteration.             |
+-------------------------------+-----------+------------------------------+
| linsol_plugin                 | OT_STRING | Linear solver plugin         |
+-------------------------------+-----------+------------------------------+
| maxDualJump                   | OT_DOUBLE | Maximum allowed jump in dual |
|                               |           | variables in linear          |
|                               |           | independence tests.          |
+-------------------------------+-----------+------------------------------+
| maxPrimalJump                 | OT_DOUBLE | Maximum allowed jump in      |
|                               |           | primal variables in nonzero  |
|                               |           | curvature tests.             |
+-------------------------------+-----------+------------------------------+
| max_schur                     | OT_INT    | Maximal number of Schur      |
|                               |           | updates [75]                 |
+-------------------------------+-----------+------------------------------+
| nWSR                          | OT_INT    | The maximum number of        |
|                               |           | working set recalculations   |
|                               |           | to be performed during the   |
|                               |           | initial homotopy. Default is |
|                               |           | 5(nx + nc)                   |
+-------------------------------+-----------+------------------------------+
| numRefinementSteps            | OT_INT    | Maximum number of iterative  |
|                               |           | refinement steps.            |
+-------------------------------+-----------+------------------------------+
| numRegularisationSteps        | OT_INT    | Maximum number of successive |
|                               |           | regularisation steps.        |
+-------------------------------+-----------+------------------------------+
| printLevel                    | OT_STRING | Defines the amount of text   |
|                               |           | output during QP solution,   |
|                               |           | see Section 5.7              |
+-------------------------------+-----------+------------------------------+
| schur                         | OT_BOOL   | Use Schur Complement         |
|                               |           | Approach [false]             |
+-------------------------------+-----------+------------------------------+
| sparse                        | OT_BOOL   | Formulate the QP using       |
|                               |           | sparse matrices. [false]     |
+-------------------------------+-----------+------------------------------+
| terminationTolerance          | OT_DOUBLE | Relative termination         |
|                               |           | tolerance to stop homotopy.  |
+-------------------------------+-----------+------------------------------+



--------------------------------------------------------------------------------

sqic
----



Interface to the SQIC solver for quadratic programming

Extra doc: https://github.com/casadi/casadi/wiki/L_21s



--------------------------------------------------------------------------------

superscs
--------



Interface to the SuperSCS solver for conic programming

Joris Gillis, 2019

Extra doc: https://github.com/casadi/casadi/wiki/L_21z

>List of available options

+----------+---------+-----------------------------------+
|    Id    |  Type   |            Description            |
+==========+=========+===================================+
| superscs | OT_DICT | Options to be passed to superscs. |
+----------+---------+-----------------------------------+



--------------------------------------------------------------------------------

ipqp
----



Solves QPs using a Mehrotra predictor-corrector interior point method

Extra doc: https://github.com/casadi/casadi/wiki/L_23c

>List of available options

+-----------------------+-----------+--------------------------------------+
|          Id           |   Type    |             Description              |
+=======================+===========+======================================+
| constr_viol_tol       | OT_DOUBLE | Constraint violation tolerance       |
|                       |           | [1e-8].                              |
+-----------------------+-----------+--------------------------------------+
| dual_inf_tol          | OT_DOUBLE | Dual feasibility violation tolerance |
|                       |           | [1e-8]                               |
+-----------------------+-----------+--------------------------------------+
| linear_solver         | OT_STRING | A custom linear solver creator       |
|                       |           | function [default: ldl]              |
+-----------------------+-----------+--------------------------------------+
| linear_solver_options | OT_DICT   | Options to be passed to the linear   |
|                       |           | solver                               |
+-----------------------+-----------+--------------------------------------+
| max_iter              | OT_INT    | Maximum number of iterations [1000]. |
+-----------------------+-----------+--------------------------------------+
| min_lam               | OT_DOUBLE | Smallest multiplier treated as       |
|                       |           | inactive for the initial active set  |
|                       |           | [0].                                 |
+-----------------------+-----------+--------------------------------------+
| print_header          | OT_BOOL   | Print header [true].                 |
+-----------------------+-----------+--------------------------------------+
| print_info            | OT_BOOL   | Print info [true].                   |
+-----------------------+-----------+--------------------------------------+
| print_iter            | OT_BOOL   | Print iterations [true].             |
+-----------------------+-----------+--------------------------------------+



--------------------------------------------------------------------------------

nlpsol
------



Solve QPs using an  Nlpsol Use the 'nlpsol' option to specify the NLP solver
 to use.

Extra doc: https://github.com/casadi/casadi/wiki/L_235

>List of available options

+----------------+-----------+---------------------------------+
|       Id       |   Type    |           Description           |
+================+===========+=================================+
| nlpsol         | OT_STRING | Name of solver.                 |
+----------------+-----------+---------------------------------+
| nlpsol_options | OT_DICT   | Options to be passed to solver. |
+----------------+-----------+---------------------------------+



--------------------------------------------------------------------------------

qrqp
----



Solve QPs using an active-set method

Extra doc: https://github.com/casadi/casadi/wiki/L_22y

>List of available options

+-----------------+-----------+--------------------------------------------+
|       Id        |   Type    |                Description                 |
+=================+===========+============================================+
| constr_viol_tol | OT_DOUBLE | Constraint violation tolerance [1e-8].     |
+-----------------+-----------+--------------------------------------------+
| dual_inf_tol    | OT_DOUBLE | Dual feasibility violation tolerance       |
|                 |           | [1e-8]                                     |
+-----------------+-----------+--------------------------------------------+
| max_iter        | OT_INT    | Maximum number of iterations [1000].       |
+-----------------+-----------+--------------------------------------------+
| min_lam         | OT_DOUBLE | Smallest multiplier treated as inactive    |
|                 |           | for the initial active set [0].            |
+-----------------+-----------+--------------------------------------------+
| print_header    | OT_BOOL   | Print header [true].                       |
+-----------------+-----------+--------------------------------------------+
| print_info      | OT_BOOL   | Print info [true].                         |
+-----------------+-----------+--------------------------------------------+
| print_iter      | OT_BOOL   | Print iterations [true].                   |
+-----------------+-----------+--------------------------------------------+
| print_lincomb   | OT_BOOL   | Print dependant linear combinations of     |
|                 |           | constraints [false]. Printed numbers are   |
|                 |           | 0-based indices into the vector of [simple |
|                 |           | bounds;linear bounds]                      |
+-----------------+-----------+--------------------------------------------+

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_21n

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L43

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L43-L46

";

%feature("docstring") casadi::IndexRecution::conic_debug "

[INTERNAL] 
Generate native code in the interfaced language for debugging

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L54

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L54-L58

>  void casadi::conic_debug(const Function &f, std::ostream &file)
------------------------------------------------------------------------
[INTERNAL] 
Generate native code in the interfaced language for debugging

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L54

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L54-L58

";

";

%feature("docstring") casadi::IndexRecution::conic_in "

[INTERNAL] 
Get QP solver input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1eg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L72

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L72-L89

>  std::string casadi::conic_in(casadi_int ind)
------------------------------------------------------------------------
[INTERNAL] 
Get QP solver input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1eg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L72

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L72-L89

";

";

%feature("docstring") casadi::IndexRecution::conic_out "

[INTERNAL] 
Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1eh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L91

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L91-L100

>  std::string casadi::conic_out(casadi_int ind)
------------------------------------------------------------------------
[INTERNAL] 
Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1eh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L91

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L91-L100

";

";

%feature("docstring") casadi::IndexRecution::conic_n_in "

[INTERNAL] 
Get the number of QP solver inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ei

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L102

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L102-L104

";

%feature("docstring") casadi::IndexRecution::conic_n_out "

[INTERNAL] 
Get the number of QP solver outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ej

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L106

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L106-L108

";

%feature("docstring") casadi::IndexRecution::qpsol "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::conic_options "

[INTERNAL] 
Get all options for a plugin.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ek

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L542

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L542-L544

";

%feature("docstring") casadi::IndexRecution::conic_option_type "

[INTERNAL] 
Get type info for a particular option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1el

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L546

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L546-L548

";

%feature("docstring") casadi::IndexRecution::conic_option_info "

[INTERNAL] 
Get documentation for a particular option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1em

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L550

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L550-L552

";

%feature("docstring") casadi::IndexRecution::has_dple "

[INTERNAL] 
Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L31

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L31-L33

";

%feature("docstring") casadi::IndexRecution::load_dple "

[INTERNAL] 
Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L35

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L35-L37

";

%feature("docstring") casadi::IndexRecution::doc_dple "

[INTERNAL] 
Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L39

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L39-L41

";

%feature("docstring") casadi::IndexRecution::dplesol "

[INTERNAL]

Discrete periodic Lyapunov Equation solver Given matrices  $A_k$ and 
symmetric  $V_k, k = 0..K-1$

::

  A_k in R^(n x n)
  V_k in R^n
  

provides all of  $P_k$ that satisfy:

::

  P_0 = A_(K-1)*P_(K-1)*A_(K-1)' + V_k
  P_k+1 = A_k*P_k*A_k' + V_k  for k = 1..K-1
  

General information

>List of available options

+------------------+-----------------+------------------+------------------+
|        Id        |      Type       |   Description    |     Used in      |
+==================+=================+==================+==================+
| ad_weight        | OT_DOUBLE       | Weighting factor | casadi::Function |
|                  |                 | for derivative   | Internal         |
|                  |                 | calculation.When |                  |
|                  |                 | there is an      |                  |
|                  |                 | option of either |                  |
|                  |                 | using forward or |                  |
|                  |                 | reverse mode     |                  |
|                  |                 | directional      |                  |
|                  |                 | derivatives, the |                  |
|                  |                 | condition ad_wei |                  |
|                  |                 | ght*nf<=(1-      |                  |
|                  |                 | ad_weight)*na is |                  |
|                  |                 | used where nf    |                  |
|                  |                 | and na are       |                  |
|                  |                 | estimates of the |                  |
|                  |                 | number of        |                  |
|                  |                 | forward/reverse  |                  |
|                  |                 | mode directional |                  |
|                  |                 | derivatives      |                  |
|                  |                 | needed. By       |                  |
|                  |                 | default,         |                  |
|                  |                 | ad_weight is     |                  |
|                  |                 | calculated       |                  |
|                  |                 | automatically,   |                  |
|                  |                 | but this can be  |                  |
|                  |                 | overridden by    |                  |
|                  |                 | setting this     |                  |
|                  |                 | option. In       |                  |
|                  |                 | particular, 0    |                  |
|                  |                 | means forcing    |                  |
|                  |                 | forward mode and |                  |
|                  |                 | 1 forcing        |                  |
|                  |                 | reverse mode.    |                  |
|                  |                 | Leave unset for  |                  |
|                  |                 | (class specific) |                  |
|                  |                 | heuristics.      |                  |
+------------------+-----------------+------------------+------------------+
| ad_weight_sp     | OT_DOUBLE       | Weighting factor | casadi::Function |
|                  |                 | for sparsity     | Internal         |
|                  |                 | pattern          |                  |
|                  |                 | calculation calc |                  |
|                  |                 | ulation.Override |                  |
|                  |                 | s default        |                  |
|                  |                 | behavior. Set to |                  |
|                  |                 | 0 and 1 to force |                  |
|                  |                 | forward and      |                  |
|                  |                 | reverse mode     |                  |
|                  |                 | respectively.    |                  |
|                  |                 | Cf. option       |                  |
|                  |                 | \"ad_weight\".     |                  |
|                  |                 | When set to -1,  |                  |
|                  |                 | sparsity is      |                  |
|                  |                 | completely       |                  |
|                  |                 | ignored and      |                  |
|                  |                 | dense matrices   |                  |
|                  |                 | are used.        |                  |
+------------------+-----------------+------------------+------------------+
| always_inline    | OT_BOOL         | Force inlining.  | casadi::Function |
|                  |                 |                  | Internal         |
+------------------+-----------------+------------------+------------------+
| cache            | OT_DICT         | Prepopulate the  | casadi::Function |
|                  |                 | function cache.  | Internal         |
|                  |                 | Default: empty   |                  |
+------------------+-----------------+------------------+------------------+
| compiler         | OT_STRING       | Just-in-time     | casadi::Function |
|                  |                 | compiler plugin  | Internal         |
|                  |                 | to be used.      |                  |
+------------------+-----------------+------------------+------------------+
| const_dim        | OT_BOOL         | Assume constant  | casadi::Dple     |
|                  |                 | dimension of P   |                  |
+------------------+-----------------+------------------+------------------+
| custom_jacobian  | OT_FUNCTION     | Override         | casadi::Function |
|                  |                 | CasADi's AD. Use | Internal         |
|                  |                 | together with    |                  |
|                  |                 | 'jac_penalty':   |                  |
|                  |                 | 0. Note: Highly  |                  |
|                  |                 | experimental.    |                  |
|                  |                 | Syntax may break |                  |
|                  |                 | often.           |                  |
+------------------+-----------------+------------------+------------------+
| der_options      | OT_DICT         | Default options  | casadi::Function |
|                  |                 | to be used to    | Internal         |
|                  |                 | populate         |                  |
|                  |                 | forward_options, |                  |
|                  |                 | reverse_options, |                  |
|                  |                 | and              |                  |
|                  |                 | jacobian_options |                  |
|                  |                 | before those     |                  |
|                  |                 | options are      |                  |
|                  |                 | merged in.       |                  |
+------------------+-----------------+------------------+------------------+
| derivative_of    | OT_FUNCTION     | The function is  | casadi::Function |
|                  |                 | a derivative of  | Internal         |
|                  |                 | another          |                  |
|                  |                 | function. The    |                  |
|                  |                 | type of          |                  |
|                  |                 | derivative       |                  |
|                  |                 | (directional     |                  |
|                  |                 | derivative,      |                  |
|                  |                 | Jacobian) is     |                  |
|                  |                 | inferred from    |                  |
|                  |                 | the function     |                  |
|                  |                 | name.            |                  |
+------------------+-----------------+------------------+------------------+
| dump             | OT_BOOL         | Dump function to | casadi::Function |
|                  |                 | file upon first  | Internal         |
|                  |                 | evaluation.      |                  |
|                  |                 | [false]          |                  |
+------------------+-----------------+------------------+------------------+
| dump_dir         | OT_STRING       | Directory to     | casadi::Function |
|                  |                 | dump             | Internal         |
|                  |                 | inputs/outputs   |                  |
|                  |                 | to. Make sure    |                  |
|                  |                 | the directory    |                  |
|                  |                 | exists [.]       |                  |
+------------------+-----------------+------------------+------------------+
| dump_format      | OT_STRING       | Choose file      | casadi::Function |
|                  |                 | format to dump   | Internal         |
|                  |                 | matrices. See    |                  |
|                  |                 | DM.from_file     |                  |
|                  |                 | [mtx]            |                  |
+------------------+-----------------+------------------+------------------+
| dump_in          | OT_BOOL         | Dump numerical   | casadi::Function |
|                  |                 | values of inputs | Internal         |
|                  |                 | to file          |                  |
|                  |                 | (readable with   |                  |
|                  |                 | DM.from_file )   |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| dump_out         | OT_BOOL         | Dump numerical   | casadi::Function |
|                  |                 | values of        | Internal         |
|                  |                 | outputs to file  |                  |
|                  |                 | (readable with   |                  |
|                  |                 | DM.from_file )   |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| enable_fd        | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation by   |                  |
|                  |                 | finite           |                  |
|                  |                 | differencing.    |                  |
|                  |                 | [default:        |                  |
|                  |                 | false]]          |                  |
+------------------+-----------------+------------------+------------------+
| enable_forward   | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | Jacobian-times-  |                  |
|                  |                 | vector products  |                  |
|                  |                 | - typically      |                  |
|                  |                 | using forward    |                  |
|                  |                 | mode AD - if     |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| enable_jacobian  | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | Jacobians of all |                  |
|                  |                 | differentiable   |                  |
|                  |                 | outputs with     |                  |
|                  |                 | respect to all   |                  |
|                  |                 | differentiable   |                  |
|                  |                 | inputs - if      |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| enable_reverse   | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | transposed       |                  |
|                  |                 | Jacobian-times-  |                  |
|                  |                 | vector products  |                  |
|                  |                 | - typically      |                  |
|                  |                 | using reverse    |                  |
|                  |                 | mode AD - if     |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| eps_unstable     | OT_DOUBLE       | A margin for     | casadi::Dple     |
|                  |                 | unstability      |                  |
|                  |                 | detection        |                  |
+------------------+-----------------+------------------+------------------+
| error_on_fail    | OT_BOOL         | Throw exceptions | casadi::ProtoFun |
|                  |                 | when function    | ction            |
|                  |                 | evaluation fails |                  |
|                  |                 | (default true).  |                  |
+------------------+-----------------+------------------+------------------+
| error_unstable   | OT_BOOL         | Throw an         | casadi::Dple     |
|                  |                 | exception when   |                  |
|                  |                 | it is detected   |                  |
|                  |                 | that             |                  |
|                  |                 | Product(A_i,     |                  |
|                  |                 | i=N..1)has       |                  |
|                  |                 | eigenvalues      |                  |
|                  |                 | greater than     |                  |
|                  |                 | 1-eps_unstable   |                  |
+------------------+-----------------+------------------+------------------+
| fd_method        | OT_STRING       | Method for       | casadi::Function |
|                  |                 | finite           | Internal         |
|                  |                 | differencing     |                  |
|                  |                 | [default         |                  |
|                  |                 | 'central']       |                  |
+------------------+-----------------+------------------+------------------+
| fd_options       | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to the    | Internal         |
|                  |                 | finite           |                  |
|                  |                 | difference       |                  |
|                  |                 | instance         |                  |
+------------------+-----------------+------------------+------------------+
| forward_options  | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | forward mode     |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| gather_stats     | OT_BOOL         | Deprecated       | casadi::Function |
|                  |                 | option           | Internal         |
|                  |                 | (ignored):       |                  |
|                  |                 | Statistics are   |                  |
|                  |                 | now always       |                  |
|                  |                 | collected.       |                  |
+------------------+-----------------+------------------+------------------+
| input_scheme     | OT_STRINGVECTOR | Deprecated       | casadi::Function |
|                  |                 | option (ignored) | Internal         |
+------------------+-----------------+------------------+------------------+
| inputs_check     | OT_BOOL         | Throw exceptions | casadi::Function |
|                  |                 | when the         | Internal         |
|                  |                 | numerical values |                  |
|                  |                 | of the inputs    |                  |
|                  |                 | don't make sense |                  |
+------------------+-----------------+------------------+------------------+
| is_diff_in       | OT_BOOLVECTOR   | Indicate for     | casadi::Function |
|                  |                 | each input if it | Internal         |
|                  |                 | should be        |                  |
|                  |                 | differentiable.  |                  |
+------------------+-----------------+------------------+------------------+
| is_diff_out      | OT_BOOLVECTOR   | Indicate for     | casadi::Function |
|                  |                 | each output if   | Internal         |
|                  |                 | it should be     |                  |
|                  |                 | differentiable.  |                  |
+------------------+-----------------+------------------+------------------+
| jac_penalty      | OT_DOUBLE       | When requested   | casadi::Function |
|                  |                 | for a number of  | Internal         |
|                  |                 | forward/reverse  |                  |
|                  |                 | directions, it   |                  |
|                  |                 | may be cheaper   |                  |
|                  |                 | to compute first |                  |
|                  |                 | the full         |                  |
|                  |                 | jacobian and     |                  |
|                  |                 | then multiply    |                  |
|                  |                 | with seeds,      |                  |
|                  |                 | rather than      |                  |
|                  |                 | obtain the       |                  |
|                  |                 | requested        |                  |
|                  |                 | directions in a  |                  |
|                  |                 | straightforward  |                  |
|                  |                 | manner. Casadi   |                  |
|                  |                 | uses a heuristic |                  |
|                  |                 | to decide which  |                  |
|                  |                 | is cheaper. A    |                  |
|                  |                 | high value of    |                  |
|                  |                 | 'jac_penalty'    |                  |
|                  |                 | makes it less    |                  |
|                  |                 | likely for the   |                  |
|                  |                 | heurstic to      |                  |
|                  |                 | chose the full   |                  |
|                  |                 | Jacobian         |                  |
|                  |                 | strategy. The    |                  |
|                  |                 | special value -1 |                  |
|                  |                 | indicates never  |                  |
|                  |                 | to use the full  |                  |
|                  |                 | Jacobian         |                  |
|                  |                 | strategy         |                  |
+------------------+-----------------+------------------+------------------+
| jacobian_options | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | Jacobian         |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| jit              | OT_BOOL         | Use just-in-time | casadi::Function |
|                  |                 | compiler to      | Internal         |
|                  |                 | speed up the     |                  |
|                  |                 | evaluation       |                  |
+------------------+-----------------+------------------+------------------+
| jit_cleanup      | OT_BOOL         | Cleanup up the   | casadi::Function |
|                  |                 | temporary source | Internal         |
|                  |                 | file that jit    |                  |
|                  |                 | creates.         |                  |
|                  |                 | Default: true    |                  |
+------------------+-----------------+------------------+------------------+
| jit_name         | OT_STRING       | The file name    | casadi::Function |
|                  |                 | used to write    | Internal         |
|                  |                 | out code. The    |                  |
|                  |                 | actual file      |                  |
|                  |                 | names used       |                  |
|                  |                 | depend on 'jit_t |                  |
|                  |                 | emp_suffix' and  |                  |
|                  |                 | include          |                  |
|                  |                 | extensions.      |                  |
|                  |                 | Default:         |                  |
|                  |                 | 'jit_tmp'        |                  |
+------------------+-----------------+------------------+------------------+
| jit_options      | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to the    | Internal         |
|                  |                 | jit compiler.    |                  |
+------------------+-----------------+------------------+------------------+
| jit_serialize    | OT_STRING       | Specify          | casadi::Function |
|                  |                 | behaviour when   | Internal         |
|                  |                 | serializing a    |                  |
|                  |                 | jitted function: |                  |
|                  |                 | SOURCE|link|embe |                  |
|                  |                 | d.               |                  |
+------------------+-----------------+------------------+------------------+
| jit_temp_suffix  | OT_BOOL         | Use a temporary  | casadi::Function |
|                  |                 | (seemingly       | Internal         |
|                  |                 | random) filename |                  |
|                  |                 | suffix for       |                  |
|                  |                 | generated code   |                  |
|                  |                 | and libraries.   |                  |
|                  |                 | This is desired  |                  |
|                  |                 | for thread-      |                  |
|                  |                 | safety. This     |                  |
|                  |                 | behaviour may    |                  |
|                  |                 | defeat caching   |                  |
|                  |                 | compiler         |                  |
|                  |                 | wrappers.        |                  |
|                  |                 | Default: true    |                  |
+------------------+-----------------+------------------+------------------+
| max_io           | OT_INT          | Acceptable       | casadi::Function |
|                  |                 | number of inputs | Internal         |
|                  |                 | and outputs.     |                  |
|                  |                 | Warn if          |                  |
|                  |                 | exceeded.        |                  |
+------------------+-----------------+------------------+------------------+
| max_num_dir      | OT_INT          | Specify the      | casadi::Function |
|                  |                 | maximum number   | Internal         |
|                  |                 | of directions    |                  |
|                  |                 | for derivative   |                  |
|                  |                 | functions.       |                  |
|                  |                 | Overrules the    |                  |
|                  |                 | builtin optimize |                  |
|                  |                 | d_num_dir.       |                  |
+------------------+-----------------+------------------+------------------+
| never_inline     | OT_BOOL         | Forbid inlining. | casadi::Function |
|                  |                 |                  | Internal         |
+------------------+-----------------+------------------+------------------+
| output_scheme    | OT_STRINGVECTOR | Deprecated       | casadi::Function |
|                  |                 | option (ignored) | Internal         |
+------------------+-----------------+------------------+------------------+
| pos_def          | OT_BOOL         | Assume P         | casadi::Dple     |
|                  |                 | positive         |                  |
|                  |                 | definite         |                  |
+------------------+-----------------+------------------+------------------+
| post_expand      | OT_BOOL         | After            | casadi::Function |
|                  |                 | construction,    | Internal         |
|                  |                 | expand this      |                  |
|                  |                 | Function .       |                  |
|                  |                 | Default: False   |                  |
+------------------+-----------------+------------------+------------------+
| post_expand_opti | OT_DICT         | Options to be    | casadi::Function |
| ons              |                 | passed to post-  | Internal         |
|                  |                 | construction     |                  |
|                  |                 | expansion.       |                  |
|                  |                 | Default: empty   |                  |
+------------------+-----------------+------------------+------------------+
| print_in         | OT_BOOL         | Print numerical  | casadi::Function |
|                  |                 | values of inputs | Internal         |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| print_out        | OT_BOOL         | Print numerical  | casadi::Function |
|                  |                 | values of        | Internal         |
|                  |                 | outputs          |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| print_time       | OT_BOOL         | print            | casadi::ProtoFun |
|                  |                 | information      | ction            |
|                  |                 | about execution  |                  |
|                  |                 | time. Implies    |                  |
|                  |                 | record_time.     |                  |
+------------------+-----------------+------------------+------------------+
| record_time      | OT_BOOL         | record           | casadi::ProtoFun |
|                  |                 | information      | ction            |
|                  |                 | about execution  |                  |
|                  |                 | time, for        |                  |
|                  |                 | retrieval with   |                  |
|                  |                 | stats().         |                  |
+------------------+-----------------+------------------+------------------+
| regularity_check | OT_BOOL         | Throw exceptions | casadi::ProtoFun |
|                  |                 | when NaN or Inf  | ction            |
|                  |                 | appears during   |                  |
|                  |                 | evaluation       |                  |
+------------------+-----------------+------------------+------------------+
| reverse_options  | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | reverse mode     |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| user_data        | OT_VOIDPTR      | A user-defined   | casadi::Function |
|                  |                 | field that can   | Internal         |
|                  |                 | be used to       |                  |
|                  |                 | identify the     |                  |
|                  |                 | function or pass |                  |
|                  |                 | additional       |                  |
|                  |                 | information      |                  |
+------------------+-----------------+------------------+------------------+
| verbose          | OT_BOOL         | Verbose          | casadi::ProtoFun |
|                  |                 | evaluation  for  | ction            |
|                  |                 | debugging        |                  |
+------------------+-----------------+------------------+------------------+

>Input scheme: casadi::DpleInput (DPLE_NUM_IN = 2)

+-----------+-------+------------------------------------------------------+
| Full name | Short |                     Description                      |
+===========+=======+======================================================+
| DPLE_A    | a     | A matrices (horzcat when const_dim, diagcat          |
|           |       | otherwise) [a].                                      |
+-----------+-------+------------------------------------------------------+
| DPLE_V    | v     | V matrices (horzcat when const_dim, diagcat          |
|           |       | otherwise) [v].                                      |
+-----------+-------+------------------------------------------------------+

>Output scheme: casadi::DpleOutput (DPLE_NUM_OUT = 1)

+-----------+-------+------------------------------------------------------+
| Full name | Short |                     Description                      |
+===========+=======+======================================================+
| DPLE_P    | p     | Lyapunov matrix (horzcat when const_dim, diagcat     |
|           |       | otherwise) (Cholesky of P if pos_def) [p].           |
+-----------+-------+------------------------------------------------------+

List of plugins
- slicot

Note: some of the plugins in this list might not be available on your 

system.  Also, there might be extra plugins available to you that are 
not 
listed here. You can obtain their documentation with  
Dple.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

slicot
------



An efficient solver for Discrete Periodic Lyapunov Equations using 
SLICOT

Uses Periodic Schur Decomposition ('psd') and does not assume positive
 
definiteness. Based on Periodic Lyapunov equations: some applications
 and 
new algorithms. Int. J. Control, vol. 67, pp. 69-87, 1997.

Overview of the method: J. Gillis Practical Methods for Approximate 
Robust 
Periodic Optimal Control ofNonlinear Mechanical Systems, PhD 
Thesis, 
KULeuven, 2015

Extra doc: https://github.com/casadi/casadi/wiki/L_22j

>List of available options

+-----------------------+-----------+--------------------------------------+
|          Id           |   Type    |             Description              |
+=======================+===========+======================================+
| linear_solver         | OT_STRING | User-defined linear solver class.    |
|                       |           | Needed for sensitivities.            |
+-----------------------+-----------+--------------------------------------+
| linear_solver_options | OT_DICT   | Options to be passed to the linear   |
|                       |           | solver.                              |
+-----------------------+-----------+--------------------------------------+
| psd_num_zero          | OT_DOUBLE | Numerical zero used in Periodic      |
|                       |           | Schur decomposition with slicot.This |
|                       |           | option is needed when your systems   |
|                       |           | has Floquet multiplierszero or close |
|                       |           | to zero                              |
+-----------------------+-----------+--------------------------------------+

Joris Gillis

Extra doc: https://github.com/casadi/casadi/wiki/L_21o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L97

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L97-L100

>  Function casadi::dplesol(const std::string &name, const std::string &solver, const SpDict &st, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL]

Discrete periodic Lyapunov Equation solver Given matrices  $A_k$ and 
symmetric  $V_k, k = 0..K-1$

::

  A_k in R^(n x n)
  V_k in R^n
  

provides all of  $P_k$ that satisfy:

::

  P_0 = A_(K-1)*P_(K-1)*A_(K-1)' + V_k
  P_k+1 = A_k*P_k*A_k' + V_k  for k = 1..K-1
  

General information

>List of available options

+------------------+-----------------+------------------+------------------+
|        Id        |      Type       |   Description    |     Used in      |
+==================+=================+==================+==================+
| ad_weight        | OT_DOUBLE       | Weighting factor | casadi::Function |
|                  |                 | for derivative   | Internal         |
|                  |                 | calculation.When |                  |
|                  |                 | there is an      |                  |
|                  |                 | option of either |                  |
|                  |                 | using forward or |                  |
|                  |                 | reverse mode     |                  |
|                  |                 | directional      |                  |
|                  |                 | derivatives, the |                  |
|                  |                 | condition ad_wei |                  |
|                  |                 | ght*nf<=(1-      |                  |
|                  |                 | ad_weight)*na is |                  |
|                  |                 | used where nf    |                  |
|                  |                 | and na are       |                  |
|                  |                 | estimates of the |                  |
|                  |                 | number of        |                  |
|                  |                 | forward/reverse  |                  |
|                  |                 | mode directional |                  |
|                  |                 | derivatives      |                  |
|                  |                 | needed. By       |                  |
|                  |                 | default,         |                  |
|                  |                 | ad_weight is     |                  |
|                  |                 | calculated       |                  |
|                  |                 | automatically,   |                  |
|                  |                 | but this can be  |                  |
|                  |                 | overridden by    |                  |
|                  |                 | setting this     |                  |
|                  |                 | option. In       |                  |
|                  |                 | particular, 0    |                  |
|                  |                 | means forcing    |                  |
|                  |                 | forward mode and |                  |
|                  |                 | 1 forcing        |                  |
|                  |                 | reverse mode.    |                  |
|                  |                 | Leave unset for  |                  |
|                  |                 | (class specific) |                  |
|                  |                 | heuristics.      |                  |
+------------------+-----------------+------------------+------------------+
| ad_weight_sp     | OT_DOUBLE       | Weighting factor | casadi::Function |
|                  |                 | for sparsity     | Internal         |
|                  |                 | pattern          |                  |
|                  |                 | calculation calc |                  |
|                  |                 | ulation.Override |                  |
|                  |                 | s default        |                  |
|                  |                 | behavior. Set to |                  |
|                  |                 | 0 and 1 to force |                  |
|                  |                 | forward and      |                  |
|                  |                 | reverse mode     |                  |
|                  |                 | respectively.    |                  |
|                  |                 | Cf. option       |                  |
|                  |                 | \"ad_weight\".     |                  |
|                  |                 | When set to -1,  |                  |
|                  |                 | sparsity is      |                  |
|                  |                 | completely       |                  |
|                  |                 | ignored and      |                  |
|                  |                 | dense matrices   |                  |
|                  |                 | are used.        |                  |
+------------------+-----------------+------------------+------------------+
| always_inline    | OT_BOOL         | Force inlining.  | casadi::Function |
|                  |                 |                  | Internal         |
+------------------+-----------------+------------------+------------------+
| cache            | OT_DICT         | Prepopulate the  | casadi::Function |
|                  |                 | function cache.  | Internal         |
|                  |                 | Default: empty   |                  |
+------------------+-----------------+------------------+------------------+
| compiler         | OT_STRING       | Just-in-time     | casadi::Function |
|                  |                 | compiler plugin  | Internal         |
|                  |                 | to be used.      |                  |
+------------------+-----------------+------------------+------------------+
| const_dim        | OT_BOOL         | Assume constant  | casadi::Dple     |
|                  |                 | dimension of P   |                  |
+------------------+-----------------+------------------+------------------+
| custom_jacobian  | OT_FUNCTION     | Override         | casadi::Function |
|                  |                 | CasADi's AD. Use | Internal         |
|                  |                 | together with    |                  |
|                  |                 | 'jac_penalty':   |                  |
|                  |                 | 0. Note: Highly  |                  |
|                  |                 | experimental.    |                  |
|                  |                 | Syntax may break |                  |
|                  |                 | often.           |                  |
+------------------+-----------------+------------------+------------------+
| der_options      | OT_DICT         | Default options  | casadi::Function |
|                  |                 | to be used to    | Internal         |
|                  |                 | populate         |                  |
|                  |                 | forward_options, |                  |
|                  |                 | reverse_options, |                  |
|                  |                 | and              |                  |
|                  |                 | jacobian_options |                  |
|                  |                 | before those     |                  |
|                  |                 | options are      |                  |
|                  |                 | merged in.       |                  |
+------------------+-----------------+------------------+------------------+
| derivative_of    | OT_FUNCTION     | The function is  | casadi::Function |
|                  |                 | a derivative of  | Internal         |
|                  |                 | another          |                  |
|                  |                 | function. The    |                  |
|                  |                 | type of          |                  |
|                  |                 | derivative       |                  |
|                  |                 | (directional     |                  |
|                  |                 | derivative,      |                  |
|                  |                 | Jacobian) is     |                  |
|                  |                 | inferred from    |                  |
|                  |                 | the function     |                  |
|                  |                 | name.            |                  |
+------------------+-----------------+------------------+------------------+
| dump             | OT_BOOL         | Dump function to | casadi::Function |
|                  |                 | file upon first  | Internal         |
|                  |                 | evaluation.      |                  |
|                  |                 | [false]          |                  |
+------------------+-----------------+------------------+------------------+
| dump_dir         | OT_STRING       | Directory to     | casadi::Function |
|                  |                 | dump             | Internal         |
|                  |                 | inputs/outputs   |                  |
|                  |                 | to. Make sure    |                  |
|                  |                 | the directory    |                  |
|                  |                 | exists [.]       |                  |
+------------------+-----------------+------------------+------------------+
| dump_format      | OT_STRING       | Choose file      | casadi::Function |
|                  |                 | format to dump   | Internal         |
|                  |                 | matrices. See    |                  |
|                  |                 | DM.from_file     |                  |
|                  |                 | [mtx]            |                  |
+------------------+-----------------+------------------+------------------+
| dump_in          | OT_BOOL         | Dump numerical   | casadi::Function |
|                  |                 | values of inputs | Internal         |
|                  |                 | to file          |                  |
|                  |                 | (readable with   |                  |
|                  |                 | DM.from_file )   |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| dump_out         | OT_BOOL         | Dump numerical   | casadi::Function |
|                  |                 | values of        | Internal         |
|                  |                 | outputs to file  |                  |
|                  |                 | (readable with   |                  |
|                  |                 | DM.from_file )   |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| enable_fd        | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation by   |                  |
|                  |                 | finite           |                  |
|                  |                 | differencing.    |                  |
|                  |                 | [default:        |                  |
|                  |                 | false]]          |                  |
+------------------+-----------------+------------------+------------------+
| enable_forward   | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | Jacobian-times-  |                  |
|                  |                 | vector products  |                  |
|                  |                 | - typically      |                  |
|                  |                 | using forward    |                  |
|                  |                 | mode AD - if     |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| enable_jacobian  | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | Jacobians of all |                  |
|                  |                 | differentiable   |                  |
|                  |                 | outputs with     |                  |
|                  |                 | respect to all   |                  |
|                  |                 | differentiable   |                  |
|                  |                 | inputs - if      |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| enable_reverse   | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | transposed       |                  |
|                  |                 | Jacobian-times-  |                  |
|                  |                 | vector products  |                  |
|                  |                 | - typically      |                  |
|                  |                 | using reverse    |                  |
|                  |                 | mode AD - if     |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| eps_unstable     | OT_DOUBLE       | A margin for     | casadi::Dple     |
|                  |                 | unstability      |                  |
|                  |                 | detection        |                  |
+------------------+-----------------+------------------+------------------+
| error_on_fail    | OT_BOOL         | Throw exceptions | casadi::ProtoFun |
|                  |                 | when function    | ction            |
|                  |                 | evaluation fails |                  |
|                  |                 | (default true).  |                  |
+------------------+-----------------+------------------+------------------+
| error_unstable   | OT_BOOL         | Throw an         | casadi::Dple     |
|                  |                 | exception when   |                  |
|                  |                 | it is detected   |                  |
|                  |                 | that             |                  |
|                  |                 | Product(A_i,     |                  |
|                  |                 | i=N..1)has       |                  |
|                  |                 | eigenvalues      |                  |
|                  |                 | greater than     |                  |
|                  |                 | 1-eps_unstable   |                  |
+------------------+-----------------+------------------+------------------+
| fd_method        | OT_STRING       | Method for       | casadi::Function |
|                  |                 | finite           | Internal         |
|                  |                 | differencing     |                  |
|                  |                 | [default         |                  |
|                  |                 | 'central']       |                  |
+------------------+-----------------+------------------+------------------+
| fd_options       | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to the    | Internal         |
|                  |                 | finite           |                  |
|                  |                 | difference       |                  |
|                  |                 | instance         |                  |
+------------------+-----------------+------------------+------------------+
| forward_options  | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | forward mode     |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| gather_stats     | OT_BOOL         | Deprecated       | casadi::Function |
|                  |                 | option           | Internal         |
|                  |                 | (ignored):       |                  |
|                  |                 | Statistics are   |                  |
|                  |                 | now always       |                  |
|                  |                 | collected.       |                  |
+------------------+-----------------+------------------+------------------+
| input_scheme     | OT_STRINGVECTOR | Deprecated       | casadi::Function |
|                  |                 | option (ignored) | Internal         |
+------------------+-----------------+------------------+------------------+
| inputs_check     | OT_BOOL         | Throw exceptions | casadi::Function |
|                  |                 | when the         | Internal         |
|                  |                 | numerical values |                  |
|                  |                 | of the inputs    |                  |
|                  |                 | don't make sense |                  |
+------------------+-----------------+------------------+------------------+
| is_diff_in       | OT_BOOLVECTOR   | Indicate for     | casadi::Function |
|                  |                 | each input if it | Internal         |
|                  |                 | should be        |                  |
|                  |                 | differentiable.  |                  |
+------------------+-----------------+------------------+------------------+
| is_diff_out      | OT_BOOLVECTOR   | Indicate for     | casadi::Function |
|                  |                 | each output if   | Internal         |
|                  |                 | it should be     |                  |
|                  |                 | differentiable.  |                  |
+------------------+-----------------+------------------+------------------+
| jac_penalty      | OT_DOUBLE       | When requested   | casadi::Function |
|                  |                 | for a number of  | Internal         |
|                  |                 | forward/reverse  |                  |
|                  |                 | directions, it   |                  |
|                  |                 | may be cheaper   |                  |
|                  |                 | to compute first |                  |
|                  |                 | the full         |                  |
|                  |                 | jacobian and     |                  |
|                  |                 | then multiply    |                  |
|                  |                 | with seeds,      |                  |
|                  |                 | rather than      |                  |
|                  |                 | obtain the       |                  |
|                  |                 | requested        |                  |
|                  |                 | directions in a  |                  |
|                  |                 | straightforward  |                  |
|                  |                 | manner. Casadi   |                  |
|                  |                 | uses a heuristic |                  |
|                  |                 | to decide which  |                  |
|                  |                 | is cheaper. A    |                  |
|                  |                 | high value of    |                  |
|                  |                 | 'jac_penalty'    |                  |
|                  |                 | makes it less    |                  |
|                  |                 | likely for the   |                  |
|                  |                 | heurstic to      |                  |
|                  |                 | chose the full   |                  |
|                  |                 | Jacobian         |                  |
|                  |                 | strategy. The    |                  |
|                  |                 | special value -1 |                  |
|                  |                 | indicates never  |                  |
|                  |                 | to use the full  |                  |
|                  |                 | Jacobian         |                  |
|                  |                 | strategy         |                  |
+------------------+-----------------+------------------+------------------+
| jacobian_options | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | Jacobian         |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| jit              | OT_BOOL         | Use just-in-time | casadi::Function |
|                  |                 | compiler to      | Internal         |
|                  |                 | speed up the     |                  |
|                  |                 | evaluation       |                  |
+------------------+-----------------+------------------+------------------+
| jit_cleanup      | OT_BOOL         | Cleanup up the   | casadi::Function |
|                  |                 | temporary source | Internal         |
|                  |                 | file that jit    |                  |
|                  |                 | creates.         |                  |
|                  |                 | Default: true    |                  |
+------------------+-----------------+------------------+------------------+
| jit_name         | OT_STRING       | The file name    | casadi::Function |
|                  |                 | used to write    | Internal         |
|                  |                 | out code. The    |                  |
|                  |                 | actual file      |                  |
|                  |                 | names used       |                  |
|                  |                 | depend on 'jit_t |                  |
|                  |                 | emp_suffix' and  |                  |
|                  |                 | include          |                  |
|                  |                 | extensions.      |                  |
|                  |                 | Default:         |                  |
|                  |                 | 'jit_tmp'        |                  |
+------------------+-----------------+------------------+------------------+
| jit_options      | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to the    | Internal         |
|                  |                 | jit compiler.    |                  |
+------------------+-----------------+------------------+------------------+
| jit_serialize    | OT_STRING       | Specify          | casadi::Function |
|                  |                 | behaviour when   | Internal         |
|                  |                 | serializing a    |                  |
|                  |                 | jitted function: |                  |
|                  |                 | SOURCE|link|embe |                  |
|                  |                 | d.               |                  |
+------------------+-----------------+------------------+------------------+
| jit_temp_suffix  | OT_BOOL         | Use a temporary  | casadi::Function |
|                  |                 | (seemingly       | Internal         |
|                  |                 | random) filename |                  |
|                  |                 | suffix for       |                  |
|                  |                 | generated code   |                  |
|                  |                 | and libraries.   |                  |
|                  |                 | This is desired  |                  |
|                  |                 | for thread-      |                  |
|                  |                 | safety. This     |                  |
|                  |                 | behaviour may    |                  |
|                  |                 | defeat caching   |                  |
|                  |                 | compiler         |                  |
|                  |                 | wrappers.        |                  |
|                  |                 | Default: true    |                  |
+------------------+-----------------+------------------+------------------+
| max_io           | OT_INT          | Acceptable       | casadi::Function |
|                  |                 | number of inputs | Internal         |
|                  |                 | and outputs.     |                  |
|                  |                 | Warn if          |                  |
|                  |                 | exceeded.        |                  |
+------------------+-----------------+------------------+------------------+
| max_num_dir      | OT_INT          | Specify the      | casadi::Function |
|                  |                 | maximum number   | Internal         |
|                  |                 | of directions    |                  |
|                  |                 | for derivative   |                  |
|                  |                 | functions.       |                  |
|                  |                 | Overrules the    |                  |
|                  |                 | builtin optimize |                  |
|                  |                 | d_num_dir.       |                  |
+------------------+-----------------+------------------+------------------+
| never_inline     | OT_BOOL         | Forbid inlining. | casadi::Function |
|                  |                 |                  | Internal         |
+------------------+-----------------+------------------+------------------+
| output_scheme    | OT_STRINGVECTOR | Deprecated       | casadi::Function |
|                  |                 | option (ignored) | Internal         |
+------------------+-----------------+------------------+------------------+
| pos_def          | OT_BOOL         | Assume P         | casadi::Dple     |
|                  |                 | positive         |                  |
|                  |                 | definite         |                  |
+------------------+-----------------+------------------+------------------+
| post_expand      | OT_BOOL         | After            | casadi::Function |
|                  |                 | construction,    | Internal         |
|                  |                 | expand this      |                  |
|                  |                 | Function .       |                  |
|                  |                 | Default: False   |                  |
+------------------+-----------------+------------------+------------------+
| post_expand_opti | OT_DICT         | Options to be    | casadi::Function |
| ons              |                 | passed to post-  | Internal         |
|                  |                 | construction     |                  |
|                  |                 | expansion.       |                  |
|                  |                 | Default: empty   |                  |
+------------------+-----------------+------------------+------------------+
| print_in         | OT_BOOL         | Print numerical  | casadi::Function |
|                  |                 | values of inputs | Internal         |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| print_out        | OT_BOOL         | Print numerical  | casadi::Function |
|                  |                 | values of        | Internal         |
|                  |                 | outputs          |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| print_time       | OT_BOOL         | print            | casadi::ProtoFun |
|                  |                 | information      | ction            |
|                  |                 | about execution  |                  |
|                  |                 | time. Implies    |                  |
|                  |                 | record_time.     |                  |
+------------------+-----------------+------------------+------------------+
| record_time      | OT_BOOL         | record           | casadi::ProtoFun |
|                  |                 | information      | ction            |
|                  |                 | about execution  |                  |
|                  |                 | time, for        |                  |
|                  |                 | retrieval with   |                  |
|                  |                 | stats().         |                  |
+------------------+-----------------+------------------+------------------+
| regularity_check | OT_BOOL         | Throw exceptions | casadi::ProtoFun |
|                  |                 | when NaN or Inf  | ction            |
|                  |                 | appears during   |                  |
|                  |                 | evaluation       |                  |
+------------------+-----------------+------------------+------------------+
| reverse_options  | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | reverse mode     |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| user_data        | OT_VOIDPTR      | A user-defined   | casadi::Function |
|                  |                 | field that can   | Internal         |
|                  |                 | be used to       |                  |
|                  |                 | identify the     |                  |
|                  |                 | function or pass |                  |
|                  |                 | additional       |                  |
|                  |                 | information      |                  |
+------------------+-----------------+------------------+------------------+
| verbose          | OT_BOOL         | Verbose          | casadi::ProtoFun |
|                  |                 | evaluation  for  | ction            |
|                  |                 | debugging        |                  |
+------------------+-----------------+------------------+------------------+

>Input scheme: casadi::DpleInput (DPLE_NUM_IN = 2)

+-----------+-------+------------------------------------------------------+
| Full name | Short |                     Description                      |
+===========+=======+======================================================+
| DPLE_A    | a     | A matrices (horzcat when const_dim, diagcat          |
|           |       | otherwise) [a].                                      |
+-----------+-------+------------------------------------------------------+
| DPLE_V    | v     | V matrices (horzcat when const_dim, diagcat          |
|           |       | otherwise) [v].                                      |
+-----------+-------+------------------------------------------------------+

>Output scheme: casadi::DpleOutput (DPLE_NUM_OUT = 1)

+-----------+-------+------------------------------------------------------+
| Full name | Short |                     Description                      |
+===========+=======+======================================================+
| DPLE_P    | p     | Lyapunov matrix (horzcat when const_dim, diagcat     |
|           |       | otherwise) (Cholesky of P if pos_def) [p].           |
+-----------+-------+------------------------------------------------------+

List of plugins
- slicot

Note: some of the plugins in this list might not be available on your 

system.  Also, there might be extra plugins available to you that are 
not 
listed here. You can obtain their documentation with  
Dple.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

slicot
------



An efficient solver for Discrete Periodic Lyapunov Equations using 
SLICOT

Uses Periodic Schur Decomposition ('psd') and does not assume positive
 
definiteness. Based on Periodic Lyapunov equations: some applications
 and 
new algorithms. Int. J. Control, vol. 67, pp. 69-87, 1997.

Overview of the method: J. Gillis Practical Methods for Approximate 
Robust 
Periodic Optimal Control ofNonlinear Mechanical Systems, PhD 
Thesis, 
KULeuven, 2015

Extra doc: https://github.com/casadi/casadi/wiki/L_22j

>List of available options

+-----------------------+-----------+--------------------------------------+
|          Id           |   Type    |             Description              |
+=======================+===========+======================================+
| linear_solver         | OT_STRING | User-defined linear solver class.    |
|                       |           | Needed for sensitivities.            |
+-----------------------+-----------+--------------------------------------+
| linear_solver_options | OT_DICT   | Options to be passed to the linear   |
|                       |           | solver.                              |
+-----------------------+-----------+--------------------------------------+
| psd_num_zero          | OT_DOUBLE | Numerical zero used in Periodic      |
|                       |           | Schur decomposition with slicot.This |
|                       |           | option is needed when your systems   |
|                       |           | has Floquet multiplierszero or close |
|                       |           | to zero                              |
+-----------------------+-----------+--------------------------------------+

Joris Gillis

Extra doc: https://github.com/casadi/casadi/wiki/L_21o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L97

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L97-L100

";

";

%feature("docstring") casadi::IndexRecution::dple_in "

[INTERNAL] 
Get DPLE input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ne

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L114

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L114-L121

>  std::string casadi::dple_in(casadi_int ind)
------------------------------------------------------------------------
[INTERNAL] 
Get DPLE input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ne

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L114

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L114-L121

";

";

%feature("docstring") casadi::IndexRecution::dple_out "

[INTERNAL] 
Get DPLE output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1nf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L123

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L123-L129

>  std::string casadi::dple_out(casadi_int ind)
------------------------------------------------------------------------
[INTERNAL] 
Get DPLE output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1nf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L123

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L123-L129

";

";

%feature("docstring") casadi::IndexRecution::dple_n_in "

[INTERNAL] 
Get the number of QP solver inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ng

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L131

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L131-L133

";

%feature("docstring") casadi::IndexRecution::dple_n_out "

[INTERNAL] 
Get the number of QP solver outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1nh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L135-L137

";

%feature("docstring") casadi::IndexRecution::trim_path "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::message_prefix "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::has_expm "

[INTERNAL] 
Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.hpp#L32

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.cpp#L32-L34

";

%feature("docstring") casadi::IndexRecution::load_expm "

[INTERNAL] 
Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.hpp#L36

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.cpp#L36-L38

";

%feature("docstring") casadi::IndexRecution::doc_expm "

[INTERNAL] 
Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.hpp#L40

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.cpp#L40-L42

";

%feature("docstring") casadi::IndexRecution::expmsol "

[INTERNAL]

Performs a matrix exponentiation expm(A)
General information

>List of available options

+------------------+-----------------+------------------+------------------+
|        Id        |      Type       |   Description    |     Used in      |
+==================+=================+==================+==================+
| ad_weight        | OT_DOUBLE       | Weighting factor | casadi::Function |
|                  |                 | for derivative   | Internal         |
|                  |                 | calculation.When |                  |
|                  |                 | there is an      |                  |
|                  |                 | option of either |                  |
|                  |                 | using forward or |                  |
|                  |                 | reverse mode     |                  |
|                  |                 | directional      |                  |
|                  |                 | derivatives, the |                  |
|                  |                 | condition ad_wei |                  |
|                  |                 | ght*nf<=(1-      |                  |
|                  |                 | ad_weight)*na is |                  |
|                  |                 | used where nf    |                  |
|                  |                 | and na are       |                  |
|                  |                 | estimates of the |                  |
|                  |                 | number of        |                  |
|                  |                 | forward/reverse  |                  |
|                  |                 | mode directional |                  |
|                  |                 | derivatives      |                  |
|                  |                 | needed. By       |                  |
|                  |                 | default,         |                  |
|                  |                 | ad_weight is     |                  |
|                  |                 | calculated       |                  |
|                  |                 | automatically,   |                  |
|                  |                 | but this can be  |                  |
|                  |                 | overridden by    |                  |
|                  |                 | setting this     |                  |
|                  |                 | option. In       |                  |
|                  |                 | particular, 0    |                  |
|                  |                 | means forcing    |                  |
|                  |                 | forward mode and |                  |
|                  |                 | 1 forcing        |                  |
|                  |                 | reverse mode.    |                  |
|                  |                 | Leave unset for  |                  |
|                  |                 | (class specific) |                  |
|                  |                 | heuristics.      |                  |
+------------------+-----------------+------------------+------------------+
| ad_weight_sp     | OT_DOUBLE       | Weighting factor | casadi::Function |
|                  |                 | for sparsity     | Internal         |
|                  |                 | pattern          |                  |
|                  |                 | calculation calc |                  |
|                  |                 | ulation.Override |                  |
|                  |                 | s default        |                  |
|                  |                 | behavior. Set to |                  |
|                  |                 | 0 and 1 to force |                  |
|                  |                 | forward and      |                  |
|                  |                 | reverse mode     |                  |
|                  |                 | respectively.    |                  |
|                  |                 | Cf. option       |                  |
|                  |                 | \"ad_weight\".     |                  |
|                  |                 | When set to -1,  |                  |
|                  |                 | sparsity is      |                  |
|                  |                 | completely       |                  |
|                  |                 | ignored and      |                  |
|                  |                 | dense matrices   |                  |
|                  |                 | are used.        |                  |
+------------------+-----------------+------------------+------------------+
| always_inline    | OT_BOOL         | Force inlining.  | casadi::Function |
|                  |                 |                  | Internal         |
+------------------+-----------------+------------------+------------------+
| cache            | OT_DICT         | Prepopulate the  | casadi::Function |
|                  |                 | function cache.  | Internal         |
|                  |                 | Default: empty   |                  |
+------------------+-----------------+------------------+------------------+
| compiler         | OT_STRING       | Just-in-time     | casadi::Function |
|                  |                 | compiler plugin  | Internal         |
|                  |                 | to be used.      |                  |
+------------------+-----------------+------------------+------------------+
| const_A          | OT_BOOL         | Assume A is      | casadi::Expm     |
|                  |                 | constant.        |                  |
|                  |                 | Default: false.  |                  |
+------------------+-----------------+------------------+------------------+
| custom_jacobian  | OT_FUNCTION     | Override         | casadi::Function |
|                  |                 | CasADi's AD. Use | Internal         |
|                  |                 | together with    |                  |
|                  |                 | 'jac_penalty':   |                  |
|                  |                 | 0. Note: Highly  |                  |
|                  |                 | experimental.    |                  |
|                  |                 | Syntax may break |                  |
|                  |                 | often.           |                  |
+------------------+-----------------+------------------+------------------+
| der_options      | OT_DICT         | Default options  | casadi::Function |
|                  |                 | to be used to    | Internal         |
|                  |                 | populate         |                  |
|                  |                 | forward_options, |                  |
|                  |                 | reverse_options, |                  |
|                  |                 | and              |                  |
|                  |                 | jacobian_options |                  |
|                  |                 | before those     |                  |
|                  |                 | options are      |                  |
|                  |                 | merged in.       |                  |
+------------------+-----------------+------------------+------------------+
| derivative_of    | OT_FUNCTION     | The function is  | casadi::Function |
|                  |                 | a derivative of  | Internal         |
|                  |                 | another          |                  |
|                  |                 | function. The    |                  |
|                  |                 | type of          |                  |
|                  |                 | derivative       |                  |
|                  |                 | (directional     |                  |
|                  |                 | derivative,      |                  |
|                  |                 | Jacobian) is     |                  |
|                  |                 | inferred from    |                  |
|                  |                 | the function     |                  |
|                  |                 | name.            |                  |
+------------------+-----------------+------------------+------------------+
| dump             | OT_BOOL         | Dump function to | casadi::Function |
|                  |                 | file upon first  | Internal         |
|                  |                 | evaluation.      |                  |
|                  |                 | [false]          |                  |
+------------------+-----------------+------------------+------------------+
| dump_dir         | OT_STRING       | Directory to     | casadi::Function |
|                  |                 | dump             | Internal         |
|                  |                 | inputs/outputs   |                  |
|                  |                 | to. Make sure    |                  |
|                  |                 | the directory    |                  |
|                  |                 | exists [.]       |                  |
+------------------+-----------------+------------------+------------------+
| dump_format      | OT_STRING       | Choose file      | casadi::Function |
|                  |                 | format to dump   | Internal         |
|                  |                 | matrices. See    |                  |
|                  |                 | DM.from_file     |                  |
|                  |                 | [mtx]            |                  |
+------------------+-----------------+------------------+------------------+
| dump_in          | OT_BOOL         | Dump numerical   | casadi::Function |
|                  |                 | values of inputs | Internal         |
|                  |                 | to file          |                  |
|                  |                 | (readable with   |                  |
|                  |                 | DM.from_file )   |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| dump_out         | OT_BOOL         | Dump numerical   | casadi::Function |
|                  |                 | values of        | Internal         |
|                  |                 | outputs to file  |                  |
|                  |                 | (readable with   |                  |
|                  |                 | DM.from_file )   |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| enable_fd        | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation by   |                  |
|                  |                 | finite           |                  |
|                  |                 | differencing.    |                  |
|                  |                 | [default:        |                  |
|                  |                 | false]]          |                  |
+------------------+-----------------+------------------+------------------+
| enable_forward   | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | Jacobian-times-  |                  |
|                  |                 | vector products  |                  |
|                  |                 | - typically      |                  |
|                  |                 | using forward    |                  |
|                  |                 | mode AD - if     |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| enable_jacobian  | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | Jacobians of all |                  |
|                  |                 | differentiable   |                  |
|                  |                 | outputs with     |                  |
|                  |                 | respect to all   |                  |
|                  |                 | differentiable   |                  |
|                  |                 | inputs - if      |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| enable_reverse   | OT_BOOL         | Enable           | casadi::Function |
|                  |                 | derivative       | Internal         |
|                  |                 | calculation      |                  |
|                  |                 | using generated  |                  |
|                  |                 | functions for    |                  |
|                  |                 | transposed       |                  |
|                  |                 | Jacobian-times-  |                  |
|                  |                 | vector products  |                  |
|                  |                 | - typically      |                  |
|                  |                 | using reverse    |                  |
|                  |                 | mode AD - if     |                  |
|                  |                 | available.       |                  |
|                  |                 | [default: true]  |                  |
+------------------+-----------------+------------------+------------------+
| error_on_fail    | OT_BOOL         | Throw exceptions | casadi::ProtoFun |
|                  |                 | when function    | ction            |
|                  |                 | evaluation fails |                  |
|                  |                 | (default true).  |                  |
+------------------+-----------------+------------------+------------------+
| fd_method        | OT_STRING       | Method for       | casadi::Function |
|                  |                 | finite           | Internal         |
|                  |                 | differencing     |                  |
|                  |                 | [default         |                  |
|                  |                 | 'central']       |                  |
+------------------+-----------------+------------------+------------------+
| fd_options       | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to the    | Internal         |
|                  |                 | finite           |                  |
|                  |                 | difference       |                  |
|                  |                 | instance         |                  |
+------------------+-----------------+------------------+------------------+
| forward_options  | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | forward mode     |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| gather_stats     | OT_BOOL         | Deprecated       | casadi::Function |
|                  |                 | option           | Internal         |
|                  |                 | (ignored):       |                  |
|                  |                 | Statistics are   |                  |
|                  |                 | now always       |                  |
|                  |                 | collected.       |                  |
+------------------+-----------------+------------------+------------------+
| input_scheme     | OT_STRINGVECTOR | Deprecated       | casadi::Function |
|                  |                 | option (ignored) | Internal         |
+------------------+-----------------+------------------+------------------+
| inputs_check     | OT_BOOL         | Throw exceptions | casadi::Function |
|                  |                 | when the         | Internal         |
|                  |                 | numerical values |                  |
|                  |                 | of the inputs    |                  |
|                  |                 | don't make sense |                  |
+------------------+-----------------+------------------+------------------+
| is_diff_in       | OT_BOOLVECTOR   | Indicate for     | casadi::Function |
|                  |                 | each input if it | Internal         |
|                  |                 | should be        |                  |
|                  |                 | differentiable.  |                  |
+------------------+-----------------+------------------+------------------+
| is_diff_out      | OT_BOOLVECTOR   | Indicate for     | casadi::Function |
|                  |                 | each output if   | Internal         |
|                  |                 | it should be     |                  |
|                  |                 | differentiable.  |                  |
+------------------+-----------------+------------------+------------------+
| jac_penalty      | OT_DOUBLE       | When requested   | casadi::Function |
|                  |                 | for a number of  | Internal         |
|                  |                 | forward/reverse  |                  |
|                  |                 | directions, it   |                  |
|                  |                 | may be cheaper   |                  |
|                  |                 | to compute first |                  |
|                  |                 | the full         |                  |
|                  |                 | jacobian and     |                  |
|                  |                 | then multiply    |                  |
|                  |                 | with seeds,      |                  |
|                  |                 | rather than      |                  |
|                  |                 | obtain the       |                  |
|                  |                 | requested        |                  |
|                  |                 | directions in a  |                  |
|                  |                 | straightforward  |                  |
|                  |                 | manner. Casadi   |                  |
|                  |                 | uses a heuristic |                  |
|                  |                 | to decide which  |                  |
|                  |                 | is cheaper. A    |                  |
|                  |                 | high value of    |                  |
|                  |                 | 'jac_penalty'    |                  |
|                  |                 | makes it less    |                  |
|                  |                 | likely for the   |                  |
|                  |                 | heurstic to      |                  |
|                  |                 | chose the full   |                  |
|                  |                 | Jacobian         |                  |
|                  |                 | strategy. The    |                  |
|                  |                 | special value -1 |                  |
|                  |                 | indicates never  |                  |
|                  |                 | to use the full  |                  |
|                  |                 | Jacobian         |                  |
|                  |                 | strategy         |                  |
+------------------+-----------------+------------------+------------------+
| jacobian_options | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | Jacobian         |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| jit              | OT_BOOL         | Use just-in-time | casadi::Function |
|                  |                 | compiler to      | Internal         |
|                  |                 | speed up the     |                  |
|                  |                 | evaluation       |                  |
+------------------+-----------------+------------------+------------------+
| jit_cleanup      | OT_BOOL         | Cleanup up the   | casadi::Function |
|                  |                 | temporary source | Internal         |
|                  |                 | file that jit    |                  |
|                  |                 | creates.         |                  |
|                  |                 | Default: true    |                  |
+------------------+-----------------+------------------+------------------+
| jit_name         | OT_STRING       | The file name    | casadi::Function |
|                  |                 | used to write    | Internal         |
|                  |                 | out code. The    |                  |
|                  |                 | actual file      |                  |
|                  |                 | names used       |                  |
|                  |                 | depend on 'jit_t |                  |
|                  |                 | emp_suffix' and  |                  |
|                  |                 | include          |                  |
|                  |                 | extensions.      |                  |
|                  |                 | Default:         |                  |
|                  |                 | 'jit_tmp'        |                  |
+------------------+-----------------+------------------+------------------+
| jit_options      | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to the    | Internal         |
|                  |                 | jit compiler.    |                  |
+------------------+-----------------+------------------+------------------+
| jit_serialize    | OT_STRING       | Specify          | casadi::Function |
|                  |                 | behaviour when   | Internal         |
|                  |                 | serializing a    |                  |
|                  |                 | jitted function: |                  |
|                  |                 | SOURCE|link|embe |                  |
|                  |                 | d.               |                  |
+------------------+-----------------+------------------+------------------+
| jit_temp_suffix  | OT_BOOL         | Use a temporary  | casadi::Function |
|                  |                 | (seemingly       | Internal         |
|                  |                 | random) filename |                  |
|                  |                 | suffix for       |                  |
|                  |                 | generated code   |                  |
|                  |                 | and libraries.   |                  |
|                  |                 | This is desired  |                  |
|                  |                 | for thread-      |                  |
|                  |                 | safety. This     |                  |
|                  |                 | behaviour may    |                  |
|                  |                 | defeat caching   |                  |
|                  |                 | compiler         |                  |
|                  |                 | wrappers.        |                  |
|                  |                 | Default: true    |                  |
+------------------+-----------------+------------------+------------------+
| max_io           | OT_INT          | Acceptable       | casadi::Function |
|                  |                 | number of inputs | Internal         |
|                  |                 | and outputs.     |                  |
|                  |                 | Warn if          |                  |
|                  |                 | exceeded.        |                  |
+------------------+-----------------+------------------+------------------+
| max_num_dir      | OT_INT          | Specify the      | casadi::Function |
|                  |                 | maximum number   | Internal         |
|                  |                 | of directions    |                  |
|                  |                 | for derivative   |                  |
|                  |                 | functions.       |                  |
|                  |                 | Overrules the    |                  |
|                  |                 | builtin optimize |                  |
|                  |                 | d_num_dir.       |                  |
+------------------+-----------------+------------------+------------------+
| never_inline     | OT_BOOL         | Forbid inlining. | casadi::Function |
|                  |                 |                  | Internal         |
+------------------+-----------------+------------------+------------------+
| output_scheme    | OT_STRINGVECTOR | Deprecated       | casadi::Function |
|                  |                 | option (ignored) | Internal         |
+------------------+-----------------+------------------+------------------+
| post_expand      | OT_BOOL         | After            | casadi::Function |
|                  |                 | construction,    | Internal         |
|                  |                 | expand this      |                  |
|                  |                 | Function .       |                  |
|                  |                 | Default: False   |                  |
+------------------+-----------------+------------------+------------------+
| post_expand_opti | OT_DICT         | Options to be    | casadi::Function |
| ons              |                 | passed to post-  | Internal         |
|                  |                 | construction     |                  |
|                  |                 | expansion.       |                  |
|                  |                 | Default: empty   |                  |
+------------------+-----------------+------------------+------------------+
| print_in         | OT_BOOL         | Print numerical  | casadi::Function |
|                  |                 | values of inputs | Internal         |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| print_out        | OT_BOOL         | Print numerical  | casadi::Function |
|                  |                 | values of        | Internal         |
|                  |                 | outputs          |                  |
|                  |                 | [default: false] |                  |
+------------------+-----------------+------------------+------------------+
| print_time       | OT_BOOL         | print            | casadi::ProtoFun |
|                  |                 | information      | ction            |
|                  |                 | about execution  |                  |
|                  |                 | time. Implies    |                  |
|                  |                 | record_time.     |                  |
+------------------+-----------------+------------------+------------------+
| record_time      | OT_BOOL         | record           | casadi::ProtoFun |
|                  |                 | information      | ction            |
|                  |                 | about execution  |                  |
|                  |                 | time, for        |                  |
|                  |                 | retrieval with   |                  |
|                  |                 | stats().         |                  |
+------------------+-----------------+------------------+------------------+
| regularity_check | OT_BOOL         | Throw exceptions | casadi::ProtoFun |
|                  |                 | when NaN or Inf  | ction            |
|                  |                 | appears during   |                  |
|                  |                 | evaluation       |                  |
+------------------+-----------------+------------------+------------------+
| reverse_options  | OT_DICT         | Options to be    | casadi::Function |
|                  |                 | passed to a      | Internal         |
|                  |                 | reverse mode     |                  |
|                  |                 | constructor      |                  |
+------------------+-----------------+------------------+------------------+
| user_data        | OT_VOIDPTR      | A user-defined   | casadi::Function |
|                  |                 | field that can   | Internal         |
|                  |                 | be used to       |                  |
|                  |                 | identify the     |                  |
|                  |                 | function or pass |                  |
|                  |                 | additional       |                  |
|                  |                 | information      |                  |
+------------------+-----------------+------------------+------------------+
| verbose          | OT_BOOL         | Verbose          | casadi::ProtoFun |
|                  |                 | evaluation  for  | ction            |
|                  |                 | debugging        |                  |
+------------------+-----------------+------------------+------------------+

List of plugins
- slicot

Note: some of the plugins in this list might not be available on your 

system.  Also, there might be extra plugins available to you that are 
not 
listed here. You can obtain their documentation with  
Expm.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

slicot
------



Extra doc: https://github.com/casadi/casadi/wiki/L_22l

Joris Gillis

Extra doc: https://github.com/casadi/casadi/wiki/L_21l

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.hpp#L44

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.cpp#L44-L47

";

%feature("docstring") casadi::IndexRecution::expm_n_in "

[INTERNAL] 
Get the number of expm solver inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_rs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.hpp#L49

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.cpp#L49-L51

";

%feature("docstring") casadi::IndexRecution::expm_n_out "

[INTERNAL] 
Get the number of expm solver outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_rt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.hpp#L53

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.cpp#L53-L55

";

%feature("docstring") casadi::IndexRecution::external "

[INTERNAL] 
Load an external function from a shared library.

Parameters:
-----------

name: 
Name as in the label assigned to a CasADi  Function object: 
Function(name,...,...) Will be used to look up 
symbols/functions named eg. 
<name>_eval Use  nm (linux/osx) or  depends.exe (win) to check which symbols
 are present in your shared library

bin_name: 
File name of the shared library

Extra doc: https://github.com/casadi/casadi/wiki/L_i1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/external.hpp#L51

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/external.cpp#L51-L54

>  Function casadi::external(const std::string &name, const std::string &bin_name, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Load an external function from a shared library.

Parameters:
-----------

name: 
Name as in the label assigned to a CasADi  Function object: 
Function(name,...,...) Will be used to look up 
symbols/functions named eg. 
<name>_eval Use  nm (linux/osx) or  depends.exe (win) to check which symbols
 are present in your shared library

bin_name: 
File name of the shared library

Extra doc: https://github.com/casadi/casadi/wiki/L_i1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/external.hpp#L51

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/external.cpp#L51-L54

";

";

%feature("docstring") casadi::IndexRecution::_function_buffer_eval "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::index_interp1d "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::combine "

[INTERNAL] 
Combine two dicts. First has priority.

Extra doc: https://github.com/casadi/casadi/wiki/L_17t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L593

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L593-L599

";

%feature("docstring") casadi::IndexRecution::update_dict "

[INTERNAL] 
Update the target dictionary in place with source elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_17u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L601

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L601-L614

";

%feature("docstring") casadi::IndexRecution::get_from_dict "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::extract_from_dict "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::extract_from_dict_inplace "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::is_slice "

[INTERNAL] 
Check if an index vector can be represented more efficiently as 
a 
slice.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L170

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L170-L200

>  bool casadi::is_slice(const std::vector< casadi_int > &v, bool ind1=false)
------------------------------------------------------------------------
[INTERNAL] 
Check if an index vector can be represented more efficiently as a 
slice.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L170

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L170-L200

";

";

%feature("docstring") casadi::IndexRecution::to_slice "

[INTERNAL] 
Construct from an index vector (requires is_slice(v) to be true)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L152

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L152-L168

>  Slice casadi::to_slice(const std::vector< casadi_int > &v, bool ind1=false)
------------------------------------------------------------------------
[INTERNAL] 
Construct from an index vector (requires is_slice(v) to be true)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L152

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L152-L168

";

";

%feature("docstring") casadi::IndexRecution::text2type "

[INTERNAL] 
Convert to a type.

Extra doc: https://github.com/casadi/casadi/wiki/L_15y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L39

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L39-L44

";

%feature("docstring") casadi::IndexRecution::text2vector "

[INTERNAL] 
Get entry as a vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_15z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L50

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L50-L56

";

%feature("docstring") casadi::IndexRecution::text2set "

[INTERNAL] 
Get entry as a set.

Extra doc: https://github.com/casadi/casadi/wiki/L_160

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L62

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L62-L67

";

%feature("docstring") casadi::IndexRecution::simpleRK "

[INTERNAL] 
Construct an explicit Runge-Kutta integrator.

The constructed function has three inputs, corresponding to initial 
state 
(x0), parameter (p) and integration time (h) and one output, 
corresponding 
to final state (xf).

Parameters:
-----------

f: 
ODE function with two inputs (x and p) and one output (xdot)

N: 
Number of integrator steps

order: 
Order of interpolating polynomials

Extra doc: https://github.com/casadi/casadi/wiki/L_1sr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L128

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L128-L187

";

%feature("docstring") casadi::IndexRecution::collocation_interpolators "

[INTERNAL] 
Obtain collocation interpolating matrices.

A collocation method poses a polynomial Pi that interpolates exactly 

through an initial state (0,X_0) and helper states at collocation 
points 
(tau_j,X:collPoint(j)).

This function computes the linear mapping between dPi/dt and 
coefficients 
Z=[X_0 X:collPoints].

Parameters:
-----------

tau: 
location of collocation points, as obtained from collocation_points

C: 
interpolating coefficients to obtain derivatives. Length: order+1, 

order+1



::

dPi/dt @Z_j = (1/h) Sum_i C[j][i]*Z_i,



with h the length of the integration interval.

Parameters:
-----------

D: 
interpolating coefficients to obtain end state. Length: order+1



::

Pi @X_f = Sum_i D[i]*Z_i



Extra doc: https://github.com/casadi/casadi/wiki/L_1sp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L189-L231

";

%feature("docstring") casadi::IndexRecution::collocation_coeff "

[INTERNAL] 
Obtain collocation interpolating matrices.

A collocation method poses a polynomial Pi that interpolates exactly 

through an initial state (0,X_0) and helper states at collocation 
points 
(tau_j,Xc_j) with j=1..degree.

This function computes the linear mapping between dPi/dt and 
coefficients 
Z=[X_0 Xc].

Parameters:
-----------

tau: 
location of collocation points (length: degree), as obtained from 

collocation_points

C: 
interpolating coefficients to obtain derivatives. Size: (degree+1)-by-

degree

You may find the slopes of Pi at the collocation points as

::

dPi/dt @ Xc = (1/h) Z*C,



with h the length of the integration interval.

Parameters:
-----------

D: 
interpolating coefficients to obtain end state. Size: (degree+1)-by-1

You may find the end point of Pi as

::

Pi @X_f = Z*D



Parameters:
-----------

B: 
quadrature coefficients Size: degree-by-1

Given quadrature righ-hand-sides 'quad' evaluated at the collocation 

points, you may find the integrated quadratures as

::

q = quad*B*h



Extra doc: https://github.com/casadi/casadi/wiki/L_1sq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L233-L287

";

%feature("docstring") casadi::IndexRecution::simpleIRK "

[INTERNAL] 
Construct an implicit Runge-Kutta integrator using a collocation
 
scheme.

The constructed function has three inputs, corresponding to initial 
state 
(x0), parameter (p) and integration time (h) and one output, 
corresponding 
to final state (xf).

Parameters:
-----------

f: 
ODE function with two inputs (x and p) and one output (xdot)

N: 
Number of integrator steps

order: 
Order of interpolating polynomials

scheme: 
 Collocation scheme, as excepted by collocationPoints function.

solver: 
Solver plugin

solver_options: 
Options to be passed to the solver plugin

Extra doc: https://github.com/casadi/casadi/wiki/L_1ss

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L289

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L289-L354

";

%feature("docstring") casadi::IndexRecution::simpleIntegrator "

[INTERNAL] 
Simplified wrapper for the  Integrator class.

Extra doc: https://github.com/casadi/casadi/wiki/L_1st

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L356

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L356-L395

";

%feature("docstring") casadi::IndexRecution::has_integrator "

[INTERNAL] 
Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L97

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L97-L99

";

%feature("docstring") casadi::IndexRecution::load_integrator "

[INTERNAL] 
Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L101

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L101-L103

";

%feature("docstring") casadi::IndexRecution::doc_integrator "

[INTERNAL] 
Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L105

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L105-L107

";

%feature("docstring") casadi::IndexRecution::integrator "

[INTERNAL]

>  Function casadi::integrator(const std::string &name, const std::string &solver, const MXDict &dae, const Dict &opts)

>  Function casadi::integrator(const std::string &name, const std::string &solver, const Function &dae, const Dict &opts)

>  Function casadi::integrator(const std::string &name, const std::string &solver, const SXDict &dae, double t0, const std::vector< double > &tout, const Dict &opts)

>  Function casadi::integrator(const std::string &name, const std::string &solver, const MXDict &dae, double t0, const std::vector< double > &tout, const Dict &opts)

>  Function casadi::integrator(const std::string &name, const std::string &solver, const Function &dae, double t0, const std::vector< double > &tout, const Dict &opts)

>  Function casadi::integrator(const std::string &name, const std::string &solver, const SXDict &dae, double t0, double tf, const Dict &opts)

>  Function casadi::integrator(const std::string &name, const std::string &solver, const MXDict &dae, double t0, double tf, const Dict &opts)

>  Function casadi::integrator(const std::string &name, const std::string &solver, const Function &dae, double t0, double tf, const Dict &opts)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::IndexRecution::integrator_in "

[INTERNAL] 
Get integrator input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_7d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L171

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L171-L183

>  std::string casadi::integrator_in(casadi_int ind)
------------------------------------------------------------------------
[INTERNAL] 
Get integrator input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_7d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L171

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L171-L183

";

";

%feature("docstring") casadi::IndexRecution::integrator_out "

[INTERNAL] 
Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_7e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L185

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L185-L197

>  std::string casadi::integrator_out(casadi_int ind)
------------------------------------------------------------------------
[INTERNAL] 
Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_7e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L185

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L185-L197

";

";

%feature("docstring") casadi::IndexRecution::integrator_n_in "

[INTERNAL] 
Get the number of integrator inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_7f

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L199

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L199-L201

";

%feature("docstring") casadi::IndexRecution::integrator_n_out "

[INTERNAL] 
Get the number of integrator outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_7g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L203

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L203-L205

";

%feature("docstring") casadi::IndexRecution::dyn_in "

[INTERNAL] 
Get simulator input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_25r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L215

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L215-L217

>  std::string casadi::dyn_in(casadi_int ind)
------------------------------------------------------------------------
[INTERNAL] 
Get simulator input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_25r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L215

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L215-L217

";

";

%feature("docstring") casadi::IndexRecution::dyn_out "

[INTERNAL] 
Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_25s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L219

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L219-L221

>  std::string casadi::dyn_out(casadi_int ind)
------------------------------------------------------------------------
[INTERNAL] 
Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_25s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L219

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L219-L221

";

";

%feature("docstring") casadi::IndexRecution::dyn_n_in "

[INTERNAL] 
Get the number of simulator inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_25t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L223

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L223-L225

";

%feature("docstring") casadi::IndexRecution::dyn_n_out "

[INTERNAL] 
Get the number of simulator outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_25u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L227

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L227-L229

";

%feature("docstring") casadi::IndexRecution::has_interpolant "

[INTERNAL] 
Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.hpp#L34

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.cpp#L34-L36

";

%feature("docstring") casadi::IndexRecution::load_interpolant "

[INTERNAL] 
Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.hpp#L38

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.cpp#L38-L40

";

%feature("docstring") casadi::IndexRecution::doc_interpolant "

[INTERNAL] 
Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.hpp#L42

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.cpp#L42-L44

";

%feature("docstring") casadi::IndexRecution::interpolant "

[INTERNAL] 
Parametric variant of interpolant.

The resulting function will have additional arguments for the grid and
 
coefficients

By default, derivatives wrt the coefficients are not supported (zero).
 Some
 interpolant plugins may support the  inline=true which enables correct 
derivatives

Extra doc: https://github.com/casadi/casadi/wiki/L_1p4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.hpp#L205

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.cpp#L205-L213

>  Function casadi::interpolant(const std::string &name, const std::string &solver, const std::vector< casadi_int > &grid_dims, casadi_int m=1, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Parametric variant of interpolant.

The resulting function will have additional arguments for the grid and
 
coefficients

By default, derivatives wrt the coefficients are not supported (zero).
 Some
 interpolant plugins may support the  inline=true which enables correct 
derivatives

Extra doc: https://github.com/casadi/casadi/wiki/L_1p4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.hpp#L205

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.cpp#L205-L213

";

";

%feature("docstring") casadi::IndexRecution::has_linsol "

[INTERNAL] 
Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L205

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L205-L207

";

%feature("docstring") casadi::IndexRecution::load_linsol "

[INTERNAL] 
Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L209

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L209-L211

";

%feature("docstring") casadi::IndexRecution::doc_linsol "

[INTERNAL] 
Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L213

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L213-L215

";

%feature("docstring") casadi::IndexRecution::detect_simple_bounds "

[INTERNAL]

>  void casadi::detect_simple_bounds(const MX &x, const MX &p, const MX &g, const MX &lbg, const MX &ubg, std::vector< casadi_int > &gi, MX &lbx, MX &ubx, Function &lam_forward, Function &lam_backward)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::IndexRecution::check_sos "

[INTERNAL] 
Check sos structure and generate defaults.

Extra doc: https://github.com/casadi/casadi/wiki/L_1sx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_tools.hpp#L79

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_tools.hpp#L79-L120

";

%feature("docstring") casadi::IndexRecution::has_nlpsol "

[INTERNAL] 
Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L34

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L34-L36

";

%feature("docstring") casadi::IndexRecution::load_nlpsol "

[INTERNAL] 
Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L38

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L38-L40

";

%feature("docstring") casadi::IndexRecution::doc_nlpsol "

[INTERNAL] 
Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L42

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L42-L44

";

%feature("docstring") casadi::IndexRecution::nlpsol "

[INTERNAL]

>  Function casadi::nlpsol(const std::string &name, const std::string &solver, const MXDict &nlp, const Dict &opts)

>  Function casadi::nlpsol(const std::string &name, const std::string &solver, const NlpBuilder &nl, const Dict &opts)

>  Function casadi::nlpsol(const std::string &name, const std::string &solver, const std::string &fname, const Dict &opts)

>  Function casadi::nlpsol(const std::string &name, const std::string &solver, const Importer &compiler, const Dict &opts)

>  Function casadi::nlpsol(const std::string &name, const std::string &solver, const Function &nlp, const Dict &opts)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::IndexRecution::nlpsol_in "

[INTERNAL] 
Get NLP solver input scheme name by index.

>Input scheme: casadi::NlpsolInput (NLPSOL_NUM_IN = 8)

+---------------+--------+-------------------------------------------------+
|   Full name   | Short  |                   Description                   |
+===============+========+=================================================+
| NLPSOL_X0     | x0     | Decision variables, initial guess (nx x 1)      |
+---------------+--------+-------------------------------------------------+
| NLPSOL_P      | p      | Value of fixed parameters (np x 1)              |
+---------------+--------+-------------------------------------------------+
| NLPSOL_LBX    | lbx    | Decision variables lower bound (nx x 1),        |
|               |        | default -inf.                                   |
+---------------+--------+-------------------------------------------------+
| NLPSOL_UBX    | ubx    | Decision variables upper bound (nx x 1),        |
|               |        | default +inf.                                   |
+---------------+--------+-------------------------------------------------+
| NLPSOL_LBG    | lbg    | Constraints lower bound (ng x 1), default -inf. |
+---------------+--------+-------------------------------------------------+
| NLPSOL_UBG    | ubg    | Constraints upper bound (ng x 1), default +inf. |
+---------------+--------+-------------------------------------------------+
| NLPSOL_LAM_X0 | lam_x0 | Lagrange multipliers for bounds on X, initial   |
|               |        | guess (nx x 1)                                  |
+---------------+--------+-------------------------------------------------+
| NLPSOL_LAM_G0 | lam_g0 | Lagrange multipliers for bounds on G, initial   |
|               |        | guess (ng x 1)                                  |
+---------------+--------+-------------------------------------------------+

Extra doc: https://github.com/casadi/casadi/wiki/L_1t0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L228

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L228-L241

>  std::string casadi::nlpsol_in(casadi_int ind)
------------------------------------------------------------------------
[INTERNAL] 
Get NLP solver input scheme name by index.

>Input scheme: casadi::NlpsolInput (NLPSOL_NUM_IN = 8)

+---------------+--------+-------------------------------------------------+
|   Full name   | Short  |                   Description                   |
+===============+========+=================================================+
| NLPSOL_X0     | x0     | Decision variables, initial guess (nx x 1)      |
+---------------+--------+-------------------------------------------------+
| NLPSOL_P      | p      | Value of fixed parameters (np x 1)              |
+---------------+--------+-------------------------------------------------+
| NLPSOL_LBX    | lbx    | Decision variables lower bound (nx x 1),        |
|               |        | default -inf.                                   |
+---------------+--------+-------------------------------------------------+
| NLPSOL_UBX    | ubx    | Decision variables upper bound (nx x 1),        |
|               |        | default +inf.                                   |
+---------------+--------+-------------------------------------------------+
| NLPSOL_LBG    | lbg    | Constraints lower bound (ng x 1), default -inf. |
+---------------+--------+-------------------------------------------------+
| NLPSOL_UBG    | ubg    | Constraints upper bound (ng x 1), default +inf. |
+---------------+--------+-------------------------------------------------+
| NLPSOL_LAM_X0 | lam_x0 | Lagrange multipliers for bounds on X, initial   |
|               |        | guess (nx x 1)                                  |
+---------------+--------+-------------------------------------------------+
| NLPSOL_LAM_G0 | lam_g0 | Lagrange multipliers for bounds on G, initial   |
|               |        | guess (ng x 1)                                  |
+---------------+--------+-------------------------------------------------+

Extra doc: https://github.com/casadi/casadi/wiki/L_1t0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L228

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L228-L241

";

";

%feature("docstring") casadi::IndexRecution::nlpsol_out "

[INTERNAL] 
Get output scheme name by index.

>Output scheme: casadi::NlpsolOutput (NLPSOL_NUM_OUT = 6)

+--------------+-------+---------------------------------------------------+
|  Full name   | Short |                    Description                    |
+==============+=======+===================================================+
| NLPSOL_X     | x     | Decision variables at the optimal solution (nx x  |
|              |       | 1)                                                |
+--------------+-------+---------------------------------------------------+
| NLPSOL_F     | f     | Cost function value at the optimal solution (1 x  |
|              |       | 1)                                                |
+--------------+-------+---------------------------------------------------+
| NLPSOL_G     | g     | Constraints function at the optimal solution (ng  |
|              |       | x 1)                                              |
+--------------+-------+---------------------------------------------------+
| NLPSOL_LAM_X | lam_x | Lagrange multipliers for bounds on X at the       |
|              |       | solution (nx x 1)                                 |
+--------------+-------+---------------------------------------------------+
| NLPSOL_LAM_G | lam_g | Lagrange multipliers for bounds on G at the       |
|              |       | solution (ng x 1)                                 |
+--------------+-------+---------------------------------------------------+
| NLPSOL_LAM_P | lam_p | Lagrange multipliers for bounds on P at the       |
|              |       | solution (np x 1)                                 |
+--------------+-------+---------------------------------------------------+

Extra doc: https://github.com/casadi/casadi/wiki/L_1t1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L243

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L243-L254

>  std::string casadi::nlpsol_out(casadi_int ind)
------------------------------------------------------------------------
[INTERNAL] 
Get output scheme name by index.

>Output scheme: casadi::NlpsolOutput (NLPSOL_NUM_OUT = 6)

+--------------+-------+---------------------------------------------------+
|  Full name   | Short |                    Description                    |
+==============+=======+===================================================+
| NLPSOL_X     | x     | Decision variables at the optimal solution (nx x  |
|              |       | 1)                                                |
+--------------+-------+---------------------------------------------------+
| NLPSOL_F     | f     | Cost function value at the optimal solution (1 x  |
|              |       | 1)                                                |
+--------------+-------+---------------------------------------------------+
| NLPSOL_G     | g     | Constraints function at the optimal solution (ng  |
|              |       | x 1)                                              |
+--------------+-------+---------------------------------------------------+
| NLPSOL_LAM_X | lam_x | Lagrange multipliers for bounds on X at the       |
|              |       | solution (nx x 1)                                 |
+--------------+-------+---------------------------------------------------+
| NLPSOL_LAM_G | lam_g | Lagrange multipliers for bounds on G at the       |
|              |       | solution (ng x 1)                                 |
+--------------+-------+---------------------------------------------------+
| NLPSOL_LAM_P | lam_p | Lagrange multipliers for bounds on P at the       |
|              |       | solution (np x 1)                                 |
+--------------+-------+---------------------------------------------------+

Extra doc: https://github.com/casadi/casadi/wiki/L_1t1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L243

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L243-L254

";

";

%feature("docstring") casadi::IndexRecution::nlpsol_n_in "

[INTERNAL] 
Number of NLP solver inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L256

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L256-L258

";

%feature("docstring") casadi::IndexRecution::nlpsol_n_out "

[INTERNAL] 
Number of NLP solver outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L260

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L260-L262

";

%feature("docstring") casadi::IndexRecution::nlpsol_options "

[INTERNAL] 
Get all options for a plugin.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L900

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L900-L902

";

%feature("docstring") casadi::IndexRecution::nlpsol_option_type "

[INTERNAL] 
Get type info for a particular option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L904

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L904-L906

";

%feature("docstring") casadi::IndexRecution::nlpsol_option_info "

[INTERNAL] 
Get documentation for a particular option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L908

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L908-L910

";

%feature("docstring") casadi::IndexRecution::rootfinder_in "

[INTERNAL] 
Get rootfinder input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L47

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L47-L54

>  std::string casadi::rootfinder_in(casadi_int ind)
------------------------------------------------------------------------
[INTERNAL] 
Get rootfinder input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L47

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L47-L54

";

";

%feature("docstring") casadi::IndexRecution::rootfinder_out "

[INTERNAL] 
Get rootfinder output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L56

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L56-L62

>  std::string casadi::rootfinder_out(casadi_int ind)
------------------------------------------------------------------------
[INTERNAL] 
Get rootfinder output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L56

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L56-L62

";

";

%feature("docstring") casadi::IndexRecution::rootfinder_n_in "

[INTERNAL] 
Number of rootfinder inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L64

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L64-L66

";

%feature("docstring") casadi::IndexRecution::rootfinder_n_out "

[INTERNAL] 
Number of rootfinder outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L68

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L68-L70

";

%feature("docstring") casadi::IndexRecution::rootfinder_options "

[INTERNAL] 
Get all options for a plugin.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L72

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L72-L74

";

%feature("docstring") casadi::IndexRecution::rootfinder_option_type "

[INTERNAL] 
Get type info for a particular option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L76

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L76-L78

";

%feature("docstring") casadi::IndexRecution::rootfinder_option_info "

[INTERNAL] 
Get documentation for a particular option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L80

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L80-L82

";

%feature("docstring") casadi::IndexRecution::has_rootfinder "

[INTERNAL] 
Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L84

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L84-L86

";

%feature("docstring") casadi::IndexRecution::load_rootfinder "

[INTERNAL] 
Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L88

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L88-L90

";

%feature("docstring") casadi::IndexRecution::doc_rootfinder "

[INTERNAL] 
Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L92

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L92-L94

";

%feature("docstring") casadi::IndexRecution::rootfinder "

[INTERNAL]

>  Function casadi::rootfinder(const std::string &name, const std::string &solver, const MXDict &rfp, const Dict &opts)

>  Function casadi::rootfinder(const std::string &name, const std::string &solver, const Function &f, const Dict &opts)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::IndexRecution::einstein_process "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::Contraction "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::einstein_eval "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::shared_cast "

[INTERNAL] 
Typecast a shared object to a base class to a shared object to a
 
derived class,.

cf. dynamic_cast (const)

Extra doc: https://github.com/casadi/casadi/wiki/L_b7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L259

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L259-L262

>  const B casadi::shared_cast(const SharedObject &A)
------------------------------------------------------------------------
[INTERNAL] 
Typecast a shared object to a base class to a shared object to a 
derived class,.

cf. dynamic_cast (const)

Extra doc: https://github.com/casadi/casadi/wiki/L_b7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L259

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L259-L262

";

";

%feature("docstring") casadi::IndexRecution::is_slice2 "

[INTERNAL] 
Check if an index vector can be represented more efficiently as 
two 
nested slices.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L202

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L202-L253

";

%feature("docstring") casadi::IndexRecution::to_slice2 "

[INTERNAL] 
Construct nested slices from an index vector (requires 
is_slice2(v) to
 be true)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L255

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L255-L289

";

%feature("docstring") casadi::IndexRecution::matrixName< SXElem > "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::slicot_mb03vd "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::slicot_mb03vy "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::slicot_mb03wd "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::slicot_mb05nd "

[INTERNAL] ";


// File: namespacecasadi_1_1IndexRecution.xml


// File: namespaceproxsuite_1_1proxqp.xml


// File: namespacestd.xml
%feature("docstring") std::mul_overflows "
[INTERNAL] ";


// File: namespacestd_1_1chrono.xml


// File: chapter1.xml


// File: chapter2.xml


// File: chapter3.xml


// File: chapter4.xml


// File: chapter5.xml


// File: chapter6.xml

