
// File: index.xml

// File: classcasadi_1_1Blocksqp.xml
%feature("docstring") casadi::Blocksqp "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1BonMinMessageHandler.xml
%feature("docstring") casadi::BonMinMessageHandler "

Helper class to direct messages to  uout()

IPOPT has the concept of a Journal/Journalist BONMIN and CBC do not.

";


// File: classcasadi_1_1BSplineInterpolant.xml
%feature("docstring") casadi::BSplineInterpolant "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Callback.xml
%feature("docstring") casadi::Callback "

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

%feature("docstring") casadi::Callback::has_jacobian "

Return Jacobian of all input elements with respect to all output 
elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_oh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L168

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L94-L96

";

%feature("docstring") casadi::Callback::get_jacobian "

Return Jacobian of all input elements with respect to all output 
elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_oh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L169

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L99-L104

";

%feature("docstring") casadi::Callback::has_forward "

Return function that calculates forward derivatives.

forward(nfwd) returns a cached instance if available, and calls   Function 
get_forward(casadi_int nfwd) if no cached version is available.

Extra doc: https://github.com/casadi/casadi/wiki/L_oi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L183

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L114-L116

";

%feature("docstring") casadi::Callback::get_forward "

Return function that calculates forward derivatives.

forward(nfwd) returns a cached instance if available, and calls   Function 
get_forward(casadi_int nfwd) if no cached version is available.

Extra doc: https://github.com/casadi/casadi/wiki/L_oi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L184

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L107-L112

";

%feature("docstring") casadi::Callback::has_reverse "

Return function that calculates adjoint derivatives.

reverse(nadj) returns a cached instance if available, and calls   Function 
get_reverse(casadi_int nadj) if no cached version is available.

Extra doc: https://github.com/casadi/casadi/wiki/L_oj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L198

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L126-L128

";

%feature("docstring") casadi::Callback::get_reverse "

Return function that calculates adjoint derivatives.

reverse(nadj) returns a cached instance if available, and calls   Function 
get_reverse(casadi_int nadj) if no cached version is available.

Extra doc: https://github.com/casadi/casadi/wiki/L_oj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L199

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L119-L124

";

%feature("docstring") casadi::Callback::has_jac_sparsity "

Return sparsity of Jacobian of all input elements.

with respect to all output elements

Extra doc: https://github.com/casadi/casadi/wiki/L_ok

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L211

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L211-L211

";

%feature("docstring") casadi::Callback::get_jac_sparsity "

Return sparsity of Jacobian of all input elements.

with respect to all output elements

Extra doc: https://github.com/casadi/casadi/wiki/L_ok

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L212

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L212-L213

";

%feature("docstring") casadi::Callback::jit "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L260-L272

>  Function casadi::Function::jit(const std::string &name, const std::string &body, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const std::vector< Sparsity > &sparsity_in, const std::vector< Sparsity > &sparsity_out, const Dict &opts=Dict())
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L260-L272

";

";

%feature("docstring") casadi::Callback::expand "

Expand a function to SX.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L207

>  Function casadi::Function::expand(const std::string &name, const Dict &opts=Dict()) const
------------------------------------------------------------------------

Expand a function to SX.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L207

";

";

%feature("docstring") casadi::Callback::size1_in "

Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240-L240

>  casadi_int casadi::Function::size1_in(const std::string &iname) const
------------------------------------------------------------------------

Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240-L240

";

";

%feature("docstring") casadi::Callback::size2_in "

Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242-L242

>  casadi_int casadi::Function::size2_in(const std::string &iname) const
------------------------------------------------------------------------

Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242-L242

";

";

%feature("docstring") casadi::Callback::size_in "

Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244-L246

>  std::pair<casadi_int, casadi_int> casadi::Function::size_in(const std::string &iname) const
------------------------------------------------------------------------

Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244-L246

";

";

%feature("docstring") casadi::Callback::size1_out "

Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254-L254

>  casadi_int casadi::Function::size1_out(const std::string &oname) const
------------------------------------------------------------------------

Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254-L254

";

";

%feature("docstring") casadi::Callback::size2_out "

Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256-L256

>  casadi_int casadi::Function::size2_out(const std::string &oname) const
------------------------------------------------------------------------

Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256-L256

";

";

%feature("docstring") casadi::Callback::size_out "

Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258-L260

>  std::pair<casadi_int, casadi_int> casadi::Function::size_out(const std::string &oname) const
------------------------------------------------------------------------

Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258-L260

";

";

%feature("docstring") casadi::Callback::nnz_in "

Get number of input nonzeros.

For a particular input or for all of the inputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L271

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L271-L271

>  casadi_int casadi::Function::nnz_in(const std::string &iname) const
------------------------------------------------------------------------

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

Get number of output nonzeros.

For a particular output or for all of the outputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L282

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L282-L282

>  casadi_int casadi::Function::nnz_out(const std::string &oname) const
------------------------------------------------------------------------

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

Get number of input elements.

For a particular input or for all of the inputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1ve

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L293-L293

>  casadi_int casadi::Function::numel_in(const std::string &iname) const
------------------------------------------------------------------------

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

Get number of output elements.

For a particular output or for all of the outputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L304

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L304-L304

>  casadi_int casadi::Function::numel_out(const std::string &oname) const
------------------------------------------------------------------------

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

Get sparsity of a given input.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L373

>  const Sparsity& casadi::Function::sparsity_in(const std::string &iname) const
------------------------------------------------------------------------

Get sparsity of a given input.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L373

";

";

%feature("docstring") casadi::Callback::sparsity_out "

Get sparsity of a given output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L381

>  const Sparsity& casadi::Function::sparsity_out(const std::string &iname) const
------------------------------------------------------------------------

Get sparsity of a given output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L381

";

";

%feature("docstring") casadi::Callback::is_diff_in "

Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L390

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1005-L1011

>  std::vector< bool > casadi::Function::is_diff_in() const
------------------------------------------------------------------------

Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L390

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1005-L1011

";

";

%feature("docstring") casadi::Callback::is_diff_out "

Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L391

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1013-L1019

>  std::vector< bool > casadi::Function::is_diff_out() const
------------------------------------------------------------------------

Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L391

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1013-L1019

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

Evaluate the function symbolically or numerically.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L524

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1407-L1414

>  void casadi::Function::call(const MXDict &arg, MXDict &res, bool always_inline=false, bool never_inline=false) const
------------------------------------------------------------------------

Evaluate the function symbolically or numerically.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L524

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1407-L1414

";

";

%feature("docstring") casadi::Callback::call_gen "

[INTERNAL] 
Call using a map.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1129

>  void casadi::Function::call_gen(const std::map< std::string, M > &arg, std::map< std::string, M > &res, bool always_inline, bool never_inline) const
------------------------------------------------------------------------
[INTERNAL] 
Call using a map.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1129

";

";

%feature("docstring") casadi::Callback::buf_in "

[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L564

>  std::vector<const double*> casadi::Function::buf_in(MapArg arg) const
------------------------------------------------------------------------
[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L564

";

";

%feature("docstring") casadi::Callback::buf_out "

[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L568

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L412-L425

>  vector< double * > casadi::Function::buf_out(MPrRes res) const
------------------------------------------------------------------------
[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L568

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L412-L425

";

";

%feature("docstring") casadi::Callback::mapaccum "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L482-L484

>  Function casadi::Function::mapaccum(casadi_int N, const Dict &opts=Dict()) const
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L482-L484

";

";

%feature("docstring") casadi::Callback::fold "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L475-L481

";

%feature("docstring") casadi::Callback::map "

";

";

%feature("docstring") casadi::Callback::generate_in "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1174-L1182

>  std::vector< DM > casadi::Function::generate_in(const std::string &fname)
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1174-L1182

";

";

%feature("docstring") casadi::Callback::generate_out "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1184-L1192

>  std::vector< DM > casadi::Function::generate_out(const std::string &fname)
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1184-L1192

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1194-L1197

>  void casadi::Function::export_code(const std::string &lang, std::ostream &stream, const Dict &options=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Export function in specific language.

Only allowed for (a subset of) SX/MX Functions

Extra doc: https://github.com/casadi/casadi/wiki/L_1wz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L904

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1194-L1197

";

";

%feature("docstring") casadi::Callback::serialize "

[INTERNAL] 
Serialize.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L893

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1211-L1215

>  std::string casadi::Function::serialize(const Dict &opts=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Serialize.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L893

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1211-L1215

";

";

%feature("docstring") casadi::Callback::save "

Save  Function to a file.

See: 
 load

Extra doc: https://github.com/casadi/casadi/wiki/L_240

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L900

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1206-L1209

";

%feature("docstring") casadi::Callback::sx_in "

Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L944

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1503-L1509

>  const vector< SX > casadi::Function::sx_in() const
------------------------------------------------------------------------

Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L944

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1503-L1509

";

";

%feature("docstring") casadi::Callback::mx_in "

Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L949

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1527-L1529

>  const vector< MX > casadi::Function::mx_in() const
------------------------------------------------------------------------

Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L949

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1527-L1529

";

";

%feature("docstring") casadi::Callback::sx_out "

Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L962

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1511-L1517

>  const vector< SX > casadi::Function::sx_out() const
------------------------------------------------------------------------

Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L962

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1511-L1517

";

";

%feature("docstring") casadi::Callback::mx_out "

Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L967

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1531-L1533

>  const vector< MX > casadi::Function::mx_out() const
------------------------------------------------------------------------

Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L967

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1531-L1533

";

";

%feature("docstring") casadi::Callback::nz_from_in "

Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L974

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1535-L1537

";

%feature("docstring") casadi::Callback::nz_from_out "

Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L975

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1539-L1541

";

%feature("docstring") casadi::Callback::nz_to_in "

Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L976

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1543-L1545

";

%feature("docstring") casadi::Callback::nz_to_out "

Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L977

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1547-L1549

";

%feature("docstring") casadi::Callback::convert_in "

Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L996

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1587-L1589

>  std::vector< MX > casadi::Function::convert_in(const MXDict &arg) const
------------------------------------------------------------------------

Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L996

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1587-L1589

";

";

%feature("docstring") casadi::Callback::convert_out "

Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L998

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1595-L1597

>  std::vector< MX > casadi::Function::convert_out(const MXDict &arg) const
------------------------------------------------------------------------

Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L998

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1595-L1597

";

";

%feature("docstring") casadi::Callback::has_spfwd "

Is the class able to propagate seeds through the algorithm?

Extra doc: https://github.com/casadi/casadi/wiki/L_1xl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1078

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1619-L1621

";

%feature("docstring") casadi::Callback::has_sprev "

Is the class able to propagate seeds through the algorithm?

Extra doc: https://github.com/casadi/casadi/wiki/L_1xl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1079

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1623-L1625

";

%feature("docstring") casadi::Callback::Callback "

Copy constructor (throws an error)

Extra doc: https://github.com/casadi/casadi/wiki/L_o3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L64

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L34-L36

>  casadi::Callback::Callback(const Callback &obj)
------------------------------------------------------------------------

Copy constructor (throws an error)

Extra doc: https://github.com/casadi/casadi/wiki/L_o3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L64

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L34-L36

";

";

%feature("docstring") casadi::Callback::~Callback "

Destructor.

Extra doc: https://github.com/casadi/casadi/wiki/L_o4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L69

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L46-L53

";

%feature("docstring") casadi::Callback::construct "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L38-L44

";

%feature("docstring") casadi::Callback::init "

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

Evaluate numerically, using temporary matrices and work vectors.

This signature is not thread-safe. For guaranteed thread-safety, use  
eval_buffer

Extra doc: https://github.com/casadi/casadi/wiki/L_o8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L106

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L62-L64

";

%feature("docstring") casadi::Callback::eval_buffer "

A copy-free low level interface.

In Python, you will be passed two tuples of memoryview objects

Extra doc: https://github.com/casadi/casadi/wiki/L_o9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L113

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L55-L58

";

%feature("docstring") casadi::Callback::get_n_in "

Get the number of inputs.

This function is called during construction.

Extra doc: https://github.com/casadi/casadi/wiki/L_oa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L122

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L66-L68

";

%feature("docstring") casadi::Callback::get_n_out "

Get the number of outputs.

This function is called during construction.

Extra doc: https://github.com/casadi/casadi/wiki/L_ob

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L129

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L70-L72

";

%feature("docstring") casadi::Callback::get_sparsity_in "

Get the sparsity of an input.

This function is called during construction.

Extra doc: https://github.com/casadi/casadi/wiki/L_oc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L136

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L74-L76

";

%feature("docstring") casadi::Callback::get_sparsity_out "

Get the sparsity of an output.

This function is called during construction.

Extra doc: https://github.com/casadi/casadi/wiki/L_od

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L143

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L78-L80

";

%feature("docstring") casadi::Callback::get_name_in "

Get the name of an input.

This function is called during construction.

Extra doc: https://github.com/casadi/casadi/wiki/L_oe

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L82-L84

";

%feature("docstring") casadi::Callback::get_name_out "

Get the name of an output.

This function is called during construction.

Extra doc: https://github.com/casadi/casadi/wiki/L_of

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L86-L88

";

%feature("docstring") casadi::Callback::uses_output "

Do the derivative functions need nondifferentiated outputs?

Extra doc: https://github.com/casadi/casadi/wiki/L_og

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.hpp#L162

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/callback.cpp#L90-L92

";

%feature("docstring") casadi::Callback::n_in "

Get the number of function inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L228

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L785-L787

";

%feature("docstring") casadi::Callback::n_out "

Get the number of function outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L789-L791

";

%feature("docstring") casadi::Callback::name_in "

Get input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L320

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L941-L947

>  const string & casadi::Function::name_in(casadi_int ind) const
------------------------------------------------------------------------

Get input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L320

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L941-L947

";

";

%feature("docstring") casadi::Callback::name_out "

Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L325

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L949-L955

>  const string & casadi::Function::name_out(casadi_int ind) const
------------------------------------------------------------------------

Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L325

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L949-L955

";

";

%feature("docstring") casadi::Callback::index_in "

Find the index for a string describing a particular entry of an input 

scheme.

example: schemeEntry(\"x_opt\") -> returns NLPSOL_X if 
FunctionInternal 
adheres to SCHEME_NLPINput

Extra doc: https://github.com/casadi/casadi/wiki/L_1vk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L333

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L925-L931

";

%feature("docstring") casadi::Callback::index_out "

Find the index for a string describing a particular entry of an output
 
scheme.

example: schemeEntry(\"x_opt\") -> returns NLPSOL_X if 
FunctionInternal 
adheres to SCHEME_NLPINput

Extra doc: https://github.com/casadi/casadi/wiki/L_1vl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L341

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L933-L939

";

%feature("docstring") casadi::Callback::default_in "

Get default input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L346

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1416-L1418

";

%feature("docstring") casadi::Callback::max_in "

Get largest input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L351

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1420-L1422

";

%feature("docstring") casadi::Callback::min_in "

Get smallest input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L356

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1424-L1426

";

%feature("docstring") casadi::Callback::nominal_in "

Get nominal input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L361

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1428-L1430

";

%feature("docstring") casadi::Callback::nominal_out "

Get nominal output value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L366

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1432-L1434

";

%feature("docstring") casadi::Callback::oracle "

Get oracle.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L407

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1816-L1822

";

%feature("docstring") casadi::Callback::wrap "

Wrap in an  Function instance consisting of only one  MX call.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L412

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1824-L1826

";

%feature("docstring") casadi::Callback::wrap_as_needed "

Wrap in a  Function with options.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L417

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1828-L1830

";

%feature("docstring") casadi::Callback::which_depends "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1745-L1751

";

%feature("docstring") casadi::Callback::print_dimensions "

Print dimensions of inputs and outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L432

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1093-L1095

";

%feature("docstring") casadi::Callback::print_options "

Print options to a stream.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L437

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1097-L1099

";

%feature("docstring") casadi::Callback::print_option "

Print all information there is to know about a certain option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L442

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1101-L1103

";

%feature("docstring") casadi::Callback::has_option "

Does a particular option exist.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L447

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1105-L1112

";

%feature("docstring") casadi::Callback::change_option "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1114-L1124

";

%feature("docstring") casadi::Callback::jacobian_old "

[DEPRECATED] Replaced by  Function::factory.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L466

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L854-L860

";

%feature("docstring") casadi::Callback::hessian_old "

[DEPRECATED] Replaced by  Function::factory.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L471

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L862-L870

";

%feature("docstring") casadi::Callback::jacobian "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L882-L888

";

%feature("docstring") casadi::Callback::rev "

[INTERNAL] 
Propagate sparsity backward with temporary memory allocation.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L633

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L458-L473

>  int casadi::Function::rev(std::vector< bvec_t * > arg, std::vector< bvec_t * > res) const
------------------------------------------------------------------------
[INTERNAL] 
Propagate sparsity backward with temporary memory allocation.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L633

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L458-L473

";

";

%feature("docstring") casadi::Callback::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization: 
Type of parallelization used: unroll|serial|openmp

Extra doc: https://github.com/casadi/casadi/wiki/L_1wh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L642

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L722-L729

";

%feature("docstring") casadi::Callback::slice "

returns a new function with a selection of inputs/outputs of the 
original

Extra doc: https://github.com/casadi/casadi/wiki/L_1wl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L754

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L713-L720

";

%feature("docstring") casadi::Callback::forward "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1077-L1083

";

%feature("docstring") casadi::Callback::reverse "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1085-L1091

";

%feature("docstring") casadi::Callback::jac_sparsity "

Get, if necessary generate, the sparsity of a single Jacobian block.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L830

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L909-L915

>  Sparsity casadi::Function::jac_sparsity(casadi_int oind, casadi_int iind, bool compact=false) const
------------------------------------------------------------------------

Get, if necessary generate, the sparsity of a single Jacobian block.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L830

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L909-L915

";

";

%feature("docstring") casadi::Callback::generate "

Export / Generate C code for the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L840

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1130-L1132

>  std::string casadi::Function::generate(const Dict &opts=Dict()) const
------------------------------------------------------------------------

Export / Generate C code for the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L840

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1130-L1132

";

";

%feature("docstring") casadi::Callback::generate_dependencies "

Export / Generate C code for the dependency function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ww

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L845

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1140-L1142

";

%feature("docstring") casadi::Callback::stats "

Get all statistics obtained at the end of the last evaluate call.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L932

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L894-L896

";

%feature("docstring") casadi::Callback::has_free "

Does the function have free variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1004

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1627-L1629

";

%feature("docstring") casadi::Callback::get_free "

Get free variables as a string.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1009

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1126-L1128

";

%feature("docstring") casadi::Callback::free_sx "

Get all the free variables of the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1014

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1603-L1609

";

%feature("docstring") casadi::Callback::free_mx "

Get all the free variables of the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1019

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1611-L1617

";

%feature("docstring") casadi::Callback::generate_lifted "

Extract the functions needed for the Lifted  Newton method.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1024

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1631-L1637

";

%feature("docstring") casadi::Callback::n_nodes "

Number of nodes in the algorithm.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1030

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1695-L1701

";

%feature("docstring") casadi::Callback::n_instructions "

Number of instruction in the algorithm (SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xe

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1035

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1639-L1645

";

%feature("docstring") casadi::Callback::instruction_id "

Identifier index of the instruction (SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1040

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1663-L1669

";

%feature("docstring") casadi::Callback::instruction_input "

Locations in the work vector for the inputs of the instruction.

(SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1047

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1671-L1677

";

%feature("docstring") casadi::Callback::instruction_constant "

Get the floating point output argument of an instruction (SXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1052

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1679-L1685

";

%feature("docstring") casadi::Callback::instruction_output "

Location in the work vector for the output of the instruction.

(SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1059

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1687-L1693

";

%feature("docstring") casadi::Callback::instruction_MX "

Get the  MX node corresponding to an instruction (MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1064

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1647-L1653

";

%feature("docstring") casadi::Callback::instructions_sx "

Get the SX node corresponding to all instructions (SXFunction)

Note: input and output instructions have no SX representation. This 
method 
returns nan for those instructions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1072

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1655-L1661

";

%feature("docstring") casadi::Callback::sz_arg "

Get required length of arg field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1085

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1025-L1025

";

%feature("docstring") casadi::Callback::sz_res "

Get required length of res field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1090

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1027-L1027

";

%feature("docstring") casadi::Callback::sz_iw "

Get required length of iw field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1095

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1029-L1029

";

%feature("docstring") casadi::Callback::sz_w "

Get required length of w field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1031-L1031

";

%feature("docstring") casadi::Callback::sz_work "

[INTERNAL] 
Get number of temporary variables needed.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1106

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1021-L1023

";

%feature("docstring") casadi::Callback::set_work "

[INTERNAL] 
Set the (persistent) work vectors.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1111

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1050-L1057

";

%feature("docstring") casadi::Callback::set_temp "

[INTERNAL] 
Set the (temporary) work vectors.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1117

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1059-L1066

";

%feature("docstring") casadi::Callback::setup "

[INTERNAL] 
Set the (persistent and temporary) work vectors.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1123

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1068-L1075

";

%feature("docstring") casadi::Callback::name "

Name of the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1137

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1244-L1250

";

%feature("docstring") casadi::Callback::is_a "

Check if the function is of a particular type.

Optionally check if name matches one of the base classes (default 
true)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1144

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1599-L1601

";

%feature("docstring") casadi::Callback::assert_size_in "

Assert that an input dimension is equal so some given value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1715-L1721

";

%feature("docstring") casadi::Callback::assert_size_out "

Assert that an output dimension is equal so some given value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1192

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1723-L1728

";

%feature("docstring") casadi::Callback::checkout "

Checkout a memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1703-L1705

";

%feature("docstring") casadi::Callback::release "

Release a memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1198

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1707-L1709

";

%feature("docstring") casadi::Callback::memory "

[INTERNAL] 
Get memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1202

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1711-L1713

";

%feature("docstring") casadi::Callback::get_function "

Get a dependency function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1213

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1761-L1767

>  Function casadi::Function::get_function(const std::string &name) const
------------------------------------------------------------------------

Get a dependency function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1213

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1761-L1767

";

";

%feature("docstring") casadi::Callback::has_function "

Check if a particular dependency exists.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1218

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1769-L1776

";

%feature("docstring") casadi::Callback::find "

Get a specific function embedded in the expression graphs.

Parameters:
-----------

max_depth: 
Maximum depth - a negative number indicates no maximum

name: 
Name of function needed

Extra doc: https://github.com/casadi/casadi/wiki/L_1y7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1796-L1813

>  Function casadi::Function::find(casadi_int max_depth, const std::string &name) const
------------------------------------------------------------------------

Get a specific function embedded in the expression graphs.

Parameters:
-----------

max_depth: 
Maximum depth - a negative number indicates no maximum

name: 
Name of function needed

Extra doc: https://github.com/casadi/casadi/wiki/L_1y7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1796-L1813

";

";

%feature("docstring") casadi::Callback::info "

Obtain information about function

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1236

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1842-L1844

";

%feature("docstring") casadi::Callback::class_name "

Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L100-L102

";

%feature("docstring") casadi::Callback::disp "

Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L104-L110

";

%feature("docstring") casadi::Callback::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::Callback::is_null "

Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L74-L76

";

%feature("docstring") casadi::Callback::__hash__ "

Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L130-L132

";


// File: classcasadi_1_1casadi__limits.xml
%feature("docstring") casadi::casadi_limits "

casadi_limits class

The following class, which acts as a complements to the standard 

numeric_limits class, allows specifying certain properties of scalar 

objects. The template can be specialized for e.g. symbolic scalars 
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

Form message string.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L68

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L68-L68

>  casadi::CasadiException::CasadiException(const std::string &msg)
------------------------------------------------------------------------

Form message string.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L68

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L68-L68

";

";

%feature("docstring") casadi::CasadiException::~CasadiException "

throw ()
Destructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L71

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L71-L71

";

%feature("docstring") casadi::CasadiException::what "

throw ()
Display error.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L74

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L74-L76

";


// File: classcasadi_1_1CasadiHandler.xml



// File: classcasadi_1_1CasadiMeta.xml
%feature("docstring") casadi::CasadiMeta "

Collects global CasADi meta information.

Joris Gillis

Extra doc: https://github.com/casadi/casadi/wiki/L_23k

C++ includes: casadi_meta.hpp
";


// File: classcasadi_1_1ClangCompiler.xml
%feature("docstring") casadi::ClangCompiler "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1CodeGenerator.xml
%feature("docstring") casadi::CodeGenerator "

Helper class for C code generation.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_ru

C++ includes: code_generator.hpp
";

%feature("docstring") casadi::CodeGenerator::CodeGenerator "

Constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L46

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L36-L169

";

%feature("docstring") casadi::CodeGenerator::add "

Add a function (name generated)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L49

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L277-L302

";

%feature("docstring") casadi::CodeGenerator::dump "

[INTERNAL] 
Generate a file, return code as string.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L57

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L304-L308

>  string casadi::CodeGenerator::dump()
------------------------------------------------------------------------
[INTERNAL] 
Generate a file, return code as string.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L57

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L304-L308

";

";

%feature("docstring") casadi::CodeGenerator::generate "

Generate file(s)

The \"prefix\" argument will be prepended to the generated files and 
may be
 a directory or a file prefix. returns the filename

Extra doc: https://github.com/casadi/casadi/wiki/L_rv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L66

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L384-L429

";

%feature("docstring") casadi::CodeGenerator::add_include "

Add an include file optionally using a relative path \"...\" instead 
of an 
absolute path <...>

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L69

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L719-L739

";

%feature("docstring") casadi::CodeGenerator::add_dependency "

[INTERNAL] 
Add a function dependency.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L74

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L196-L275

";

%feature("docstring") casadi::CodeGenerator::add_external "

[INTERNAL] 
Add an external function declaration.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L77

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L763-L765

";

%feature("docstring") casadi::CodeGenerator::shorthand "

[INTERNAL] 
Add/get a shorthand.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L83

>  std::string casadi::CodeGenerator::shorthand(const std::string &name, bool allow_adding=true)
------------------------------------------------------------------------
[INTERNAL] 
Add/get a shorthand.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L83

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
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L788-L790

";

%feature("docstring") casadi::CodeGenerator::get_constant "

[INTERNAL] 
Get or add an integer constant.

Extra doc: https://github.com/casadi/casadi/wiki/L_ry

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L104

>  casadi_int casadi::CodeGenerator::get_constant(const std::vector< casadi_int > &v, bool allow_adding=false)
------------------------------------------------------------------------
[INTERNAL] 
Get or add an integer constant.

Extra doc: https://github.com/casadi/casadi/wiki/L_ry

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L104

";

";

%feature("docstring") casadi::CodeGenerator::constant "

[INTERNAL]

>  string casadi::CodeGenerator::constant(casadi_int v)
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::CodeGenerator::constant_copy "

[INTERNAL] 
Represent an array constant; adding it when new.

Extra doc: https://github.com/casadi/casadi/wiki/L_s0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L114

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L861-L874

";

%feature("docstring") casadi::CodeGenerator::define_rom_double "

[INTERNAL] 
Allocate file scope double read-only memory.

Extra doc: https://github.com/casadi/casadi/wiki/L_s2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L124

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L504-L509

";

%feature("docstring") casadi::CodeGenerator::rom_double "

[INTERNAL] 
Access file scope double read-only memory.

Extra doc: https://github.com/casadi/casadi/wiki/L_s3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L129

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L511-L516

";

%feature("docstring") casadi::CodeGenerator::define_rom_integer "

[INTERNAL] 
Allocate file scope integer read-only memory.

Extra doc: https://github.com/casadi/casadi/wiki/L_s4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L134

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L518-L523

";

%feature("docstring") casadi::CodeGenerator::rom_integer "

[INTERNAL] 
Access file scope integer read-only memory.

Extra doc: https://github.com/casadi/casadi/wiki/L_s5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L139

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L525-L530

";

%feature("docstring") casadi::CodeGenerator::print_formatted "

[INTERNAL] 
Print without newline characters.

Extra doc: https://github.com/casadi/casadi/wiki/L_s8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L156

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1624-L1648

";

%feature("docstring") casadi::CodeGenerator::flush "

[INTERNAL] 
Flush the buffer to a stream of choice.

Extra doc: https://github.com/casadi/casadi/wiki/L_sa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L171

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1671-L1674

";

%feature("docstring") casadi::CodeGenerator::local "

[INTERNAL] 
Declare a local variable.

Extra doc: https://github.com/casadi/casadi/wiki/L_sb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L176

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1676-L1688

";

%feature("docstring") casadi::CodeGenerator::scope_enter "

[INTERNAL] 
Enter a local scope.

Extra doc: https://github.com/casadi/casadi/wiki/L_sc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L181

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L171-L174

";

%feature("docstring") casadi::CodeGenerator::scope_exit "

[INTERNAL] 
Exit a local scope.

Extra doc: https://github.com/casadi/casadi/wiki/L_sd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L186

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L176-L194

";

%feature("docstring") casadi::CodeGenerator::sx_work "

[INTERNAL] 
Declare a work vector element.

Extra doc: https://github.com/casadi/casadi/wiki/L_se

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L191

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1690-L1701

";

%feature("docstring") casadi::CodeGenerator::init_local "

[INTERNAL] 
Specify the default value for a local variable.

Extra doc: https://github.com/casadi/casadi/wiki/L_sf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L196

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1703-L1706

";

%feature("docstring") casadi::CodeGenerator::indent "

[INTERNAL] 
Increase indentation.

Extra doc: https://github.com/casadi/casadi/wiki/L_sg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L201-L201

";

%feature("docstring") casadi::CodeGenerator::unindent "

[INTERNAL] 
Decrease indentation.

Extra doc: https://github.com/casadi/casadi/wiki/L_sh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L206

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L206-L206

";

%feature("docstring") casadi::CodeGenerator::avoid_stack "

[INTERNAL] 
Avoid stack?

Extra doc: https://github.com/casadi/casadi/wiki/L_si

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L211

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L211-L211

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
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L228

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1709-L1804

";

%feature("docstring") casadi::CodeGenerator::dot "

[INTERNAL] 
Codegen inner product.

Extra doc: https://github.com/casadi/casadi/wiki/L_sm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L235

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1416-L1422

";

%feature("docstring") casadi::CodeGenerator::mv "

[INTERNAL] 
Codegen dense matrix-vector multiplication.

Extra doc: https://github.com/casadi/casadi/wiki/L_so

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L246

>  std::string casadi::CodeGenerator::mv(const std::string &x, casadi_int nrow_x, casadi_int ncol_x, const std::string &y, const std::string &z, bool tr)
------------------------------------------------------------------------
[INTERNAL] 
Codegen dense matrix-vector multiplication.

Extra doc: https://github.com/casadi/casadi/wiki/L_so

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L246

";

";

%feature("docstring") casadi::CodeGenerator::axpy "

[INTERNAL] 
Codegen axpy: y += a*x.

Extra doc: https://github.com/casadi/casadi/wiki/L_sp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L252

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1568-L1572

";

%feature("docstring") casadi::CodeGenerator::scal "

[INTERNAL] 
Codegen axpy: x *= alpha.

Extra doc: https://github.com/casadi/casadi/wiki/L_sq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1574-L1577

";

%feature("docstring") casadi::CodeGenerator::mtimes "

[INTERNAL] 
Codegen sparse matrix-matrix multiplication.

Extra doc: https://github.com/casadi/casadi/wiki/L_sr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L263

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1593-L1600

";

%feature("docstring") casadi::CodeGenerator::trilsolve "

[INTERNAL] 
Codegen lower triangular solve.

Extra doc: https://github.com/casadi/casadi/wiki/L_ss

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L271

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1602-L1607

";

%feature("docstring") casadi::CodeGenerator::triusolve "

[INTERNAL] 
Codegen upper triangular solve.

Extra doc: https://github.com/casadi/casadi/wiki/L_st

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L277

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1609-L1614

";

%feature("docstring") casadi::CodeGenerator::bilin "

[INTERNAL] 
Codegen bilinear form.

Extra doc: https://github.com/casadi/casadi/wiki/L_su

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L283

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1424-L1430

";

%feature("docstring") casadi::CodeGenerator::rank1 "

[INTERNAL] 
Rank-1 update.

Extra doc: https://github.com/casadi/casadi/wiki/L_sv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L289

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1432-L1440

";

%feature("docstring") casadi::CodeGenerator::logsumexp "

[INTERNAL] 
\\\\brie LogSumExp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1617-L1622

";

%feature("docstring") casadi::CodeGenerator::interpn "

[INTERNAL] 
Multilinear interpolation.

Extra doc: https://github.com/casadi/casadi/wiki/L_sw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L298

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1442-L1452

";

%feature("docstring") casadi::CodeGenerator::interpn_grad "

[INTERNAL] 
Multilinear interpolation - calculate gradient.

Extra doc: https://github.com/casadi/casadi/wiki/L_sx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L307

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1454-L1464

";

%feature("docstring") casadi::CodeGenerator::trans "

[INTERNAL] 
Transpose.

Extra doc: https://github.com/casadi/casadi/wiki/L_sy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L317

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1466-L1472

";

%feature("docstring") casadi::CodeGenerator::qr "

[INTERNAL] 
QR factorization.

Extra doc: https://github.com/casadi/casadi/wiki/L_sz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L323

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1839-L1846

";

%feature("docstring") casadi::CodeGenerator::qr_solve "

[INTERNAL] 
QR solve.

Extra doc: https://github.com/casadi/casadi/wiki/L_t0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L332

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1849-L1858

";

%feature("docstring") casadi::CodeGenerator::lsqr_solve "

[INTERNAL] 
\\\\brief LSQR solve

Extra doc: https://github.com/casadi/casadi/wiki/L_t1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L341

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1861-L1866

";

%feature("docstring") casadi::CodeGenerator::ldl "

[INTERNAL] 
LDL factorization.

Extra doc: https://github.com/casadi/casadi/wiki/L_t2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L347

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1869-L1875

";

%feature("docstring") casadi::CodeGenerator::ldl_solve "

[INTERNAL] 
LDL solve.

Extra doc: https://github.com/casadi/casadi/wiki/L_t3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L355

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1878-L1884

";

%feature("docstring") casadi::CodeGenerator::fmax "

[INTERNAL] 
fmax

Extra doc: https://github.com/casadi/casadi/wiki/L_t4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L363

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1887-L1890

";

%feature("docstring") casadi::CodeGenerator::fmin "

[INTERNAL] 
fmin

Extra doc: https://github.com/casadi/casadi/wiki/L_t5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L368

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1893-L1896

";

%feature("docstring") casadi::CodeGenerator::mmax "

[INTERNAL] 
mmax

Extra doc: https://github.com/casadi/casadi/wiki/L_t6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L373

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1935-L1938

";

%feature("docstring") casadi::CodeGenerator::mmin "

[INTERNAL] 
mmin

Extra doc: https://github.com/casadi/casadi/wiki/L_t7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L378

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1941-L1944

";

%feature("docstring") casadi::CodeGenerator::vfmax "

[INTERNAL] 
vfmax

Extra doc: https://github.com/casadi/casadi/wiki/L_ta

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L393

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1911-L1914

>  std::string casadi::CodeGenerator::vfmax(const std::string &x, const std::string &n, const std::string &y)
------------------------------------------------------------------------
[INTERNAL] 
vfmax

Extra doc: https://github.com/casadi/casadi/wiki/L_ta

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L393

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1911-L1914

";

";

%feature("docstring") casadi::CodeGenerator::vfmin "

[INTERNAL] 
vfmin

Extra doc: https://github.com/casadi/casadi/wiki/L_tb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L398

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1917-L1920

>  std::string casadi::CodeGenerator::vfmin(const std::string &x, const std::string &n, const std::string &y)
------------------------------------------------------------------------
[INTERNAL] 
vfmin

Extra doc: https://github.com/casadi/casadi/wiki/L_tb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L398

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1917-L1920

";

";

%feature("docstring") casadi::CodeGenerator::max "

[INTERNAL] 
max

Extra doc: https://github.com/casadi/casadi/wiki/L_tc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L403

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1923-L1926

";

%feature("docstring") casadi::CodeGenerator::min "

[INTERNAL] 
min

Extra doc: https://github.com/casadi/casadi/wiki/L_td

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L408

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1929-L1932

";

%feature("docstring") casadi::CodeGenerator::norm_inf "

[INTERNAL] 
norm_inf

Extra doc: https://github.com/casadi/casadi/wiki/L_te

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L413

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1959-L1962

";

%feature("docstring") casadi::CodeGenerator::max_viol "

[INTERNAL] 
max_viol

Extra doc: https://github.com/casadi/casadi/wiki/L_tf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L418

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1947-L1950

";

%feature("docstring") casadi::CodeGenerator::sum_viol "

[INTERNAL] 
sum_viol

Extra doc: https://github.com/casadi/casadi/wiki/L_tg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L424

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1953-L1956

";

%feature("docstring") casadi::CodeGenerator::bound_consistency "

[INTERNAL] 
bound_consistency

Extra doc: https://github.com/casadi/casadi/wiki/L_th

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L430

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1990-L1995

";

%feature("docstring") casadi::CodeGenerator::lb_eig "

[INTERNAL] 
lb_eig

Extra doc: https://github.com/casadi/casadi/wiki/L_ti

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L436

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1965-L1968

";

%feature("docstring") casadi::CodeGenerator::regularize "

[INTERNAL] 
regularize

Extra doc: https://github.com/casadi/casadi/wiki/L_tj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L441

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1971-L1974

";

%feature("docstring") casadi::CodeGenerator::convexify_eval "

[INTERNAL] 
convexify

Extra doc: https://github.com/casadi/casadi/wiki/L_tk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L446

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1977-L1981

";

%feature("docstring") casadi::CodeGenerator::low "

[INTERNAL] 
low

Extra doc: https://github.com/casadi/casadi/wiki/L_tl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L452

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1984-L1987

";

%feature("docstring") casadi::CodeGenerator::declare "

[INTERNAL] 
Declare a function.

Extra doc: https://github.com/casadi/casadi/wiki/L_tm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L458

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1474-L1485

";

%feature("docstring") casadi::CodeGenerator::comment "

[INTERNAL] 
Write a comment line (ignored if not verbose)

Extra doc: https://github.com/casadi/casadi/wiki/L_tn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L463

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1806-L1810

";

%feature("docstring") casadi::CodeGenerator::add_auxiliary "

[INTERNAL] 
Add a built-in auxiliary function.

Extra doc: https://github.com/casadi/casadi/wiki/L_tp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L547

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L880-L1284

";

%feature("docstring") casadi::CodeGenerator::add_io_sparsities "

[INTERNAL] 
Add io sparsity patterns of a function.

Extra doc: https://github.com/casadi/casadi/wiki/L_tq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L552

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1813-L1836

";

%feature("docstring") casadi::CodeGenerator::work "

[INTERNAL] 
Get work vector name from index

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L557

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L640-L648

";

%feature("docstring") casadi::CodeGenerator::workel "

[INTERNAL] 
Get work vector element from index

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L560

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L650-L656

";

%feature("docstring") casadi::CodeGenerator::print_vector "

[INTERNAL] 
Print real vector to a c file.

Extra doc: https://github.com/casadi/casadi/wiki/L_ts

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L575

>  void casadi::CodeGenerator::print_vector(std::ostream &s, const std::string &name, const std::vector< double > &v)
------------------------------------------------------------------------
[INTERNAL] 
Print real vector to a c file.

Extra doc: https://github.com/casadi/casadi/wiki/L_ts

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L575

";

";

%feature("docstring") casadi::CodeGenerator::copy "

[INTERNAL] 
Create a copy operation.

Extra doc: https://github.com/casadi/casadi/wiki/L_tt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L581

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1355-L1362

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
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L590

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1406-L1414

";

%feature("docstring") casadi::CodeGenerator::clear "

[INTERNAL] 
Create a fill operation.

Extra doc: https://github.com/casadi/casadi/wiki/L_tv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L595

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1384-L1390

";

%feature("docstring") casadi::CodeGenerator::arg "

[INTERNAL] 
Refer to argument.

Extra doc: https://github.com/casadi/casadi/wiki/L_tw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L600

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1392-L1394

";

%feature("docstring") casadi::CodeGenerator::res "

[INTERNAL] 
Refer to resuly.

Extra doc: https://github.com/casadi/casadi/wiki/L_tx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L605

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1396-L1398

";

%feature("docstring") casadi::CodeGenerator::mem "

[INTERNAL] 
Access thread-local memory.

Extra doc: https://github.com/casadi/casadi/wiki/L_ty

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L610

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1400-L1404

";

%feature("docstring") casadi::CodeGenerator::project "

[INTERNAL] 
Sparse assignment.

Extra doc: https://github.com/casadi/casadi/wiki/L_tz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L615

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1488-L1500

";

%feature("docstring") casadi::CodeGenerator::tri_project "

[INTERNAL] 
Project triangular part.

Extra doc: https://github.com/casadi/casadi/wiki/L_u0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L622

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1503-L1511

";

%feature("docstring") casadi::CodeGenerator::densify "

[INTERNAL] 
Densify.

Extra doc: https://github.com/casadi/casadi/wiki/L_u1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L628

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1514-L1522

";

%feature("docstring") casadi::CodeGenerator::sparsify "

[INTERNAL] 
Sparsify.

Extra doc: https://github.com/casadi/casadi/wiki/L_u2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L634

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1525-L1533

";

%feature("docstring") casadi::CodeGenerator::to_mex "

[INTERNAL] 
Create matrix in MATLAB's MEX format.

Extra doc: https://github.com/casadi/casadi/wiki/L_u3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L640

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1286-L1291

";

%feature("docstring") casadi::CodeGenerator::from_mex "

[INTERNAL] 
Get matrix from MATLAB's MEX format.

Extra doc: https://github.com/casadi/casadi/wiki/L_u4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L645

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1293-L1304

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
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L668

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L1998-L2001

";

%feature("docstring") casadi::CodeGenerator::cache_check "

[INTERNAL] 
cache check

Extra doc: https://github.com/casadi/casadi/wiki/L_u8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.hpp#L673

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/code_generator.cpp#L2004-L2009

";


// File: classcasadi_1_1Collocation.xml
%feature("docstring") casadi::Collocation "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Conic.xml
%feature("docstring") casadi::Conic "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Context.xml
%feature("docstring") casadi::Context "

Context.

Joris Gillis

Joris Gillis

Extra doc: https://github.com/casadi/casadi/wiki/L_1142

C++ includes: context.hpp
";

%feature("docstring") casadi::Context::Context "

Constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/context.hpp#L42
";


// File: classcasadi_1_1CustomNlpsol.xml
%feature("docstring") casadi::CustomNlpsol "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1CvodesSimulator.xml
%feature("docstring") casadi::CvodesSimulator "

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
 Specify all variables of a type 
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

A symbolic representation of a differential-algebraic equations model.

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

Independent variable (usually time)

Extra doc: https://github.com/casadi/casadi/wiki/L_5e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L93

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L63-L65

";

%feature("docstring") casadi::DaeBuilder::x "

Differential states.

Extra doc: https://github.com/casadi/casadi/wiki/L_5f

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L98

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L67-L69

";

%feature("docstring") casadi::DaeBuilder::ode "

Ordinary differential equations (ODE)

Extra doc: https://github.com/casadi/casadi/wiki/L_5g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L103

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L71-L73

";

%feature("docstring") casadi::DaeBuilder::z "

Algebraic variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_5h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L108

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L75-L77

";

%feature("docstring") casadi::DaeBuilder::alg "

Algebraic equations.

Extra doc: https://github.com/casadi/casadi/wiki/L_5i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L113

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L79-L81

";

%feature("docstring") casadi::DaeBuilder::q "

Quadrature states.

Extra doc: https://github.com/casadi/casadi/wiki/L_5j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L118

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L83-L85

";

%feature("docstring") casadi::DaeBuilder::quad "

Quadrature equations.

Extra doc: https://github.com/casadi/casadi/wiki/L_5k

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L123

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L87-L89

";

%feature("docstring") casadi::DaeBuilder::y "

Output variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_5l

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L128

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L91-L93

";

%feature("docstring") casadi::DaeBuilder::ydef "

Definitions of output variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_5m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L133

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L95-L97

";

%feature("docstring") casadi::DaeBuilder::u "

Free controls.

Extra doc: https://github.com/casadi/casadi/wiki/L_5n

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L99-L101

";

%feature("docstring") casadi::DaeBuilder::p "

Parameters.

Extra doc: https://github.com/casadi/casadi/wiki/L_5o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L143

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L103-L105

";

%feature("docstring") casadi::DaeBuilder::c "

Named constants.

Extra doc: https://github.com/casadi/casadi/wiki/L_5p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L148

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L107-L109

";

%feature("docstring") casadi::DaeBuilder::cdef "

Definitions of named constants.

Extra doc: https://github.com/casadi/casadi/wiki/L_5q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L111-L113

";

%feature("docstring") casadi::DaeBuilder::d "

Dependent parameters.

Extra doc: https://github.com/casadi/casadi/wiki/L_5r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L158

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L115-L117

";

%feature("docstring") casadi::DaeBuilder::ddef "

Definitions of dependent parameters.

Interdependencies are allowed but must be non-cyclic.

Extra doc: https://github.com/casadi/casadi/wiki/L_5s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L165

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L119-L121

";

%feature("docstring") casadi::DaeBuilder::w "

Dependent variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_5t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L170

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L123-L125

";

%feature("docstring") casadi::DaeBuilder::wdef "

Dependent variables and corresponding definitions.

Interdependencies are allowed but must be non-cyclic.

Extra doc: https://github.com/casadi/casadi/wiki/L_5u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L177

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L127-L129

";

%feature("docstring") casadi::DaeBuilder::aux "

Auxiliary variables: Used e.g. to define functions.

Extra doc: https://github.com/casadi/casadi/wiki/L_5v

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L182

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L131-L133

";

%feature("docstring") casadi::DaeBuilder::init_lhs "

Initial conditions, left-hand-side.

Extra doc: https://github.com/casadi/casadi/wiki/L_5w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L187

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L135-L137

";

%feature("docstring") casadi::DaeBuilder::init_rhs "

Initial conditions, right-hand-side.

Extra doc: https://github.com/casadi/casadi/wiki/L_5x

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L192

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L139-L141

";

%feature("docstring") casadi::DaeBuilder::when_cond "

When statement: triggering condition.

Extra doc: https://github.com/casadi/casadi/wiki/L_5y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L197

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L143-L145

";

%feature("docstring") casadi::DaeBuilder::when_lhs "

When statement: left-hand-side.

Extra doc: https://github.com/casadi/casadi/wiki/L_5z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L202

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L147-L149

";

%feature("docstring") casadi::DaeBuilder::when_rhs "

When statement: right-hand-side.

Extra doc: https://github.com/casadi/casadi/wiki/L_60

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L207

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L151-L153

";

%feature("docstring") casadi::DaeBuilder::has_t "

Is there a time variable?

Extra doc: https://github.com/casadi/casadi/wiki/L_64

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L231

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L182-L184

";

%feature("docstring") casadi::DaeBuilder::nx "

Differential states.

Extra doc: https://github.com/casadi/casadi/wiki/L_65

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L236

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L186-L188

";

%feature("docstring") casadi::DaeBuilder::nz "

Algebraic variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_66

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L241

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L190-L192

";

%feature("docstring") casadi::DaeBuilder::nq "

Quadrature states.

Extra doc: https://github.com/casadi/casadi/wiki/L_67

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L246

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L194-L196

";

%feature("docstring") casadi::DaeBuilder::ny "

Output variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_68

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L251

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L198-L200

";

%feature("docstring") casadi::DaeBuilder::nu "

Free controls.

Extra doc: https://github.com/casadi/casadi/wiki/L_69

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L256

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L202-L204

";

%feature("docstring") casadi::DaeBuilder::np "

Parameters.

Extra doc: https://github.com/casadi/casadi/wiki/L_6a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L261

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L206-L208

";

%feature("docstring") casadi::DaeBuilder::nc "

Named constants.

Extra doc: https://github.com/casadi/casadi/wiki/L_6b

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L266

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L210-L212

";

%feature("docstring") casadi::DaeBuilder::nd "

Dependent parameters.

Extra doc: https://github.com/casadi/casadi/wiki/L_6c

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L271

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L214-L216

";

%feature("docstring") casadi::DaeBuilder::nw "

Dependent variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_6d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L276

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L218-L220

";

%feature("docstring") casadi::DaeBuilder::add_t "

Add an independent variable (time)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L284

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L427-L432

";

%feature("docstring") casadi::DaeBuilder::add_p "

Add a new parameter.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L287

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L434-L441

";

%feature("docstring") casadi::DaeBuilder::add_u "

Add a new control.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L290

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L443-L450

";

%feature("docstring") casadi::DaeBuilder::add_x "

Add a new differential state.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L452-L459

";

%feature("docstring") casadi::DaeBuilder::add_z "

Add a new algebraic variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L296

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L461-L468

";

%feature("docstring") casadi::DaeBuilder::add_q "

Add a new quadrature state.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L299

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L470-L477

";

%feature("docstring") casadi::DaeBuilder::add_c "

Add a new constant.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L302

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L479-L486

";

%feature("docstring") casadi::DaeBuilder::add_d "

Add a new dependent parameter.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L305

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L488-L495

";

%feature("docstring") casadi::DaeBuilder::add_w "

Add a new dependent variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L308

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L497-L504

";

%feature("docstring") casadi::DaeBuilder::add_y "

Add a new output.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L311

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L506-L513

";

%feature("docstring") casadi::DaeBuilder::set_ode "

Specify the ordinary differential equation for a state.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L314

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L532-L538

";

%feature("docstring") casadi::DaeBuilder::set_alg "

Specificy the residual equation for an algebraic variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L317

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L540-L546

";

%feature("docstring") casadi::DaeBuilder::add_aux "

Add an auxiliary variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L320

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L515-L519

";

%feature("docstring") casadi::DaeBuilder::add_init "

Add an initial equation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L323

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L521-L524

";

%feature("docstring") casadi::DaeBuilder::add_when "

Add a when statement.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L326

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L526-L530

";

%feature("docstring") casadi::DaeBuilder::sanity_check "

Check if dimensions match.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L329

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L548-L554

";

%feature("docstring") casadi::DaeBuilder::clear_in "

Clear input variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L362

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L419-L425

";

%feature("docstring") casadi::DaeBuilder::eliminate_w "

Eliminate all dependent variables.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L365

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L608-L614

";

%feature("docstring") casadi::DaeBuilder::lift "

Lift problem formulation by extracting shared subexpressions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L368

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L616-L622

";

%feature("docstring") casadi::DaeBuilder::eliminate_quad "

Eliminate quadrature states and turn them into ODE states.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L371

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L230-L236

";

%feature("docstring") casadi::DaeBuilder::sort_d "

Sort dependent parameters.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L374

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L238-L244

";

%feature("docstring") casadi::DaeBuilder::sort_w "

Sort dependent variables.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L377

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L246-L252

";

%feature("docstring") casadi::DaeBuilder::sort_z "

Sort algebraic variables.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L380

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L254-L260

";

%feature("docstring") casadi::DaeBuilder::prune "

Prune unused controls.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L383

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L262-L268

";

%feature("docstring") casadi::DaeBuilder::tear "

Identify iteration variables and residual equations using naming 

convention.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L386

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L270-L276

";

%feature("docstring") casadi::DaeBuilder::add_fun "

Add an external function.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L403

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L747-L750

>  Function casadi::DaeBuilder::add_fun(const std::string &name, const Importer &compiler, const Dict &opts=Dict())
------------------------------------------------------------------------

Add an external function.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L403

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L747-L750

";

";

%feature("docstring") casadi::DaeBuilder::has_fun "

Does a particular function already exist?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L407

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L752-L759

";

%feature("docstring") casadi::DaeBuilder::fun "

Get all functions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L413

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L791-L793

>  std::vector< Function > casadi::DaeBuilder::fun() const
------------------------------------------------------------------------

Get all functions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L413

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L791-L793

";

";

%feature("docstring") casadi::DaeBuilder::gather_fun "

Collect embedded functions from the expression graph.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L416

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L770-L789

";

%feature("docstring") casadi::DaeBuilder::parse_fmi "

Import existing problem from FMI/XML

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L423

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L423-L423

";

%feature("docstring") casadi::DaeBuilder::load_fmi_description "

Import problem description from FMI or XML.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L426

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L222-L228

";

%feature("docstring") casadi::DaeBuilder::add_lc "

Add a named linear combination of output expressions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L429

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L688-L695

";

%feature("docstring") casadi::DaeBuilder::create "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L456

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L719-L726

>  Function casadi::DaeBuilder::create(const std::string &name, const Dict &opts=Dict()) const
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L456

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L719-L726

";

";

%feature("docstring") casadi::DaeBuilder::var "

[INTERNAL] 
Get variable expressions by index.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L678

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L869-L876

>  std::vector< MX > casadi::DaeBuilder::var(const std::vector< size_t > &ind) const
------------------------------------------------------------------------
[INTERNAL] 
Get variable expressions by index.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L678

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L869-L876

";

";

%feature("docstring") casadi::DaeBuilder::beq "

Get/set the binding equation for a variable

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L474

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L591-L598

";

%feature("docstring") casadi::DaeBuilder::set_beq "

Get/set the binding equation for a variable

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L475

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L600-L606

";

%feature("docstring") casadi::DaeBuilder::value_reference "

Get/set value reference

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L480

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L624-L626

";

%feature("docstring") casadi::DaeBuilder::set_value_reference "

Get/set value reference

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L481

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L628-L630

";

%feature("docstring") casadi::DaeBuilder::description "

Get/set description

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L486

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L632-L634

";

%feature("docstring") casadi::DaeBuilder::set_description "

Get/set description

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L487

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L636-L638

";

%feature("docstring") casadi::DaeBuilder::type "

Get/set the type

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L492

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L640-L642

";

%feature("docstring") casadi::DaeBuilder::set_type "

Get/set the type

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L493

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L644-L646

";

%feature("docstring") casadi::DaeBuilder::causality "

Get/set the causality

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L498

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L648-L650

";

%feature("docstring") casadi::DaeBuilder::set_causality "

Get/set the causality

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L499

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L652-L654

";

%feature("docstring") casadi::DaeBuilder::variability "

Get/set the variability

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L504

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L656-L658

";

%feature("docstring") casadi::DaeBuilder::set_variability "

Get/set the variability

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L505

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L660-L662

";

%feature("docstring") casadi::DaeBuilder::initial "

Get/set the initial property

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L510

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L664-L666

";

%feature("docstring") casadi::DaeBuilder::set_initial "

Get/set the initial property

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L511

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L668-L670

";

%feature("docstring") casadi::DaeBuilder::unit "

Get/set the unit

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L516

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L672-L674

";

%feature("docstring") casadi::DaeBuilder::set_unit "

Get/set the unit

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L517

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L676-L678

";

%feature("docstring") casadi::DaeBuilder::display_unit "

Get/set the display unit

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L522

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L680-L682

";

%feature("docstring") casadi::DaeBuilder::set_display_unit "

Get/set the display unit

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L523

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L684-L686

";

%feature("docstring") casadi::DaeBuilder::variable "

[INTERNAL] 
Access a variable by index

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L662

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L839-L845

>  const Variable & casadi::DaeBuilder::variable(size_t ind) const
------------------------------------------------------------------------
[INTERNAL] 
Access a variable by index

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L662

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L839-L845

";

";

%feature("docstring") casadi::DaeBuilder::type_name "

Readable name of the class.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L74

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L74-L74

";

%feature("docstring") casadi::DaeBuilder::DaeBuilder "

Construct a  DaeBuilder instance.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L80

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L54-L57

>  casadi::DaeBuilder::DaeBuilder(const std::string &name, const std::string &path=\"\", const Dict &opts=Dict())
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L694

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L906-L915

>  std::vector< std::string > casadi::DaeBuilder::name(const std::vector< size_t > &ind) const
------------------------------------------------------------------------
[INTERNAL] 
Get variable names by indices.

Extra doc: https://github.com/casadi/casadi/wiki/L_6i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L694

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L906-L915

";

";

%feature("docstring") casadi::DaeBuilder::outputs "

Model structure: outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_61

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L213

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L155-L162

";

%feature("docstring") casadi::DaeBuilder::derivatives "

Model structure: derivatives.

Extra doc: https://github.com/casadi/casadi/wiki/L_62

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L218

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L164-L171

";

%feature("docstring") casadi::DaeBuilder::initial_unknowns "

Model structure: initial unknowns.

Extra doc: https://github.com/casadi/casadi/wiki/L_63

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L223

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L173-L180

";

%feature("docstring") casadi::DaeBuilder::dependent_fun "

Construct a function for evaluating dependent parameters.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L459

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L804-L813

";

%feature("docstring") casadi::DaeBuilder::der "

Get the time derivative of an expression, single variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L530

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L565-L578

>  std::string casadi::DaeBuilder::der(const std::string &name) const
------------------------------------------------------------------------

Get the time derivative of an expression, single variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L530

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L565-L578

";

";

%feature("docstring") casadi::DaeBuilder::attribute "

Get an attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L577

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L926-L934

>  std::vector< double > casadi::DaeBuilder::attribute(const std::string &a, const std::vector< std::string > &name) const
------------------------------------------------------------------------

Get an attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L577

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L926-L934

";

";

%feature("docstring") casadi::DaeBuilder::set_attribute "

Set an attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L580

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L944-L951

>  void casadi::DaeBuilder::set_attribute(const std::string &a, const std::vector< std::string > &name, const std::vector< double > &val)
------------------------------------------------------------------------

Set an attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L580

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L944-L951

";

";

%feature("docstring") casadi::DaeBuilder::min "

Get the lower bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L584

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L962-L969

>  std::vector< double > casadi::DaeBuilder::min(const std::vector< std::string > &name) const
------------------------------------------------------------------------

Get the lower bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L584

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L962-L969

";

";

%feature("docstring") casadi::DaeBuilder::set_min "

Set the lower bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L587

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L979-L985

>  void casadi::DaeBuilder::set_min(const std::vector< std::string > &name, const std::vector< double > &val)
------------------------------------------------------------------------

Set the lower bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L587

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L979-L985

";

";

%feature("docstring") casadi::DaeBuilder::max "

Get the upper bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L590

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L996-L1003

>  std::vector< double > casadi::DaeBuilder::max(const std::vector< std::string > &name) const
------------------------------------------------------------------------

Get the upper bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L590

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L996-L1003

";

";

%feature("docstring") casadi::DaeBuilder::set_max "

Set the upper bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L593

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1013-L1019

>  void casadi::DaeBuilder::set_max(const std::vector< std::string > &name, const std::vector< double > &val)
------------------------------------------------------------------------

Set the upper bound.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L593

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1013-L1019

";

";

%feature("docstring") casadi::DaeBuilder::nominal "

Get the nominal value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L596

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1030-L1037

>  std::vector< double > casadi::DaeBuilder::nominal(const std::vector< std::string > &name) const
------------------------------------------------------------------------

Get the nominal value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L596

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1030-L1037

";

";

%feature("docstring") casadi::DaeBuilder::set_nominal "

Set the nominal value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L599

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1047-L1053

>  void casadi::DaeBuilder::set_nominal(const std::vector< std::string > &name, const std::vector< double > &val)
------------------------------------------------------------------------

Set the nominal value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L599

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1047-L1053

";

";

%feature("docstring") casadi::DaeBuilder::start "

Get the start attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L602

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1064-L1071

>  std::vector< double > casadi::DaeBuilder::start(const std::vector< std::string > &name) const
------------------------------------------------------------------------

Get the start attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L602

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1064-L1071

";

";

%feature("docstring") casadi::DaeBuilder::set_start "

Set the start attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L605

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1081-L1087

>  void casadi::DaeBuilder::set_start(const std::vector< std::string > &name, const std::vector< double > &val)
------------------------------------------------------------------------

Set the start attribute.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L605

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1081-L1087

";

";

%feature("docstring") casadi::DaeBuilder::set "

Set the current value (string)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L611

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1121-L1128

>  void casadi::DaeBuilder::set(const std::vector< std::string > &name, const std::vector< std::string > &val)
------------------------------------------------------------------------

Set the current value (string)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L611

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1121-L1128

";

";

%feature("docstring") casadi::DaeBuilder::get "

Evaluate the values for a set of variables at the initial time.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L614

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1134-L1149

>  std::vector< GenericType > casadi::DaeBuilder::get(const std::vector< std::string > &name) const
------------------------------------------------------------------------

Evaluate the values for a set of variables at the initial time.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L614

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1134-L1149

";

";

%feature("docstring") casadi::DaeBuilder::add_variable "

[INTERNAL] 
Add a variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L651

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L296-L302

>  size_t casadi::DaeBuilder::add_variable(const std::string &name, const Variable &var)
------------------------------------------------------------------------
[INTERNAL] 
Add a variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L651

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L296-L302

";

";

%feature("docstring") casadi::DaeBuilder::add_variable_new "

Add a new variable from symbolic expressions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L632

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L331-L335

>  size_t casadi::DaeBuilder::add_variable_new(const MX &new_v)
------------------------------------------------------------------------

Add a new variable from symbolic expressions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L632

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L331-L335

";

";

%feature("docstring") casadi::DaeBuilder::has_variable "

Check if a particular variable exists.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L635

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L278-L285

";

%feature("docstring") casadi::DaeBuilder::all_variables "

Get a list of all variables.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L638

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L287-L294

";

%feature("docstring") casadi::DaeBuilder::oracle "

Get the (cached) oracle, SX or  MX.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L641

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L795-L802

";

%feature("docstring") casadi::DaeBuilder::jac_sparsity "

Get Jacobian sparsity.

Extra doc: https://github.com/casadi/casadi/wiki/L_6g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L646

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L1151-L1159

";

%feature("docstring") casadi::DaeBuilder::find "

[INTERNAL] 
Get indices of variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L684

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L887-L894

>  std::vector< size_t > casadi::DaeBuilder::find(const std::vector< std::string > &name) const
------------------------------------------------------------------------
[INTERNAL] 
Get indices of variable.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.hpp#L684

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dae_builder.cpp#L887-L894

";

";

%feature("docstring") casadi::DaeBuilder::class_name "

Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L100-L102

";

%feature("docstring") casadi::DaeBuilder::disp "

Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L104-L110

";

%feature("docstring") casadi::DaeBuilder::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::DaeBuilder::is_null "

Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L74-L76

";

%feature("docstring") casadi::DaeBuilder::__hash__ "

Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L130-L132

";


// File: classcasadi_1_1DeserializerBase.xml
%feature("docstring") casadi::DeserializerBase "

C++ includes: serializer.hpp
";

%feature("docstring") casadi::DeserializerBase::DeserializerBase "

[INTERNAL] ";


// File: classcasadi_1_1DeserializingStream.xml
%feature("docstring") casadi::DeserializingStream "

Helper class for Serialization.

Joris Gillis

Extra doc: https://github.com/casadi/casadi/wiki/L_ak

C++ includes: serializing_stream.hpp
";

%feature("docstring") casadi::DeserializingStream::DeserializingStream "

";

";

%feature("docstring") casadi::DeserializingStream::unpack "

";

";


// File: classcasadi_1_1DllLibrary.xml
%feature("docstring") casadi::DllLibrary "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1DotWriter.xml
%feature("docstring") casadi::DotWriter "

C++ includes: dot_writer.hpp
";

%feature("docstring") casadi::DotWriter::DotWriter "

Constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dot_writer.hpp#L41

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dot_writer.cpp#L30-L34

";


// File: classcasadi_1_1Dple.xml
%feature("docstring") casadi::Dple "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Expm.xml
%feature("docstring") casadi::Expm "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1FastNewton.xml
%feature("docstring") casadi::FastNewton "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1FileDeserializer.xml
%feature("docstring") casadi::FileDeserializer "

C++ includes: serializer.hpp
";

%feature("docstring") casadi::FileDeserializer::FileDeserializer "

Advanced deserialization of CasADi objects.

See: 
 FileSerializer

Extra doc: https://github.com/casadi/casadi/wiki/L_7t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.hpp#L250

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.cpp#L107-L113

";


// File: classcasadi_1_1FileSerializer.xml
%feature("docstring") casadi::FileSerializer "

C++ includes: serializer.hpp
";

%feature("docstring") casadi::FileSerializer::FileSerializer "

Advanced serialization of CasADi objects.

See: 
 StringSerializer,  FileDeserializer

Extra doc: https://github.com/casadi/casadi/wiki/L_7q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.hpp#L221

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.cpp#L42-L50

";


// File: classcasadi_1_1FixedStepIntegrator.xml
%feature("docstring") casadi::FixedStepIntegrator "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Function.xml
%feature("docstring") casadi::Function "

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
|                  |                 | ght*nf<=(1-ad_we |                  |
|                  |                 | ight)*na is used |                  |
|                  |                 | where nf and na  |                  |
|                  |                 | are estimates of |                  |
|                  |                 | the number of    |                  |
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
| print_time       | OT_BOOL         | print            | casadi::Function |
|                  |                 | information      | Internal         |
|                  |                 | about execution  |                  |
|                  |                 | time. Implies    |                  |
|                  |                 | record_time.     |                  |
+------------------+-----------------+------------------+------------------+
| record_time      | OT_BOOL         | record           | casadi::Function |
|                  |                 | information      | Internal         |
|                  |                 | about execution  |                  |
|                  |                 | time, for        |                  |
|                  |                 | retrieval with   |                  |
|                  |                 | stats() .        |                  |
+------------------+-----------------+------------------+------------------+
| regularity_check | OT_BOOL         | Throw exceptions | casadi::Function |
|                  |                 | when NaN or Inf  | Internal         |
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
| verbose          | OT_BOOL         | Verbose          | casadi::Function |
|                  |                 | evaluation  for  | Internal         |
|                  |                 | debugging        |                  |
+------------------+-----------------+------------------+------------------+

C++ includes: function.hpp
";

%feature("docstring") casadi::Function::Function "

Construct from a file.

Extra doc: https://github.com/casadi/casadi/wiki/L_1uz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L75

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L92-L94

>  casadi::Function::Function(const std::string &fname)
------------------------------------------------------------------------

Construct from a file.

Extra doc: https://github.com/casadi/casadi/wiki/L_1uz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L75

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L92-L94

";

";

%feature("docstring") casadi::Function::jit "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L260-L272

>  Function casadi::Function::jit(const std::string &name, const std::string &body, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const std::vector< Sparsity > &sparsity_in, const std::vector< Sparsity > &sparsity_out, const Dict &opts=Dict())
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L260-L272

";

";

%feature("docstring") casadi::Function::expand "

Expand a function to SX.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L207

>  Function casadi::Function::expand(const std::string &name, const Dict &opts=Dict()) const
------------------------------------------------------------------------

Expand a function to SX.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L207

";

";

%feature("docstring") casadi::Function::size1_in "

Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240-L240

>  casadi_int casadi::Function::size1_in(const std::string &iname) const
------------------------------------------------------------------------

Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L240-L240

";

";

%feature("docstring") casadi::Function::size2_in "

Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242-L242

>  casadi_int casadi::Function::size2_in(const std::string &iname) const
------------------------------------------------------------------------

Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L242-L242

";

";

%feature("docstring") casadi::Function::size_in "

Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244-L246

>  std::pair<casadi_int, casadi_int> casadi::Function::size_in(const std::string &iname) const
------------------------------------------------------------------------

Get input dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1va

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L244-L246

";

";

%feature("docstring") casadi::Function::size1_out "

Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254-L254

>  casadi_int casadi::Function::size1_out(const std::string &oname) const
------------------------------------------------------------------------

Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L254-L254

";

";

%feature("docstring") casadi::Function::size2_out "

Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256-L256

>  casadi_int casadi::Function::size2_out(const std::string &oname) const
------------------------------------------------------------------------

Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L256-L256

";

";

%feature("docstring") casadi::Function::size_out "

Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258-L260

>  std::pair<casadi_int, casadi_int> casadi::Function::size_out(const std::string &oname) const
------------------------------------------------------------------------

Get output dimension.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L258-L260

";

";

%feature("docstring") casadi::Function::nnz_in "

Get number of input nonzeros.

For a particular input or for all of the inputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L271

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L271-L271

>  casadi_int casadi::Function::nnz_in(const std::string &iname) const
------------------------------------------------------------------------

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

Get number of output nonzeros.

For a particular output or for all of the outputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L282

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L282-L282

>  casadi_int casadi::Function::nnz_out(const std::string &oname) const
------------------------------------------------------------------------

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

Get number of input elements.

For a particular input or for all of the inputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1ve

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L293-L293

>  casadi_int casadi::Function::numel_in(const std::string &iname) const
------------------------------------------------------------------------

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

Get number of output elements.

For a particular output or for all of the outputs

Extra doc: https://github.com/casadi/casadi/wiki/L_1vf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L304

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L304-L304

>  casadi_int casadi::Function::numel_out(const std::string &oname) const
------------------------------------------------------------------------

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

Get sparsity of a given input.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L373

>  const Sparsity& casadi::Function::sparsity_in(const std::string &iname) const
------------------------------------------------------------------------

Get sparsity of a given input.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L373

";

";

%feature("docstring") casadi::Function::sparsity_out "

Get sparsity of a given output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L381

>  const Sparsity& casadi::Function::sparsity_out(const std::string &iname) const
------------------------------------------------------------------------

Get sparsity of a given output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L381

";

";

%feature("docstring") casadi::Function::is_diff_in "

Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L390

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1005-L1011

>  std::vector< bool > casadi::Function::is_diff_in() const
------------------------------------------------------------------------

Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L390

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1005-L1011

";

";

%feature("docstring") casadi::Function::is_diff_out "

Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L391

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1013-L1019

>  std::vector< bool > casadi::Function::is_diff_out() const
------------------------------------------------------------------------

Get differentiability of inputs/output.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L391

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1013-L1019

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

Evaluate the function symbolically or numerically.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L524

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1407-L1414

>  void casadi::Function::call(const MXDict &arg, MXDict &res, bool always_inline=false, bool never_inline=false) const
------------------------------------------------------------------------

Evaluate the function symbolically or numerically.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L524

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1407-L1414

";

";

%feature("docstring") casadi::Function::call_gen "

[INTERNAL] 
Call using a map.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1129

>  void casadi::Function::call_gen(const std::map< std::string, M > &arg, std::map< std::string, M > &res, bool always_inline, bool never_inline) const
------------------------------------------------------------------------
[INTERNAL] 
Call using a map.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1129

";

";

%feature("docstring") casadi::Function::buf_in "

[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L564

>  std::vector<const double*> casadi::Function::buf_in(MapArg arg) const
------------------------------------------------------------------------
[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L564

";

";

%feature("docstring") casadi::Function::buf_out "

[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L568

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L412-L425

>  vector< double * > casadi::Function::buf_out(MPrRes res) const
------------------------------------------------------------------------
[INTERNAL] 
Supported arguments for numerical evaluation and converters.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L568

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L412-L425

";

";

%feature("docstring") casadi::Function::mapaccum "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L482-L484

>  Function casadi::Function::mapaccum(casadi_int N, const Dict &opts=Dict()) const
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L482-L484

";

";

%feature("docstring") casadi::Function::fold "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L475-L481

";

%feature("docstring") casadi::Function::map "

";

";

%feature("docstring") casadi::Function::generate_in "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1174-L1182

>  std::vector< DM > casadi::Function::generate_in(const std::string &fname)
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1174-L1182

";

";

%feature("docstring") casadi::Function::generate_out "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1184-L1192

>  std::vector< DM > casadi::Function::generate_out(const std::string &fname)
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1184-L1192

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1194-L1197

>  void casadi::Function::export_code(const std::string &lang, std::ostream &stream, const Dict &options=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Export function in specific language.

Only allowed for (a subset of) SX/MX Functions

Extra doc: https://github.com/casadi/casadi/wiki/L_1wz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L904

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1194-L1197

";

";

%feature("docstring") casadi::Function::serialize "

[INTERNAL] 
Serialize.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L893

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1211-L1215

>  std::string casadi::Function::serialize(const Dict &opts=Dict()) const
------------------------------------------------------------------------
[INTERNAL] 
Serialize.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L893

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1211-L1215

";

";

%feature("docstring") casadi::Function::save "

Save  Function to a file.

See: 
 load

Extra doc: https://github.com/casadi/casadi/wiki/L_240

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L900

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1206-L1209

";

%feature("docstring") casadi::Function::sx_in "

Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L944

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1503-L1509

>  const vector< SX > casadi::Function::sx_in() const
------------------------------------------------------------------------

Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L944

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1503-L1509

";

";

%feature("docstring") casadi::Function::mx_in "

Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L949

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1527-L1529

>  const vector< MX > casadi::Function::mx_in() const
------------------------------------------------------------------------

Get symbolic primitives equivalent to the input expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L949

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1527-L1529

";

";

%feature("docstring") casadi::Function::sx_out "

Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L962

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1511-L1517

>  const vector< SX > casadi::Function::sx_out() const
------------------------------------------------------------------------

Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L962

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1511-L1517

";

";

%feature("docstring") casadi::Function::mx_out "

Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L967

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1531-L1533

>  const vector< MX > casadi::Function::mx_out() const
------------------------------------------------------------------------

Get symbolic primitives equivalent to the output expressions.

There is no guarantee that subsequent calls return unique answers

Extra doc: https://github.com/casadi/casadi/wiki/L_1x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L967

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1531-L1533

";

";

%feature("docstring") casadi::Function::nz_from_in "

Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L974

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1535-L1537

";

%feature("docstring") casadi::Function::nz_from_out "

Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L975

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1539-L1541

";

%feature("docstring") casadi::Function::nz_to_in "

Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L976

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1543-L1545

";

%feature("docstring") casadi::Function::nz_to_out "

Convert from/to flat vector of input/output nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L977

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1547-L1549

";

%feature("docstring") casadi::Function::convert_in "

Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L996

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1587-L1589

>  std::vector< MX > casadi::Function::convert_in(const MXDict &arg) const
------------------------------------------------------------------------

Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L996

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1587-L1589

";

";

%feature("docstring") casadi::Function::convert_out "

Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L998

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1595-L1597

>  std::vector< MX > casadi::Function::convert_out(const MXDict &arg) const
------------------------------------------------------------------------

Convert from/to input/output lists/map.

Will raise an error when an unknown key is used or a list has 
incorrect 
size. Does not perform sparsity checking.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L998

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1595-L1597

";

";

%feature("docstring") casadi::Function::has_spfwd "

Is the class able to propagate seeds through the algorithm?

Extra doc: https://github.com/casadi/casadi/wiki/L_1xl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1078

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1619-L1621

";

%feature("docstring") casadi::Function::has_sprev "

Is the class able to propagate seeds through the algorithm?

Extra doc: https://github.com/casadi/casadi/wiki/L_1xl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1079

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1623-L1625

";

%feature("docstring") casadi::Function::~Function "

Destructor.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L61-L62

";

%feature("docstring") casadi::Function::n_in "

Get the number of function inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L228

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L785-L787

";

%feature("docstring") casadi::Function::n_out "

Get the number of function outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1v9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L789-L791

";

%feature("docstring") casadi::Function::name_in "

Get input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L320

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L941-L947

>  const string & casadi::Function::name_in(casadi_int ind) const
------------------------------------------------------------------------

Get input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L320

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L941-L947

";

";

%feature("docstring") casadi::Function::name_out "

Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L325

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L949-L955

>  const string & casadi::Function::name_out(casadi_int ind) const
------------------------------------------------------------------------

Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L325

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L949-L955

";

";

%feature("docstring") casadi::Function::index_in "

Find the index for a string describing a particular entry of an input 

scheme.

example: schemeEntry(\"x_opt\") -> returns NLPSOL_X if 
FunctionInternal 
adheres to SCHEME_NLPINput

Extra doc: https://github.com/casadi/casadi/wiki/L_1vk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L333

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L925-L931

";

%feature("docstring") casadi::Function::index_out "

Find the index for a string describing a particular entry of an output
 
scheme.

example: schemeEntry(\"x_opt\") -> returns NLPSOL_X if 
FunctionInternal 
adheres to SCHEME_NLPINput

Extra doc: https://github.com/casadi/casadi/wiki/L_1vl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L341

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L933-L939

";

%feature("docstring") casadi::Function::default_in "

Get default input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L346

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1416-L1418

";

%feature("docstring") casadi::Function::max_in "

Get largest input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L351

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1420-L1422

";

%feature("docstring") casadi::Function::min_in "

Get smallest input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L356

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1424-L1426

";

%feature("docstring") casadi::Function::nominal_in "

Get nominal input value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L361

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1428-L1430

";

%feature("docstring") casadi::Function::nominal_out "

Get nominal output value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L366

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1432-L1434

";

%feature("docstring") casadi::Function::oracle "

Get oracle.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L407

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1816-L1822

";

%feature("docstring") casadi::Function::wrap "

Wrap in an  Function instance consisting of only one  MX call.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L412

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1824-L1826

";

%feature("docstring") casadi::Function::wrap_as_needed "

Wrap in a  Function with options.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L417

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1828-L1830

";

%feature("docstring") casadi::Function::which_depends "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1745-L1751

";

%feature("docstring") casadi::Function::print_dimensions "

Print dimensions of inputs and outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L432

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1093-L1095

";

%feature("docstring") casadi::Function::print_options "

Print options to a stream.

Extra doc: https://github.com/casadi/casadi/wiki/L_1vz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L437

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1097-L1099

";

%feature("docstring") casadi::Function::print_option "

Print all information there is to know about a certain option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L442

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1101-L1103

";

%feature("docstring") casadi::Function::has_option "

Does a particular option exist.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L447

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1105-L1112

";

%feature("docstring") casadi::Function::change_option "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1114-L1124

";

%feature("docstring") casadi::Function::uses_output "

Do the derivative functions need nondifferentiated outputs?

Extra doc: https://github.com/casadi/casadi/wiki/L_1w3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L460

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L849-L851

";

%feature("docstring") casadi::Function::jacobian_old "

[DEPRECATED] Replaced by  Function::factory.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L466

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L854-L860

";

%feature("docstring") casadi::Function::hessian_old "

[DEPRECATED] Replaced by  Function::factory.

Extra doc: https://github.com/casadi/casadi/wiki/L_1w5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L471

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L862-L870

";

%feature("docstring") casadi::Function::jacobian "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L882-L888

";

%feature("docstring") casadi::Function::rev "

[INTERNAL] 
Propagate sparsity backward with temporary memory allocation.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L633

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L458-L473

>  int casadi::Function::rev(std::vector< bvec_t * > arg, std::vector< bvec_t * > res) const
------------------------------------------------------------------------
[INTERNAL] 
Propagate sparsity backward with temporary memory allocation.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L633

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L458-L473

";

";

%feature("docstring") casadi::Function::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization: 
Type of parallelization used: unroll|serial|openmp

Extra doc: https://github.com/casadi/casadi/wiki/L_1wh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L642

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L722-L729

";

%feature("docstring") casadi::Function::slice "

returns a new function with a selection of inputs/outputs of the 
original

Extra doc: https://github.com/casadi/casadi/wiki/L_1wl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L754

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L713-L720

";

%feature("docstring") casadi::Function::forward "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1077-L1083

";

%feature("docstring") casadi::Function::reverse "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1085-L1091

";

%feature("docstring") casadi::Function::jac_sparsity "

Get, if necessary generate, the sparsity of a single Jacobian block.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L830

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L909-L915

>  Sparsity casadi::Function::jac_sparsity(casadi_int oind, casadi_int iind, bool compact=false) const
------------------------------------------------------------------------

Get, if necessary generate, the sparsity of a single Jacobian block.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L830

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L909-L915

";

";

%feature("docstring") casadi::Function::generate "

Export / Generate C code for the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L840

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1130-L1132

>  std::string casadi::Function::generate(const Dict &opts=Dict()) const
------------------------------------------------------------------------

Export / Generate C code for the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1wv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L840

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1130-L1132

";

";

%feature("docstring") casadi::Function::generate_dependencies "

Export / Generate C code for the dependency function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ww

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L845

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1140-L1142

";

%feature("docstring") casadi::Function::stats "

Get all statistics obtained at the end of the last evaluate call.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L932

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L894-L896

";

%feature("docstring") casadi::Function::has_free "

Does the function have free variables.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1004

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1627-L1629

";

%feature("docstring") casadi::Function::get_free "

Get free variables as a string.

Extra doc: https://github.com/casadi/casadi/wiki/L_1x9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1009

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1126-L1128

";

%feature("docstring") casadi::Function::free_sx "

Get all the free variables of the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1014

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1603-L1609

";

%feature("docstring") casadi::Function::free_mx "

Get all the free variables of the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1019

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1611-L1617

";

%feature("docstring") casadi::Function::generate_lifted "

Extract the functions needed for the Lifted  Newton method.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1024

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1631-L1637

";

%feature("docstring") casadi::Function::n_nodes "

Number of nodes in the algorithm.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1030

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1695-L1701

";

%feature("docstring") casadi::Function::n_instructions "

Number of instruction in the algorithm (SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xe

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1035

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1639-L1645

";

%feature("docstring") casadi::Function::instruction_id "

Identifier index of the instruction (SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1040

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1663-L1669

";

%feature("docstring") casadi::Function::instruction_input "

Locations in the work vector for the inputs of the instruction.

(SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1047

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1671-L1677

";

%feature("docstring") casadi::Function::instruction_constant "

Get the floating point output argument of an instruction (SXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1052

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1679-L1685

";

%feature("docstring") casadi::Function::instruction_output "

Location in the work vector for the output of the instruction.

(SXFunction/MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1059

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1687-L1693

";

%feature("docstring") casadi::Function::instruction_MX "

Get the  MX node corresponding to an instruction (MXFunction)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1064

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1647-L1653

";

%feature("docstring") casadi::Function::instructions_sx "

Get the SX node corresponding to all instructions (SXFunction)

Note: input and output instructions have no SX representation. This 
method 
returns nan for those instructions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1072

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1655-L1661

";

%feature("docstring") casadi::Function::sz_arg "

Get required length of arg field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1085

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1025-L1025

";

%feature("docstring") casadi::Function::sz_res "

Get required length of res field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1090

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1027-L1027

";

%feature("docstring") casadi::Function::sz_iw "

Get required length of iw field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1095

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1029-L1029

";

%feature("docstring") casadi::Function::sz_w "

Get required length of w field.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1031-L1031

";

%feature("docstring") casadi::Function::sz_work "

[INTERNAL] 
Get number of temporary variables needed.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1106

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1021-L1023

";

%feature("docstring") casadi::Function::set_work "

[INTERNAL] 
Set the (persistent) work vectors.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1111

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1050-L1057

";

%feature("docstring") casadi::Function::set_temp "

[INTERNAL] 
Set the (temporary) work vectors.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1117

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1059-L1066

";

%feature("docstring") casadi::Function::setup "

[INTERNAL] 
Set the (persistent and temporary) work vectors.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1123

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1068-L1075

";

%feature("docstring") casadi::Function::name "

Name of the function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1xv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1137

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1244-L1250

";

%feature("docstring") casadi::Function::is_a "

Check if the function is of a particular type.

Optionally check if name matches one of the base classes (default 
true)

Extra doc: https://github.com/casadi/casadi/wiki/L_1xw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1144

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1599-L1601

";

%feature("docstring") casadi::Function::assert_size_in "

Assert that an input dimension is equal so some given value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1715-L1721

";

%feature("docstring") casadi::Function::assert_size_out "

Assert that an output dimension is equal so some given value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1192

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1723-L1728

";

%feature("docstring") casadi::Function::checkout "

Checkout a memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1703-L1705

";

%feature("docstring") casadi::Function::release "

Release a memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1198

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1707-L1709

";

%feature("docstring") casadi::Function::memory "

[INTERNAL] 
Get memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1202

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1711-L1713

";

%feature("docstring") casadi::Function::get_function "

Get a dependency function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1213

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1761-L1767

>  Function casadi::Function::get_function(const std::string &name) const
------------------------------------------------------------------------

Get a dependency function.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1213

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1761-L1767

";

";

%feature("docstring") casadi::Function::has_function "

Check if a particular dependency exists.

Extra doc: https://github.com/casadi/casadi/wiki/L_1y5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1218

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1769-L1776

";

%feature("docstring") casadi::Function::find "

Get a specific function embedded in the expression graphs.

Parameters:
-----------

max_depth: 
Maximum depth - a negative number indicates no maximum

name: 
Name of function needed

Extra doc: https://github.com/casadi/casadi/wiki/L_1y7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1796-L1813

>  Function casadi::Function::find(casadi_int max_depth, const std::string &name) const
------------------------------------------------------------------------

Get a specific function embedded in the expression graphs.

Parameters:
-----------

max_depth: 
Maximum depth - a negative number indicates no maximum

name: 
Name of function needed

Extra doc: https://github.com/casadi/casadi/wiki/L_1y7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1796-L1813

";

";

%feature("docstring") casadi::Function::info "

Obtain information about function

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1236

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1842-L1844

";

%feature("docstring") casadi::Function::class_name "

Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L100-L102

";

%feature("docstring") casadi::Function::disp "

Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L104-L110

";

%feature("docstring") casadi::Function::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::Function::is_null "

Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L74-L76

";

%feature("docstring") casadi::Function::__hash__ "

Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L130-L132

";


// File: classcasadi_1_1FunctionBuffer.xml
%feature("docstring") casadi::FunctionBuffer "

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

Set input buffer for input i.

mem.set_arg(0, memoryview(a))

Note that CasADi uses 'fortran' order: column-by-column

Extra doc: https://github.com/casadi/casadi/wiki/L_1yb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1314

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1882-L1887

";

%feature("docstring") casadi::FunctionBuffer::set_res "

Set output buffer for ouput i.

mem.set_res(0, memoryview(a))

Note that CasADi uses 'fortran' order: column-by-column

Extra doc: https://github.com/casadi/casadi/wiki/L_1yc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1323

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1888-L1893

";

%feature("docstring") casadi::FunctionBuffer::ret "

Get last return value.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.hpp#L1325

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/function.cpp#L1901-L1903

";


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
Remainder after division: (x,y) -> mod(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L598

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L598-L600

>  ExType casadi::GenericExpression::mod(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Remainder after division: (x,y) -> mod(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L598

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L598-L600

";

";

%feature("docstring") casadi::GenericExpressionCommon::fmod "

[INTERNAL] 
Remainder after division: (x,y) -> mod(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L601

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L601-L603

";

%feature("docstring") casadi::GenericExpressionCommon::atan2 "

[INTERNAL] 
Two argument arc tangent: (x,y) -> atan2(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L613

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L613-L615

>  ExType casadi::GenericExpression::atan2(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Two argument arc tangent: (x,y) -> atan2(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L613

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L613-L615

";

";

%feature("docstring") casadi::GenericExpressionCommon::if_else_zero "

[INTERNAL] 
Conditional assignment: (x,y) -> x ? y : 0.

Extra doc: https://github.com/casadi/casadi/wiki/L_ps

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L625

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L625-L627

>  ExType casadi::GenericExpression::if_else_zero(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Conditional assignment: (x,y) -> x ? y : 0.

Extra doc: https://github.com/casadi/casadi/wiki/L_ps

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L625

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L625-L627

";

";

%feature("docstring") casadi::GenericExpressionCommon::fmin "

[INTERNAL] 
Smallest of two values: (x,y) -> min(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L637

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L637-L639

>  ExType casadi::GenericExpression::fmin(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Smallest of two values: (x,y) -> min(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L637

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L637-L639

";

";

%feature("docstring") casadi::GenericExpressionCommon::fmax "

[INTERNAL] 
Largest of two values: (x,y) -> max(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L649

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L649-L651

>  ExType casadi::GenericExpression::fmax(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Largest of two values: (x,y) -> max(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L649

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L649-L651

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L665

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L665-L667

";

%feature("docstring") casadi::GenericExpressionCommon::copysign "

[INTERNAL] 
Copy sign

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L675

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L675-L677

>  ExType casadi::GenericExpression::copysign(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Copy sign

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L675

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L675-L677

";

";

%feature("docstring") casadi::GenericExpressionCommon::constpow "

[INTERNAL] 
Elementwise power with const power

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L685

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L685-L687

>  ExType casadi::GenericExpression::constpow(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Elementwise power with const power

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L685

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L685-L687

";

";

%feature("docstring") casadi::GenericExpressionCommon::printme "

[INTERNAL] 
Debug printing

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L695

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L695-L697

>  ExType casadi::GenericExpression::printme(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Debug printing

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L695

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L695-L697

";

";

%feature("docstring") casadi::GenericExpressionCommon::hypot "

[INTERNAL] 
Precision variant for 2 norm: (x,y) -> sqrt(x^2+y^2)

Extra doc: https://github.com/casadi/casadi/wiki/L_pw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L707

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L707-L709

>  ExType casadi::GenericExpression::hypot(const ExType &x, const ExType &y)
------------------------------------------------------------------------
[INTERNAL] 
Precision variant for 2 norm: (x,y) -> sqrt(x^2+y^2)

Extra doc: https://github.com/casadi/casadi/wiki/L_pw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L707

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L707-L709

";

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

Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L194

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L194-L194

";

%feature("docstring") casadi::GenericMatrixCommon::get_colind "

Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L195-L195

";

%feature("docstring") casadi::GenericMatrixCommon::row "

Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200-L200

>  casadi_int casadi::GenericMatrix< MatType >::row(casadi_int el) const
------------------------------------------------------------------------

Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200-L200

";

";

%feature("docstring") casadi::GenericMatrixCommon::colind "

Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201-L201

>  casadi_int casadi::GenericMatrix< MatType >::colind(casadi_int col) const
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L424

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L424-L426

>  MatType casadi::GenericMatrix::sumsqr(const MatType &x)
------------------------------------------------------------------------
[INTERNAL] 
Calculate sum of squares: sum_ij X_ij^2.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L424

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L424-L426

";

";

%feature("docstring") casadi::GenericMatrixCommon::linspace "

[INTERNAL] 
Matlab's  linspace command.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L455

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L455-L457

>  MatType casadi::GenericMatrix::linspace(const MatType &a, const MatType &b, casadi_int nsteps)
------------------------------------------------------------------------
[INTERNAL] 
Matlab's  linspace command.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L455

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L455-L457

";

";

%feature("docstring") casadi::GenericMatrixCommon::cross "

[INTERNAL] 
Matlab's  cross command.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L462

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L462-L464

>  MatType casadi::GenericMatrix::cross(const MatType &a, const MatType &b, casadi_int dim=-1)
------------------------------------------------------------------------
[INTERNAL] 
Matlab's  cross command.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L462

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L462-L464

";

";

%feature("docstring") casadi::GenericMatrixCommon::skew "

[INTERNAL] 
Generate a skew symmetric matrix from a 3-vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L469

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L469-L471

>  MatType casadi::GenericMatrix::skew(const MatType &a)
------------------------------------------------------------------------
[INTERNAL] 
Generate a skew symmetric matrix from a 3-vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L469

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L469-L471

";

";

%feature("docstring") casadi::GenericMatrixCommon::inv_skew "

[INTERNAL] 
Generate the 3-vector progenitor of a skew symmetric matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L476

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L476-L478

>  MatType casadi::GenericMatrix::inv_skew(const MatType &a)
------------------------------------------------------------------------
[INTERNAL] 
Generate the 3-vector progenitor of a skew symmetric matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L476

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L476-L478

";

";

%feature("docstring") casadi::GenericMatrixCommon::tril2symm "

[INTERNAL] 
Convert a lower triangular matrix to a symmetric one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L514

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L514-L514

>  MatType casadi::GenericMatrix::tril2symm(const MatType &a)
------------------------------------------------------------------------
[INTERNAL] 
Convert a lower triangular matrix to a symmetric one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L514

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L514-L514

";

";

%feature("docstring") casadi::GenericMatrixCommon::triu2symm "

[INTERNAL] 
Convert a upper triangular matrix to a symmetric one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L519

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L519-L519

>  MatType casadi::GenericMatrix::triu2symm(const MatType &a)
------------------------------------------------------------------------
[INTERNAL] 
Convert a upper triangular matrix to a symmetric one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L519

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L519-L519

";

";

%feature("docstring") casadi::GenericMatrixCommon::repsum "

[INTERNAL] 
Given a repeated matrix, computes the sum of repeated parts.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L966

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L966-L968

>  MatType casadi::GenericMatrix::repsum(const MatType &A, casadi_int n, casadi_int m=1)
------------------------------------------------------------------------
[INTERNAL] 
Given a repeated matrix, computes the sum of repeated parts.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L966

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L966-L968

";

";

%feature("docstring") casadi::GenericMatrixCommon::diff "

[INTERNAL] 
Returns difference (n-th order) along given axis (MATLAB 
convention)

Extra doc: https://github.com/casadi/casadi/wiki/L_1c8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L544

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L544-L546

>  MatType casadi::GenericMatrix::diff(const MatType &x, casadi_int n=1, casadi_int axis=-1)
------------------------------------------------------------------------
[INTERNAL] 
Returns difference (n-th order) along given axis (MATLAB convention)

Extra doc: https://github.com/casadi/casadi/wiki/L_1c8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L544

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L544-L546

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L881

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L881-L883

>  bool casadi::GenericMatrix::is_linear(const MatType &expr, const MatType &var)
------------------------------------------------------------------------
[INTERNAL] 
Is expr linear in var?

False negatives are possible (an expression may not be recognised as 
linear
 while it really is), false positives not.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L881

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L881-L883

";

";

%feature("docstring") casadi::GenericMatrixCommon::is_quadratic "

[INTERNAL] 
Is expr quadratic in var?

False negatives are possible (an expression may not be recognised as 

quadratic while it really is), false positives not.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L892

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L892-L894

>  bool casadi::GenericMatrix::is_quadratic(const MatType &expr, const MatType &var)
------------------------------------------------------------------------
[INTERNAL] 
Is expr quadratic in var?

False negatives are possible (an expression may not be recognised as 

quadratic while it really is), false positives not.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L892

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L892-L894

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L906

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L906-L909

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L906

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L906-L909

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L919

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L919-L922

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L919

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L919-L922

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

%feature("docstring") casadi::GenericMatrixCommon::bilin "

Calculate bilinear form x^T A y.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L401

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L401-L403

>  MatType casadi::GenericMatrix< MatType >::bilin(const MatType &A, const MatType &x, const MatType &y)
------------------------------------------------------------------------

Calculate bilinear form x^T A y.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L401

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L401-L403

";

";

%feature("docstring") casadi::GenericMatrixCommon::rank1 "

Make a rank-1 update to a matrix A.

Calculates A + 1/2 * alpha * x*y'

Extra doc: https://github.com/casadi/casadi/wiki/L_1bp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L413

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L413-L416

>  MatType casadi::GenericMatrix< MatType >::rank1(const MatType &A, const MatType &alpha, const MatType &x, const MatType &y)
------------------------------------------------------------------------

Make a rank-1 update to a matrix A.

Calculates A + 1/2 * alpha * x*y'

Extra doc: https://github.com/casadi/casadi/wiki/L_1bp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L413

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L413-L416

";

";

%feature("docstring") casadi::GenericMatrixCommon::hessian "

[INTERNAL] 
Hessian and (optionally) gradient.

Extra doc: https://github.com/casadi/casadi/wiki/L_23z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L860

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L860-L863

>  MatType casadi::GenericMatrix::hessian(const MatType &ex, const MatType &arg, MatType &output_g, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Hessian and (optionally) gradient.

Extra doc: https://github.com/casadi/casadi/wiki/L_23z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L860

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L860-L863

";

";

%feature("docstring") casadi::GenericMatrixCommon::mmin "

[INTERNAL] 
Smallest element in a matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L974

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L974-L976

";

%feature("docstring") casadi::GenericMatrixCommon::mmax "

[INTERNAL] 
Largest element in a matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L983

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L983-L985

";

%feature("docstring") casadi::GenericMatrixCommon::jtimes "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L827

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L827-L830

>  MatType casadi::GenericMatrix< MatType >::jtimes(const MatType &ex, const MatType &arg, const MatType &v, bool tr=false, const Dict &opts=Dict())
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L827

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L827-L830

";

";

%feature("docstring") casadi::GenericMatrixCommon::gradient "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L807

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L807-L809

>  MatType casadi::GenericMatrix< MatType >::gradient(const MatType &ex, const MatType &arg, const Dict &opts=Dict())
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L807

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L807-L809

";

";

%feature("docstring") casadi::GenericMatrixCommon::tangent "

Calculate the tangent of an expression.

Extra doc: https://github.com/casadi/casadi/wiki/L_23y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L814

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L814-L816

>  MatType casadi::GenericMatrix< MatType >::tangent(const MatType &ex, const MatType &arg, const Dict &opts=Dict())
------------------------------------------------------------------------

Calculate the tangent of an expression.

Extra doc: https://github.com/casadi/casadi/wiki/L_23y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L814

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L814-L816

";

";

%feature("docstring") casadi::GenericMatrixCommon::linearize "

Linearize an expression.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L735

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L735-L738

>  MatType casadi::GenericMatrix< MatType >::linearize(const MatType &f, const MatType &x, const MatType &x0, const Dict &opts=Dict())
------------------------------------------------------------------------

Linearize an expression.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L735

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L735-L738

";

";

%feature("docstring") casadi::GenericMatrixCommon::mpower "

Matrix power x^n.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L319

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L319-L321

>  MatType casadi::GenericMatrix< MatType >::mpower(const MatType &x, const MatType &n)
------------------------------------------------------------------------

Matrix power x^n.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bi

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L319

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L319-L321

";

";

%feature("docstring") casadi::GenericMatrixCommon::soc "

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

Create a vector of length r of vectors of length p.

with nrow-by-ncol symbolic primitives

Extra doc: https://github.com/casadi/casadi/wiki/L_1dg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1060

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1060-L1062

>  static std::vector<std::vector<MatType> > casadi::GenericMatrix< MatType >::sym(const std::string &name, casadi_int nrow, casadi_int ncol, casadi_int p, casadi_int r)
------------------------------------------------------------------------

Create a vector of length r of vectors of length p.

with nrow-by-ncol symbolic primitives

Extra doc: https://github.com/casadi/casadi/wiki/L_1dg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1060

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1060-L1062

";

";

%feature("docstring") casadi::GenericMatrixCommon::zeros "

Create a dense matrix or a matrix with specified sparsity with all 
entries 
zero.

Extra doc: https://github.com/casadi/casadi/wiki/L_1dh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1073

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1073-L1075

>  static MatType casadi::GenericMatrix< MatType >::zeros(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------

Create a dense matrix or a matrix with specified sparsity with all 
entries 
zero.

Extra doc: https://github.com/casadi/casadi/wiki/L_1dh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1073

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1073-L1075

";

";

%feature("docstring") casadi::GenericMatrixCommon::ones "

Create a dense matrix or a matrix with specified sparsity with all 
entries 
one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1di

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1086

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1086-L1088

>  static MatType casadi::GenericMatrix< MatType >::ones(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------

Create a dense matrix or a matrix with specified sparsity with all 
entries 
one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1di

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1086

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1086-L1088

";

";

%feature("docstring") casadi::GenericMatrixCommon::nnz "

Get the number of (structural) non-zero elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1an

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L84

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1105-L1107

";

%feature("docstring") casadi::GenericMatrixCommon::nnz_lower "

Get the number of non-zeros in the lower triangular half.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ao

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L89

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1110-L1112

";

%feature("docstring") casadi::GenericMatrixCommon::nnz_upper "

Get the number of non-zeros in the upper triangular half.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ap

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L94

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1115-L1117

";

%feature("docstring") casadi::GenericMatrixCommon::nnz_diag "

Get get the number of non-zeros on the diagonal.

Extra doc: https://github.com/casadi/casadi/wiki/L_1aq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L99

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1120-L1122

";

%feature("docstring") casadi::GenericMatrixCommon::numel "

Get the number of elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ar

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L104

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1125-L1127

";

%feature("docstring") casadi::GenericMatrixCommon::size1 "

Get the first dimension (i.e. number of rows)

Extra doc: https://github.com/casadi/casadi/wiki/L_1as

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L109

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1130-L1132

";

%feature("docstring") casadi::GenericMatrixCommon::rows "

Get the number of rows, Octave-style syntax.

Extra doc: https://github.com/casadi/casadi/wiki/L_1at

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L114

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L114-L114

";

%feature("docstring") casadi::GenericMatrixCommon::size2 "

Get the second dimension (i.e. number of columns)

Extra doc: https://github.com/casadi/casadi/wiki/L_1au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L119

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1135-L1137

";

%feature("docstring") casadi::GenericMatrixCommon::columns "

Get the number of columns, Octave-style syntax.

Extra doc: https://github.com/casadi/casadi/wiki/L_1av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L124

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L124-L124

";

%feature("docstring") casadi::GenericMatrixCommon::dim "

Get string representation of dimensions.

The representation is e.g. \"4x5\" or \"4x5,10nz\"

Extra doc: https://github.com/casadi/casadi/wiki/L_1aw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L131

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1150-L1152

";

%feature("docstring") casadi::GenericMatrixCommon::size "

Get the size along a particular dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ay

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L141

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1145-L1147

>  casadi_int casadi::GenericMatrix< MatType >::size(casadi_int axis) const
------------------------------------------------------------------------

Get the size along a particular dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ay

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L141

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1145-L1147

";

";

%feature("docstring") casadi::GenericMatrixCommon::is_empty "

Check if the sparsity is empty, i.e. if one of the dimensions is zero.

(or optionally both dimensions)

Extra doc: https://github.com/casadi/casadi/wiki/L_1az

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L148

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L148-L148

";

%feature("docstring") casadi::GenericMatrixCommon::is_dense "

Check if the matrix expression is dense.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L153-L153

";

%feature("docstring") casadi::GenericMatrixCommon::is_scalar "

Check if the matrix expression is scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L158

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1155-L1157

";

%feature("docstring") casadi::GenericMatrixCommon::is_square "

Check if the matrix expression is square.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L163

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L163-L163

";

%feature("docstring") casadi::GenericMatrixCommon::is_vector "

Check if the matrix is a row or column vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L168

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L168-L168

";

%feature("docstring") casadi::GenericMatrixCommon::is_row "

Check if the matrix is a row vector (i.e.  size1()==1)

Extra doc: https://github.com/casadi/casadi/wiki/L_1b4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L173

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L173-L173

";

%feature("docstring") casadi::GenericMatrixCommon::is_column "

Check if the matrix is a column vector (i.e.  size2()==1)

Extra doc: https://github.com/casadi/casadi/wiki/L_1b5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L178

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L178-L178

";

%feature("docstring") casadi::GenericMatrixCommon::is_triu "

Check if the matrix is upper triangular.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L183

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L183-L183

";

%feature("docstring") casadi::GenericMatrixCommon::is_tril "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100-L1102

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L446

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L446-L449

>  MatType casadi::GenericMatrix::logsumexp(const MatType &x, const MatType &margin)
------------------------------------------------------------------------
[INTERNAL] 
Scaled version of logsumexp.

Scaled such that max(x) <= logsumexp(x, margin) <= max(x)+margin

Extra doc: https://github.com/casadi/casadi/wiki/L_1bs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L446

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L446-L449

";

";

%feature("docstring") casadi::GenericMatrixCommon::det "

[INTERNAL] 
 Matrix determinant (experimental)

Extra doc: https://github.com/casadi/casadi/wiki/L_1bx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L483

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L483-L483

";

%feature("docstring") casadi::GenericMatrixCommon::inv_minor "

[INTERNAL] 
 Matrix inverse (experimental)

Extra doc: https://github.com/casadi/casadi/wiki/L_1by

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L488

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L488-L488

";

%feature("docstring") casadi::GenericMatrixCommon::inv "

[INTERNAL] 
 Matrix inverse.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L500

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L500-L504

>  MatType casadi::GenericMatrix::inv(const MatType &A, const std::string &lsolver, const Dict &options=Dict())
------------------------------------------------------------------------
[INTERNAL] 
 Matrix inverse.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L500

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L500-L504

";

";

%feature("docstring") casadi::GenericMatrixCommon::trace "

[INTERNAL] 
 Matrix trace.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L509

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L509-L509

";

%feature("docstring") casadi::GenericMatrixCommon::norm_fro "

[INTERNAL] 
Frobenius norm.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L524

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L524-L524

";

%feature("docstring") casadi::GenericMatrixCommon::norm_2 "

[INTERNAL] 
2-norm

Extra doc: https://github.com/casadi/casadi/wiki/L_1c5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L529

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L529-L529

";

%feature("docstring") casadi::GenericMatrixCommon::norm_1 "

[INTERNAL] 
1-norm

Extra doc: https://github.com/casadi/casadi/wiki/L_1c6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L534

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L534-L534

";

%feature("docstring") casadi::GenericMatrixCommon::norm_inf "

[INTERNAL] 
Infinity-norm.

Extra doc: https://github.com/casadi/casadi/wiki/L_1c7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L539

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L539-L539

";

%feature("docstring") casadi::GenericMatrixCommon::cumsum "

[INTERNAL] 
Returns cumulative sum along given axis (MATLAB convention)

Extra doc: https://github.com/casadi/casadi/wiki/L_1c9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L551

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L551-L553

";

%feature("docstring") casadi::GenericMatrixCommon::dot "

[INTERNAL] 
Inner product of two matrices.

with x and y matrices of the same dimension

Extra doc: https://github.com/casadi/casadi/wiki/L_1ca

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L560

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L560-L562

";

%feature("docstring") casadi::GenericMatrixCommon::nullspace "

[INTERNAL] 
Computes the nullspace of a matrix A.

Finds Z m-by-(m-n) such that AZ = 0 with A n-by-m with m > n

Assumes A is full rank

Inspired by Numerical Methods in Scientific Computing by Ake Bjorck

Extra doc: https://github.com/casadi/casadi/wiki/L_1cb

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L574

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L574-L576

";

%feature("docstring") casadi::GenericMatrixCommon::polyval "

[INTERNAL] 
Evaluate a polynomial with coefficients p in x.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L581

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L581-L583

";

%feature("docstring") casadi::GenericMatrixCommon::diag "

[INTERNAL] 
Get the diagonal of a matrix or construct a diagonal.

When the input is square, the diagonal elements are returned. If the 
input 
is vector-like, a diagonal matrix is constructed with it.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L591

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L591-L593

";

%feature("docstring") casadi::GenericMatrixCommon::unite "

[INTERNAL] 
Unite two matrices no overlapping sparsity.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ce

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L598

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L598-L600

";

%feature("docstring") casadi::GenericMatrixCommon::densify "

[INTERNAL] 
Make the matrix dense and assign nonzeros to a value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L612

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L612-L614

>  MatType casadi::GenericMatrix::densify(const MatType &x, const MatType &val)
------------------------------------------------------------------------
[INTERNAL] 
Make the matrix dense and assign nonzeros to a value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L612

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L612-L614

";

";

%feature("docstring") casadi::GenericMatrixCommon::project "

[INTERNAL] 
Create a new matrix with a given sparsity pattern but with the.

nonzeros taken from an existing matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_1ch

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L621

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L621-L624

";

%feature("docstring") casadi::GenericMatrixCommon::if_else "

[INTERNAL] 
Branching on  MX nodes.

Ternary operator, \"cond ? if_true : if_false\"

Extra doc: https://github.com/casadi/casadi/wiki/L_1ci

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L631

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L631-L634

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L642

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L642-L645

";

%feature("docstring") casadi::GenericMatrixCommon::depends_on "

[INTERNAL] 
Check if expression depends on the argument.

The argument must be symbolic

Extra doc: https://github.com/casadi/casadi/wiki/L_1ck

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L652

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L652-L654

";

%feature("docstring") casadi::GenericMatrixCommon::substitute "

[INTERNAL] 
Substitute variable var with expression expr in multiple 
expressions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L668

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L668-L671

>  std::vector<MatType> casadi::GenericMatrix::substitute(const std::vector< MatType > &ex, const std::vector< MatType > &v, const std::vector< MatType > &vdef)
------------------------------------------------------------------------
[INTERNAL] 
Substitute variable var with expression expr in multiple expressions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L668

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L668-L671

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L680

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L680-L684

";

%feature("docstring") casadi::GenericMatrixCommon::cse "

[INTERNAL] 
Common subexpression elimination.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L697

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L697-L699

>  std::vector<MatType> casadi::GenericMatrix::cse(const std::vector< MatType > &e)
------------------------------------------------------------------------
[INTERNAL] 
Common subexpression elimination.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L697

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L697-L699

";

";

%feature("docstring") casadi::GenericMatrixCommon::solve "

[INTERNAL] 
 Solve a system of equations: A*x = b.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L726

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L726-L730

>  MatType casadi::GenericMatrix::solve(const MatType &A, const MatType &b, const std::string &lsolver, const Dict &dict=Dict())
------------------------------------------------------------------------
[INTERNAL] 
 Solve a system of equations: A*x = b.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L726

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L726-L730

";

";

%feature("docstring") casadi::GenericMatrixCommon::pinv "

[INTERNAL] 
Computes the Moore-Penrose pseudo-inverse.

If the matrix A is fat (size1>size2), mul(A, pinv(A)) is unity. If the
 
matrix A is slender (size2<size1), mul(pinv(A), A) is unity.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L762

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L762-L765

>  MatType casadi::GenericMatrix::pinv(const MatType &A, const std::string &lsolver, const Dict &dict=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Computes the Moore-Penrose pseudo-inverse.

If the matrix A is fat (size1>size2), mul(A, pinv(A)) is unity. If the
 
matrix A is slender (size2<size1), mul(pinv(A), A) is unity.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L762

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L762-L765

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L777

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L777-L779

";

%feature("docstring") casadi::GenericMatrixCommon::expm "

[INTERNAL] 
Calculate  Matrix exponential.

Extra doc: https://github.com/casadi/casadi/wiki/L_23w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L785

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L785-L787

";

%feature("docstring") casadi::GenericMatrixCommon::jacobian "

[INTERNAL] 
Calculate Jacobian.

Sparse matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_1cv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L794

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L794-L797

";

%feature("docstring") casadi::GenericMatrixCommon::forward "

[INTERNAL] 
Forward directional derivative.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L836

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L836-L840

";

%feature("docstring") casadi::GenericMatrixCommon::reverse "

[INTERNAL] 
Reverse directional derivative.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L846

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L846-L850

";

%feature("docstring") casadi::GenericMatrixCommon::which_depends "

[INTERNAL] 
Find out which variables enter with some order.

Extra doc: https://github.com/casadi/casadi/wiki/L_1cz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L869

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L869-L872

";

%feature("docstring") casadi::GenericMatrixCommon::n_nodes "

[INTERNAL] 
Count number of nodes

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L925

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L925-L927

";

%feature("docstring") casadi::GenericMatrixCommon::simplify "

[INTERNAL] 
Simplify an expression.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L930

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L930-L932

";

%feature("docstring") casadi::GenericMatrixCommon::print_operator "

[INTERNAL] 
Get a string representation for a binary MatType, using custom 

arguments.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L938

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L938-L940

";

%feature("docstring") casadi::GenericMatrixCommon::extract "

[INTERNAL] 
Introduce intermediate variables for selected nodes in a graph.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L945

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L945-L950

";

%feature("docstring") casadi::GenericMatrixCommon::shared "

[INTERNAL] 
Extract shared subexpressions from an set of expressions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1d6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L955

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L955-L961

";


// File: classcasadi_1_1GenericType.xml
%feature("docstring") casadi::GenericType "

Generic data type, can hold different types such as bool, casadi_int, 

string etc.

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L136-L138

";

%feature("docstring") casadi::GenericType::is_int "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L162

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L140-L142

";

%feature("docstring") casadi::GenericType::is_double "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L163

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L144-L146

";

%feature("docstring") casadi::GenericType::is_string "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L164

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L148-L150

";

%feature("docstring") casadi::GenericType::is_empty_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L165

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L152-L159

";

%feature("docstring") casadi::GenericType::is_int_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L166

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L161-L163

";

%feature("docstring") casadi::GenericType::is_int_vector_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L167

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L169-L171

";

%feature("docstring") casadi::GenericType::is_double_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L168

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L173-L175

";

%feature("docstring") casadi::GenericType::is_double_vector_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L169

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L177-L179

";

%feature("docstring") casadi::GenericType::is_bool_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L170

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L165-L167

";

%feature("docstring") casadi::GenericType::is_string_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L171

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L181-L183

";

%feature("docstring") casadi::GenericType::is_dict "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L172

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L197-L199

";

%feature("docstring") casadi::GenericType::is_function "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L173

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L185-L187

";

%feature("docstring") casadi::GenericType::is_function_vector "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L174

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L189-L191

";

%feature("docstring") casadi::GenericType::is_void_pointer "

[INTERNAL] 
Check if a particular type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L175

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L193-L195

";

%feature("docstring") casadi::GenericType::as_bool "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L182

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L264-L267

";

%feature("docstring") casadi::GenericType::as_int "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L183

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L269-L272

";

%feature("docstring") casadi::GenericType::as_double "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L184

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L274-L277

";

%feature("docstring") casadi::GenericType::as_string "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L185

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L279-L282

";

%feature("docstring") casadi::GenericType::as_int_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L186

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L284-L287

";

%feature("docstring") casadi::GenericType::as_bool_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L187

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L289-L292

";

%feature("docstring") casadi::GenericType::as_int_vector_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L188

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L294-L297

";

%feature("docstring") casadi::GenericType::as_double_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L299-L302

";

%feature("docstring") casadi::GenericType::as_double_vector_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L190

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L304-L307

";

%feature("docstring") casadi::GenericType::as_string_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L191

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L309-L312

";

%feature("docstring") casadi::GenericType::as_dict "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L192

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L314-L317

";

%feature("docstring") casadi::GenericType::as_function "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L193

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L319-L322

";

%feature("docstring") casadi::GenericType::as_function_vector "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L194

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L324-L327

";

%feature("docstring") casadi::GenericType::as_void_pointer "

[INTERNAL] 
Cast to the internal type.

Extra doc: https://github.com/casadi/casadi/wiki/L_17q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L329-L332

";

%feature("docstring") casadi::GenericType::to_bool "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L334-L343

";

%feature("docstring") casadi::GenericType::to_int "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L345-L354

";

%feature("docstring") casadi::GenericType::to_double "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L202

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L356-L363

";

%feature("docstring") casadi::GenericType::to_string "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L203

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L365-L368

";

%feature("docstring") casadi::GenericType::to_int_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L204

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L375-L378

";

%feature("docstring") casadi::GenericType::to_bool_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L205

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L386-L395

";

%feature("docstring") casadi::GenericType::to_int_vector_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L206

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L397-L400

";

%feature("docstring") casadi::GenericType::to_double_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L207

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L402-L410

";

%feature("docstring") casadi::GenericType::to_double_vector_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L208

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L412-L423

";

%feature("docstring") casadi::GenericType::to_string_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L209

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L425-L441

";

%feature("docstring") casadi::GenericType::to_dict "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L210

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L443-L446

";

%feature("docstring") casadi::GenericType::to_function "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L211

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L448-L451

";

%feature("docstring") casadi::GenericType::to_function_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L212

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L453-L456

";

%feature("docstring") casadi::GenericType::to_void_pointer "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L213

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L458-L466

";

%feature("docstring") casadi::GenericType::to_int_type_vector "

[INTERNAL] 
Convert to a type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L214

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L370-L373

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

Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_17r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L225

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L549-L552

";

%feature("docstring") casadi::GenericType::class_name "

[INTERNAL] 
Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L100-L102

";

%feature("docstring") casadi::GenericType::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L104-L110

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
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L74-L76

";

%feature("docstring") casadi::GenericType::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L130-L132

";


// File: classcasadi_1_1GlobalOptions.xml
%feature("docstring") casadi::GlobalOptions "

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

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1ImplicitToNlp.xml
%feature("docstring") casadi::ImplicitToNlp "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Importer.xml
%feature("docstring") casadi::Importer "

Importer.

Just-in-time compilation of code
General informationList of plugins
- clang

- shell

Note: some of the plugins in this list might not be available on your 

system. Also, there might be extra plugins available to you that are 
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

Importer factory.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L97

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L35-L46

>  casadi::Importer::Importer(const std::string &name, const std::string &compiler, const Dict &opts=Dict())
------------------------------------------------------------------------

Importer factory.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L97

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L35-L46

";

";

%feature("docstring") casadi::Importer::plugin_name "

Query plugin name.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L118

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L72-L74

";

%feature("docstring") casadi::Importer::get_function "

[INTERNAL] 
Get a function pointer for numerical evaluation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L125

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L80-L82

";

%feature("docstring") casadi::Importer::has_meta "

Does a meta entry exist?

Extra doc: https://github.com/casadi/casadi/wiki/L_165

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L145

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L84-L86

";

%feature("docstring") casadi::Importer::get_meta "

Get entry as a text.

Extra doc: https://github.com/casadi/casadi/wiki/L_166

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L88-L90

";

%feature("docstring") casadi::Importer::inlined "

Check if a function is inlined.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L92-L94

";

%feature("docstring") casadi::Importer::body "

Get the function body, if inlined.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L156

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L96-L98

";

%feature("docstring") casadi::Importer::library "

Get library name.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L159

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L100-L102

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

Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_16c

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.hpp#L211

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/importer.cpp#L104-L106

";

%feature("docstring") casadi::Importer::class_name "

Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L100-L102

";

%feature("docstring") casadi::Importer::disp "

Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L104-L110

";

%feature("docstring") casadi::Importer::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::Importer::is_null "

Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L74-L76

";

%feature("docstring") casadi::Importer::__hash__ "

Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L130-L132

";


// File: classcasadi_1_1IncrementalSerializer.xml



// File: classcasadi_1_1Input.xml


// File: classcasadi_1_1Integrator.xml
%feature("docstring") casadi::Integrator "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Interpolant.xml
%feature("docstring") casadi::Interpolant "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1IOInstruction.xml


// File: classcasadi_1_1Ipqp.xml
%feature("docstring") casadi::Ipqp "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1KeyboardInterruptException.xml
%feature("docstring") casadi::KeyboardInterruptException "

C++ includes: exception.hpp
";

%feature("docstring") 
casadi::KeyboardInterruptException::KeyboardInterruptException "

Default constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L85

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L85-L85

";

%feature("docstring") 
casadi::KeyboardInterruptException::~KeyboardInterruptException "

throw ()
Destructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L87

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L87-L87

";

%feature("docstring") casadi::KeyboardInterruptException::what "

throw ()
Display error.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L74

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/exception.hpp#L74-L76

";


// File: classcasadi_1_1LapackLu.xml
%feature("docstring") casadi::LapackLu "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1LapackQr.xml
%feature("docstring") casadi::LapackQr "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1LinearInterpolant.xml
%feature("docstring") casadi::LinearInterpolant "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Linsol.xml
%feature("docstring") casadi::Linsol "

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

system. Also, there might be extra plugins available to you that are 
not 
listed here. You can obtain their documentation with   
Linsol.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

csparsecholesky
---------------



Linsol with CSparseCholesky Interface

Extra doc: https://github.com/casadi/casadi/wiki/L_21u

Linsol with CSparseCholesky Interface

Extra doc: https://github.com/casadi/casadi/wiki/L_22s



--------------------------------------------------------------------------------

csparse
-------



Linsol with CSparse Interface

Extra doc: https://github.com/casadi/casadi/wiki/L_21t

Linsol with CSparse Interface

Extra doc: https://github.com/casadi/casadi/wiki/L_22r



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
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L188-L195

>  int casadi::Linsol::solve(const double *A, double *x, casadi_int nrhs=1, bool tr=false, int mem=0) const
------------------------------------------------------------------------
[INTERNAL] 
Low-level API

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L136

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L188-L195

";

";

%feature("docstring") casadi::Linsol::sfact "

[INTERNAL] 
Symbolic factorization of the linear system, e.g. selecting 
pivots.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L103

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L105-L108

>  void casadi::Linsol::sfact(const DM &A) const
------------------------------------------------------------------------
[INTERNAL] 
Symbolic factorization of the linear system, e.g. selecting pivots.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L103

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L105-L108

";

";

%feature("docstring") casadi::Linsol::nfact "

[INTERNAL] 
Numeric factorization of the linear system.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L106

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L127-L130

>  void casadi::Linsol::nfact(const DM &A) const
------------------------------------------------------------------------
[INTERNAL] 
Numeric factorization of the linear system.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L106

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L127-L130

";

";

%feature("docstring") casadi::Linsol::neig "

Number of negative eigenvalues.

Not available for all solvers

Extra doc: https://github.com/casadi/casadi/wiki/L_1kk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L119

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L166-L171

>  casadi_int casadi::Linsol::neig(const DM &A) const
------------------------------------------------------------------------

Number of negative eigenvalues.

Not available for all solvers

Extra doc: https://github.com/casadi/casadi/wiki/L_1kk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L119

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L166-L171

";

";

%feature("docstring") casadi::Linsol::rank "

Matrix rank.

Not available for all solvers

Extra doc: https://github.com/casadi/casadi/wiki/L_1kl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L126

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L177-L182

>  casadi_int casadi::Linsol::rank(const DM &A) const
------------------------------------------------------------------------

Matrix rank.

Not available for all solvers

Extra doc: https://github.com/casadi/casadi/wiki/L_1kl

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L126

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L177-L182

";

";

%feature("docstring") casadi::Linsol::Linsol "

Constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L66

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L35-L39

>  casadi::Linsol::Linsol(const std::string &name, const std::string &solver, const Sparsity &sp, const Dict &opts=Dict())
------------------------------------------------------------------------

Constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L66

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L35-L39

";

";

%feature("docstring") casadi::Linsol::plugin_name "

Query plugin name.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L97

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L65-L67

";

%feature("docstring") casadi::Linsol::sparsity "

Get linear system sparsity.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L69-L71

";

%feature("docstring") casadi::Linsol::stats "

Get all statistics obtained at the end of the last evaluate call.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L129

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L218-L220

";

%feature("docstring") casadi::Linsol::checkout "

[INTERNAL] 
Checkout a memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L142

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L197-L199

";

%feature("docstring") casadi::Linsol::release "

[INTERNAL] 
Release a memory object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L145

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L201-L203

";

%feature("docstring") casadi::Linsol::serialize "

[INTERNAL] 
Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_1km

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.cpp#L222-L225

";

%feature("docstring") casadi::Linsol::class_name "

Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L100-L102

";

%feature("docstring") casadi::Linsol::disp "

Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L104-L110

";

%feature("docstring") casadi::Linsol::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::Linsol::is_null "

Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L74-L76

";

%feature("docstring") casadi::Linsol::__hash__ "

Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L130-L132

";


// File: classcasadi_1_1LinsolCall.xml


// File: classcasadi_1_1LinsolLdl.xml
%feature("docstring") casadi::LinsolLdl "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1LinsolQr.xml
%feature("docstring") casadi::LinsolQr "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Logger.xml
%feature("docstring") casadi::Logger "

Keeps track of logging output to screen and/or files.

All printout from CasADi routines should go through this files.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_23u

C++ includes: casadi_logger.hpp
";


// File: classcasadi_1_1Matrix.xml
%feature("docstring") casadi::MatrixCommon "

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

C++ includes: casadi_common.hpp
";

%feature("docstring") casadi::MatrixCommon::get "

Get a submatrix, two arguments

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L239

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L106-L126

>  void casadi::Matrix< Scalar >::get(Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &rr, const Matrix< casadi_int > &cc) const
------------------------------------------------------------------------

Get a submatrix, two arguments

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L239

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L106-L126

";

";

%feature("docstring") casadi::MatrixCommon::set "

Set a submatrix, two arguments

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L255

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L213-L277

>  void casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &rr, const Matrix< casadi_int > &cc)
------------------------------------------------------------------------

Set a submatrix, two arguments

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L255

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L213-L277

";

";

%feature("docstring") casadi::MatrixCommon::get_nz "

Get a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L262

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L401-L427

>  void casadi::Matrix< Scalar >::get_nz(Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &k) const
------------------------------------------------------------------------

Get a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L262

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L401-L427

";

";

%feature("docstring") casadi::MatrixCommon::set_nz "

Set a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L268

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L442-L484

>  void casadi::Matrix< Scalar >::set_nz(const Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &k)
------------------------------------------------------------------------

Set a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L268

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L442-L484

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

%feature("docstring") casadi::MatrixCommon::hessian "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L306

>  static Matrix<Scalar> casadi::Matrix< Scalar >::hessian(const Matrix< Scalar > &f, const Matrix< Scalar > &x, Matrix< Scalar > &g, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L306

";

";

%feature("docstring") casadi::MatrixCommon::substitute "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L313

>  static std::vector<Matrix<Scalar> > casadi::Matrix< Scalar >::substitute(const std::vector< Matrix< Scalar > > &ex, const std::vector< Matrix< Scalar > > &v, const std::vector< Matrix< Scalar > > &vdef)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L313

";

";

%feature("docstring") casadi::MatrixCommon::substitute_inplace "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L316

";

%feature("docstring") casadi::MatrixCommon::pinv "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L321

>  static Matrix<Scalar> casadi::Matrix< Scalar >::pinv(const Matrix< Scalar > &A, const std::string &lsolver, const Dict &opts)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L321

";

";

%feature("docstring") casadi::MatrixCommon::expm_const "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L323

";

%feature("docstring") casadi::MatrixCommon::expm "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L324

";

%feature("docstring") casadi::MatrixCommon::solve "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L326

>  static Matrix<Scalar> casadi::Matrix< Scalar >::solve(const Matrix< Scalar > &A, const Matrix< Scalar > &b, const std::string &lsolver, const Dict &opts)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L326

";

";

%feature("docstring") casadi::MatrixCommon::inv "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L329

>  static Matrix<Scalar> casadi::Matrix< Scalar >::inv(const Matrix< Scalar > &A, const std::string &lsolver, const Dict &opts)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L329

";

";

%feature("docstring") casadi::MatrixCommon::n_nodes "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L332

";

%feature("docstring") casadi::MatrixCommon::print_operator "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L333

";

%feature("docstring") casadi::MatrixCommon::extract "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L335

";

%feature("docstring") casadi::MatrixCommon::shared "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L337

";

%feature("docstring") casadi::MatrixCommon::_bilin "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L342

";

%feature("docstring") casadi::MatrixCommon::_rank1 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L345

";

%feature("docstring") casadi::MatrixCommon::if_else "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L349

";

%feature("docstring") casadi::MatrixCommon::conditional "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L353

";

%feature("docstring") casadi::MatrixCommon::depends_on "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L357

";

%feature("docstring") casadi::MatrixCommon::mrdivide "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L358

";

%feature("docstring") casadi::MatrixCommon::mldivide "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L359

";

%feature("docstring") casadi::MatrixCommon::symvar "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L360

";

%feature("docstring") casadi::MatrixCommon::det "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L361

";

%feature("docstring") casadi::MatrixCommon::inv_minor "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L362

";

%feature("docstring") casadi::MatrixCommon::trace "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L363

";

%feature("docstring") casadi::MatrixCommon::norm_1 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L364

";

%feature("docstring") casadi::MatrixCommon::norm_2 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L365

";

%feature("docstring") casadi::MatrixCommon::norm_fro "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L366

";

%feature("docstring") casadi::MatrixCommon::norm_inf "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L367

";

%feature("docstring") casadi::MatrixCommon::sum2 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L368

";

%feature("docstring") casadi::MatrixCommon::sum1 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L369

";

%feature("docstring") casadi::MatrixCommon::dot "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L370

";

%feature("docstring") casadi::MatrixCommon::nullspace "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L371

";

%feature("docstring") casadi::MatrixCommon::diag "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L372

";

%feature("docstring") casadi::MatrixCommon::unite "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L373

";

%feature("docstring") casadi::MatrixCommon::project "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L374

";

%feature("docstring") casadi::MatrixCommon::polyval "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L376

";

%feature("docstring") casadi::MatrixCommon::densify "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L378

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L487-L489

>  Matrix< Scalar > casadi::Matrix< Scalar >::densify(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L378

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L487-L489

";

";

%feature("docstring") casadi::MatrixCommon::einstein "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L386

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L561-L568

>  Matrix< Scalar > casadi::Matrix< Scalar >::einstein(const Matrix< Scalar > &A, const Matrix< Scalar > &B, const std::vector< casadi_int > &dim_a, const std::vector< casadi_int > &dim_b, const std::vector< casadi_int > &dim_c, const std::vector< casadi_int > &a, const std::vector< casadi_int > &b, const std::vector< casadi_int > &c)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L386

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L561-L568

";

";

%feature("docstring") casadi::MatrixCommon::cumsum "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L391

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L522-L533

";

%feature("docstring") casadi::MatrixCommon::_logsumexp "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L392

";

%feature("docstring") casadi::MatrixCommon::cse "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L393

";

%feature("docstring") casadi::MatrixCommon::blockcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L398

";

%feature("docstring") casadi::MatrixCommon::horzcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L399

";

%feature("docstring") casadi::MatrixCommon::horzsplit "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L401

";

%feature("docstring") casadi::MatrixCommon::vertcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L403

";

%feature("docstring") casadi::MatrixCommon::vertsplit "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L405

";

%feature("docstring") casadi::MatrixCommon::diagsplit "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L408

";

%feature("docstring") casadi::MatrixCommon::reshape "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L412

>  static Matrix<Scalar> casadi::Matrix< Scalar >::reshape(const Matrix< Scalar > &x, const Sparsity &sp)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L412

";

";

%feature("docstring") casadi::MatrixCommon::kron "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L413

";

%feature("docstring") casadi::MatrixCommon::mtimes "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L414

";

%feature("docstring") casadi::MatrixCommon::mac "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L415

";

%feature("docstring") casadi::MatrixCommon::sparsify "

[INTERNAL] 
Make a matrix sparse by removing numerical zeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_191

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L611

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L611-L613

>  Matrix<Scalar> casadi::Matrix::sparsify(const Matrix< Scalar > &A, double tol=0)
------------------------------------------------------------------------
[INTERNAL] 
Make a matrix sparse by removing numerical zeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_191

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L611

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L611-L613

";

";

%feature("docstring") casadi::MatrixCommon::expand "

[INTERNAL] 
Expand the expression as a weighted sum (with constant weights)

Extra doc: https://github.com/casadi/casadi/wiki/L_192

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L618

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L618-L621

>  void casadi::Matrix::expand(const Matrix< Scalar > &ex, Matrix< Scalar > &weights, Matrix< Scalar > &terms)
------------------------------------------------------------------------
[INTERNAL] 
Expand the expression as a weighted sum (with constant weights)

Extra doc: https://github.com/casadi/casadi/wiki/L_192

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L618

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L618-L621

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L633

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L633-L637

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L633

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L633-L637

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L650

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L650-L653

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L650

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L650-L653

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L666

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L666-L668

>  Matrix<Scalar> casadi::Matrix::heaviside(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
Heaviside function.

\\\\[ \\\\begin {cases} H(x) = 0 & x<0 \\\\\\\\ H(x) = 1/2 & x=0 
\\\\\\\\ 
H(x) = 1 & x>0 \\\\\\\\ \\\\end {cases} \\\\]

Extra doc: https://github.com/casadi/casadi/wiki/L_195

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L666

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L666-L668

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L683

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L683-L685

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L683

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L683-L685

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L698

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L698-L700

>  Matrix<Scalar> casadi::Matrix::triangle(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
triangle function

\\\\[ \\\\begin {cases} \\\\Lambda(x) = 0 & |x| >= 1 \\\\\\\\ 
\\\\Lambda(x)
 = 1-|x| & |x| < 1 \\\\end {cases} \\\\]

Extra doc: https://github.com/casadi/casadi/wiki/L_23o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L698

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L698-L700

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L715

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L715-L717

>  Matrix<Scalar> casadi::Matrix::ramp(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
ramp function

\\\\[ \\\\begin {cases} R(x) = 0 & x <= 1 \\\\\\\\ R(x) = x & x > 1 

\\\\\\\\ \\\\end {cases} \\\\]

Also called: slope function

Extra doc: https://github.com/casadi/casadi/wiki/L_23p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L715

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L715-L717

";

";

%feature("docstring") casadi::MatrixCommon::gauss_quadrature "

[INTERNAL] 
Integrate f from a to b using Gaussian quadrature with n points.

Extra doc: https://github.com/casadi/casadi/wiki/L_196

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L730

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L730-L734

>  Matrix<Scalar> casadi::Matrix::gauss_quadrature(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Matrix< Scalar > &a, const Matrix< Scalar > &b, casadi_int order, const Matrix< Scalar > &w)
------------------------------------------------------------------------
[INTERNAL] 
Integrate f from a to b using Gaussian quadrature with n points.

Extra doc: https://github.com/casadi/casadi/wiki/L_196

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L730

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L730-L734

";

";

%feature("docstring") casadi::MatrixCommon::forward "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L442

";

%feature("docstring") casadi::MatrixCommon::reverse "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L447

";

%feature("docstring") casadi::MatrixCommon::which_depends "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L451

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L756

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L756-L758

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L756

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L756-L758

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L799

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L799-L803

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L799

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L799-L803

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L811

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L811-L814

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L811

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L811-L814

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L822

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L822-L824

>  Matrix<Scalar> casadi::Matrix::poly_roots(const Matrix< Scalar > &p)
------------------------------------------------------------------------
[INTERNAL] 
Attempts to find the roots of a polynomial.

This will only work for polynomials up to order 3 It is assumed that 
the 
roots are real.

Extra doc: https://github.com/casadi/casadi/wiki/L_198

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L822

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L822-L824

";

";

%feature("docstring") casadi::MatrixCommon::eig_symbolic "

[INTERNAL] 
Attempts to find the eigenvalues of a symbolic matrix.

This will only work for up to 3x3 matrices

Extra doc: https://github.com/casadi/casadi/wiki/L_199

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L831

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L831-L833

>  Matrix<Scalar> casadi::Matrix::eig_symbolic(const Matrix< Scalar > &m)
------------------------------------------------------------------------
[INTERNAL] 
Attempts to find the eigenvalues of a symbolic matrix.

This will only work for up to 3x3 matrices

Extra doc: https://github.com/casadi/casadi/wiki/L_199

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L831

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L831-L833

";

";

%feature("docstring") casadi::MatrixCommon::evalf "

[INTERNAL] 
Evaluates the expression numerically.

An error is raised when the expression contains symbols

Extra doc: https://github.com/casadi/casadi/wiki/L_19a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L841

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L841-L843

>  Matrix<double> casadi::Matrix::evalf(const Matrix< Scalar > &expr)
------------------------------------------------------------------------
[INTERNAL] 
Evaluates the expression numerically.

An error is raised when the expression contains symbols

Extra doc: https://github.com/casadi/casadi/wiki/L_19a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L841

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L841-L843

";

";

%feature("docstring") casadi::MatrixCommon::qr_sparse "

[INTERNAL] 
Sparse direct QR factorization.

See T. Davis: Direct Methods for Sparse Linear Systems

Extra doc: https://github.com/casadi/casadi/wiki/L_18t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L537

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L537-L541

>  void casadi::Matrix::qr_sparse(const Matrix< Scalar > &A, Matrix< Scalar > &V, Matrix< Scalar > &R, Matrix< Scalar > &beta, std::vector< casadi_int > &prinv, std::vector< casadi_int > &pc, bool amd=true)
------------------------------------------------------------------------
[INTERNAL] 
Sparse direct QR factorization.

See T. Davis: Direct Methods for Sparse Linear Systems

Extra doc: https://github.com/casadi/casadi/wiki/L_18t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L537

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L537-L541

";

";

%feature("docstring") casadi::MatrixCommon::qr_solve "

[INTERNAL] 
 Solve using a sparse QR factorization.

Extra doc: https://github.com/casadi/casadi/wiki/L_18u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L547

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L547-L552

>  Matrix<Scalar> casadi::Matrix::qr_solve(const Matrix< Scalar > &b, const Matrix< Scalar > &v, const Matrix< Scalar > &r, const Matrix< Scalar > &beta, const std::vector< casadi_int > &prinv, const std::vector< casadi_int > &pc, bool tr=false)
------------------------------------------------------------------------
[INTERNAL] 
 Solve using a sparse QR factorization.

Extra doc: https://github.com/casadi/casadi/wiki/L_18u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L547

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L547-L552

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L528

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L528-L530

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L528

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L528-L530

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L571

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L571-L574

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L571

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L571-L574

";

";

%feature("docstring") casadi::MatrixCommon::ldl_solve "

[INTERNAL] 
 Solve using a sparse LDL^T factorization.

Extra doc: https://github.com/casadi/casadi/wiki/L_18x

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L580

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L580-L583

>  Matrix<Scalar> casadi::Matrix::ldl_solve(const Matrix< Scalar > &b, const Matrix< Scalar > &D, const Matrix< Scalar > &LT, const std::vector< casadi_int > &p)
------------------------------------------------------------------------
[INTERNAL] 
 Solve using a sparse LDL^T factorization.

Extra doc: https://github.com/casadi/casadi/wiki/L_18x

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L580

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L580-L583

";

";

%feature("docstring") casadi::MatrixCommon::all "

[INTERNAL] 
Returns true only if every element in the matrix is true.

Extra doc: https://github.com/casadi/casadi/wiki/L_18z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L595

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L595-L597

>  Matrix<Scalar> casadi::Matrix::all(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
Returns true only if every element in the matrix is true.

Extra doc: https://github.com/casadi/casadi/wiki/L_18z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L595

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L595-L597

";

";

%feature("docstring") casadi::MatrixCommon::any "

[INTERNAL] 
Returns true only if any element in the matrix is true.

Extra doc: https://github.com/casadi/casadi/wiki/L_18y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L588

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L588-L590

>  Matrix<Scalar> casadi::Matrix::any(const Matrix< Scalar > &x)
------------------------------------------------------------------------
[INTERNAL] 
Returns true only if any element in the matrix is true.

Extra doc: https://github.com/casadi/casadi/wiki/L_18y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L588

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L588-L590

";

";

%feature("docstring") casadi::MatrixCommon::adj "

[INTERNAL] 
 Matrix adjoint.

Extra doc: https://github.com/casadi/casadi/wiki/L_18p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L502

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L502-L504

>  Matrix<Scalar> casadi::Matrix::adj(const Matrix< Scalar > &A)
------------------------------------------------------------------------
[INTERNAL] 
 Matrix adjoint.

Extra doc: https://github.com/casadi/casadi/wiki/L_18p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L502

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L502-L504

";

";

%feature("docstring") casadi::MatrixCommon::minor "

[INTERNAL] 
Get the (i,j) minor matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_18q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L509

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L509-L511

>  Matrix<Scalar> casadi::Matrix::minor(const Matrix< Scalar > &x, casadi_int i, casadi_int j)
------------------------------------------------------------------------
[INTERNAL] 
Get the (i,j) minor matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_18q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L509

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L509-L511

";

";

%feature("docstring") casadi::MatrixCommon::cofactor "

[INTERNAL] 
Get the (i,j) cofactor matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_18r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L516

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L516-L518

>  Matrix<Scalar> casadi::Matrix::cofactor(const Matrix< Scalar > &x, casadi_int i, casadi_int j)
------------------------------------------------------------------------
[INTERNAL] 
Get the (i,j) cofactor matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_18r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L516

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L516-L518

";

";

%feature("docstring") casadi::MatrixCommon::chol "

[INTERNAL] 
Obtain a Cholesky factorisation of a matrix.

Performs and LDL transformation [L,D] = ldl(A) and returns 
diag(sqrt(D))*L'

Extra doc: https://github.com/casadi/casadi/wiki/L_18v

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L560

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L560-L562

>  Matrix<Scalar> casadi::Matrix::chol(const Matrix< Scalar > &A)
------------------------------------------------------------------------
[INTERNAL] 
Obtain a Cholesky factorisation of a matrix.

Performs and LDL transformation [L,D] = ldl(A) and returns 
diag(sqrt(D))*L'

Extra doc: https://github.com/casadi/casadi/wiki/L_18v

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L560

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L560-L562

";

";

%feature("docstring") casadi::MatrixCommon::norm_inf_mul "

[INTERNAL] 
Inf-norm of a Matrix-Matrix product.

Extra doc: https://github.com/casadi/casadi/wiki/L_190

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L603

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L603-L605

>  Matrix<Scalar> casadi::Matrix::norm_inf_mul(const Matrix< Scalar > &x, const Matrix< Scalar > &y)
------------------------------------------------------------------------
[INTERNAL] 
Inf-norm of a Matrix-Matrix product.

Extra doc: https://github.com/casadi/casadi/wiki/L_190

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L603

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L603-L605

";

";

%feature("docstring") casadi::MatrixCommon::diagcat "

[INTERNAL] 
Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L484

";

%feature("docstring") casadi::MatrixCommon::nonzeros "

[INTERNAL] 
Access the non-zero elements

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L958

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L958-L958

>  const std::vector<Scalar>& casadi::Matrix< Scalar >::nonzeros() const
------------------------------------------------------------------------
[INTERNAL] 
Access the non-zero elements

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L958

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L958-L958

";

";

%feature("docstring") casadi::MatrixCommon::ptr "

[INTERNAL] 
Get a pointer to the data

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L964

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L964-L964

>  const Scalar* casadi::Matrix< Scalar >::ptr() const
------------------------------------------------------------------------
[INTERNAL] 
Get a pointer to the data

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L964

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L964-L964

";

";

%feature("docstring") casadi::MatrixCommon::get_ptr "

[INTERNAL] 
Get a pointer to the data

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L966

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L966-L966

>  const friend Scalar* casadi::Matrix::get_ptr(const Matrix< Scalar > &v)
------------------------------------------------------------------------
[INTERNAL] 
Get a pointer to the data

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L966

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L966-L966

";

";

%feature("docstring") casadi::MatrixCommon::triplet "

Construct a sparse matrix from triplet form.

Default matrix size is max(col) x max(row)

Extra doc: https://github.com/casadi/casadi/wiki/L_23t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L991

>  static Matrix<Scalar> casadi::Matrix< Scalar >::triplet(const std::vector< casadi_int > &row, const std::vector< casadi_int > &col, const Matrix< Scalar > &d, const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------

Construct a sparse matrix from triplet form.

Default matrix size is max(col) x max(row)

Extra doc: https://github.com/casadi/casadi/wiki/L_23t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L991

";

";

%feature("docstring") casadi::MatrixCommon::inf "

create a matrix with all inf

Extra doc: https://github.com/casadi/casadi/wiki/L_19k

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1003

>  static Matrix<Scalar> casadi::Matrix< Scalar >::inf(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------

create a matrix with all inf

Extra doc: https://github.com/casadi/casadi/wiki/L_19k

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1003

";

";

%feature("docstring") casadi::MatrixCommon::nan "

create a matrix with all nan

Extra doc: https://github.com/casadi/casadi/wiki/L_19l

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1012

>  static Matrix<Scalar> casadi::Matrix< Scalar >::nan(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------

create a matrix with all nan

Extra doc: https://github.com/casadi/casadi/wiki/L_19l

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1012

";

";

%feature("docstring") casadi::MatrixCommon::rand "

Create a matrix with uniformly distributed random numbers.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ab

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1201-L1203

>  static Matrix<Scalar> casadi::Matrix< Scalar >::rand(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------

Create a matrix with uniformly distributed random numbers.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ab

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1201-L1203

";

";

%feature("docstring") casadi::MatrixCommon::MatrixCommon "

[INTERNAL] 
Sparse matrix with a given sparsity and non-zero elements.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1284

>  casadi::Matrix< Scalar >::Matrix(const Sparsity &sp, const std::vector< Scalar > &d, bool dummy)
------------------------------------------------------------------------
[INTERNAL] 
Sparse matrix with a given sparsity and non-zero elements.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1284

";

";

%feature("docstring") casadi::MatrixCommon::scalar "

[INTERNAL] 
Convert to scalar type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L180

";

%feature("docstring") casadi::MatrixCommon::has_nz "

Returns true if the matrix has a non-zero at location rr, cc.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L219

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L219-L219

";

%feature("docstring") casadi::MatrixCommon::__nonzero__ "

Returns the truth value of a  Matrix.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L65-L71

";

%feature("docstring") casadi::MatrixCommon::T "

Transpose the matrix.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L492

";

%feature("docstring") casadi::MatrixCommon::print_split "

Get strings corresponding to the nonzeros and the interdependencies.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L871

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L681-L685

";

%feature("docstring") casadi::MatrixCommon::disp "

Print a representation of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L875

";

%feature("docstring") casadi::MatrixCommon::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L878

";

%feature("docstring") casadi::MatrixCommon::print_scalar "

Print scalar.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L881

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L594-L618

";

%feature("docstring") casadi::MatrixCommon::print_vector "

Print vector-style.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L884

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L621-L623

";

%feature("docstring") casadi::MatrixCommon::print_dense "

Print dense matrix-stype.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L887

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L671-L673

";

%feature("docstring") casadi::MatrixCommon::print_sparse "

Print sparse matrix style.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L890

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_impl.hpp#L676-L678

";

%feature("docstring") casadi::MatrixCommon::erase "

Erase a submatrix (leaving structural zeros in its place)

Erase elements of a matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_19g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L935

>  void casadi::Matrix< Scalar >::erase(const std::vector< casadi_int > &rr, bool ind1=false)
------------------------------------------------------------------------

Erase a submatrix (leaving structural zeros in its place)

Erase elements of a matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_19g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L935

";

";

%feature("docstring") casadi::MatrixCommon::remove "

Remove columns and rows.

Remove/delete rows and/or columns of a matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_19h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L942

";

%feature("docstring") casadi::MatrixCommon::enlarge "

Enlarge matrix.

Make the matrix larger by inserting empty rows and columns, keeping 
the 
existing non-zeros

Extra doc: https://github.com/casadi/casadi/wiki/L_19i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L950

";

%feature("docstring") casadi::MatrixCommon::sparsity "

[INTERNAL] 
Const access the sparsity - reference to data member.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L970

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L970-L970

";

%feature("docstring") casadi::MatrixCommon::get_sparsity "

Get an owning reference to the sparsity pattern.

Extra doc: https://github.com/casadi/casadi/wiki/L_19j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L977

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L977-L977

";

%feature("docstring") casadi::MatrixCommon::element_hash "

Returns a number that is unique for a given symbolic scalar.

Only defined if symbolic scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_19n

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1025

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L127-L129

";

%feature("docstring") casadi::MatrixCommon::is_regular "

Checks if expression does not contain NaN or Inf.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1028

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L101-L114

";

%feature("docstring") casadi::MatrixCommon::is_smooth "

Check if smooth.

Extra doc: https://github.com/casadi/casadi/wiki/L_19o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1033

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L117-L124

";

%feature("docstring") casadi::MatrixCommon::is_leaf "

Check if SX is a leaf of the SX graph.

Only defined if symbolic scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_19p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1040

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L132-L134

";

%feature("docstring") casadi::MatrixCommon::is_commutative "

Check whether a binary SX is commutative.

Only defined if symbolic scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_19q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1047

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L137-L139

";

%feature("docstring") casadi::MatrixCommon::is_symbolic "

Check if symbolic (Dense)

Sparse matrices invariable return false

Extra doc: https://github.com/casadi/casadi/wiki/L_19r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1054

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L151-L157

";

%feature("docstring") casadi::MatrixCommon::is_valid_input "

Check if matrix can be used to define function inputs.

Sparse matrices can return true if all non-zero elements are symbolic

Extra doc: https://github.com/casadi/casadi/wiki/L_19s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1061

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L142-L148

";

%feature("docstring") casadi::MatrixCommon::is_constant "

Check if the matrix is constant (note that false negative answers are 

possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_19v

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1083

";

%feature("docstring") casadi::MatrixCommon::is_integer "

Check if the matrix is integer-valued.

(note that false negative answers are possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_19w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1090

";

%feature("docstring") casadi::MatrixCommon::is_zero "

check if the matrix is 0 (note that false negative answers are 
possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_19x

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1095

";

%feature("docstring") casadi::MatrixCommon::is_one "

check if the matrix is 1 (note that false negative answers are 
possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_19y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1100

";

%feature("docstring") casadi::MatrixCommon::is_minus_one "

check if the matrix is -1 (note that false negative answers are 
possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_19z

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1105

";

%feature("docstring") casadi::MatrixCommon::is_eye "

check if the matrix is an identity matrix (note that false negative 
answers

are possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_1a0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1112

";

%feature("docstring") casadi::MatrixCommon::op "

Get operation type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1115

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L160-L162

";

%feature("docstring") casadi::MatrixCommon::is_op "

Is it a certain operation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1118

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L165-L167

";

%feature("docstring") casadi::MatrixCommon::has_zeros "

Check if the matrix has any zero entries which are not structural 
zeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1a1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1123

";

%feature("docstring") casadi::MatrixCommon::get_nonzeros "

Get all nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1a4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1140

>  std::vector<A> casadi::Matrix< Scalar >::get_nonzeros() const
------------------------------------------------------------------------

Get all nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1a4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1140

";

";

%feature("docstring") casadi::MatrixCommon::get_elements "

Get all elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1a3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1133

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1133-L1133

";

%feature("docstring") casadi::MatrixCommon::name "

Get name (only if symbolic scalar)

Extra doc: https://github.com/casadi/casadi/wiki/L_1a8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1164

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L189-L191

";

%feature("docstring") casadi::MatrixCommon::dep "

Get expressions of the children of the expression.

Only defined if symbolic scalar. Wraps  SXElem SXElem::dep(casadi_int ch=0) 
const.

Extra doc: https://github.com/casadi/casadi/wiki/L_1a9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1172

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L194-L196

";

%feature("docstring") casadi::MatrixCommon::n_dep "

Get the number of dependencies of a binary  SXElem.

Only defined if symbolic scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_1aa

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1179

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_instantiator.cpp#L199-L201

";

%feature("docstring") casadi::MatrixCommon::export_code "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1220

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dm_instantiator.cpp#L96-L194

";

%feature("docstring") casadi::MatrixCommon::info "

Obtain information about sparsity

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1224

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dm_instantiator.cpp#L197-L199

";

%feature("docstring") casadi::MatrixCommon::serialize "

Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ah

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1250

>  void casadi::Matrix< Scalar >::serialize(SerializingStream &s) const
------------------------------------------------------------------------

Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ah

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1250

";

";

%feature("docstring") casadi::MatrixCommon::to_file "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L1270

";

%feature("docstring") casadi::MatrixCommon::nnz "

[INTERNAL] 
Expose base class functions.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1105-L1107

";

%feature("docstring") casadi::MatrixCommon::nnz_lower "

[INTERNAL] 
Get the number of non-zeros in the lower triangular half.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ao

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L190

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1110-L1112

";

%feature("docstring") casadi::MatrixCommon::nnz_upper "

[INTERNAL] 
Get the number of non-zeros in the upper triangular half.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ap

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L191

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1115-L1117

";

%feature("docstring") casadi::MatrixCommon::numel "

[INTERNAL] 
Get the number of elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ar

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L192

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1125-L1127

";

%feature("docstring") casadi::MatrixCommon::size1 "

[INTERNAL] 
Get the first dimension (i.e. number of rows)

Extra doc: https://github.com/casadi/casadi/wiki/L_1as

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L193

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1130-L1132

";

%feature("docstring") casadi::MatrixCommon::size2 "

[INTERNAL] 
Get the second dimension (i.e. number of columns)

Extra doc: https://github.com/casadi/casadi/wiki/L_1au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L194

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1135-L1137

";

%feature("docstring") casadi::MatrixCommon::size "

[INTERNAL] 
Get the size along a particular dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ay

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1145-L1147

>  casadi_int casadi::GenericMatrix< MatType >::size(casadi_int axis) const
------------------------------------------------------------------------
[INTERNAL] 
Get the size along a particular dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ay

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1145-L1147

";

";

%feature("docstring") casadi::MatrixCommon::is_empty "

[INTERNAL] 
Check if the sparsity is empty, i.e. if one of the dimensions is
 zero.

(or optionally both dimensions)

Extra doc: https://github.com/casadi/casadi/wiki/L_1az

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L196

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L148-L148

";

%feature("docstring") casadi::MatrixCommon::is_scalar "

[INTERNAL] 
Check if the matrix expression is scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L197

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1155-L1157

";

%feature("docstring") casadi::MatrixCommon::is_dense "

[INTERNAL] 
Check if the matrix expression is dense.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L198

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L153-L153

";

%feature("docstring") casadi::MatrixCommon::is_vector "

[INTERNAL] 
Check if the matrix is a row or column vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L199

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L168-L168

";

%feature("docstring") casadi::MatrixCommon::is_row "

[INTERNAL] 
Check if the matrix is a row vector (i.e. size1()==1)

Extra doc: https://github.com/casadi/casadi/wiki/L_1b4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L173-L173

";

%feature("docstring") casadi::MatrixCommon::is_column "

[INTERNAL] 
Check if the matrix is a column vector (i.e. size2()==1)

Extra doc: https://github.com/casadi/casadi/wiki/L_1b5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L178-L178

";

%feature("docstring") casadi::MatrixCommon::is_tril "

[INTERNAL] 
Check if the matrix is lower triangular.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L202

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L188-L188

";

%feature("docstring") casadi::MatrixCommon::is_triu "

[INTERNAL] 
Check if the matrix is upper triangular.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L203

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L183-L183

";

%feature("docstring") casadi::MatrixCommon::colind "

[INTERNAL] ";

%feature("docstring") casadi::MatrixCommon::row "

[INTERNAL] ";

%feature("docstring") casadi::MatrixCommon::dim "

[INTERNAL] 
Get string representation of dimensions.

The representation is e.g. \"4x5\" or \"4x5,10nz\"

Extra doc: https://github.com/casadi/casadi/wiki/L_1aw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L206

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1150-L1152

";

%feature("docstring") casadi::MatrixCommon::nz "

[INTERNAL] 
Access vector nonzero or slice of nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L210

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L258-L260

>  NonZeros<MatType, K> casadi::GenericMatrix< MatType >::nz(const K &k)
------------------------------------------------------------------------
[INTERNAL] 
Access vector nonzero or slice of nonzeros.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L210

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

%feature("docstring") casadi::MX::printme "

";

";

%feature("docstring") casadi::MX::get_row "

Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L194

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L194-L194

";

%feature("docstring") casadi::MX::get_colind "

Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L195

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L195-L195

";

%feature("docstring") casadi::MX::row "

Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200-L200

>  casadi_int casadi::GenericMatrix< MX  >::row(casadi_int el) const
------------------------------------------------------------------------

Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L200-L200

";

";

%feature("docstring") casadi::MX::colind "

Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201-L201

>  casadi_int casadi::GenericMatrix< MX  >::colind(casadi_int col) const
------------------------------------------------------------------------

Get the sparsity pattern. See the Sparsity class for details.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L201-L201

";

";

%feature("docstring") casadi::MX::interp1d "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L213

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1238-L1294

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
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1188-L1198

";

%feature("docstring") casadi::MX::cross "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L227

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1201-L1232

";

%feature("docstring") casadi::MX::skew "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L228

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1297-L1305

";

%feature("docstring") casadi::MX::inv_skew "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L229

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1308-L1313

";

%feature("docstring") casadi::MX::tril2symm "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L230

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1317-L1323

";

%feature("docstring") casadi::MX::triu2symm "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L231

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1341-L1347

";

%feature("docstring") casadi::MX::diff "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1551-L1576

";

%feature("docstring") casadi::MX::is_linear "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L235

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1520-L1522

";

%feature("docstring") casadi::MX::is_quadratic "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L236

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1525-L1527

";

%feature("docstring") casadi::MX::quadratic_coeff "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L237

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1530-L1538

";

%feature("docstring") casadi::MX::linear_coeff "

[INTERNAL] 
Functions called by friend functions defined here.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L239

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1541-L1548

";

%feature("docstring") casadi::MX::bilin "

Calculate bilinear form x^T A y.

Extra doc: https://github.com/casadi/casadi/wiki/L_1bo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L404

";

%feature("docstring") casadi::MX::rank1 "

Make a rank-1 update to a matrix A.

Calculates A + 1/2 * alpha * x*y'

Extra doc: https://github.com/casadi/casadi/wiki/L_1bp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L417

";

%feature("docstring") casadi::MX::mpower "

Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L996

";

%feature("docstring") casadi::MX::soc "

Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L997

";

%feature("docstring") casadi::MX::linearize "

Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L994

";

%feature("docstring") casadi::MX::gradient "

Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L992

";

%feature("docstring") casadi::MX::tangent "

Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L993

";

%feature("docstring") casadi::MX::jtimes "

Functions called by friend functions defined here

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L990

";

%feature("docstring") casadi::MX::sym "

Create a vector of length r of vectors of length p.

with nrow-by-ncol symbolic primitives

Extra doc: https://github.com/casadi/casadi/wiki/L_1dg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1060

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1060-L1062

>  static std::vector<std::vector<MX > > casadi::GenericMatrix< MX  >::sym(const std::string &name, casadi_int nrow, casadi_int ncol, casadi_int p, casadi_int r)
------------------------------------------------------------------------

Create a vector of length r of vectors of length p.

with nrow-by-ncol symbolic primitives

Extra doc: https://github.com/casadi/casadi/wiki/L_1dg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1060

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1060-L1062

";

";

%feature("docstring") casadi::MX::zeros "

Create a dense matrix or a matrix with specified sparsity with all 
entries 
zero.

Extra doc: https://github.com/casadi/casadi/wiki/L_1dh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1073

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1073-L1075

>  static MX  casadi::GenericMatrix< MX  >::zeros(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------

Create a dense matrix or a matrix with specified sparsity with all 
entries 
zero.

Extra doc: https://github.com/casadi/casadi/wiki/L_1dh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1073

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1073-L1075

";

";

%feature("docstring") casadi::MX::ones "

Create a dense matrix or a matrix with specified sparsity with all 
entries 
one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1di

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1086

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1086-L1088

>  static MX  casadi::GenericMatrix< MX  >::ones(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------

Create a dense matrix or a matrix with specified sparsity with all 
entries 
one.

Extra doc: https://github.com/casadi/casadi/wiki/L_1di

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1086

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1086-L1088

";

";

%feature("docstring") casadi::MX::binary "

Create nodes by their ID.

Extra doc: https://github.com/casadi/casadi/wiki/L_r1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L398

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L508-L533

";

%feature("docstring") casadi::MX::unary "

Create nodes by their ID.

Extra doc: https://github.com/casadi/casadi/wiki/L_r1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L399

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L535-L537

";

%feature("docstring") casadi::MX::inf "

create a matrix with all inf

Extra doc: https://github.com/casadi/casadi/wiki/L_r2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L408

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L555-L557

>  MX casadi::MX::inf(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------

create a matrix with all inf

Extra doc: https://github.com/casadi/casadi/wiki/L_r2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L408

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L555-L557

";

";

%feature("docstring") casadi::MX::nan "

create a matrix with all nan

Extra doc: https://github.com/casadi/casadi/wiki/L_r3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L417

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L567-L569

>  MX casadi::MX::nan(const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------

create a matrix with all nan

Extra doc: https://github.com/casadi/casadi/wiki/L_r3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L417

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L567-L569

";

";

%feature("docstring") casadi::MX::get "

[INTERNAL] 
Get a const pointer to the node.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L427

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L539-L541

>  MXNode * casadi::MX::get() const
------------------------------------------------------------------------
[INTERNAL] 
Get a const pointer to the node.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L427

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L539-L541

";

";

%feature("docstring") casadi::MX::set "

Set a submatrix, two arguments

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L475

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L239-L297

>  void casadi::MX::set(const MX &m, bool ind1, const Matrix< casadi_int > &rr, const Matrix< casadi_int > &cc)
------------------------------------------------------------------------

Set a submatrix, two arguments

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L475

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L239-L297

";

";

%feature("docstring") casadi::MX::get_nz "

Get a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L488

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L428-L431

>  void casadi::MX::get_nz(MX &m, bool ind1, const MX &inner, const MX &outer) const
------------------------------------------------------------------------

Get a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L488

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L428-L431

";

";

%feature("docstring") casadi::MX::set_nz "

Set a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L496

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L496-L496

>  void casadi::MX::set_nz(const MX &m, bool ind1, casadi_int kk)
------------------------------------------------------------------------

Set a set of nonzeros

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L496

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L496-L496

";

";

%feature("docstring") casadi::MX::einstein "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L647-L653

>  MX casadi::MX::einstein(const MX &A, const MX &B, const std::vector< casadi_int > &dim_a, const std::vector< casadi_int > &dim_b, const std::vector< casadi_int > &dim_c, const std::vector< casadi_int > &a, const std::vector< casadi_int > &b, const std::vector< casadi_int > &c)
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L647-L653

";

";

%feature("docstring") casadi::MX::is_equal "

[INTERNAL] 
Functions called by friend functions defined for  
GenericExpression

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L531

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L811-L813

";

%feature("docstring") casadi::MX::mmin "

[INTERNAL] 
Functions called by friend functions defined for  
GenericExpression

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L532

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L815-L817

";

%feature("docstring") casadi::MX::mmax "

[INTERNAL] 
Functions called by friend functions defined for  
GenericExpression

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L533

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L819-L821

";

%feature("docstring") casadi::MX::horzcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L538

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L994-L1028

";

%feature("docstring") casadi::MX::diagcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L539

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1030-L1053

";

%feature("docstring") casadi::MX::vertcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L540

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1055-L1094

";

%feature("docstring") casadi::MX::horzsplit "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L541

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1096-L1111

";

%feature("docstring") casadi::MX::diagsplit "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L542

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1113-L1128

";

%feature("docstring") casadi::MX::vertsplit "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L544

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1130-L1151

";

%feature("docstring") casadi::MX::blockcat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L545

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1153-L1172

";

%feature("docstring") casadi::MX::mtimes "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L546

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L629-L637

";

%feature("docstring") casadi::MX::mac "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L547

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L669-L690

";

%feature("docstring") casadi::MX::reshape "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L549

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1206-L1212

>  MX casadi::MX::reshape(const MX &x, const Sparsity &sp)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L549

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1206-L1212

";

";

%feature("docstring") casadi::MX::kron "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L550

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1830-L1843

";

%feature("docstring") casadi::MX::repmat "

[INTERNAL] 
Functions called by friend functions defined for  
SparsityInterface

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L551

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1845-L1857

";

%feature("docstring") casadi::MX::jacobian "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L556

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1685-L1694

";

%feature("docstring") casadi::MX::hessian "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L558

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1701-L1710

>  MX casadi::MX::hessian(const MX &f, const MX &x, MX &g, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L558

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1701-L1710

";

";

%feature("docstring") casadi::MX::forward "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L560

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1713-L1740

";

%feature("docstring") casadi::MX::reverse "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L565

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1743-L1772

";

%feature("docstring") casadi::MX::which_depends "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L569

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1774-L1776

";

%feature("docstring") casadi::MX::substitute "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L572

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1369-L1389

>  std::vector< MX > casadi::MX::substitute(const std::vector< MX > &ex, const std::vector< MX > &v, const std::vector< MX > &vdef)
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L572

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1369-L1389

";

";

%feature("docstring") casadi::MX::substitute_inplace "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L575

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1343-L1363

";

%feature("docstring") casadi::MX::solve "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L579

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1876-L1879

>  MX casadi::MX::solve(const MX &a, const MX &b, const std::string &lsolver, const Dict &dict=Dict())
------------------------------------------------------------------------
[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L579

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1876-L1879

";

";

%feature("docstring") casadi::MX::inv_minor "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L581

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1786-L1788

";

%feature("docstring") casadi::MX::inv_node "

Inverse node.

Extra doc: https://github.com/casadi/casadi/wiki/L_re

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L789

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L789-L791

>  MX casadi::MX::inv_node(const MX &x)
------------------------------------------------------------------------

Inverse node.

Extra doc: https://github.com/casadi/casadi/wiki/L_re

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L789

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L789-L791

";

";

%feature("docstring") casadi::MX::inv "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L583

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1790-L1792

";

%feature("docstring") casadi::MX::pinv "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L584

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1881-L1887

";

%feature("docstring") casadi::MX::expm_const "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L586

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1889-L1894

";

%feature("docstring") casadi::MX::expm "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L587

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1896-L1899

";

%feature("docstring") casadi::MX::n_nodes "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L588

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1316-L1319

";

%feature("docstring") casadi::MX::print_operator "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L589

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1339-L1341

";

%feature("docstring") casadi::MX::extract "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L590

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1505-L1676

";

%feature("docstring") casadi::MX::shared "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L592

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1678-L1683

";

%feature("docstring") casadi::MX::if_else "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L594

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1214-L1234

";

%feature("docstring") casadi::MX::conditional "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L596

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1236-L1268

";

%feature("docstring") casadi::MX::depends_on "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L598

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1907-L1923

";

%feature("docstring") casadi::MX::simplify "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L599

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1194-L1196

";

%feature("docstring") casadi::MX::dot "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L600

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L692-L694

";

%feature("docstring") casadi::MX::mrdivide "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L601

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L721-L724

";

%feature("docstring") casadi::MX::mldivide "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L602

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L726-L729

";

%feature("docstring") casadi::MX::norm_2 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L603

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1174-L1180

";

%feature("docstring") casadi::MX::norm_fro "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L604

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1182-L1184

";

%feature("docstring") casadi::MX::norm_1 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L605

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1186-L1188

";

%feature("docstring") casadi::MX::norm_inf "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L606

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1190-L1192

";

%feature("docstring") casadi::MX::unite "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L607

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1270-L1294

";

%feature("docstring") casadi::MX::trace "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L608

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1296-L1303

";

%feature("docstring") casadi::MX::diag "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L609

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1305-L1314

";

%feature("docstring") casadi::MX::sum2 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L610

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1321-L1323

";

%feature("docstring") casadi::MX::sum1 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L611

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1325-L1327

";

%feature("docstring") casadi::MX::polyval "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L612

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1329-L1337

";

%feature("docstring") casadi::MX::det "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L613

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1778-L1780

";

%feature("docstring") casadi::MX::symvar "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L614

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1794-L1797

";

%feature("docstring") casadi::MX::nullspace "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L615

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1901-L1905

";

%feature("docstring") casadi::MX::repsum "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L616

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1859-L1861

";

%feature("docstring") casadi::MX::densify "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L617

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L867-L878

";

%feature("docstring") casadi::MX::_bilin "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L618

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2158-L2160

";

%feature("docstring") casadi::MX::_rank1 "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L619

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2162-L2164

";

%feature("docstring") casadi::MX::project "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L620

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L850-L865

";

%feature("docstring") casadi::MX::cumsum "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L621

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L655-L667

";

%feature("docstring") casadi::MX::_logsumexp "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L622

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2166-L2168

";

%feature("docstring") casadi::MX::cse "

[INTERNAL] 
Functions called by friend functions defined for  GenericMatrix

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L623

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L1996-L2062

";

%feature("docstring") casadi::MX::find "

Find first nonzero, returned as row index.

If failed, returns the number of rows

Extra doc: https://github.com/casadi/casadi/wiki/L_r7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L686

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L686-L688

>  MX casadi::MX::find(const MX &x)
------------------------------------------------------------------------

Find first nonzero, returned as row index.

If failed, returns the number of rows

Extra doc: https://github.com/casadi/casadi/wiki/L_r7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L686

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L686-L688

";

";

%feature("docstring") casadi::MX::low "

Find first nonzero.

If failed, returns the number of rows

Extra doc: https://github.com/casadi/casadi/wiki/L_r8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L695

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L695-L697

>  MX casadi::MX::low(const MX &v, const MX &p, const Dict &options=Dict())
------------------------------------------------------------------------

Find first nonzero.

If failed, returns the number of rows

Extra doc: https://github.com/casadi/casadi/wiki/L_r8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L695

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L695-L697

";

";

%feature("docstring") casadi::MX::graph_substitute "

Substitute multiple expressions in graph.

Substitute variable var with expression expr in multiple expressions, 

preserving nodes

Extra doc: https://github.com/casadi/casadi/wiki/L_ra

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L716

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L716-L720

>  std::vector<MX> casadi::MX::graph_substitute(const std::vector< MX > &ex, const std::vector< MX > &v, const std::vector< MX > &vdef)
------------------------------------------------------------------------

Substitute multiple expressions in graph.

Substitute variable var with expression expr in multiple expressions, 

preserving nodes

Extra doc: https://github.com/casadi/casadi/wiki/L_ra

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L716

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L716-L720

";

";

%feature("docstring") casadi::MX::matrix_expand "

Expand  MX graph to SXFunction call.

Expand the given expression e, optionally supplying expressions 
contained 
in it at which expansion should stop.

Extra doc: https://github.com/casadi/casadi/wiki/L_rc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L741

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L741-L745

>  std::vector<MX> casadi::MX::matrix_expand(const std::vector< MX > &e, const std::vector< MX > &boundary=std::vector< MX >(), const Dict &options=Dict())
------------------------------------------------------------------------

Expand  MX graph to SXFunction call.

Expand the given expression e, optionally supplying expressions 
contained 
in it at which expansion should stop.

Extra doc: https://github.com/casadi/casadi/wiki/L_rc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L741

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L741-L745

";

";

%feature("docstring") casadi::MX::lift "

Lift the expression.

Experimental feature

Extra doc: https://github.com/casadi/casadi/wiki/L_rd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L782

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L782-L784

>  MX casadi::MX::lift(const MX &x, const MX &x_guess)
------------------------------------------------------------------------

Lift the expression.

Experimental feature

Extra doc: https://github.com/casadi/casadi/wiki/L_rd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L782

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L782-L784

";

";

%feature("docstring") casadi::MX::evalf "

Evaluates the expression numerically.

An error is raised when the expression contains symbols

Extra doc: https://github.com/casadi/casadi/wiki/L_rf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L798

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L798-L800

>  DM casadi::MX::evalf(const MX &expr)
------------------------------------------------------------------------

Evaluates the expression numerically.

An error is raised when the expression contains symbols

Extra doc: https://github.com/casadi/casadi/wiki/L_rf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L798

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L798-L800

";

";

%feature("docstring") casadi::MX::bspline "

";

";

%feature("docstring") casadi::MX::convexify "

";

";

%feature("docstring") casadi::MX::ad_forward "

[INTERNAL] 
Called from MXFunction.

Extra doc: https://github.com/casadi/casadi/wiki/L_ro

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L863

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2179-L2186

";

%feature("docstring") casadi::MX::ad_reverse "

[INTERNAL] 
Called from MXFunction.

Extra doc: https://github.com/casadi/casadi/wiki/L_ro

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L865

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2188-L2195

";

%feature("docstring") casadi::MX::MX "

[INTERNAL] 
Construct constant matrix with a given sparsity and values.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L870

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L115-L117

>  casadi::MX::MX(const Sparsity &sp, double val, bool dummy)
------------------------------------------------------------------------
[INTERNAL] 
Construct constant matrix with a given sparsity and values.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L870

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L115-L117

";

";

%feature("docstring") casadi::MX::sparsity "

[INTERNAL] 
Get the sparsity pattern.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L207

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100-L1102

>  const Sparsity & casadi::GenericMatrix< MX  >::sparsity() const
------------------------------------------------------------------------
[INTERNAL] 
Get the sparsity pattern.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b9

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L207

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1100-L1102

";

";

%feature("docstring") casadi::MX::__nonzero__ "

Returns the truth value of an  MX expression.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L184

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L138-L140

";

%feature("docstring") casadi::MX::get_sparsity "

Get an owning reference to the sparsity pattern.

Extra doc: https://github.com/casadi/casadi/wiki/L_qd

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L189

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L189-L189

";

%feature("docstring") casadi::MX::erase "

Erase a submatrix (leaving structural zeros in its place)

Erase elements of a matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_qf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L204

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L605-L617

>  void casadi::MX::erase(const std::vector< casadi_int > &rr, bool ind1=false)
------------------------------------------------------------------------

Erase a submatrix (leaving structural zeros in its place)

Erase elements of a matrix

Extra doc: https://github.com/casadi/casadi/wiki/L_qf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L204

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L605-L617

";

";

%feature("docstring") casadi::MX::enlarge "

Enlarge matrix.

Make the matrix larger by inserting empty rows and columns, keeping 
the 
existing non-zeros

Extra doc: https://github.com/casadi/casadi/wiki/L_qg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L211

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L619-L627

";

%feature("docstring") casadi::MX::dep "

Get the nth dependency as  MX.

Extra doc: https://github.com/casadi/casadi/wiki/L_qj

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L236

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L731-L733

";

%feature("docstring") casadi::MX::n_out "

Number of outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_qk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L241

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L842-L844

";

%feature("docstring") casadi::MX::get_output "

Get an output.

Extra doc: https://github.com/casadi/casadi/wiki/L_ql

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L246

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L846-L848

";

%feature("docstring") casadi::MX::n_dep "

Get the number of dependencies of a binary  SXElem.

Extra doc: https://github.com/casadi/casadi/wiki/L_qm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L251

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L735-L737

";

%feature("docstring") casadi::MX::name "

Get the name.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L254

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L739-L741

";

%feature("docstring") casadi::MX::is_symbolic "

Check if symbolic.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L263

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L743-L745

";

%feature("docstring") casadi::MX::is_constant "

Check if constant.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L266

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L747-L749

";

%feature("docstring") casadi::MX::is_call "

Check if evaluation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L269

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L751-L753

";

%feature("docstring") casadi::MX::which_function "

Get function - only valid when  is_call() is true.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L272

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L755-L757

";

%feature("docstring") casadi::MX::is_output "

Check if evaluation output.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L275

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L759-L761

";

%feature("docstring") casadi::MX::which_output "

Get the index of evaluation output - only valid when  is_output() is true.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L278

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L763-L765

";

%feature("docstring") casadi::MX::is_op "

Is it a certain operation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L281

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L767-L769

";

%feature("docstring") casadi::MX::is_multiplication "

Check if multiplication.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L284

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L771-L773

";

%feature("docstring") casadi::MX::is_commutative "

Check if commutative operation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L287

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L823-L828

";

%feature("docstring") casadi::MX::is_norm "

Check if norm.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L290

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L775-L777

";

%feature("docstring") casadi::MX::is_valid_input "

Check if matrix can be used to define function inputs.

Valid inputs for MXFunctions are combinations of Reshape, 
concatenations 
and SymbolicMX

Extra doc: https://github.com/casadi/casadi/wiki/L_qn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L297

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L898-L900

";

%feature("docstring") casadi::MX::n_primitives "

Get the number of primitives for MXFunction inputs/outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_qo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L302

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L902-L904

";

%feature("docstring") casadi::MX::primitives "

Get primitives.

Extra doc: https://github.com/casadi/casadi/wiki/L_qp

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L307

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L906-L912

";

%feature("docstring") casadi::MX::split_primitives "

Split up an expression along symbolic primitives.

Extra doc: https://github.com/casadi/casadi/wiki/L_qq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L312

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L914-L920

";

%feature("docstring") casadi::MX::join_primitives "

Join an expression along symbolic primitives.

Extra doc: https://github.com/casadi/casadi/wiki/L_qr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L317

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L922-L928

";

%feature("docstring") casadi::MX::is_eye "

check if identity

Extra doc: https://github.com/casadi/casadi/wiki/L_qu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L339

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L938-L940

";

%feature("docstring") casadi::MX::is_zero "

check if zero (note that false negative answers are possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_qv

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L344

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L942-L948

";

%feature("docstring") casadi::MX::is_one "

check if zero (note that false negative answers are possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_qw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L349

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L950-L952

";

%feature("docstring") casadi::MX::is_minus_one "

check if zero (note that false negative answers are possible)

Extra doc: https://github.com/casadi/casadi/wiki/L_qx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L354

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L954-L956

";

%feature("docstring") casadi::MX::is_transpose "

Is the expression a transpose?

Extra doc: https://github.com/casadi/casadi/wiki/L_qy

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L359

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L958-L960

";

%feature("docstring") casadi::MX::is_regular "

Checks if expression does not contain NaN or Inf.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L362

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L962-L968

";

%feature("docstring") casadi::MX::is_binary "

Is binary operation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L365

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L787-L789

";

%feature("docstring") casadi::MX::is_unary "

Is unary operation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L368

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L791-L793

";

%feature("docstring") casadi::MX::op "

Get operation type.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L371

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L795-L797

";

%feature("docstring") casadi::MX::info "

Obtain information about node

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L374

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L799-L801

";

%feature("docstring") casadi::MX::serialize "

Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_qz

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L379

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L803-L805

";

%feature("docstring") casadi::MX::attachAssert "

returns itself, but with an assertion attached

If y does not evaluate to 1, a runtime error is raised

Extra doc: https://github.com/casadi/casadi/wiki/L_rg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L810

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L700-L705

";

%feature("docstring") casadi::MX::monitor "

Monitor an expression.

Returns itself, but with the side effect of printing the nonzeros 
along 
with a comment

Extra doc: https://github.com/casadi/casadi/wiki/L_rh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L817

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L707-L709

";

%feature("docstring") casadi::MX::T "

Transpose the matrix.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L820

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L970-L972

";

%feature("docstring") casadi::MX::mapping "

Get an IM representation of a GetNonzeros or SetNonzeros node.

Extra doc: https://github.com/casadi/casadi/wiki/L_ri

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L825

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L830-L832

";

%feature("docstring") casadi::MX::eval_mx "

Evaluate the  MX node with new symbolic dependencies.

Extra doc: https://github.com/casadi/casadi/wiki/L_rn

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.hpp#L856

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/mx.cpp#L2171-L2177

";

%feature("docstring") casadi::MX::nnz "

Get the number of (structural) non-zero elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1an

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L84

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1105-L1107

";

%feature("docstring") casadi::MX::nnz_lower "

Get the number of non-zeros in the lower triangular half.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ao

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L89

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1110-L1112

";

%feature("docstring") casadi::MX::nnz_upper "

Get the number of non-zeros in the upper triangular half.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ap

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L94

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1115-L1117

";

%feature("docstring") casadi::MX::nnz_diag "

Get get the number of non-zeros on the diagonal.

Extra doc: https://github.com/casadi/casadi/wiki/L_1aq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L99

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1120-L1122

";

%feature("docstring") casadi::MX::numel "

Get the number of elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ar

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L104

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1125-L1127

";

%feature("docstring") casadi::MX::size1 "

Get the first dimension (i.e. number of rows)

Extra doc: https://github.com/casadi/casadi/wiki/L_1as

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L109

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1130-L1132

";

%feature("docstring") casadi::MX::rows "

Get the number of rows, Octave-style syntax.

Extra doc: https://github.com/casadi/casadi/wiki/L_1at

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L114

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L114-L114

";

%feature("docstring") casadi::MX::size2 "

Get the second dimension (i.e. number of columns)

Extra doc: https://github.com/casadi/casadi/wiki/L_1au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L119

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1135-L1137

";

%feature("docstring") casadi::MX::columns "

Get the number of columns, Octave-style syntax.

Extra doc: https://github.com/casadi/casadi/wiki/L_1av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L124

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L124-L124

";

%feature("docstring") casadi::MX::dim "

Get string representation of dimensions.

The representation is e.g. \"4x5\" or \"4x5,10nz\"

Extra doc: https://github.com/casadi/casadi/wiki/L_1aw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L131

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1150-L1152

";

%feature("docstring") casadi::MX::size "

Get the size along a particular dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ay

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L141

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1145-L1147

>  casadi_int casadi::GenericMatrix< MX  >::size(casadi_int axis) const
------------------------------------------------------------------------

Get the size along a particular dimensions.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ay

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L141

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1145-L1147

";

";

%feature("docstring") casadi::MX::is_empty "

Check if the sparsity is empty, i.e. if one of the dimensions is zero.

(or optionally both dimensions)

Extra doc: https://github.com/casadi/casadi/wiki/L_1az

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L148

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L148-L148

";

%feature("docstring") casadi::MX::is_dense "

Check if the matrix expression is dense.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L153-L153

";

%feature("docstring") casadi::MX::is_scalar "

Check if the matrix expression is scalar.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L158

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L1155-L1157

";

%feature("docstring") casadi::MX::is_square "

Check if the matrix expression is square.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L163

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L163-L163

";

%feature("docstring") casadi::MX::is_vector "

Check if the matrix is a row or column vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L168

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L168-L168

";

%feature("docstring") casadi::MX::is_row "

Check if the matrix is a row vector (i.e.  size1()==1)

Extra doc: https://github.com/casadi/casadi/wiki/L_1b4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L173

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L173-L173

";

%feature("docstring") casadi::MX::is_column "

Check if the matrix is a column vector (i.e.  size2()==1)

Extra doc: https://github.com/casadi/casadi/wiki/L_1b5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L178

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L178-L178

";

%feature("docstring") casadi::MX::is_triu "

Check if the matrix is upper triangular.

Extra doc: https://github.com/casadi/casadi/wiki/L_1b6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L183

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_matrix.hpp#L183-L183

";

%feature("docstring") casadi::MX::is_tril "

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

Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L100-L102

";

%feature("docstring") casadi::MX::disp "

Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L104-L110

";

%feature("docstring") casadi::MX::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::MX::is_null "

Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L74-L76

";

%feature("docstring") casadi::MX::__hash__ "

Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L130-L132

";


// File: classcasadi_1_1Newton.xml
%feature("docstring") casadi::Newton "

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

A symbolic NLP representation.

Joel Andersson

::

  Extra doc: https://github.com/casadi/casadi/wiki/L_1e2 
  



C++ includes: nlp_builder.hpp
";

%feature("docstring") casadi::NlpBuilder::import_nl "

Import an .nl file.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.hpp#L74

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.cpp#L33-L36

";

%feature("docstring") casadi::NlpBuilder::type_name "

Readable name of the class.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.hpp#L77

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.hpp#L77-L77

";

%feature("docstring") casadi::NlpBuilder::disp "

Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.hpp#L80

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.cpp#L38-L46

";

%feature("docstring") casadi::NlpBuilder::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.hpp#L83

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_builder.hpp#L83-L87

";


// File: classcasadi_1_1Nlpsol.xml
%feature("docstring") casadi::Nlpsol "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1NonZeros.xml
%feature("docstring") casadi::NonZeros "

Access to a set of nonzeros.

NonZeros class for  Matrix NonZeros is the return type for operator[] of the
  Matrix class, it allows access to the value as well as changing the parent
 
object 
Joel Andersson

C++ includes: nonzeros.hpp
";

%feature("docstring") casadi::NonZeros::NonZeros "

Default copy constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nonzeros.hpp#L46

>  casadi::NonZeros< M, K >::NonZeros(const NonZeros< M, K > &y)=default
------------------------------------------------------------------------

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

Clear constraints.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L166

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L104-L110

>  void casadi::Opti::subject_to()
------------------------------------------------------------------------

Clear constraints.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L166

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L104-L110

";

";

%feature("docstring") casadi::Opti::set_initial "

Set initial guess for decision variables

::

  * opti.set_initial(x, 2)
  * opti.set_initial(10*x(1), 2)
  * 



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L190

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L130-L136

>  void casadi::Opti::set_initial(const std::vector< MX > &assignments)
------------------------------------------------------------------------

Set initial guess for decision variables

::

  * opti.set_initial(x, 2)
  * opti.set_initial(10*x(1), 2)
  * 



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L190

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L130-L136

";

";

%feature("docstring") casadi::Opti::set_value "

Set value of parameter.

Each parameter must be given a value before 'solve' can be called

Extra doc: https://github.com/casadi/casadi/wiki/L_1d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L147-L153

>  void casadi::Opti::set_value(const std::vector< MX > &assignments)
------------------------------------------------------------------------

Set value of parameter.

Each parameter must be given a value before 'solve' can be called

Extra doc: https://github.com/casadi/casadi/wiki/L_1d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L147-L153

";

";

%feature("docstring") casadi::Opti::value "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L188-L194

>  DM casadi::Opti::value(const SX &x, const std::vector< MX > &values=std::vector< MX >()) const
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L188-L194

";

";

%feature("docstring") casadi::Opti::callback_class "

Helper methods for callback()

Do not use directly.

Extra doc: https://github.com/casadi/casadi/wiki/L_1p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L381

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L373-L385

>  void casadi::Opti::callback_class()
------------------------------------------------------------------------

Helper methods for callback()

Do not use directly.

Extra doc: https://github.com/casadi/casadi/wiki/L_1p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L381

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L373-L385

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
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L53-L59

";

%feature("docstring") casadi::Opti::parameter "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L76-L82

";

%feature("docstring") casadi::Opti::minimize "

Set objective.

Objective must be a scalar. Default objective: 0 When method is called
 
multiple times, the last call takes effect

Extra doc: https://github.com/casadi/casadi/wiki/L_1a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L133

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L84-L90

";

%feature("docstring") casadi::Opti::solver "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L113-L121

";

%feature("docstring") casadi::Opti::solve "

Crunch the numbers; solve the problem.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L204

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L155-L161

";

%feature("docstring") casadi::Opti::solve_limited "

Crunch the numbers; solve the problem.

Allows the solver to return without error when an iteration or time 
limit 
is reached

Extra doc: https://github.com/casadi/casadi/wiki/L_1e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L212

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L163-L169

";

%feature("docstring") casadi::Opti::stats "

Get statistics.

nlpsol stats are passed as-is. No stability can be guaranteed about 
this 
part of the API

Extra doc: https://github.com/casadi/casadi/wiki/L_1f

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L234

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L196-L202

";

%feature("docstring") casadi::Opti::return_status "

Get return status of solver.



::

     passed as-is from nlpsol
  

No stability can be guaranteed about this part of the API

Extra doc: https://github.com/casadi/casadi/wiki/L_1g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L204-L210

";

%feature("docstring") casadi::Opti::initial "

get assignment expressions for initial values

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L245

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L212-L218

";

%feature("docstring") casadi::Opti::value_variables "

get assignment expressions for latest values

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L248

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L220-L226

";

%feature("docstring") casadi::Opti::dual "

get the dual variable

m must be a constraint expression. The returned value is still a 
symbolic 
expression. Use  value on it to obtain the numerical value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L236-L242

";

%feature("docstring") casadi::Opti::nx "

Number of (scalarised) decision variables.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L261

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L244-L250

";

%feature("docstring") casadi::Opti::np "

Number of (scalarised) parameters.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L264

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L252-L258

";

%feature("docstring") casadi::Opti::ng "

Number of (scalarised) constraints.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L267

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L260-L266

";

%feature("docstring") casadi::Opti::x "

Get all (scalarised) decision variables as a symbolic column vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L270

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L268-L274

";

%feature("docstring") casadi::Opti::p "

Get all (scalarised) parameters as a symbolic column vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L273

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L276-L282

";

%feature("docstring") casadi::Opti::g "

Get all (scalarised) constraint expressions as a column vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L276

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L284-L290

";

%feature("docstring") casadi::Opti::f "

Get objective expression.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L279

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L292-L298

";

%feature("docstring") casadi::Opti::lbg "

Get all (scalarised) bounds on constraints as a column vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L282

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L300-L306

";

%feature("docstring") casadi::Opti::lam_g "

Get all (scalarised) dual variables as a symbolic column vector.

Useful for obtaining the Lagrange Hessian:

::

  * sol.value(hessian(opti.f+opti.lam_g'*opti.g,opti.x)) % MATLAB
  * sol.value(hessian(opti.f+dot(opti.lam_g,opti.g),opti.x)[0]) # Python
  * 



Extra doc: https://github.com/casadi/casadi/wiki/L_1i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L294

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L317-L323

";

%feature("docstring") casadi::Opti::to_function "

";

";

%feature("docstring") casadi::Opti::debug "

Get a copy with advanced functionality.

You get access to more methods, but you have no guarantees about API 

stability

The copy is effectively a deep copy: Updating the state of the copy 
does 
not update the original.

Extra doc: https://github.com/casadi/casadi/wiki/L_1l

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L334

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L741-L743

";

%feature("docstring") casadi::Opti::advanced "

Get a copy with advanced functionality.

You get access to more methods, but you have no guarantees about API 

stability

The copy is effectively a deep copy: Updating the state of the copy 
does 
not update the original.

Extra doc: https://github.com/casadi/casadi/wiki/L_1m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L344

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L744-L746

";

%feature("docstring") casadi::Opti::copy "

Get a copy of the.

The copy is effectively a deep copy: Updating the state of the copy 
does 
not update the original.

Extra doc: https://github.com/casadi/casadi/wiki/L_1n

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L352

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L747-L749

";

%feature("docstring") casadi::Opti::update_user_dict "

";

";

%feature("docstring") casadi::Opti::user_dict "

Get user data.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L363

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L400-L406

";

%feature("docstring") casadi::Opti::type_name "

Readable name of the class.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L366

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L366-L366

";

%feature("docstring") casadi::Opti::disp "

Print representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L369

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L644-L664

";

%feature("docstring") casadi::Opti::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L372

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L666-L670

";

%feature("docstring") casadi::Opti::~Opti "

[INTERNAL] 
Destructor.

Extra doc: https://github.com/casadi/casadi/wiki/L_1q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L390

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L390-L390

";

%feature("docstring") casadi::Opti::class_name "

Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L100-L102

";

%feature("docstring") casadi::Opti::is_null "

Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L74-L76

";

%feature("docstring") casadi::Opti::__hash__ "

Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L130-L132

";


// File: classcasadi_1_1OptiAdvanced.xml
%feature("docstring") casadi::OptiAdvanced "

C++ includes: optistack.hpp
";

%feature("docstring") casadi::OptiAdvanced::subject_to "

Clear constraints.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L166

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L104-L110

>  void casadi::Opti::subject_to()
------------------------------------------------------------------------

Clear constraints.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L166

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L104-L110

";

";

%feature("docstring") casadi::OptiAdvanced::set_initial "

Set initial guess for decision variables

::

  * opti.set_initial(x, 2)
  * opti.set_initial(10*x(1), 2)
  * 



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L190

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L130-L136

>  void casadi::Opti::set_initial(const std::vector< MX > &assignments)
------------------------------------------------------------------------

Set initial guess for decision variables

::

  * opti.set_initial(x, 2)
  * opti.set_initial(10*x(1), 2)
  * 



Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L190

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L130-L136

";

";

%feature("docstring") casadi::OptiAdvanced::set_value "

Set value of parameter.

Each parameter must be given a value before 'solve' can be called

Extra doc: https://github.com/casadi/casadi/wiki/L_1d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L147-L153

>  void casadi::Opti::set_value(const std::vector< MX > &assignments)
------------------------------------------------------------------------

Set value of parameter.

Each parameter must be given a value before 'solve' can be called

Extra doc: https://github.com/casadi/casadi/wiki/L_1d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L200

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L147-L153

";

";

%feature("docstring") casadi::OptiAdvanced::value "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L188-L194

>  DM casadi::Opti::value(const SX &x, const std::vector< MX > &values=std::vector< MX >()) const
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L188-L194

";

";

%feature("docstring") casadi::OptiAdvanced::callback_class "

Helper methods for callback()

Do not use directly.

Extra doc: https://github.com/casadi/casadi/wiki/L_1p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L381

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L373-L385

>  void casadi::Opti::callback_class()
------------------------------------------------------------------------

Helper methods for callback()

Do not use directly.

Extra doc: https://github.com/casadi/casadi/wiki/L_1p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L381

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L373-L385

";

";

%feature("docstring") casadi::OptiAdvanced::symvar "

Get symbols present in expression.

Returned vector is ordered according to the order of  variable()/parameter()
 calls used to create the variables

Extra doc: https://github.com/casadi/casadi/wiki/L_1u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L497

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L440-L446

>  std::vector< MX > casadi::OptiAdvanced::symvar(const MX &expr, VariableType type) const
------------------------------------------------------------------------

Get symbols present in expression.

Returned vector is ordered according to the order of  variable()/parameter()
 calls used to create the variables

Extra doc: https://github.com/casadi/casadi/wiki/L_1u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L497

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L440-L446

";

";

%feature("docstring") casadi::OptiAdvanced::~OptiAdvanced "

Destructor.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L479

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L479-L479

";

%feature("docstring") casadi::OptiAdvanced::casadi_solver "

Get the underlying CasADi solver of the  Opti stack.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L483

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L408-L414

";

%feature("docstring") casadi::OptiAdvanced::is_parametric "

return true if expression is only dependant on  Opti parameters, not 
variables

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L486

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L416-L422

";

%feature("docstring") casadi::OptiAdvanced::canon_expr "

Interpret an expression (for internal use only)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L501

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L448-L454

";

%feature("docstring") casadi::OptiAdvanced::get_meta "

Get meta-data of symbol (for internal use only)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L504

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L456-L462

";

%feature("docstring") casadi::OptiAdvanced::get_meta_con "

Get meta-data of symbol (for internal use only)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L507

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L464-L470

";

%feature("docstring") casadi::OptiAdvanced::set_meta "

Set meta-data of an expression.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L510

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L472-L478

";

%feature("docstring") casadi::OptiAdvanced::set_meta_con "

Set meta-data of an expression.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L513

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L480-L486

";

%feature("docstring") casadi::OptiAdvanced::bake "

Fix the structure of the optimization problem.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L544

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L672-L678

";

%feature("docstring") casadi::OptiAdvanced::variable "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L53-L59

";

%feature("docstring") casadi::OptiAdvanced::parameter "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L76-L82

";

%feature("docstring") casadi::OptiAdvanced::minimize "

Set objective.

Objective must be a scalar. Default objective: 0 When method is called
 
multiple times, the last call takes effect

Extra doc: https://github.com/casadi/casadi/wiki/L_1a

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L133

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L84-L90

";

%feature("docstring") casadi::OptiAdvanced::solver "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L113-L121

";

%feature("docstring") casadi::OptiAdvanced::solve "

Crunch the numbers; solve the problem.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L204

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L155-L161

";

%feature("docstring") casadi::OptiAdvanced::solve_limited "

Crunch the numbers; solve the problem.

Allows the solver to return without error when an iteration or time 
limit 
is reached

Extra doc: https://github.com/casadi/casadi/wiki/L_1e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L212

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L163-L169

";

%feature("docstring") casadi::OptiAdvanced::stats "

Get statistics.

nlpsol stats are passed as-is. No stability can be guaranteed about 
this 
part of the API

Extra doc: https://github.com/casadi/casadi/wiki/L_1f

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L234

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L196-L202

";

%feature("docstring") casadi::OptiAdvanced::return_status "

Get return status of solver.



::

     passed as-is from nlpsol
  

No stability can be guaranteed about this part of the API

Extra doc: https://github.com/casadi/casadi/wiki/L_1g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L242

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L204-L210

";

%feature("docstring") casadi::OptiAdvanced::initial "

get assignment expressions for initial values

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L245

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L212-L218

";

%feature("docstring") casadi::OptiAdvanced::value_variables "

get assignment expressions for latest values

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L248

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L220-L226

";

%feature("docstring") casadi::OptiAdvanced::dual "

get the dual variable

m must be a constraint expression. The returned value is still a 
symbolic 
expression. Use  value on it to obtain the numerical value.

Extra doc: https://github.com/casadi/casadi/wiki/L_1h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L258

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L236-L242

";

%feature("docstring") casadi::OptiAdvanced::nx "

Number of (scalarised) decision variables.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L261

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L244-L250

";

%feature("docstring") casadi::OptiAdvanced::np "

Number of (scalarised) parameters.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L264

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L252-L258

";

%feature("docstring") casadi::OptiAdvanced::ng "

Number of (scalarised) constraints.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L267

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L260-L266

";

%feature("docstring") casadi::OptiAdvanced::x "

Get all (scalarised) decision variables as a symbolic column vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L270

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L268-L274

";

%feature("docstring") casadi::OptiAdvanced::p "

Get all (scalarised) parameters as a symbolic column vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L273

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L276-L282

";

%feature("docstring") casadi::OptiAdvanced::g "

Get all (scalarised) constraint expressions as a column vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L276

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L284-L290

";

%feature("docstring") casadi::OptiAdvanced::f "

Get objective expression.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L279

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L292-L298

";

%feature("docstring") casadi::OptiAdvanced::lbg "

Get all (scalarised) bounds on constraints as a column vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L282

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L300-L306

";

%feature("docstring") casadi::OptiAdvanced::lam_g "

Get all (scalarised) dual variables as a symbolic column vector.

Useful for obtaining the Lagrange Hessian:

::

  * sol.value(hessian(opti.f+opti.lam_g'*opti.g,opti.x)) % MATLAB
  * sol.value(hessian(opti.f+dot(opti.lam_g,opti.g),opti.x)[0]) # Python
  * 



Extra doc: https://github.com/casadi/casadi/wiki/L_1i

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L294

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L317-L323

";

%feature("docstring") casadi::OptiAdvanced::to_function "

";

";

%feature("docstring") casadi::OptiAdvanced::debug "

Get a copy with advanced functionality.

You get access to more methods, but you have no guarantees about API 

stability

The copy is effectively a deep copy: Updating the state of the copy 
does 
not update the original.

Extra doc: https://github.com/casadi/casadi/wiki/L_1l

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L334

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L741-L743

";

%feature("docstring") casadi::OptiAdvanced::advanced "

Get a copy with advanced functionality.

You get access to more methods, but you have no guarantees about API 

stability

The copy is effectively a deep copy: Updating the state of the copy 
does 
not update the original.

Extra doc: https://github.com/casadi/casadi/wiki/L_1m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L344

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L744-L746

";

%feature("docstring") casadi::OptiAdvanced::copy "

Get a copy of the.

The copy is effectively a deep copy: Updating the state of the copy 
does 
not update the original.

Extra doc: https://github.com/casadi/casadi/wiki/L_1n

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L352

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L747-L749

";

%feature("docstring") casadi::OptiAdvanced::update_user_dict "

";

";

%feature("docstring") casadi::OptiAdvanced::user_dict "

Get user data.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L363

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L400-L406

";

%feature("docstring") casadi::OptiAdvanced::type_name "

Readable name of the class.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L366

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L366-L366

";

%feature("docstring") casadi::OptiAdvanced::disp "

Print representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L369

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L644-L664

";

%feature("docstring") casadi::OptiAdvanced::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L372

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L666-L670

";

%feature("docstring") casadi::OptiAdvanced::class_name "

Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L100-L102

";

%feature("docstring") casadi::OptiAdvanced::is_null "

Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L74-L76

";

%feature("docstring") casadi::OptiAdvanced::__hash__ "

Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L130-L132

";


// File: classcasadi_1_1OptiCallback.xml
%feature("docstring") casadi::OptiCallback "

C++ includes: optistack.hpp
";


// File: classcasadi_1_1OptiSol.xml
%feature("docstring") casadi::OptiSol "

A simplified interface for NLP modeling/solving.

This class offers a view with solution retrieval facilities The API is
 
guaranteed to be stable.

Joris Gillis, Erik Lambrechts

Extra doc: https://github.com/casadi/casadi/wiki/L_1v

C++ includes: optistack.hpp
";

%feature("docstring") casadi::OptiSol::value "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L594

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L768-L770

>  DM casadi::OptiSol::value(const SX &x, const std::vector< MX > &values=std::vector< MX >()) const
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L594

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L768-L770

";

";

%feature("docstring") casadi::OptiSol::value_variables "

get assignment expressions for the optimal solution

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L598

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L772-L774

";

%feature("docstring") casadi::OptiSol::stats "

Get statistics.

nlpsol stats are passed as-is. No stability can be guaranteed about 
this 
part of the API

Extra doc: https://github.com/casadi/casadi/wiki/L_1w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.hpp#L607

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/optistack.cpp#L780-L782

";


// File: classcasadi_1_1Output.xml


// File: classcasadi_1_1ParsedFile.xml
%feature("docstring") casadi::ParsedFile "

A parsed file.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_1153

C++ includes: casadi_file.hpp
";

%feature("docstring") casadi::ParsedFile::ParsedFile "

Construct from a file.

Extra doc: https://github.com/casadi/casadi/wiki/L_1156

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L62

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.cpp#L37-L39

>  casadi::ParsedFile::ParsedFile(const std::vector< std::string > &lines, int offset=0)
------------------------------------------------------------------------

Construct from a file.

Extra doc: https://github.com/casadi/casadi/wiki/L_1156

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L62

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.cpp#L37-L39

";

";

%feature("docstring") casadi::ParsedFile::parse "

Parse a list of strings.

Extra doc: https://github.com/casadi/casadi/wiki/L_1158

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L72

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.cpp#L63-L107

>  void casadi::ParsedFile::parse(const std::vector< std::string > &lines, int offset)
------------------------------------------------------------------------

Parse a list of strings.

Extra doc: https://github.com/casadi/casadi/wiki/L_1158

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L72

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.cpp#L63-L107

";

";

%feature("docstring") casadi::ParsedFile::print "

Print parsed file.

Extra doc: https://github.com/casadi/casadi/wiki/L_1159

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L77

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.cpp#L109-L115

";

%feature("docstring") casadi::ParsedFile::has "

Does an entry exist?

Extra doc: https://github.com/casadi/casadi/wiki/L_1160

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L82

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.cpp#L123-L126

";

%feature("docstring") casadi::ParsedFile::to_text "

Get entry as a text.

Extra doc: https://github.com/casadi/casadi/wiki/L_1161

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L87

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.cpp#L117-L121

";

%feature("docstring") casadi::ParsedFile::to "

Convert to a type.

Extra doc: https://github.com/casadi/casadi/wiki/L_1162

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L100

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L100-L105

";

%feature("docstring") casadi::ParsedFile::to_string "

Get entry as a string.

Extra doc: https://github.com/casadi/casadi/wiki/L_1163

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L110

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L110-L112

";

%feature("docstring") casadi::ParsedFile::to_vector "

Get entry as a vector.

Extra doc: https://github.com/casadi/casadi/wiki/L_1164

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L118

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L118-L124

";

%feature("docstring") casadi::ParsedFile::to_set "

Get entry as a set.

Extra doc: https://github.com/casadi/casadi/wiki/L_1165

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L130

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L130-L134

";

%feature("docstring") casadi::ParsedFile::to_int "

Get entry as an integer.

Extra doc: https://github.com/casadi/casadi/wiki/L_1166

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L139

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_file.hpp#L139-L139

";


// File: classcasadi_1_1Polynomial.xml
%feature("docstring") casadi::Polynomial "

Helper class for differentiating and integrating polynomials.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_8y

C++ includes: polynomial.hpp
";

%feature("docstring") casadi::Polynomial::Polynomial "

Construct from a vector of polynomial coefficients.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L56

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L56-L56

>  casadi::Polynomial::Polynomial(const std::vector< T > &coeff)
------------------------------------------------------------------------

Construct from a vector of polynomial coefficients.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L56

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L56-L56

";

";

%feature("docstring") casadi::Polynomial::coeff "

Coefficients of the polynomial.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L71

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L71-L71

";

%feature("docstring") casadi::Polynomial::degree "

Degree of the polynomial.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L74

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.cpp#L73-L75

";

%feature("docstring") casadi::Polynomial::scalar "

Get scalar value (error if  degree()!=0)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L77

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.cpp#L77-L80

";

%feature("docstring") casadi::Polynomial::derivative "

Create a new polynomial for the derivative.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L80

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.cpp#L142-L148

";

%feature("docstring") casadi::Polynomial::anti_derivative "

Create a new polynomial for the anti-derivative (primitive function)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L83

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.cpp#L150-L157

";

%feature("docstring") casadi::Polynomial::trim "

Remove excess zeros.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L86

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.cpp#L134-L140

";

%feature("docstring") casadi::Polynomial::type_name "

Readable name of the class.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L89

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L89-L89

";

%feature("docstring") casadi::Polynomial::disp "

Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.hpp#L92

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/polynomial.cpp#L57-L71

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


// File: classcasadi_1_1PrintableObject.xml
%feature("docstring") casadi::PrintableObject "

Base class for objects that have a natural string representation.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_967

C++ includes: printable_object.hpp
";

%feature("docstring") casadi::PrintableObject::getDescription "

Return a string with a description (for SWIG)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/printable_object.hpp#L51

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/printable_object.hpp#L51-L55

";

%feature("docstring") casadi::PrintableObject::getRepresentation "

Return a string with a representation (for SWIG)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/printable_object.hpp#L58

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/printable_object.hpp#L58-L62

";

%feature("docstring") casadi::PrintableObject::str "

[INTERNAL] 
Return a string with a description of the object, cf. 
str(Object) in 
Python.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/printable_object.hpp#L79

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/printable_object.hpp#L79-L81

";

%feature("docstring") casadi::PrintableObject::repr "

[INTERNAL] 
Return a string with a representation of the object, cf. 
repr(Object) 
in Python.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/printable_object.hpp#L84

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/printable_object.hpp#L84-L86

";


// File: classcasadi_1_1QpToNlp.xml
%feature("docstring") casadi::QpToNlp "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Qrqp.xml
%feature("docstring") casadi::Qrqp "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Qrsqp.xml
%feature("docstring") casadi::Qrsqp "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Rootfinder.xml
%feature("docstring") casadi::Rootfinder "

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

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1SerializerBase.xml
%feature("docstring") casadi::SerializerBase "

C++ includes: serializer.hpp
";

%feature("docstring") casadi::SerializerBase::SerializerBase "

[INTERNAL] ";


// File: classcasadi_1_1SerializingStream.xml
%feature("docstring") casadi::SerializingStream "

Helper class for Serialization.

Joris Gillis

Extra doc: https://github.com/casadi/casadi/wiki/L_ao

C++ includes: serializing_stream.hpp
";

%feature("docstring") casadi::SerializingStream::SerializingStream "

";

";

%feature("docstring") casadi::SerializingStream::pack "

";

";


// File: classcasadi_1_1SharedObject.xml
%feature("docstring") casadi::SharedObject "

SharedObject implements a reference counting framework similar for efficient
 and.

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
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L38-L41

>  casadi::SharedObject::SharedObject(const SharedObject &ref)
------------------------------------------------------------------------
[INTERNAL] 
Copy constructor (shallow copy)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L96

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L38-L41

";

";

%feature("docstring") casadi::SharedObject::~SharedObject "

[INTERNAL] 
Destructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L99

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L43-L45

";

%feature("docstring") casadi::SharedObject::class_name "

Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L100-L102

";

%feature("docstring") casadi::SharedObject::disp "

Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L104-L110

";

%feature("docstring") casadi::SharedObject::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::SharedObject::is_null "

Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L74-L76

";

%feature("docstring") casadi::SharedObject::__hash__ "

Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L130-L132

";


// File: classcasadi_1_1ShellCompiler.xml
%feature("docstring") casadi::ShellCompiler "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Simulator.xml
%feature("docstring") casadi::Simulator "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Slice.xml
%feature("docstring") casadi::Slice "

Class representing a  Slice.

Note that Python or Octave do not need to use this class. They can 
just use
 slicing utility from the host language ( M[0:6] in Python, 
M(1:7) )

Extra doc: https://github.com/casadi/casadi/wiki/L_13

C++ includes: slice.hpp
";

%feature("docstring") casadi::Slice::Slice "

";

";

%feature("docstring") casadi::Slice::all "

Get a vector of indices (nested slice)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L75

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L130-L138

>  std::vector< casadi_int > casadi::Slice::all(const Slice &outer, casadi_int len) const
------------------------------------------------------------------------

Get a vector of indices (nested slice)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L75

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L130-L138

";

";

%feature("docstring") casadi::Slice::size "

Get number of elements.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L78

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L105-L110

";

%feature("docstring") casadi::Slice::is_empty "

Check if slice is empty.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L81

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L112-L114

";

%feature("docstring") casadi::Slice::is_scalar "

Is the slice a scalar.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L84

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L140-L145

";

%feature("docstring") casadi::Slice::scalar "

Get scalar (if is_scalar)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L87

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L147-L151

";

%feature("docstring") casadi::Slice::apply "

Apply concrete length.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L98

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L67-L89

";

%feature("docstring") casadi::Slice::type_name "

Get name of the class.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L107

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L107-L107

";

%feature("docstring") casadi::Slice::disp "

Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L110

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L116-L128

";

%feature("docstring") casadi::Slice::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L113

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L113-L117

";

%feature("docstring") casadi::Slice::info "

Obtain information

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L120

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L120-L122

";

%feature("docstring") casadi::Slice::serialize "

Serialize an object.

Extra doc: https://github.com/casadi/casadi/wiki/L_14

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L127

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L293-L297

";


// File: classcasadi_1_1SlicotDple.xml
%feature("docstring") casadi::SlicotDple "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L458

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L458-L461

>  MatType casadi::SparsityInterface::horzcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w, const MatType &v, const MatType &u)
------------------------------------------------------------------------
[INTERNAL] 
Concatenate horizontally, six matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_4e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L458

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L458-L461

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::vertcat "

[INTERNAL] 
Concatenate vertically, six matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_4j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L496

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L496-L499

>  MatType casadi::SparsityInterface::vertcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w, const MatType &v, const MatType &u)
------------------------------------------------------------------------
[INTERNAL] 
Concatenate vertically, six matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_4j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L496

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L496-L499

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::horzsplit "

[INTERNAL] 
split horizontally, retaining fixed-sized groups of columns

Parameters:
-----------

incr: 
Size of each group of columns

horzcat(horzsplit(x, ...)) = x

Extra doc: https://github.com/casadi/casadi/wiki/L_3h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L130

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L130-L132

>  std::vector<MatType > casadi::SparsityInterface::horzsplit(const MatType &x, casadi_int incr=1)
------------------------------------------------------------------------
[INTERNAL] 
split horizontally, retaining fixed-sized groups of columns

Parameters:
-----------

incr: 
Size of each group of columns

horzcat(horzsplit(x, ...)) = x

Extra doc: https://github.com/casadi/casadi/wiki/L_3h

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L130

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L130-L132

";

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

::

  [SX(a_0), SX(a_1), SX(a_2), SX(a_3)]
  





::

>>> print vertsplit(SX.sym(\"a\",4),2)

::

  [SX([a_0, a_1]), SX([a_2, a_3])]
  



If the number of rows is not a multiple of  incr, the last entry returned 
will have a size smaller than  incr.



::

>>> print vertsplit(DM([0,1,2,3,4]),2)

::

  [DM([0, 1]), DM([2, 3]), DM(4)]
  



Extra doc: https://github.com/casadi/casadi/wiki/L_3k

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L182

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L182-L184

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

::

  [SX(a_0), SX(a_1), SX(a_2), SX(a_3)]
  





::

>>> print vertsplit(SX.sym(\"a\",4),2)

::

  [SX([a_0, a_1]), SX([a_2, a_3])]
  



If the number of rows is not a multiple of  incr, the last entry returned 
will have a size smaller than  incr.



::

>>> print vertsplit(DM([0,1,2,3,4]),2)

::

  [DM([0, 1]), DM([2, 3]), DM(4)]
  



Extra doc: https://github.com/casadi/casadi/wiki/L_3k

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L182

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L182-L184

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::offset "

[INTERNAL] 
Helper function, get offsets corresponding to a vector of 
matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_3j

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L150-L152

";

%feature("docstring") casadi::SparsityInterfaceCommon::blockcat "

[INTERNAL] 
Construct a matrix from 4 blocks.

Extra doc: https://github.com/casadi/casadi/wiki/L_3m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L197

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L197-L199

>  MatType casadi::SparsityInterface::blockcat(const MatType &A, const MatType &B, const MatType &C, const MatType &D)
------------------------------------------------------------------------
[INTERNAL] 
Construct a matrix from 4 blocks.

Extra doc: https://github.com/casadi/casadi/wiki/L_3m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L197

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L197-L199

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
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L225

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L225-L227

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
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L225

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L225-L227

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::diagcat "

[INTERNAL] 
Concatenate along diagonal, six matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_4o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L534

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L534-L537

>  MatType casadi::SparsityInterface::diagcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w, const MatType &v, const MatType &u)
------------------------------------------------------------------------
[INTERNAL] 
Concatenate along diagonal, six matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_4o

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L534

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L534-L537

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
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L287

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L287-L289

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
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L287

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L287-L289

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::veccat "

[INTERNAL] 
concatenate vertically while vectorizing all arguments with vec

Extra doc: https://github.com/casadi/casadi/wiki/L_3u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L294

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L294-L296

";

%feature("docstring") casadi::SparsityInterfaceCommon::mtimes "

[INTERNAL] 
 Matrix product of n matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_3w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L308

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L308-L310

>  MatType casadi::SparsityInterface::mtimes(const std::vector< MatType > &args)
------------------------------------------------------------------------
[INTERNAL] 
 Matrix product of n matrices.

Extra doc: https://github.com/casadi/casadi/wiki/L_3w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L308

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L308-L310

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
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L321

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L321-L323

";

%feature("docstring") casadi::SparsityInterfaceCommon::transpose "

[INTERNAL] 
Transpose.

Extra doc: https://github.com/casadi/casadi/wiki/L_3y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L328

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L328-L330

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
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L349

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L349-L351

";

%feature("docstring") casadi::SparsityInterfaceCommon::reshape "

[INTERNAL] 
Reshape the matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_42

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L370

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L370-L372

>  MatType casadi::SparsityInterface::reshape(const MatType &x, const Sparsity &sp)
------------------------------------------------------------------------
[INTERNAL] 
Reshape the matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_42

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L370

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L370-L372

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::sprank "

[INTERNAL] 
Obtain the structural rank of a sparsity-pattern.

Extra doc: https://github.com/casadi/casadi/wiki/L_43

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L377

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L377-L379

";

%feature("docstring") casadi::SparsityInterfaceCommon::norm_0_mul "

[INTERNAL] 
0-norm (nonzero count) of a Matrix-matrix product

Extra doc: https://github.com/casadi/casadi/wiki/L_44

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L384

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L384-L386

";

%feature("docstring") casadi::SparsityInterfaceCommon::triu "

[INTERNAL] 
Get the upper triangular part of a matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_45

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L391

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L391-L393

";

%feature("docstring") casadi::SparsityInterfaceCommon::tril "

[INTERNAL] 
Get the lower triangular part of a matrix.

Extra doc: https://github.com/casadi/casadi/wiki/L_46

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L398

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L398-L400

";

%feature("docstring") casadi::SparsityInterfaceCommon::kron "

[INTERNAL] 
Kronecker tensor product.

Creates a block matrix in which each element (i, j) is a_ij*b

Extra doc: https://github.com/casadi/casadi/wiki/L_47

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L407

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L407-L409

";

%feature("docstring") casadi::SparsityInterfaceCommon::repmat "

[INTERNAL] 
Repeat matrix A n times vertically and m times horizontally.

Extra doc: https://github.com/casadi/casadi/wiki/L_49

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L421

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L421-L423

>  MatType casadi::SparsityInterface::repmat(const MatType &A, const std::pair< casadi_int, casadi_int > &rc)
------------------------------------------------------------------------
[INTERNAL] 
Repeat matrix A n times vertically and m times horizontally.

Extra doc: https://github.com/casadi/casadi/wiki/L_49

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L421

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L421-L423

";

";

%feature("docstring") casadi::SparsityInterfaceCommon::sum1 "

[INTERNAL] 
Return a row-wise summation of elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_4p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L542

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L542-L542

";

%feature("docstring") casadi::SparsityInterfaceCommon::sum2 "

[INTERNAL] 
Return a column-wise summation of elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_4q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L547

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sparsity_interface.hpp#L547-L547

";


// File: classcasadi_1_1Sqpmethod.xml
%feature("docstring") casadi::Sqpmethod "

Diagrams
--------



C++ includes: e0_diagram.hpp
";


// File: classcasadi_1_1Logger_1_1Stream.xml
%feature("docstring") casadi::Logger::Stream "

C++ includes: casadi_logger.hpp
";


// File: classcasadi_1_1Logger_1_1Streambuf.xml
%feature("docstring") casadi::Logger::Streambuf "

C++ includes: casadi_logger.hpp
";


// File: classcasadi_1_1StringDeserializer.xml
%feature("docstring") casadi::StringDeserializer "

C++ includes: serializer.hpp
";

%feature("docstring") casadi::StringDeserializer::StringDeserializer "

Advanced deserialization of CasADi objects.

See: 
 StringDeserializer

Extra doc: https://github.com/casadi/casadi/wiki/L_7r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.hpp#L233

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.cpp#L115-L118

";

%feature("docstring") casadi::StringDeserializer::decode "

Sets the string to deserialize objects from.

Extra doc: https://github.com/casadi/casadi/wiki/L_7s

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.hpp#L240

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.cpp#L92-L97

";


// File: classcasadi_1_1StringSerializer.xml
%feature("docstring") casadi::StringSerializer "

C++ includes: serializer.hpp
";

%feature("docstring") casadi::StringSerializer::StringSerializer "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.cpp#L38-L40

";

%feature("docstring") casadi::StringSerializer::encode "

Returns a string that holds the serialized objects.

As a side effect, this method clears the internal buffer

Extra doc: https://github.com/casadi/casadi/wiki/L_7p

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.hpp#L211

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/serializer.cpp#L86-L91

";


// File: classcasadi_1_1SubIndex.xml
%feature("docstring") casadi::SubIndex "

SubIndex class for  Matrix Same as the above class but for single argument 
return for operator()
 
Joel Andersson

C++ includes: submatrix.hpp
";

%feature("docstring") casadi::SubIndex::SubIndex "

Default copy constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/submatrix.hpp#L113

>  casadi::SubIndex< M, I >::SubIndex(const SubIndex< M, I > &y)=default
------------------------------------------------------------------------

Default copy constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/submatrix.hpp#L113

";

";


// File: classcasadi_1_1SubMatrix.xml
%feature("docstring") casadi::SubMatrix "

SubMatrix class for  Matrix SubMatrix is the return type for operator() of 
the  Matrix class, it allows access to the value as well as changing the 
parent 
object 
Joel Andersson

C++ includes: submatrix.hpp
";

%feature("docstring") casadi::SubMatrix::SubMatrix "

Default copy constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/submatrix.hpp#L53

>  casadi::SubMatrix< M, I, J >::SubMatrix(const SubMatrix< M, I, J > &y)=default
------------------------------------------------------------------------

Default copy constructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/submatrix.hpp#L53

";

";


// File: classcasadi_1_1SundialsSimulator.xml
%feature("docstring") casadi::SundialsSimulator "

Diagrams
--------



C++ includes: e0_diagram.hpp
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
Remainder after division: (x,y) -> mod(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pq

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L595

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L595-L597

";

%feature("docstring") casadi::SXElem::atan2 "

[INTERNAL] 
Two argument arc tangent: (x,y) -> atan2(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pr

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L610

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L610-L612

";

%feature("docstring") casadi::SXElem::if_else_zero "

[INTERNAL] 
Conditional assignment: (x,y) -> x ? y : 0.

Extra doc: https://github.com/casadi/casadi/wiki/L_ps

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L622

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L622-L624

";

%feature("docstring") casadi::SXElem::fmin "

[INTERNAL] 
Smallest of two values: (x,y) -> min(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L634

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L634-L636

";

%feature("docstring") casadi::SXElem::fmax "

[INTERNAL] 
Largest of two values: (x,y) -> max(x,y)

Extra doc: https://github.com/casadi/casadi/wiki/L_pu

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L646

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L646-L648

";

%feature("docstring") casadi::SXElem::copysign "

[INTERNAL] 
Copy sign

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L672

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L672-L674

";

%feature("docstring") casadi::SXElem::constpow "

[INTERNAL] 
Elementwise power with const power

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L682

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L682-L684

";

%feature("docstring") casadi::SXElem::printme "

[INTERNAL] 
Debug printing

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L692

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L692-L694

";

%feature("docstring") casadi::SXElem::hypot "

[INTERNAL] 
Precision variant for 2 norm: (x,y) -> sqrt(x^2+y^2)

Extra doc: https://github.com/casadi/casadi/wiki/L_pw

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L704

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_expression.hpp#L704-L706

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
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L60-L63

>  casadi::SXElem::SXElem(const SXElem &scalar)
------------------------------------------------------------------------
[INTERNAL] 
Copy constructor.

Extra doc: https://github.com/casadi/casadi/wiki/L_10m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L107

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L60-L63

";

";

%feature("docstring") casadi::SXElem::~SXElem "

[INTERNAL] 
Destructor.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L110

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L88-L90

";

%feature("docstring") casadi::SXElem::disp "

[INTERNAL] 
Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L129

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L134-L136

";

%feature("docstring") casadi::SXElem::__nonzero__ "

[INTERNAL] 
Check the truth value of this node.

Introduced to catch bool(x) situations in python

Extra doc: https://github.com/casadi/casadi/wiki/L_10q

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L156

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L151-L154

";

%feature("docstring") casadi::SXElem::is_leaf "

[INTERNAL] 
check if this  SXElem is a leaf of the SX graph

An  SXElem qualifies as leaf when it has no dependencies.

Extra doc: https://github.com/casadi/casadi/wiki/L_10r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L163

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L436-L439

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
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L441-L444

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
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L612-L618

";

%feature("docstring") casadi::SXElem::is_nonnegative "

[INTERNAL] 
Check if a value is always nonnegative (false negatives are 
allowed)

Extra doc: https://github.com/casadi/casadi/wiki/L_10t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L188

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L509-L517

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
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L156-L158

";

%feature("docstring") casadi::SXElem::n_dep "

[INTERNAL] 
Get the number of dependencies of a binary  SXElem.

Extra doc: https://github.com/casadi/casadi/wiki/L_10v

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L205

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L532-L534

";

%feature("docstring") casadi::SXElem::__hash__ "

[INTERNAL] 
Returns a number that is unique for a given SXNode.

If the  SXElem does not point to any node, 0 is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_10w

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L212

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L536-L538

";

%feature("docstring") casadi::SXElem::inv "

[INTERNAL] 
Element-wise inverse.

Extra doc: https://github.com/casadi/casadi/wiki/L_10y

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.hpp#L222

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L160-L166

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
https://github.com/casadi/casadi/blob/develop/casadi/core/sx_elem.cpp#L620-L622

";


// File: classcasadi_1_1SymbolicQr.xml
%feature("docstring") casadi::SymbolicQr "

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

Weak reference type.

A weak reference to a  SharedObject
Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_ax

C++ includes: shared_object.hpp
";

%feature("docstring") casadi::WeakRef::WeakRef "

Construct from a shared object (also implicit type conversion)

Extra doc: https://github.com/casadi/casadi/wiki/L_az

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L193

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L158-L160

>  casadi::WeakRef::WeakRef(SharedObject shared)
------------------------------------------------------------------------

Construct from a shared object (also implicit type conversion)

Extra doc: https://github.com/casadi/casadi/wiki/L_az

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L193

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L158-L160

";

";

%feature("docstring") casadi::WeakRef::shared "

Get a shared (owning) reference.

Extra doc: https://github.com/casadi/casadi/wiki/L_b0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L198

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L142-L148

";

%feature("docstring") casadi::WeakRef::alive "

Check if alive.

Extra doc: https://github.com/casadi/casadi/wiki/L_b1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L203

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L138-L140

";

%feature("docstring") casadi::WeakRef::class_name "

Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L100-L102

";

%feature("docstring") casadi::WeakRef::disp "

Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L104-L110

";

%feature("docstring") casadi::WeakRef::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::WeakRef::is_null "

Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L74-L76

";

%feature("docstring") casadi::WeakRef::__hash__ "

Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L130-L132

";


// File: classcasadi_1_1XmlFile.xml
%feature("docstring") casadi::XmlFile "

XML parser.

Can be used for parsing XML files into CasADi data structures.

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_7k

C++ includes: xml_file.hpp
";

%feature("docstring") casadi::XmlFile::parse "

[INTERNAL] ";

%feature("docstring") casadi::XmlFile::class_name "

Get class name.

Extra doc: https://github.com/casadi/casadi/wiki/L_au

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L100-L102

";

%feature("docstring") casadi::XmlFile::disp "

Print a description of the object.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L135

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L104-L110

";

%feature("docstring") casadi::XmlFile::get_str "

Get string representation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L138-L142

";

%feature("docstring") casadi::XmlFile::is_null "

Is a null pointer?

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L150

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L74-L76

";

%feature("docstring") casadi::XmlFile::__hash__ "

Returns a number that is unique for a given Node.

If the Object does not point to any node, \"0\" is returned.

Extra doc: https://github.com/casadi/casadi/wiki/L_av

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.hpp#L157

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/shared_object.cpp#L130-L132

";


// File: namespace_0d371.xml


// File: namespacecasadi.xml
%feature("docstring") casadi::IndexRecution::collocation_points "
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
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L128

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L128-L130

";

%feature("docstring") casadi::IndexRecution::dae_reduce_index "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L1065

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L1065-L1067

>  SXDict casadi::dae_reduce_index(const SXDict &dae, Dict &stats, const Dict &opts)
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L1065

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L1065-L1067

";

";

%feature("docstring") casadi::IndexRecution::dae_map_semi_expl "

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

Obtain a generator  Function for producing consistent initial guesses of a 
reduced DAE.

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

Obtain a generator  Function for producing consistent initial guesses of a 
reduced DAE.

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

%feature("docstring") casadi::IndexRecution::matrixName "

Get typename.

Extra doc: https://github.com/casadi/casadi/wiki/L_18d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L56

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L56-L57

";

%feature("docstring") casadi::IndexRecution::matrixName< double > "

Get typename.

Extra doc: https://github.com/casadi/casadi/wiki/L_18d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L58

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L58-L58

";

%feature("docstring") casadi::IndexRecution::matrixName< casadi_int > "

Get typename.

Extra doc: https://github.com/casadi/casadi/wiki/L_18d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L59

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/matrix_decl.hpp#L59-L59

";

%feature("docstring") casadi::IndexRecution::nlpsol_default_in "

Default input for an NLP solver.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L153-L157

>  std::vector< double > casadi::nlpsol_default_in()
------------------------------------------------------------------------

Default input for an NLP solver.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L153-L157

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

%feature("docstring") casadi::IndexRecution::to_int "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::all "

[INTERNAL] 
Check if all arguments are true.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L78

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L78-L83

";

%feature("docstring") casadi::IndexRecution::any "

[INTERNAL] 
Check if any arguments are true.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L85

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L85-L90

";

%feature("docstring") casadi::IndexRecution::is_range "

[INTERNAL] 
Check if a vector matches a range.

Extra doc: https://github.com/casadi/casadi/wiki/L_1l7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L92

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L92-L102

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
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L134

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L134-L136

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
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L134

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L134-L136

";

";

%feature("docstring") casadi::IndexRecution::is_equally_spaced "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::complement "

Returns the list of all i in [0, size[ not found in supplied list.

The supplied vector may contain duplicates and may be non-monotonous 
The 
supplied vector will be checked for bounds The result vector is 
guaranteed 
to be monotonously increasing

Extra doc: https://github.com/casadi/casadi/wiki/L_1lf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L138

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L138-L154

";

%feature("docstring") casadi::IndexRecution::lookupvector "

";

";

%feature("docstring") casadi::IndexRecution::tensor_permute_mapping "

[INTERNAL] 
Computes a mapping for a (dense) tensor permutation.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L173

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L173-L238

";

%feature("docstring") casadi::IndexRecution::join "

[INTERNAL] 
Join three lists.

Extra doc: https://github.com/casadi/casadi/wiki/L_1lc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L529

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L529-L534

>  std::vector< T > casadi::join(const std::vector< T > &a, const std::vector< T > &b, const std::vector< T > &c)
------------------------------------------------------------------------
[INTERNAL] 
Join three lists.

Extra doc: https://github.com/casadi/casadi/wiki/L_1lc

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L529

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L529-L534

";

";

%feature("docstring") casadi::IndexRecution::startswith "

[INTERNAL] 
Checks if s starts with p.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L266

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L266-L272

";

%feature("docstring") casadi::IndexRecution::boolvec_not "

[INTERNAL] 
Invert all entries.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L345

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L345-L350

";

%feature("docstring") casadi::IndexRecution::boolvec_and "

[INTERNAL] 
And operation on boolean vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L352

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L352-L358

";

%feature("docstring") casadi::IndexRecution::boolvec_or "

[INTERNAL] 
Or operation on boolean vector.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L360

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.cpp#L360-L366

";

%feature("docstring") casadi::IndexRecution::boolvec_to_index "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L498

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L498-L512

";

%feature("docstring") casadi::IndexRecution::reverse "

[INTERNAL] 
Reverse a list.

Extra doc: https://github.com/casadi/casadi/wiki/L_1la

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L515

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L515-L519

";

%feature("docstring") casadi::IndexRecution::permute "

[INTERNAL] 
permute a list

Extra doc: https://github.com/casadi/casadi/wiki/L_1ld

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L537

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L537-L544

";

%feature("docstring") casadi::IndexRecution::find "

[INTERNAL] 
find nonzeros

Extra doc: https://github.com/casadi/casadi/wiki/L_1le

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L547

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L547-L553

";

%feature("docstring") casadi::IndexRecution::in_range "

Check if for each element of v holds: lower <= v_i < upper.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L601

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L601-L607

>  bool casadi::in_range(const std::vector< T > &v, casadi_int lower, casadi_int upper)
------------------------------------------------------------------------

Check if for each element of v holds: lower <= v_i < upper.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L601

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L601-L607

";

";

%feature("docstring") casadi::IndexRecution::flatten_nested_vector "

Flatten a nested std::vector tot a single flattened vector.

Contents of nested[i] ends up in 
flat[indices[i]]..flat[indices[i+1]-1]

Extra doc: https://github.com/casadi/casadi/wiki/L_1li

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L627

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L627-L642

>  void casadi::flatten_nested_vector(const std::vector< std::vector< T > > &nested, std::vector< S > &flat, std::vector< I > &indices)
------------------------------------------------------------------------

Flatten a nested std::vector tot a single flattened vector.

Contents of nested[i] ends up in 
flat[indices[i]]..flat[indices[i+1]-1]

Extra doc: https://github.com/casadi/casadi/wiki/L_1li

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L627

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L627-L642

";

";

%feature("docstring") casadi::IndexRecution::is_increasing "

Check if the vector is strictly increasing.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L651

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L651-L659

";

%feature("docstring") casadi::IndexRecution::is_decreasing "

Check if the vector is strictly decreasing.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L662

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L662-L670

";

%feature("docstring") casadi::IndexRecution::is_nonincreasing "

Check if the vector is non-increasing.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L673

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L673-L681

";

%feature("docstring") casadi::IndexRecution::is_nondecreasing "

Check if the vector is non-decreasing.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L684

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L684-L692

";

%feature("docstring") casadi::IndexRecution::is_monotone "

Check if the vector is monotone.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L695

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L695-L697

";

%feature("docstring") casadi::IndexRecution::is_strictly_monotone "

Check if the vector is strictly monotone.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L700

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L700-L702

";

%feature("docstring") casadi::IndexRecution::has_negative "

Check if the vector has negative entries.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L705

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L705-L710

";

%feature("docstring") casadi::IndexRecution::write_matlab "

Print matrix, matlab style.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L718

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L718-L723

>  void casadi::write_matlab(std::ostream &stream, const std::vector< std::vector< T > > &v)
------------------------------------------------------------------------

Print matrix, matlab style.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L718

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L718-L723

";

";

%feature("docstring") casadi::IndexRecution::read_matlab "

Read matrix, matlab style.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L746

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L746-L769

>  void casadi::read_matlab(std::ifstream &file, std::vector< std::vector< T > > &v)
------------------------------------------------------------------------

Read matrix, matlab style.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L746

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L746-L769

";

";

%feature("docstring") casadi::IndexRecution::linspace "

[INTERNAL] 
Matlab's linspace.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L772

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L772-L783

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
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L810

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L810-L835

";

%feature("docstring") casadi::IndexRecution::product "

[INTERNAL] 
product

Extra doc: https://github.com/casadi/casadi/wiki/L_1lk

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L838

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L838-L842

";

%feature("docstring") casadi::IndexRecution::sum "

[INTERNAL] 
sum

Extra doc: https://github.com/casadi/casadi/wiki/L_1ll

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L845

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L845-L849

";

%feature("docstring") casadi::IndexRecution::cumsum "

[INTERNAL] 
cumulative sum

Extra doc: https://github.com/casadi/casadi/wiki/L_1lm

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L852

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L852-L860

";

%feature("docstring") casadi::IndexRecution::diff "

[INTERNAL] 
diff

Extra doc: https://github.com/casadi/casadi/wiki/L_1ln

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L874

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L874-L881

";

%feature("docstring") casadi::IndexRecution::cumsum0 "

[INTERNAL] 
cumulative sum, starting with zero

Extra doc: https://github.com/casadi/casadi/wiki/L_1lo

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L863

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L863-L871

";

%feature("docstring") casadi::IndexRecution::is_regular "

Checks if array does not contain NaN or Inf.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L378

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/casadi_misc.hpp#L378-L384

";

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

%feature("docstring") casadi::IndexRecution::conic_debug "

Generate native code in the interfaced language for debugging

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L55

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L55-L59

>  void casadi::conic_debug(const Function &f, std::ostream &file)
------------------------------------------------------------------------

Generate native code in the interfaced language for debugging

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L55

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L55-L59

";

";

%feature("docstring") casadi::IndexRecution::conic_in "

Get QP solver input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1eg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L73

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L73-L90

>  std::string casadi::conic_in(casadi_int ind)
------------------------------------------------------------------------

Get QP solver input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1eg

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L73

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L73-L90

";

";

%feature("docstring") casadi::IndexRecution::conic_out "

Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1eh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L92

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L92-L101

>  std::string casadi::conic_out(casadi_int ind)
------------------------------------------------------------------------

Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1eh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L92

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L92-L101

";

";

%feature("docstring") casadi::IndexRecution::conic_n_in "

Get the number of QP solver inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ei

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L103

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L103-L105

";

%feature("docstring") casadi::IndexRecution::conic_n_out "

Get the number of QP solver outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ej

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L107

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L107-L109

";

%feature("docstring") casadi::IndexRecution::conic_options "

Get all options for a plugin.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ek

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L529

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L529-L531

";

%feature("docstring") casadi::IndexRecution::conic_option_type "

Get type info for a particular option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1el

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L533

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L533-L535

";

%feature("docstring") casadi::IndexRecution::conic_option_info "

Get documentation for a particular option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1em

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L537

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L537-L539

";

%feature("docstring") casadi::IndexRecution::conic "



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
|                  |                 | ght*nf<=(1-ad_we |                  |
|                  |                 | ight)*na is used |                  |
|                  |                 | where nf and na  |                  |
|                  |                 | are estimates of |                  |
|                  |                 | the number of    |                  |
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
| error_on_fail    | OT_BOOL         | When the         | casadi::Conic    |
|                  |                 | numerical        |                  |
|                  |                 | process returns  |                  |
|                  |                 | unsuccessfully,  |                  |
|                  |                 | raise an error   |                  |
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
| print_time       | OT_BOOL         | print            | casadi::Function |
|                  |                 | information      | Internal         |
|                  |                 | about execution  |                  |
|                  |                 | time. Implies    |                  |
|                  |                 | record_time.     |                  |
+------------------+-----------------+------------------+------------------+
| record_time      | OT_BOOL         | record           | casadi::Function |
|                  |                 | information      | Internal         |
|                  |                 | about execution  |                  |
|                  |                 | time, for        |                  |
|                  |                 | retrieval with   |                  |
|                  |                 | stats().         |                  |
+------------------+-----------------+------------------+------------------+
| regularity_check | OT_BOOL         | Throw exceptions | casadi::Function |
|                  |                 | when NaN or Inf  | Internal         |
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
| verbose          | OT_BOOL         | Verbose          | casadi::Function |
|                  |                 | evaluation  for  | Internal         |
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

- ecos

- gurobi

- highs

- hpmpc

- mosek

- ooqp

- osqp

- qpalm

- qpoases

- sqic

- superscs

- ipqp

- nlpsol

- qrqp

Note: some of the plugins in this list might not be available on your 

system. Also, there might be extra plugins available to you that are 
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

+---------------+-----------------------+----------------------------------+
|      Id       |         Type          |           Description            |
+===============+=======================+==================================+
| cplex         | OT_DICT               | Options to be passed to CPLEX    |
+---------------+-----------------------+----------------------------------+
| dep_check     | OT_INT                | Detect redundant constraints.    |
+---------------+-----------------------+----------------------------------+
| dump_filename | OT_STRING             | The filename to dump to.         |
+---------------+-----------------------+----------------------------------+
| dump_to_file  | OT_BOOL               | Dumps QP to file in CPLEX        |
|               |                       | format.                          |
+---------------+-----------------------+----------------------------------+
| mip_start     | OT_BOOL               | Hot start integers with x0       |
|               |                       | [Default false].                 |
+---------------+-----------------------+----------------------------------+
| qp_method     | OT_INT                | Determines which CPLEX algorithm |
|               |                       | to use.                          |
+---------------+-----------------------+----------------------------------+
| sos_groups    | OT_INTVECTORVECTOR    | Definition of SOS groups by      |
|               |                       | indices.                         |
+---------------+-----------------------+----------------------------------+
| sos_types     | OT_INTVECTOR          | Specify 1 or 2 for each SOS      |
|               |                       | group.                           |
+---------------+-----------------------+----------------------------------+
| sos_weights   | OT_DOUBLEVECTORVECTOR | Weights corresponding to SOS     |
|               |                       | entries.                         |
+---------------+-----------------------+----------------------------------+
| tol           | OT_DOUBLE             | Tolerance of solver              |
+---------------+-----------------------+----------------------------------+
| warm_start    | OT_BOOL               | Use warm start with simplex      |
|               |                       | methods (affects only the        |
|               |                       | simplex methods).                |
+---------------+-----------------------+----------------------------------+



--------------------------------------------------------------------------------

ecos
----



Interface to the ECOS Solver for quadratic programming

Extra doc: https://github.com/casadi/casadi/wiki/L_22e

>List of available options

+-------+-----------------+------------------------------------------------+
|  Id   |      Type       |                  Description                   |
+=======+=================+================================================+
| ecos  | OT_DICT         | Options to be passed to ecos.                  |
+-------+-----------------+------------------------------------------------+
| vtype | OT_STRINGVECTOR | Type of variables:                             |
|       |                 | [CONTINUOUS|binary|integer|semicont|semiint]   |
+-------+-----------------+------------------------------------------------+



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

mosek
-----



Interface to the MOSEK Solver for quadratic programming

Extra doc: https://github.com/casadi/casadi/wiki/L_21x

>List of available options

+-------+-----------------+------------------------------------------------+
|  Id   |      Type       |                  Description                   |
+=======+=================+================================================+
| mosek | OT_DICT         | Options to be passed to mosek.                 |
+-------+-----------------+------------------------------------------------+
| vtype | OT_STRINGVECTOR | Type of variables:                             |
|       |                 | [CONTINUOUS|binary|integer|semicont|semiint]   |
+-------+-----------------+------------------------------------------------+



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

qpalm
-----



Interface to the QPALM Solver for quadratic programming

Extra doc: https://github.com/casadi/casadi/wiki/L_22n

>List of available options

+-------------------+---------+--------------------------------------------+
|        Id         |  Type   |                Description                 |
+===================+=========+============================================+
| qpalm             | OT_DICT | const Options to be passed to qpalm.       |
+-------------------+---------+--------------------------------------------+
| warm_start_dual   | OT_BOOL | Use lam_a0 and lam_x0 input to warmstart   |
|                   |         | [Default: true].                           |
+-------------------+---------+--------------------------------------------+
| warm_start_primal | OT_BOOL | Use x0 input to warmstart [Default: true]. |
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
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L98
";

%feature("docstring") casadi::IndexRecution::has_conic "

Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L152
";

%feature("docstring") casadi::IndexRecution::load_conic "

Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L155
";

%feature("docstring") casadi::IndexRecution::doc_conic "

Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L158
";

%feature("docstring") casadi::IndexRecution::dplesol "



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
|                  |                 | ght*nf<=(1-ad_we |                  |
|                  |                 | ight)*na is used |                  |
|                  |                 | where nf and na  |                  |
|                  |                 | are estimates of |                  |
|                  |                 | the number of    |                  |
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
| print_time       | OT_BOOL         | print            | casadi::Function |
|                  |                 | information      | Internal         |
|                  |                 | about execution  |                  |
|                  |                 | time. Implies    |                  |
|                  |                 | record_time.     |                  |
+------------------+-----------------+------------------+------------------+
| record_time      | OT_BOOL         | record           | casadi::Function |
|                  |                 | information      | Internal         |
|                  |                 | about execution  |                  |
|                  |                 | time, for        |                  |
|                  |                 | retrieval with   |                  |
|                  |                 | stats().         |                  |
+------------------+-----------------+------------------+------------------+
| regularity_check | OT_BOOL         | Throw exceptions | casadi::Function |
|                  |                 | when NaN or Inf  | Internal         |
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
| verbose          | OT_BOOL         | Verbose          | casadi::Function |
|                  |                 | evaluation  for  | Internal         |
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

system. Also, there might be extra plugins available to you that are 
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
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L66

>  Function casadi::dplesol(const std::string &name, const std::string &solver, const SpDict &st, const Dict &opts=Dict())
------------------------------------------------------------------------



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
|                  |                 | ght*nf<=(1-ad_we |                  |
|                  |                 | ight)*na is used |                  |
|                  |                 | where nf and na  |                  |
|                  |                 | are estimates of |                  |
|                  |                 | the number of    |                  |
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
| print_time       | OT_BOOL         | print            | casadi::Function |
|                  |                 | information      | Internal         |
|                  |                 | about execution  |                  |
|                  |                 | time. Implies    |                  |
|                  |                 | record_time.     |                  |
+------------------+-----------------+------------------+------------------+
| record_time      | OT_BOOL         | record           | casadi::Function |
|                  |                 | information      | Internal         |
|                  |                 | about execution  |                  |
|                  |                 | time, for        |                  |
|                  |                 | retrieval with   |                  |
|                  |                 | stats().         |                  |
+------------------+-----------------+------------------+------------------+
| regularity_check | OT_BOOL         | Throw exceptions | casadi::Function |
|                  |                 | when NaN or Inf  | Internal         |
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
| verbose          | OT_BOOL         | Verbose          | casadi::Function |
|                  |                 | evaluation  for  | Internal         |
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

system. Also, there might be extra plugins available to you that are 
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
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L66

";

";

%feature("docstring") casadi::IndexRecution::dple_in "

Get DPLE input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ne

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L115

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L115-L122

>  std::string casadi::dple_in(casadi_int ind)
------------------------------------------------------------------------

Get DPLE input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ne

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L115

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L115-L122

";

";

%feature("docstring") casadi::IndexRecution::dple_out "

Get DPLE output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1nf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L124

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L124-L130

>  std::string casadi::dple_out(casadi_int ind)
------------------------------------------------------------------------

Get DPLE output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1nf

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L124

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L124-L130

";

";

%feature("docstring") casadi::IndexRecution::dple_n_in "

Get the number of QP solver inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1ng

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L132-L134

";

%feature("docstring") casadi::IndexRecution::dple_n_out "

Get the number of QP solver outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1nh

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L136

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.cpp#L136-L138

";

%feature("docstring") casadi::IndexRecution::has_dple "

Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L107
";

%feature("docstring") casadi::IndexRecution::load_dple "

Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L110
";

%feature("docstring") casadi::IndexRecution::doc_dple "

Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/dple.hpp#L113
";

%feature("docstring") casadi::IndexRecution::expm_n_in "

Get the number of expm solver inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_rs

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.hpp#L50

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.cpp#L50-L52

";

%feature("docstring") casadi::IndexRecution::expm_n_out "

Get the number of expm solver outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_rt

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.hpp#L54

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.cpp#L54-L56

";

%feature("docstring") casadi::IndexRecution::expmsol "



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
|                  |                 | ght*nf<=(1-ad_we |                  |
|                  |                 | ight)*na is used |                  |
|                  |                 | where nf and na  |                  |
|                  |                 | are estimates of |                  |
|                  |                 | the number of    |                  |
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
| print_time       | OT_BOOL         | print            | casadi::Function |
|                  |                 | information      | Internal         |
|                  |                 | about execution  |                  |
|                  |                 | time. Implies    |                  |
|                  |                 | record_time.     |                  |
+------------------+-----------------+------------------+------------------+
| record_time      | OT_BOOL         | record           | casadi::Function |
|                  |                 | information      | Internal         |
|                  |                 | about execution  |                  |
|                  |                 | time, for        |                  |
|                  |                 | retrieval with   |                  |
|                  |                 | stats().         |                  |
+------------------+-----------------+------------------+------------------+
| regularity_check | OT_BOOL         | Throw exceptions | casadi::Function |
|                  |                 | when NaN or Inf  | Internal         |
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
| verbose          | OT_BOOL         | Verbose          | casadi::Function |
|                  |                 | evaluation  for  | Internal         |
|                  |                 | debugging        |                  |
+------------------+-----------------+------------------+------------------+

List of plugins
- slicot

Note: some of the plugins in this list might not be available on your 

system. Also, there might be extra plugins available to you that are 
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
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.hpp#L58
";

%feature("docstring") casadi::IndexRecution::has_expm "

Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.hpp#L73
";

%feature("docstring") casadi::IndexRecution::load_expm "

Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.hpp#L76
";

%feature("docstring") casadi::IndexRecution::doc_expm "

Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/expm.hpp#L79
";

%feature("docstring") casadi::IndexRecution::external "

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

%feature("docstring") casadi::IndexRecution::combine "

[INTERNAL] 
Combine two dicts. First has priority.

Extra doc: https://github.com/casadi/casadi/wiki/L_17t

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L595

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L595-L601

";

%feature("docstring") casadi::IndexRecution::update_dict "

[INTERNAL] 
Update the target dictorionary in place with source elements.

Extra doc: https://github.com/casadi/casadi/wiki/L_17u

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.hpp#L603

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/generic_type.cpp#L603-L616

";

%feature("docstring") casadi::IndexRecution::get_from_dict "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::extract_from_dict "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::extract_from_dict_inplace "

[INTERNAL] ";

%feature("docstring") casadi::IndexRecution::is_slice "

Check if an index vector can be represented more efficiently as a 
slice.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L171

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L171-L201

>  bool casadi::is_slice(const std::vector< casadi_int > &v, bool ind1=false)
------------------------------------------------------------------------

Check if an index vector can be represented more efficiently as a 
slice.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L171

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L171-L201

";

";

%feature("docstring") casadi::IndexRecution::to_slice "

Construct from an index vector (requires is_slice(v) to be true)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L153-L169

>  Slice casadi::to_slice(const std::vector< casadi_int > &v, bool ind1=false)
------------------------------------------------------------------------

Construct from an index vector (requires is_slice(v) to be true)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L153

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L153-L169

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
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L132

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L132-L191

";

%feature("docstring") casadi::IndexRecution::collocation_interpolators "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L193

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L193-L235

";

%feature("docstring") casadi::IndexRecution::collocation_coeff "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L237

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L237-L291

";

%feature("docstring") casadi::IndexRecution::simpleIRK "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L293

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L293-L358

";

%feature("docstring") casadi::IndexRecution::simpleIntegrator "

Simplified wrapper for the  Integrator class.

Extra doc: https://github.com/casadi/casadi/wiki/L_1st

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.hpp#L360

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integration_tools.cpp#L360-L402

";

%feature("docstring") casadi::IndexRecution::integrator_in "

Get integrator input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_7d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L76

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L76-L87

>  std::string casadi::integrator_in(casadi_int ind)
------------------------------------------------------------------------

Get integrator input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_7d

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L76

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L76-L87

";

";

%feature("docstring") casadi::IndexRecution::integrator_out "

Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_7e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L89

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L89-L100

>  std::string casadi::integrator_out(casadi_int ind)
------------------------------------------------------------------------

Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_7e

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L89

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L89-L100

";

";

%feature("docstring") casadi::IndexRecution::integrator_n_in "

Get the number of integrator inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_7f

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L102

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L102-L104

";

%feature("docstring") casadi::IndexRecution::integrator_n_out "

Get the number of integrator outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_7g

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L106

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.cpp#L106-L108

";

%feature("docstring") casadi::IndexRecution::integrator "

[INTERNAL]

>  Function casadi::integrator(const std::string &name, const std::string &solver, const Function &dae, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::IndexRecution::has_integrator "

Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L96

";

%feature("docstring") casadi::IndexRecution::load_integrator "

Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L99

";

%feature("docstring") casadi::IndexRecution::doc_integrator "

Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/integrator.hpp#L102

";

%feature("docstring") casadi::IndexRecution::interpolant "

Parametric variant of interpolant.

The resulting function will have additional arguments for the grid and
 
coefficients

By default, derivatives wrt the coefficients are not supported (zero).
 Some
 interpolant plugins may support the  inline=true which enables correct 
derivatives

Extra doc: https://github.com/casadi/casadi/wiki/L_1p4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.hpp#L206

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.cpp#L206-L214

>  Function casadi::interpolant(const std::string &name, const std::string &solver, const std::vector< casadi_int > &grid_dims, casadi_int m=1, const Dict &opts=Dict())
------------------------------------------------------------------------

Parametric variant of interpolant.

The resulting function will have additional arguments for the grid and
 
coefficients

By default, derivatives wrt the coefficients are not supported (zero).
 Some
 interpolant plugins may support the  inline=true which enables correct 
derivatives

Extra doc: https://github.com/casadi/casadi/wiki/L_1p4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.hpp#L206

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.cpp#L206-L214

";

";

%feature("docstring") casadi::IndexRecution::has_interpolant "

Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.hpp#L151

";

%feature("docstring") casadi::IndexRecution::load_interpolant "

Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.hpp#L154

";

%feature("docstring") casadi::IndexRecution::doc_interpolant "

Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/interpolant.hpp#L157

";

%feature("docstring") casadi::IndexRecution::has_linsol "

Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L161
";

%feature("docstring") casadi::IndexRecution::load_linsol "

Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L164
";

%feature("docstring") casadi::IndexRecution::doc_linsol "

Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/linsol.hpp#L167
";

%feature("docstring") casadi::IndexRecution::detect_simple_bounds "

";

";

%feature("docstring") casadi::IndexRecution::check_sos "

Check sos structure and generate defaults.

Extra doc: https://github.com/casadi/casadi/wiki/L_1sx

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_tools.hpp#L79

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlp_tools.hpp#L79-L120

";

%feature("docstring") casadi::IndexRecution::nlpsol "

[INTERNAL]

>  Function casadi::nlpsol(const std::string &name, const std::string &solver, const Function &nlp, const Dict &opts=Dict())
------------------------------------------------------------------------
[INTERNAL] 
";

";

%feature("docstring") casadi::IndexRecution::nlpsol_in "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L159

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L159-L172

>  std::string casadi::nlpsol_in(casadi_int ind)
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L159

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L159-L172

";

";

%feature("docstring") casadi::IndexRecution::nlpsol_out "

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
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L174

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L174-L185

>  std::string casadi::nlpsol_out(casadi_int ind)
------------------------------------------------------------------------

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
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L174

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L174-L185

";

";

%feature("docstring") casadi::IndexRecution::nlpsol_n_in "

Number of NLP solver inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L187

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L187-L189

";

%feature("docstring") casadi::IndexRecution::nlpsol_n_out "

Number of NLP solver outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L191

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L191-L193

";

%feature("docstring") casadi::IndexRecution::nlpsol_options "

Get all options for a plugin.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L662

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L662-L664

";

%feature("docstring") casadi::IndexRecution::nlpsol_option_type "

Get type info for a particular option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L666

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L666-L668

";

%feature("docstring") casadi::IndexRecution::nlpsol_option_info "

Get documentation for a particular option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1t7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L670

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L670-L672

";

%feature("docstring") casadi::IndexRecution::has_nlpsol "

Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L158
";

%feature("docstring") casadi::IndexRecution::load_nlpsol "

Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L161
";

%feature("docstring") casadi::IndexRecution::doc_nlpsol "

Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L164
";

%feature("docstring") casadi::IndexRecution::rootfinder_in "

Get rootfinder input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L48

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L48-L55

>  std::string casadi::rootfinder_in(casadi_int ind)
------------------------------------------------------------------------

Get rootfinder input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u0

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L48

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L48-L55

";

";

%feature("docstring") casadi::IndexRecution::rootfinder_out "

Get rootfinder output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L57

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L57-L63

>  std::string casadi::rootfinder_out(casadi_int ind)
------------------------------------------------------------------------

Get rootfinder output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L57

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L57-L63

";

";

%feature("docstring") casadi::IndexRecution::rootfinder_n_in "

Number of rootfinder inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L65

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L65-L67

";

%feature("docstring") casadi::IndexRecution::rootfinder_n_out "

Number of rootfinder outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L69

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L69-L71

";

%feature("docstring") casadi::IndexRecution::rootfinder_options "

Get all options for a plugin.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L73

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L73-L75

";

%feature("docstring") casadi::IndexRecution::rootfinder_option_type "

Get type info for a particular option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L77

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L77-L79

";

%feature("docstring") casadi::IndexRecution::rootfinder_option_info "

Get documentation for a particular option.

Extra doc: https://github.com/casadi/casadi/wiki/L_1u6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L81

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.cpp#L81-L83

";

%feature("docstring") casadi::IndexRecution::rootfinder "



Create a solver for rootfinding problems Takes a function where one 
of the 
inputs is unknown and one of the outputs is a residual function
 that is 
always zero, defines a new function where the the unknown 
input has been 
replaced by a  guess for the unknown and the residual output has been 
replaced by the 
calculated value for the input.

For a function [y0, y1, ...,yi, .., yn] = F(x0, x1, ..., xj, ..., xm),
 
where xj is unknown and yi=0, defines a new function [y0, y1, ...,xj,
 .., 
yn] = G(x0, x1, ..., xj_guess, ..., xm),

xj and yi must have the same dimension and d(yi)/d(xj) must be 
invertable.

By default, the first input is unknown and the first output is the 

residual.
General information

>List of available options

+------------------+-----------------+------------------+------------------+
|        Id        |      Type       |   Description    |     Used in      |
+==================+=================+==================+==================+
| common_options   | OT_DICT         | Options for      | casadi::OracleFu |
|                  |                 | auto-generated   | nction           |
|                  |                 | functions        |                  |
+------------------+-----------------+------------------+------------------+
| constraints      | OT_INTVECTOR    | Constrain the    | casadi::Rootfind |
|                  |                 | unknowns. 0      | er               |
|                  |                 | (default): no    |                  |
|                  |                 | constraint on    |                  |
|                  |                 | ui, 1: ui >=     |                  |
|                  |                 | 0.0, -1: ui <=   |                  |
|                  |                 | 0.0, 2: ui >     |                  |
|                  |                 | 0.0, -2: ui <    |                  |
|                  |                 | 0.0.             |                  |
+------------------+-----------------+------------------+------------------+
| error_on_fail    | OT_BOOL         | When the         | casadi::Rootfind |
|                  |                 | numerical        | er               |
|                  |                 | process returns  |                  |
|                  |                 | unsuccessfully,  |                  |
|                  |                 | raise an error   |                  |
|                  |                 | (default false). |                  |
+------------------+-----------------+------------------+------------------+
| expand           | OT_BOOL         | Replace MX with  | casadi::OracleFu |
|                  |                 | SX expressions   | nction           |
|                  |                 | in problem       |                  |
|                  |                 | formulation      |                  |
|                  |                 | [false]          |                  |
+------------------+-----------------+------------------+------------------+
| implicit_input   | OT_INT          | Index of the     | casadi::Rootfind |
|                  |                 | input that       | er               |
|                  |                 | corresponds to   |                  |
|                  |                 | the actual root- |                  |
|                  |                 | finding          |                  |
+------------------+-----------------+------------------+------------------+
| implicit_output  | OT_INT          | Index of the     | casadi::Rootfind |
|                  |                 | output that      | er               |
|                  |                 | corresponds to   |                  |
|                  |                 | the actual root- |                  |
|                  |                 | finding          |                  |
+------------------+-----------------+------------------+------------------+
| jacobian_functio | OT_FUNCTION     | Function object  | casadi::Rootfind |
| n                |                 | for calculating  | er               |
|                  |                 | the Jacobian     |                  |
|                  |                 | (autogenerated   |                  |
|                  |                 | by default)      |                  |
+------------------+-----------------+------------------+------------------+
| linear_solver    | OT_STRING       | User-defined     | casadi::Rootfind |
|                  |                 | linear solver    | er               |
|                  |                 | class. Needed    |                  |
|                  |                 | for              |                  |
|                  |                 | sensitivities.   |                  |
+------------------+-----------------+------------------+------------------+
| linear_solver_op | OT_DICT         | Options to be    | casadi::Rootfind |
| tions            |                 | passed to the    | er               |
|                  |                 | linear solver.   |                  |
+------------------+-----------------+------------------+------------------+
| monitor          | OT_STRINGVECTOR | Set of user      | casadi::OracleFu |
|                  |                 | problem          | nction           |
|                  |                 | functions to be  |                  |
|                  |                 | monitored        |                  |
+------------------+-----------------+------------------+------------------+
| show_eval_warnin | OT_BOOL         | Show warnings    | casadi::OracleFu |
| gs               |                 | generated from   | nction           |
|                  |                 | function         |                  |
|                  |                 | evaluations      |                  |
|                  |                 | [true]           |                  |
+------------------+-----------------+------------------+------------------+
| specific_options | OT_DICT         | Options for      | casadi::OracleFu |
|                  |                 | specific auto-   | nction           |
|                  |                 | generated        |                  |
|                  |                 | functions,       |                  |
|                  |                 | overwriting the  |                  |
|                  |                 | defaults from    |                  |
|                  |                 | common_options.  |                  |
|                  |                 | Nested           |                  |
|                  |                 | dictionary.      |                  |
+------------------+-----------------+------------------+------------------+

>Input scheme: casadi::RootfinderInput (ROOTFINDER_NUM_IN = 2)

+---------------+-------+---------------------------------+
|   Full name   | Short |           Description           |
+===============+=======+=================================+
| ROOTFINDER_X0 | x0    | Initial guess for the solution. |
+---------------+-------+---------------------------------+
| ROOTFINDER_P  | p     | Parameters.                     |
+---------------+-------+---------------------------------+

>Output scheme: casadi::RootfinderOutput (ROOTFINDER_NUM_OUT = 1)

+--------------+-------+--------------------------------------+
|  Full name   | Short |             Description              |
+==============+=======+======================================+
| ROOTFINDER_X | x     | Solution to the system of equations. |
+--------------+-------+--------------------------------------+

List of plugins
- kinsol

- fast_newton

- nlpsol

- newton

Note: some of the plugins in this list might not be available on your 

system. Also, there might be extra plugins available to you that are 
not 
listed here. You can obtain their documentation with  
Rootfinder.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

kinsol
------



KINSOL interface from the Sundials suite

Extra doc: https://github.com/casadi/casadi/wiki/L_226

>List of available options

+---------------------------+-----------------+----------------------------+
|            Id             |      Type       |        Description         |
+===========================+=================+============================+
| abstol                    | OT_DOUBLE       | Stopping criterion         |
|                           |                 | tolerance                  |
+---------------------------+-----------------+----------------------------+
| disable_internal_warnings | OT_BOOL         | Disable KINSOL internal    |
|                           |                 | warning messages           |
+---------------------------+-----------------+----------------------------+
| exact_jacobian            | OT_BOOL         | Use exact Jacobian         |
|                           |                 | information                |
+---------------------------+-----------------+----------------------------+
| f_scale                   | OT_DOUBLEVECTOR | Equation scaling factors   |
+---------------------------+-----------------+----------------------------+
| iterative_solver          | OT_STRING       | gmres|bcgstab|tfqmr        |
+---------------------------+-----------------+----------------------------+
| linear_solver_type        | OT_STRING       | dense|banded|iterative|use |
|                           |                 | r_defined                  |
+---------------------------+-----------------+----------------------------+
| lower_bandwidth           | OT_INT          | Lower bandwidth for banded |
|                           |                 | linear solvers             |
+---------------------------+-----------------+----------------------------+
| max_iter                  | OT_INT          | Maximum number of Newton   |
|                           |                 | iterations. Putting 0 sets |
|                           |                 | the default value of       |
|                           |                 | KinSol.                    |
+---------------------------+-----------------+----------------------------+
| max_krylov                | OT_INT          | Maximum Krylov space       |
|                           |                 | dimension                  |
+---------------------------+-----------------+----------------------------+
| pretype                   | OT_STRING       | Type of preconditioner     |
+---------------------------+-----------------+----------------------------+
| print_level               | OT_INT          | Verbosity level            |
+---------------------------+-----------------+----------------------------+
| strategy                  | OT_STRING       | Globalization strategy     |
+---------------------------+-----------------+----------------------------+
| u_scale                   | OT_DOUBLEVECTOR | Variable scaling factors   |
+---------------------------+-----------------+----------------------------+
| upper_bandwidth           | OT_INT          | Upper bandwidth for banded |
|                           |                 | linear solvers             |
+---------------------------+-----------------+----------------------------+
| use_preconditioner        | OT_BOOL         | Precondition an iterative  |
|                           |                 | solver                     |
+---------------------------+-----------------+----------------------------+



--------------------------------------------------------------------------------

fast_newton
-----------



Implements simple newton iterations to solve an implicit function.

Extra doc: https://github.com/casadi/casadi/wiki/L_237

>List of available options

+------------+-----------+-------------------------------------------------+
|     Id     |   Type    |                   Description                   |
+============+===========+=================================================+
| abstol     | OT_DOUBLE | Stopping criterion tolerance on ||g||__inf)     |
+------------+-----------+-------------------------------------------------+
| abstolStep | OT_DOUBLE | Stopping criterion tolerance on step size       |
+------------+-----------+-------------------------------------------------+
| max_iter   | OT_INT    | Maximum number of Newton iterations to perform  |
|            |           | before returning.                               |
+------------+-----------+-------------------------------------------------+



--------------------------------------------------------------------------------

nlpsol
------





--------------------------------------------------------------------------------

newton
------



Implements simple newton iterations to solve an implicit function.

Extra doc: https://github.com/casadi/casadi/wiki/L_236

>List of available options

+-----------------+-----------+--------------------------------------------+
|       Id        |   Type    |                Description                 |
+=================+===========+============================================+
| abstol          | OT_DOUBLE | Stopping criterion tolerance on max(|F|)   |
+-----------------+-----------+--------------------------------------------+
| abstolStep      | OT_DOUBLE | Stopping criterion tolerance on step size  |
+-----------------+-----------+--------------------------------------------+
| line_search     | OT_BOOL   | Enable line-search (default: true)         |
+-----------------+-----------+--------------------------------------------+
| max_iter        | OT_INT    | Maximum number of Newton iterations to     |
|                 |           | perform before returning.                  |
+-----------------+-----------+--------------------------------------------+
| print_iteration | OT_BOOL   | Print information about each iteration     |
+-----------------+-----------+--------------------------------------------+

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_21r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L74

>  Function casadi::rootfinder(const std::string &name, const std::string &solver, const SXDict &rfp, const Dict &opts=Dict())
------------------------------------------------------------------------



Create a solver for rootfinding problems Takes a function where one 
of the 
inputs is unknown and one of the outputs is a residual function
 that is 
always zero, defines a new function where the the unknown 
input has been 
replaced by a  guess for the unknown and the residual output has been 
replaced by the 
calculated value for the input.

For a function [y0, y1, ...,yi, .., yn] = F(x0, x1, ..., xj, ..., xm),
 
where xj is unknown and yi=0, defines a new function [y0, y1, ...,xj,
 .., 
yn] = G(x0, x1, ..., xj_guess, ..., xm),

xj and yi must have the same dimension and d(yi)/d(xj) must be 
invertable.

By default, the first input is unknown and the first output is the 

residual.
General information

>List of available options

+------------------+-----------------+------------------+------------------+
|        Id        |      Type       |   Description    |     Used in      |
+==================+=================+==================+==================+
| common_options   | OT_DICT         | Options for      | casadi::OracleFu |
|                  |                 | auto-generated   | nction           |
|                  |                 | functions        |                  |
+------------------+-----------------+------------------+------------------+
| constraints      | OT_INTVECTOR    | Constrain the    | casadi::Rootfind |
|                  |                 | unknowns. 0      | er               |
|                  |                 | (default): no    |                  |
|                  |                 | constraint on    |                  |
|                  |                 | ui, 1: ui >=     |                  |
|                  |                 | 0.0, -1: ui <=   |                  |
|                  |                 | 0.0, 2: ui >     |                  |
|                  |                 | 0.0, -2: ui <    |                  |
|                  |                 | 0.0.             |                  |
+------------------+-----------------+------------------+------------------+
| error_on_fail    | OT_BOOL         | When the         | casadi::Rootfind |
|                  |                 | numerical        | er               |
|                  |                 | process returns  |                  |
|                  |                 | unsuccessfully,  |                  |
|                  |                 | raise an error   |                  |
|                  |                 | (default false). |                  |
+------------------+-----------------+------------------+------------------+
| expand           | OT_BOOL         | Replace MX with  | casadi::OracleFu |
|                  |                 | SX expressions   | nction           |
|                  |                 | in problem       |                  |
|                  |                 | formulation      |                  |
|                  |                 | [false]          |                  |
+------------------+-----------------+------------------+------------------+
| implicit_input   | OT_INT          | Index of the     | casadi::Rootfind |
|                  |                 | input that       | er               |
|                  |                 | corresponds to   |                  |
|                  |                 | the actual root- |                  |
|                  |                 | finding          |                  |
+------------------+-----------------+------------------+------------------+
| implicit_output  | OT_INT          | Index of the     | casadi::Rootfind |
|                  |                 | output that      | er               |
|                  |                 | corresponds to   |                  |
|                  |                 | the actual root- |                  |
|                  |                 | finding          |                  |
+------------------+-----------------+------------------+------------------+
| jacobian_functio | OT_FUNCTION     | Function object  | casadi::Rootfind |
| n                |                 | for calculating  | er               |
|                  |                 | the Jacobian     |                  |
|                  |                 | (autogenerated   |                  |
|                  |                 | by default)      |                  |
+------------------+-----------------+------------------+------------------+
| linear_solver    | OT_STRING       | User-defined     | casadi::Rootfind |
|                  |                 | linear solver    | er               |
|                  |                 | class. Needed    |                  |
|                  |                 | for              |                  |
|                  |                 | sensitivities.   |                  |
+------------------+-----------------+------------------+------------------+
| linear_solver_op | OT_DICT         | Options to be    | casadi::Rootfind |
| tions            |                 | passed to the    | er               |
|                  |                 | linear solver.   |                  |
+------------------+-----------------+------------------+------------------+
| monitor          | OT_STRINGVECTOR | Set of user      | casadi::OracleFu |
|                  |                 | problem          | nction           |
|                  |                 | functions to be  |                  |
|                  |                 | monitored        |                  |
+------------------+-----------------+------------------+------------------+
| show_eval_warnin | OT_BOOL         | Show warnings    | casadi::OracleFu |
| gs               |                 | generated from   | nction           |
|                  |                 | function         |                  |
|                  |                 | evaluations      |                  |
|                  |                 | [true]           |                  |
+------------------+-----------------+------------------+------------------+
| specific_options | OT_DICT         | Options for      | casadi::OracleFu |
|                  |                 | specific auto-   | nction           |
|                  |                 | generated        |                  |
|                  |                 | functions,       |                  |
|                  |                 | overwriting the  |                  |
|                  |                 | defaults from    |                  |
|                  |                 | common_options.  |                  |
|                  |                 | Nested           |                  |
|                  |                 | dictionary.      |                  |
+------------------+-----------------+------------------+------------------+

>Input scheme: casadi::RootfinderInput (ROOTFINDER_NUM_IN = 2)

+---------------+-------+---------------------------------+
|   Full name   | Short |           Description           |
+===============+=======+=================================+
| ROOTFINDER_X0 | x0    | Initial guess for the solution. |
+---------------+-------+---------------------------------+
| ROOTFINDER_P  | p     | Parameters.                     |
+---------------+-------+---------------------------------+

>Output scheme: casadi::RootfinderOutput (ROOTFINDER_NUM_OUT = 1)

+--------------+-------+--------------------------------------+
|  Full name   | Short |             Description              |
+==============+=======+======================================+
| ROOTFINDER_X | x     | Solution to the system of equations. |
+--------------+-------+--------------------------------------+

List of plugins
- kinsol

- fast_newton

- nlpsol

- newton

Note: some of the plugins in this list might not be available on your 

system. Also, there might be extra plugins available to you that are 
not 
listed here. You can obtain their documentation with  
Rootfinder.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

kinsol
------



KINSOL interface from the Sundials suite

Extra doc: https://github.com/casadi/casadi/wiki/L_226

>List of available options

+---------------------------+-----------------+----------------------------+
|            Id             |      Type       |        Description         |
+===========================+=================+============================+
| abstol                    | OT_DOUBLE       | Stopping criterion         |
|                           |                 | tolerance                  |
+---------------------------+-----------------+----------------------------+
| disable_internal_warnings | OT_BOOL         | Disable KINSOL internal    |
|                           |                 | warning messages           |
+---------------------------+-----------------+----------------------------+
| exact_jacobian            | OT_BOOL         | Use exact Jacobian         |
|                           |                 | information                |
+---------------------------+-----------------+----------------------------+
| f_scale                   | OT_DOUBLEVECTOR | Equation scaling factors   |
+---------------------------+-----------------+----------------------------+
| iterative_solver          | OT_STRING       | gmres|bcgstab|tfqmr        |
+---------------------------+-----------------+----------------------------+
| linear_solver_type        | OT_STRING       | dense|banded|iterative|use |
|                           |                 | r_defined                  |
+---------------------------+-----------------+----------------------------+
| lower_bandwidth           | OT_INT          | Lower bandwidth for banded |
|                           |                 | linear solvers             |
+---------------------------+-----------------+----------------------------+
| max_iter                  | OT_INT          | Maximum number of Newton   |
|                           |                 | iterations. Putting 0 sets |
|                           |                 | the default value of       |
|                           |                 | KinSol.                    |
+---------------------------+-----------------+----------------------------+
| max_krylov                | OT_INT          | Maximum Krylov space       |
|                           |                 | dimension                  |
+---------------------------+-----------------+----------------------------+
| pretype                   | OT_STRING       | Type of preconditioner     |
+---------------------------+-----------------+----------------------------+
| print_level               | OT_INT          | Verbosity level            |
+---------------------------+-----------------+----------------------------+
| strategy                  | OT_STRING       | Globalization strategy     |
+---------------------------+-----------------+----------------------------+
| u_scale                   | OT_DOUBLEVECTOR | Variable scaling factors   |
+---------------------------+-----------------+----------------------------+
| upper_bandwidth           | OT_INT          | Upper bandwidth for banded |
|                           |                 | linear solvers             |
+---------------------------+-----------------+----------------------------+
| use_preconditioner        | OT_BOOL         | Precondition an iterative  |
|                           |                 | solver                     |
+---------------------------+-----------------+----------------------------+



--------------------------------------------------------------------------------

fast_newton
-----------



Implements simple newton iterations to solve an implicit function.

Extra doc: https://github.com/casadi/casadi/wiki/L_237

>List of available options

+------------+-----------+-------------------------------------------------+
|     Id     |   Type    |                   Description                   |
+============+===========+=================================================+
| abstol     | OT_DOUBLE | Stopping criterion tolerance on ||g||__inf)     |
+------------+-----------+-------------------------------------------------+
| abstolStep | OT_DOUBLE | Stopping criterion tolerance on step size       |
+------------+-----------+-------------------------------------------------+
| max_iter   | OT_INT    | Maximum number of Newton iterations to perform  |
|            |           | before returning.                               |
+------------+-----------+-------------------------------------------------+



--------------------------------------------------------------------------------

nlpsol
------





--------------------------------------------------------------------------------

newton
------



Implements simple newton iterations to solve an implicit function.

Extra doc: https://github.com/casadi/casadi/wiki/L_236

>List of available options

+-----------------+-----------+--------------------------------------------+
|       Id        |   Type    |                Description                 |
+=================+===========+============================================+
| abstol          | OT_DOUBLE | Stopping criterion tolerance on max(|F|)   |
+-----------------+-----------+--------------------------------------------+
| abstolStep      | OT_DOUBLE | Stopping criterion tolerance on step size  |
+-----------------+-----------+--------------------------------------------+
| line_search     | OT_BOOL   | Enable line-search (default: true)         |
+-----------------+-----------+--------------------------------------------+
| max_iter        | OT_INT    | Maximum number of Newton iterations to     |
|                 |           | perform before returning.                  |
+-----------------+-----------+--------------------------------------------+
| print_iteration | OT_BOOL   | Print information about each iteration     |
+-----------------+-----------+--------------------------------------------+

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_21r

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L74

";

";

%feature("docstring") casadi::IndexRecution::has_rootfinder "

Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L128

";

%feature("docstring") casadi::IndexRecution::load_rootfinder "

Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L131

";

%feature("docstring") casadi::IndexRecution::doc_rootfinder "

Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/rootfinder.hpp#L134

";

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

%feature("docstring") casadi::IndexRecution::has_simulator "

Check if a particular plugin is available.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L53

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L53-L55

";

%feature("docstring") casadi::IndexRecution::load_simulator "

Explicitly load a plugin dynamically.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L57

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L57-L59

";

%feature("docstring") casadi::IndexRecution::doc_simulator "

Get the documentation string for a plugin.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L61

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L61-L63

";

%feature("docstring") casadi::IndexRecution::simulator "



Create an ODE/DAE simulator Solves an initial value problem (IVP) 
with the 
differential equation given as an implicit ODE coupled to an 
algebraic 
equation and a set of output equations:



::

  Time grid: [t0, t1, ..., tN]
  
  Initial conditions for forward integration
  x(t0)  = x0
  
  Sequential forward integration from t=tk to t=t{k+1}
  der(x) = fx(t, x, z, p, uk)         ODE
  0 = fz(t, x, z, p, uk)              Algebraic equations
  yk = fy(t, x, z, p, uk)             Output equations
  
  where we assume that the problem is index-1 (i.e. dfz/dz, is invertible) and furthermore that

General information

>List of available options

+------------------+-----------------+------------------+------------------+
|        Id        |      Type       |   Description    |     Used in      |
+==================+=================+==================+==================+
| common_options   | OT_DICT         | Options for      | casadi::OracleFu |
|                  |                 | auto-generated   | nction           |
|                  |                 | functions        |                  |
+------------------+-----------------+------------------+------------------+
| expand           | OT_BOOL         | Replace MX with  | casadi::OracleFu |
|                  |                 | SX expressions   | nction           |
|                  |                 | in problem       |                  |
|                  |                 | formulation      |                  |
|                  |                 | [false]          |                  |
+------------------+-----------------+------------------+------------------+
| monitor          | OT_STRINGVECTOR | Set of user      | casadi::OracleFu |
|                  |                 | problem          | nction           |
|                  |                 | functions to be  |                  |
|                  |                 | monitored        |                  |
+------------------+-----------------+------------------+------------------+
| nfwd             | OT_INT          | Number of        | casadi::Simulato |
|                  |                 | forward          | r                |
|                  |                 | sensitivities    |                  |
+------------------+-----------------+------------------+------------------+
| nondiff          | OT_BOOL         | Output nondiffer | casadi::Simulato |
|                  |                 | entiated         | r                |
+------------------+-----------------+------------------+------------------+
| print_stats      | OT_BOOL         | Print out        | casadi::Simulato |
|                  |                 | statistics after | r                |
|                  |                 | integration      |                  |
+------------------+-----------------+------------------+------------------+
| show_eval_warnin | OT_BOOL         | Show warnings    | casadi::OracleFu |
| gs               |                 | generated from   | nction           |
|                  |                 | function         |                  |
|                  |                 | evaluations      |                  |
|                  |                 | [true]           |                  |
+------------------+-----------------+------------------+------------------+
| specific_options | OT_DICT         | Options for      | casadi::OracleFu |
|                  |                 | specific auto-   | nction           |
|                  |                 | generated        |                  |
|                  |                 | functions,       |                  |
|                  |                 | overwriting the  |                  |
|                  |                 | defaults from    |                  |
|                  |                 | common_options.  |                  |
|                  |                 | Nested           |                  |
|                  |                 | dictionary.      |                  |
+------------------+-----------------+------------------+------------------+

>Input scheme: casadi::SimulatorInput (SIMULATOR_NUM_IN = 4)

+--------------+-------+---------------------------------------------------+
|  Full name   | Short |                    Description                    |
+==============+=======+===================================================+
| SIMULATOR_X0 | x0    | Differential state at the initial time.           |
+--------------+-------+---------------------------------------------------+
| SIMULATOR_U  | u     | Controls.                                         |
+--------------+-------+---------------------------------------------------+
| SIMULATOR_Z0 | z0    | Initial guess for the algebraic variable at the   |
|              |       | initial time.                                     |
+--------------+-------+---------------------------------------------------+
| SIMULATOR_P  | p     | Parameters.                                       |
+--------------+-------+---------------------------------------------------+

>Output scheme: casadi::SimulatorOutput (SIMULATOR_NUM_OUT = 3)

+-------------+-------+---------------------+
|  Full name  | Short |     Description     |
+=============+=======+=====================+
| SIMULATOR_X | x     | Differential state. |
+-------------+-------+---------------------+
| SIMULATOR_Y | y     | Outputs.            |
+-------------+-------+---------------------+
| SIMULATOR_Z | z     | Algebraic variable. |
+-------------+-------+---------------------+

List of plugins
- cvodes

Note: some of the plugins in this list might not be available on your 

system. Also, there might be extra plugins available to you that are 
not 
listed here. You can obtain their documentation with  
Simulator.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

cvodes
------



Interface to CVodes from the Sundials suite.

Extra doc: https://github.com/casadi/casadi/wiki/L_227

>List of available options

+----------------------------+-----------+---------------------------------+
|             Id             |   Type    |           Description           |
+============================+===========+=================================+
| abstol                     | OT_DOUBLE | Absolute tolerence for the IVP  |
|                            |           | solution                        |
+----------------------------+-----------+---------------------------------+
| disable_internal_warnings  | OT_BOOL   | Disable SUNDIALS internal       |
|                            |           | warning messages                |
+----------------------------+-----------+---------------------------------+
| fsens_err_con              | OT_BOOL   | include the forward             |
|                            |           | sensitivities in all error      |
|                            |           | controls                        |
+----------------------------+-----------+---------------------------------+
| interpolation_type         | OT_STRING | Type of interpolation for the   |
|                            |           | adjoint sensitivities           |
+----------------------------+-----------+---------------------------------+
| linear_multistep_method    | OT_STRING | Simulator scheme: BDF|adams     |
+----------------------------+-----------+---------------------------------+
| linear_solver              | OT_STRING | A custom linear solver creator  |
|                            |           | function [default: qr]          |
+----------------------------+-----------+---------------------------------+
| linear_solver_options      | OT_DICT   | Options to be passed to the     |
|                            |           | linear solver                   |
+----------------------------+-----------+---------------------------------+
| max_krylov                 | OT_INT    | Maximum Krylov subspace size    |
+----------------------------+-----------+---------------------------------+
| max_multistep_order        | OT_INT    | Maximum order for the           |
|                            |           | (variable-order) multistep      |
|                            |           | method                          |
+----------------------------+-----------+---------------------------------+
| max_num_steps              | OT_INT    | Maximum number of simulator     |
|                            |           | steps                           |
+----------------------------+-----------+---------------------------------+
| max_order                  | OT_DOUBLE | Maximum order                   |
+----------------------------+-----------+---------------------------------+
| max_step_size              | OT_DOUBLE | Max step size [default: 0/inf]  |
+----------------------------+-----------+---------------------------------+
| min_step_size              | OT_DOUBLE | Min step size [default: 0/0.0]  |
+----------------------------+-----------+---------------------------------+
| newton_scheme              | OT_STRING | Linear solver scheme in the     |
|                            |           | Newton method:                  |
|                            |           | DIRECT|gmres|bcgstab|tfqmr      |
+----------------------------+-----------+---------------------------------+
| nonlin_conv_coeff          | OT_DOUBLE | Coefficient in the nonlinear    |
|                            |           | convergence test                |
+----------------------------+-----------+---------------------------------+
| nonlinear_solver_iteration | OT_STRING | Nonlinear solver type:          |
|                            |           | NEWTON|functional               |
+----------------------------+-----------+---------------------------------+
| quad_err_con               | OT_BOOL   | Should the quadratures affect   |
|                            |           | the step size control           |
+----------------------------+-----------+---------------------------------+
| reltol                     | OT_DOUBLE | Relative tolerence for the IVP  |
|                            |           | solution                        |
+----------------------------+-----------+---------------------------------+
| scale_abstol               | OT_BOOL   | Scale absolute tolerance by     |
|                            |           | nominal value                   |
+----------------------------+-----------+---------------------------------+
| sensitivity_method         | OT_STRING | Sensitivity method:             |
|                            |           | SIMULTANEOUS|staggered          |
+----------------------------+-----------+---------------------------------+
| step0                      | OT_DOUBLE | initial step size [default:     |
|                            |           | 0/estimated]                    |
+----------------------------+-----------+---------------------------------+
| steps_per_checkpoint       | OT_INT    | Number of steps between two     |
|                            |           | consecutive checkpoints         |
+----------------------------+-----------+---------------------------------+
| stop_at_end                | OT_BOOL   | Stop the simulator at the end   |
|                            |           | of the interval                 |
+----------------------------+-----------+---------------------------------+
| use_preconditioner         | OT_BOOL   | Precondition the iterative      |
|                            |           | solver [default: true]          |
+----------------------------+-----------+---------------------------------+

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_21m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L75

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L75-L78

>  Function casadi::simulator(const std::string &name, const std::string &solver, const SXDict &dae, const std::vector< double > &grid, const Dict &opts=Dict())
------------------------------------------------------------------------



Create an ODE/DAE simulator Solves an initial value problem (IVP) 
with the 
differential equation given as an implicit ODE coupled to an 
algebraic 
equation and a set of output equations:



::

  Time grid: [t0, t1, ..., tN]
  
  Initial conditions for forward integration
  x(t0)  = x0
  
  Sequential forward integration from t=tk to t=t{k+1}
  der(x) = fx(t, x, z, p, uk)         ODE
  0 = fz(t, x, z, p, uk)              Algebraic equations
  yk = fy(t, x, z, p, uk)             Output equations
  
  where we assume that the problem is index-1 (i.e. dfz/dz, is invertible) and furthermore that

General information

>List of available options

+------------------+-----------------+------------------+------------------+
|        Id        |      Type       |   Description    |     Used in      |
+==================+=================+==================+==================+
| common_options   | OT_DICT         | Options for      | casadi::OracleFu |
|                  |                 | auto-generated   | nction           |
|                  |                 | functions        |                  |
+------------------+-----------------+------------------+------------------+
| expand           | OT_BOOL         | Replace MX with  | casadi::OracleFu |
|                  |                 | SX expressions   | nction           |
|                  |                 | in problem       |                  |
|                  |                 | formulation      |                  |
|                  |                 | [false]          |                  |
+------------------+-----------------+------------------+------------------+
| monitor          | OT_STRINGVECTOR | Set of user      | casadi::OracleFu |
|                  |                 | problem          | nction           |
|                  |                 | functions to be  |                  |
|                  |                 | monitored        |                  |
+------------------+-----------------+------------------+------------------+
| nfwd             | OT_INT          | Number of        | casadi::Simulato |
|                  |                 | forward          | r                |
|                  |                 | sensitivities    |                  |
+------------------+-----------------+------------------+------------------+
| nondiff          | OT_BOOL         | Output nondiffer | casadi::Simulato |
|                  |                 | entiated         | r                |
+------------------+-----------------+------------------+------------------+
| print_stats      | OT_BOOL         | Print out        | casadi::Simulato |
|                  |                 | statistics after | r                |
|                  |                 | integration      |                  |
+------------------+-----------------+------------------+------------------+
| show_eval_warnin | OT_BOOL         | Show warnings    | casadi::OracleFu |
| gs               |                 | generated from   | nction           |
|                  |                 | function         |                  |
|                  |                 | evaluations      |                  |
|                  |                 | [true]           |                  |
+------------------+-----------------+------------------+------------------+
| specific_options | OT_DICT         | Options for      | casadi::OracleFu |
|                  |                 | specific auto-   | nction           |
|                  |                 | generated        |                  |
|                  |                 | functions,       |                  |
|                  |                 | overwriting the  |                  |
|                  |                 | defaults from    |                  |
|                  |                 | common_options.  |                  |
|                  |                 | Nested           |                  |
|                  |                 | dictionary.      |                  |
+------------------+-----------------+------------------+------------------+

>Input scheme: casadi::SimulatorInput (SIMULATOR_NUM_IN = 4)

+--------------+-------+---------------------------------------------------+
|  Full name   | Short |                    Description                    |
+==============+=======+===================================================+
| SIMULATOR_X0 | x0    | Differential state at the initial time.           |
+--------------+-------+---------------------------------------------------+
| SIMULATOR_U  | u     | Controls.                                         |
+--------------+-------+---------------------------------------------------+
| SIMULATOR_Z0 | z0    | Initial guess for the algebraic variable at the   |
|              |       | initial time.                                     |
+--------------+-------+---------------------------------------------------+
| SIMULATOR_P  | p     | Parameters.                                       |
+--------------+-------+---------------------------------------------------+

>Output scheme: casadi::SimulatorOutput (SIMULATOR_NUM_OUT = 3)

+-------------+-------+---------------------+
|  Full name  | Short |     Description     |
+=============+=======+=====================+
| SIMULATOR_X | x     | Differential state. |
+-------------+-------+---------------------+
| SIMULATOR_Y | y     | Outputs.            |
+-------------+-------+---------------------+
| SIMULATOR_Z | z     | Algebraic variable. |
+-------------+-------+---------------------+

List of plugins
- cvodes

Note: some of the plugins in this list might not be available on your 

system. Also, there might be extra plugins available to you that are 
not 
listed here. You can obtain their documentation with  
Simulator.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

cvodes
------



Interface to CVodes from the Sundials suite.

Extra doc: https://github.com/casadi/casadi/wiki/L_227

>List of available options

+----------------------------+-----------+---------------------------------+
|             Id             |   Type    |           Description           |
+============================+===========+=================================+
| abstol                     | OT_DOUBLE | Absolute tolerence for the IVP  |
|                            |           | solution                        |
+----------------------------+-----------+---------------------------------+
| disable_internal_warnings  | OT_BOOL   | Disable SUNDIALS internal       |
|                            |           | warning messages                |
+----------------------------+-----------+---------------------------------+
| fsens_err_con              | OT_BOOL   | include the forward             |
|                            |           | sensitivities in all error      |
|                            |           | controls                        |
+----------------------------+-----------+---------------------------------+
| interpolation_type         | OT_STRING | Type of interpolation for the   |
|                            |           | adjoint sensitivities           |
+----------------------------+-----------+---------------------------------+
| linear_multistep_method    | OT_STRING | Simulator scheme: BDF|adams     |
+----------------------------+-----------+---------------------------------+
| linear_solver              | OT_STRING | A custom linear solver creator  |
|                            |           | function [default: qr]          |
+----------------------------+-----------+---------------------------------+
| linear_solver_options      | OT_DICT   | Options to be passed to the     |
|                            |           | linear solver                   |
+----------------------------+-----------+---------------------------------+
| max_krylov                 | OT_INT    | Maximum Krylov subspace size    |
+----------------------------+-----------+---------------------------------+
| max_multistep_order        | OT_INT    | Maximum order for the           |
|                            |           | (variable-order) multistep      |
|                            |           | method                          |
+----------------------------+-----------+---------------------------------+
| max_num_steps              | OT_INT    | Maximum number of simulator     |
|                            |           | steps                           |
+----------------------------+-----------+---------------------------------+
| max_order                  | OT_DOUBLE | Maximum order                   |
+----------------------------+-----------+---------------------------------+
| max_step_size              | OT_DOUBLE | Max step size [default: 0/inf]  |
+----------------------------+-----------+---------------------------------+
| min_step_size              | OT_DOUBLE | Min step size [default: 0/0.0]  |
+----------------------------+-----------+---------------------------------+
| newton_scheme              | OT_STRING | Linear solver scheme in the     |
|                            |           | Newton method:                  |
|                            |           | DIRECT|gmres|bcgstab|tfqmr      |
+----------------------------+-----------+---------------------------------+
| nonlin_conv_coeff          | OT_DOUBLE | Coefficient in the nonlinear    |
|                            |           | convergence test                |
+----------------------------+-----------+---------------------------------+
| nonlinear_solver_iteration | OT_STRING | Nonlinear solver type:          |
|                            |           | NEWTON|functional               |
+----------------------------+-----------+---------------------------------+
| quad_err_con               | OT_BOOL   | Should the quadratures affect   |
|                            |           | the step size control           |
+----------------------------+-----------+---------------------------------+
| reltol                     | OT_DOUBLE | Relative tolerence for the IVP  |
|                            |           | solution                        |
+----------------------------+-----------+---------------------------------+
| scale_abstol               | OT_BOOL   | Scale absolute tolerance by     |
|                            |           | nominal value                   |
+----------------------------+-----------+---------------------------------+
| sensitivity_method         | OT_STRING | Sensitivity method:             |
|                            |           | SIMULTANEOUS|staggered          |
+----------------------------+-----------+---------------------------------+
| step0                      | OT_DOUBLE | initial step size [default:     |
|                            |           | 0/estimated]                    |
+----------------------------+-----------+---------------------------------+
| steps_per_checkpoint       | OT_INT    | Number of steps between two     |
|                            |           | consecutive checkpoints         |
+----------------------------+-----------+---------------------------------+
| stop_at_end                | OT_BOOL   | Stop the simulator at the end   |
|                            |           | of the interval                 |
+----------------------------+-----------+---------------------------------+
| use_preconditioner         | OT_BOOL   | Precondition the iterative      |
|                            |           | solver [default: true]          |
+----------------------------+-----------+---------------------------------+

Joel Andersson

Extra doc: https://github.com/casadi/casadi/wiki/L_21m

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L75

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L75-L78

";

";

%feature("docstring") casadi::IndexRecution::simulator_in "

Get simulator input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L97

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L97-L106

>  std::string casadi::simulator_in(casadi_int ind)
------------------------------------------------------------------------

Get simulator input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_x5

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L97

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L97-L106

";

";

%feature("docstring") casadi::IndexRecution::simulator_out "

Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L108

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L108-L116

>  std::string casadi::simulator_out(casadi_int ind)
------------------------------------------------------------------------

Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_x6

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L108

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L108-L116

";

";

%feature("docstring") casadi::IndexRecution::dyn_in "

Get simulator input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_x1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L126

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L126-L128

>  std::string casadi::dyn_in(casadi_int ind)
------------------------------------------------------------------------

Get simulator input scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_x1

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L126

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L126-L128

";

";

%feature("docstring") casadi::IndexRecution::dyn_out "

Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_x2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L130

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L130-L132

>  std::string casadi::dyn_out(casadi_int ind)
------------------------------------------------------------------------

Get output scheme name by index.

Extra doc: https://github.com/casadi/casadi/wiki/L_x2

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L130

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.cpp#L130-L132

";

";

%feature("docstring") casadi::IndexRecution::dyn_n_in "

Get the number of simulator inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_x3

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L168

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L168-L168

";

%feature("docstring") casadi::IndexRecution::dyn_n_out "

Get the number of simulator outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_x4

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L173

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L173-L173

";

%feature("docstring") casadi::IndexRecution::simulator_n_in "

Get the number of simulator inputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_x7

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L188

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L188-L188

";

%feature("docstring") casadi::IndexRecution::simulator_n_out "

Get the number of simulator outputs.

Extra doc: https://github.com/casadi/casadi/wiki/L_x8

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L193

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/simulator.hpp#L193-L193

";

%feature("docstring") casadi::IndexRecution::is_slice2 "

Check if an index vector can be represented more efficiently as two 
nested 
slices.

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L203

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L203-L254

";

%feature("docstring") casadi::IndexRecution::to_slice2 "

Construct nested slices from an index vector (requires is_slice2(v) to
 be 
true)

Doc source: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.hpp#L256

Implementation: 
https://github.com/casadi/casadi/blob/develop/casadi/core/slice.cpp#L256-L290

";


// File: namespacecasadi_1_1IndexRecution.xml


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

