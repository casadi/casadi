
// File: index.xml

// File: classcasadi_1_1Adaptor.xml
%feature("docstring") casadi::Adaptor "[INTERNAL]  A helper class for a
Solver that delegates work to another solver.

Joris Gillis

C++ includes: adaptor.hpp ";

%feature("docstring") casadi::Adaptor::addOptions "[INTERNAL]  Add options
that are common to all Adaptor classes.

";


// File: classcasadi_1_1Assertion.xml


// File: classcasadi_1_1BinaryMX.xml


// File: classcasadi_1_1BinarySX.xml


// File: classcasadi_1_1Call.xml


// File: classcasadi_1_1Callback.xml


/*  Simple Getters & Setters  */

/*  Advanced Getters  */

/*  Option Functionality  */ %feature("docstring")
casadi::Callback::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::Callback::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::Callback::size2_in "

Get input dimension.

";

%feature("docstring") casadi::Callback::call "

Evaluate the function symbolically or numerically.

";

%feature("docstring") casadi::Callback::get_output_shape "

Specify output shape.

Specify the shape corresponding to a given output. The shape must not be
changed over the lifetime of the object

Default implementation: scalar (1,1)

";

%feature("docstring") casadi::Callback::get_n_in "

Number of input arguments.

Specify the number of input arguments that a specific instance can handle.
The number must not be changed over the lifetime of the object

Default implementation: 1

";

%feature("docstring") casadi::Callback::size_in "

Get input dimension.

";

%feature("docstring") casadi::Callback::getWorkSize "

Get the length of the work vector.

";

%feature("docstring") casadi::Callback::getOptionEnumValue "[INTERNAL]  Get
the enum value corresponding to th certain option.

";

%feature("docstring") casadi::Callback::setOptionByAllowedIndex "[INTERNAL]
Set a certain option by giving its index into the allowed values.

";

%feature("docstring") casadi::Callback::sz_arg "[INTERNAL]  Get required
length of arg field.

";

%feature("docstring") casadi::Callback::generateLiftingFunctions "

Extract the functions needed for the Lifted Newton method.

";

%feature("docstring") casadi::Callback::print "

Print a description of the object.

";

%feature("docstring") casadi::Callback::nnz_out "

Get of number of output nonzeros For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::Callback::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

";

%feature("docstring") casadi::Callback::description_in "

Get input scheme description by index.

";

%feature("docstring") casadi::Callback::free_sx "

Get all the free variables of the function.

";

%feature("docstring") casadi::Callback::name_out "

>  [str] Function.name_out() const 
------------------------------------------------------------------------

Get output scheme.

>  str Function.name_out(int ind) const 
------------------------------------------------------------------------

Get output scheme name by index.

";

%feature("docstring") casadi::Callback::nnz_in "

Get of number of input nonzeros For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::Callback::spEvaluate "[INTERNAL]  Propagate
the sparsity pattern through a set of directional.

derivatives forward or backward (for usage, see the example
propagating_sparsity.cpp)

";

%feature("docstring") casadi::Callback::getStat "

Get a single statistic obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::Callback::default_in "

Get default input value (NOTE: constant reference)

";

%feature("docstring") casadi::Callback::init "

Initialize the object This function is called after the object construction
(for the whole class hierarchy) is complete, but before the finalization
step. It is called recursively for the whole class hierarchy, starting with
the lowest level.

";

%feature("docstring") casadi::Callback::numel_out "

Get of number of output elements For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::Callback::get_n_out "

Number of output arguments.

Specify the number of output arguments that a specific instance can handle.
The number must not be changed over the lifetime of the object

Default implementation: 1

";

%feature("docstring") casadi::Callback::generate "

Export / Generate C code for the function.

";

%feature("docstring") casadi::Callback::getAtomicInput "

Get the (integer) input arguments of an atomic operation.

";

%feature("docstring") casadi::Callback::setOptionByEnumValue "[INTERNAL]
Set a certain option by giving an enum value.

";

%feature("docstring") casadi::Callback::getAtomicInputReal "

Get the floating point output argument of an atomic operation.

";

%feature("docstring") casadi::Callback::checkInputs "[INTERNAL]  Check if
the numerical values of the supplied bounds make sense.

";

%feature("docstring") casadi::Callback::derForward "

Get a function that calculates nfwd forward derivatives.

Returns a function with n_in + n_out +nfwd*n_in inputs and nfwd*n_out
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nfwd*n_in inputs correspond to forward seeds, one direction at a time The
nfwd*n_out outputs correspond to forward sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::Callback::get_n_reverse "

Return function that calculates adjoint derivatives derReverse(nadj) returns
a cached instance if available, and calls  Function getDerReverse(int nadj)
if no cached version is available.

";

%feature("docstring") casadi::Callback::type_name "

Get type name.

";

%feature("docstring") casadi::Callback::jacobian "

Generate a Jacobian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. If
compact is set to true, only the nonzeros of the input and output
expressions are considered. If symmetric is set to true, the Jacobian being
calculated is known to be symmetric (usually a Hessian), which can be
exploited by the algorithm.

The generated Jacobian has one more output than the calling function
corresponding to the Jacobian and the same number of inputs.

";

%feature("docstring") casadi::Callback::sz_iw "[INTERNAL]  Get required
length of iw field.

";

%feature("docstring") casadi::Callback::getOutput "

>  DMatrix  IOInterface< Function  >.getOutput(int oind=0) const
------------------------------------------------------------------------

Get an output by index.

Parameters:
-----------

oind:  index within the range [0..n_out()-1]

>  DMatrix  IOInterface< Function  >.getOutput(str oname) const
------------------------------------------------------------------------

Get an output by name.

Parameters:
-----------

oname:  output name. Only allowed when an output scheme is set.

>  void IOInterface< Function  >.getOutput(T val, int oind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.getOutput(T val, str oname)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::Callback::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::Callback::sx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::Callback::Callback "

>  Callback()
------------------------------------------------------------------------

Default constructor.

>  Callback(Callback obj)
------------------------------------------------------------------------

Copy constructor (throws an error)

";

%feature("docstring") casadi::Callback::~Callback "

Destructor.

";

%feature("docstring") casadi::Callback::getAtomicOutput "

Get the (integer) output argument of an atomic operation.

";

%feature("docstring") casadi::Callback::hessian "

Generate a Hessian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The generated Hessian has two more outputs than the calling function
corresponding to the Hessian and the gradients.

";

%feature("docstring") casadi::Callback::description_out "

Get output scheme description by index.

";

%feature("docstring") casadi::Callback::get_input_sparsity "

Specify input sparsity.

Specify the sparsity corresponding to a given input. The sparsity must not
be changed over the lifetime of the object

Default implementation: dense using inputShape

";

%feature("docstring") casadi::Callback::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::Callback::repr "

Print a representation of the object.

";

%feature("docstring") casadi::Callback::setFullJacobian "

Set the Jacobian of all the input nonzeros with respect to all output
nonzeros NOTE: Does not take ownership, only weak references to the Jacobian
are kept internally

";

%feature("docstring") casadi::Callback::spInit "[INTERNAL]  Reset the
sparsity propagation.

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::Callback::index_out "

Find the index for a string describing a particular entry of an output
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::Callback::construct "

Construct internal object This is the step that actually construct the
internal object, as the class constructor only creates a null pointer. It
should be called from the user constructor.

";

%feature("docstring") casadi::Callback::get_n_forward "

Return function that calculates forward derivatives derForward(nfwd) returns
a cached instance if available, and calls  Function getDerForward(int nfwd)
if no cached version is available.

";

%feature("docstring") casadi::Callback::getOptionAllowedIndex "[INTERNAL]
Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::Callback::sx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::Callback::get_forward "

Return function that calculates forward derivatives derForward(nfwd) returns
a cached instance if available, and calls  Function getDerForward(int nfwd)
if no cached version is available.

";

%feature("docstring") casadi::Callback::getStats "

Get all statistics obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::Callback::get_input_shape "

Specify input shape.

Specify the shape corresponding to a given input. The shape must not be
changed over the lifetime of the object

Default implementation: scalar (1,1)

";

%feature("docstring") casadi::Callback::dictionary "

Get the dictionary.

";

%feature("docstring") casadi::Callback::getOptionTypeName "

Get the type name of a certain option.

";

%feature("docstring") casadi::Callback::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::Callback::setOutput "

>  void IOInterface< Function  >.setOutput(T val, int oind=0)
------------------------------------------------------------------------

Set an output by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.setOutput(T val, str oname)
------------------------------------------------------------------------

Set an output by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::Callback::size1_out "

Get output dimension.

";

%feature("docstring") casadi::Callback::spCanEvaluate "[INTERNAL]  Is the
class able to propagate seeds through the algorithm?

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::Callback::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::Callback::has_jacobian "

Return Jacobian of all input elements with respect to all output elements.

";

%feature("docstring") casadi::Callback::get_output_sparsity "

Specify output sparsity.

Specify the sparsity corresponding to a given output. The sparsity must not
be changed over the lifetime of the object

Default implementation: dense using outputShape

";

%feature("docstring") casadi::Callback::finalize "

Finalize the object This function is called after the construction and init
steps are completed, but before user functions are called. It is called
recursively for the whole class hierarchy, starting with the highest level.

";

%feature("docstring") casadi::Callback::size2_out "

Get output dimension.

";

%feature("docstring") casadi::Callback::evaluate "

Evaluate.

";

%feature("docstring") casadi::Callback::mx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::Callback::is_a "

Check if the function is of a particular type Optionally check if name
matches one of the base classes (default true)

";

%feature("docstring") casadi::Callback::getRepresentation "

Return a string with a representation (for SWIG)

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
  



";

%feature("docstring") casadi::Callback::jacSparsity "

Get, if necessary generate, the sparsity of a Jacobian block

";

%feature("docstring") casadi::Callback::name_in "

>  [str] Function.name_in() const 
------------------------------------------------------------------------

Get input scheme.

>  str Function.name_in(int ind) const 
------------------------------------------------------------------------

Get input scheme name by index.

";

%feature("docstring") casadi::Callback::name "

Name of the function.

";

%feature("docstring") casadi::Callback::tangent "

Generate a tangent function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the input must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::Callback::derivative "

Get a function that calculates nfwd forward derivatives and nadj adjoint
derivatives Legacy function: Use derForward and derReverse instead.

Returns a function with (1+nfwd)*n_in+nadj*n_out inputs and (1+nfwd)*n_out +
nadj*n_in outputs. The first n_in inputs correspond to nondifferentiated
inputs. The next nfwd*n_in inputs correspond to forward seeds, one direction
at a time and the last nadj*n_out inputs correspond to adjoint seeds, one
direction at a time. The first n_out outputs correspond to nondifferentiated
outputs. The next nfwd*n_out outputs correspond to forward sensitivities,
one direction at a time and the last nadj*n_in outputs corresponds to
adjoint sensitivities, one direction at a time.

(n_in = n_in(), n_out = n_out())

";

%feature("docstring") casadi::Callback::n_in "

Get the number of function inputs.

";

%feature("docstring") casadi::Callback::sz_w "[INTERNAL]  Get required
length of w field.

";

%feature("docstring") casadi::Callback::gradient "

Generate a gradient function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the output must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::Callback::callDerivative "[INTERNAL]
Evaluate the function symbolically or numerically with directional
derivatives The first two arguments are the nondifferentiated inputs and
results of the evaluation, the next two arguments are a set of forward
directional seeds and the resulting forward directional derivatives, the
length of the vector being the number of forward directions. The next two
arguments are a set of adjoint directional seeds and the resulting adjoint
directional derivatives, the length of the vector being the number of
adjoint directions.

";

%feature("docstring") casadi::Callback::callForward "

Create call to (cached) derivative function, forward mode.

";

%feature("docstring") casadi::Callback::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::Callback::getAlgorithmSize "

Get the number of atomic operations.

";

%feature("docstring") casadi::Callback::sz_res "[INTERNAL]  Get required
length of res field.

";

%feature("docstring") casadi::Callback::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::Callback::countNodes "

Number of nodes in the algorithm.

";

%feature("docstring") casadi::Callback::get_jacobian "

Return Jacobian of all input elements with respect to all output elements.

";

%feature("docstring") casadi::Callback::printDimensions "

Print dimensions of inputs and outputs.

";

%feature("docstring") casadi::Callback::setInput "

>  void IOInterface< Function  >.setInput(T val, int iind=0)
------------------------------------------------------------------------

Set an input by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.setInput(T val, str iname)
------------------------------------------------------------------------

Set an input by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::Callback::isNull "

Is a null pointer?

";

%feature("docstring") casadi::Callback::setJacobian "

Set the Jacobian function of output oind with respect to input iind NOTE:
Does not take ownership, only weak references to the Jacobians are kept
internally

";

%feature("docstring") casadi::Callback::eval "

Evaluate numerically, temporary matrices and work vectors.

";

%feature("docstring") casadi::Callback::setDerReverse "

Set a function that calculates nadj adjoint derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::Callback::numel_in "

Get of number of input elements For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::Callback::getAtomicOperation "

Get an atomic operation operator index.

";

%feature("docstring") casadi::Callback::map "

>  [[MX] ] Function.map([[MX ] ] arg, str parallelization=\"serial\")

>  [MX] Function.map([MX ] arg, str parallelization=\"serial\")
------------------------------------------------------------------------

Evaluate symbolically in parallel (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

>  Function Function.map(str name, int N, Dict options=Dict()) const 
------------------------------------------------------------------------

Create a mapped version of this function.

Suppose the function has a signature of:

::

     f: (a, p) -> ( s )
  



The the mapaccumulated version has the signature:

::

     F: (A, P) -> (S )
  
      with
          a: horzcat([a0, a1, ..., a_(N-1)])
          p: horzcat([p0, p1, ..., p_(N-1)])
          s: horzcat([s0, s1, ..., s_(N-1)])
      and
          s0 <- f(a0, p0)
          s1 <- f(a1, p1)
          ...
          s_(N-1) <- f(a_(N-1), p_(N-1))
  



>  Function Function.map(str name, int n, [bool ] repeat_in, [bool ] repeat_out, Dict opts=Dict()) const 
------------------------------------------------------------------------

Generic map.

";

%feature("docstring") casadi::Callback::size1_in "

Get input dimension.

";

%feature("docstring") casadi::Callback::sparsity_in "

Get sparsity of a given input.

";

%feature("docstring") casadi::Callback::get_reverse "

Return function that calculates adjoint derivatives derReverse(nadj) returns
a cached instance if available, and calls  Function getDerReverse(int nadj)
if no cached version is available.

";

%feature("docstring") casadi::Callback "

Callback function functionality This class provides a public API to the
FunctionInternal class that can be subclassed by the user, who is then able
to implement the different virtual method. Note that the Function class also
provides a public API to FunctionInternal, but only allows calling, not
being called.

The user is responsible for not deleting this class for the lifetime of the
internal function object.

Joris Gillis, Joel Andersson

>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| ad_weight    | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | derivative c |              |
|              |              |              | alculation.W |              |
|              |              |              | hen there is |              |
|              |              |              | an option of |              |
|              |              |              | either using |              |
|              |              |              | forward or   |              |
|              |              |              | reverse mode |              |
|              |              |              | directional  |              |
|              |              |              | derivatives, |              |
|              |              |              | the          |              |
|              |              |              | condition ad |              |
|              |              |              | _weight*nf<= |              |
|              |              |              | (1-ad_weight |              |
|              |              |              | )*na is used |              |
|              |              |              | where nf and |              |
|              |              |              | na are       |              |
|              |              |              | estimates of |              |
|              |              |              | the number   |              |
|              |              |              | of forward/r |              |
|              |              |              | everse mode  |              |
|              |              |              | directional  |              |
|              |              |              | derivatives  |              |
|              |              |              | needed. By   |              |
|              |              |              | default,     |              |
|              |              |              | ad_weight is |              |
|              |              |              | calculated a |              |
|              |              |              | utomatically |              |
|              |              |              | , but this   |              |
|              |              |              | can be       |              |
|              |              |              | overridden   |              |
|              |              |              | by setting   |              |
|              |              |              | this option. |              |
|              |              |              | In           |              |
|              |              |              | particular,  |              |
|              |              |              | 0 means      |              |
|              |              |              | forcing      |              |
|              |              |              | forward mode |              |
|              |              |              | and 1        |              |
|              |              |              | forcing      |              |
|              |              |              | reverse      |              |
|              |              |              | mode. Leave  |              |
|              |              |              | unset for    |              |
|              |              |              | (class       |              |
|              |              |              | specific)    |              |
|              |              |              | heuristics.  |              |
+--------------+--------------+--------------+--------------+--------------+
| ad_weight_sp | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | sparsity     |              |
|              |              |              | pattern      |              |
|              |              |              | calculation  |              |
|              |              |              | calculation. |              |
|              |              |              | Overrides    |              |
|              |              |              | default      |              |
|              |              |              | behavior.    |              |
|              |              |              | Set to 0 and |              |
|              |              |              | 1 to force   |              |
|              |              |              | forward and  |              |
|              |              |              | reverse mode |              |
|              |              |              | respectively |              |
|              |              |              | . Cf. option |              |
|              |              |              | \"ad_weight\". |              |
+--------------+--------------+--------------+--------------+--------------+
| compiler     | OT_STRING    | \"clang\"      | Just-in-time | casadi::Func |
|              |              |              | compiler     | tionInternal |
|              |              |              | plugin to be |              |
|              |              |              | used.        |              |
+--------------+--------------+--------------+--------------+--------------+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode     |
|              |              |              | according to |              |
|              |              |              | a given      |              |
|              |              |              | recipe (low- |              |
|              |              |              | level)       |              |
+--------------+--------------+--------------+--------------+--------------+
| gather_stats | OT_BOOLEAN   | false        | Flag to      | casadi::Func |
|              |              |              | indicate     | tionInternal |
|              |              |              | whether      |              |
|              |              |              | statistics   |              |
|              |              |              | must be      |              |
|              |              |              | gathered     |              |
+--------------+--------------+--------------+--------------+--------------+
| input_scheme | OT_STRINGVEC | GenericType( | Custom input | casadi::Func |
|              | TOR          | )            | scheme       | tionInternal |
+--------------+--------------+--------------+--------------+--------------+
| inputs_check | OT_BOOLEAN   | true         | Throw        | casadi::Func |
|              |              |              | exceptions   | tionInternal |
|              |              |              | when the     |              |
|              |              |              | numerical    |              |
|              |              |              | values of    |              |
|              |              |              | the inputs   |              |
|              |              |              | don't make   |              |
|              |              |              | sense        |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_penalty  | OT_REAL      | 2            | When         | casadi::Func |
|              |              |              | requested    | tionInternal |
|              |              |              | for a number |              |
|              |              |              | of forward/r |              |
|              |              |              | everse       |              |
|              |              |              | directions,  |              |
|              |              |              | it may be    |              |
|              |              |              | cheaper to   |              |
|              |              |              | compute      |              |
|              |              |              | first the    |              |
|              |              |              | full         |              |
|              |              |              | jacobian and |              |
|              |              |              | then         |              |
|              |              |              | multiply     |              |
|              |              |              | with seeds,  |              |
|              |              |              | rather than  |              |
|              |              |              | obtain the   |              |
|              |              |              | requested    |              |
|              |              |              | directions   |              |
|              |              |              | in a straigh |              |
|              |              |              | tforward     |              |
|              |              |              | manner.      |              |
|              |              |              | Casadi uses  |              |
|              |              |              | a heuristic  |              |
|              |              |              | to decide    |              |
|              |              |              | which is     |              |
|              |              |              | cheaper. A   |              |
|              |              |              | high value   |              |
|              |              |              | of 'jac_pena |              |
|              |              |              | lty' makes   |              |
|              |              |              | it less      |              |
|              |              |              | likely for   |              |
|              |              |              | the heurstic |              |
|              |              |              | to chose the |              |
|              |              |              | full         |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy.    |              |
|              |              |              | The special  |              |
|              |              |              | value -1     |              |
|              |              |              | indicates    |              |
|              |              |              | never to use |              |
|              |              |              | the full     |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy     |              |
+--------------+--------------+--------------+--------------+--------------+
| jit          | OT_BOOLEAN   | false        | Use just-in- | casadi::Func |
|              |              |              | time         | tionInternal |
|              |              |              | compiler to  |              |
|              |              |              | speed up the |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| jit_options  | OT_DICT      | GenericType( | Options to   | casadi::Func |
|              |              | )            | be passed to | tionInternal |
|              |              |              | the jit      |              |
|              |              |              | compiler.    |              |
+--------------+--------------+--------------+--------------+--------------+
| monitor      | OT_STRINGVEC | GenericType( | Monitors to  | casadi::Func |
|              | TOR          | )            | be activated | tionInternal |
|              |              |              | (inputs|outp |              |
|              |              |              | uts)         |              |
+--------------+--------------+--------------+--------------+--------------+
| output_schem | OT_STRINGVEC | GenericType( | Custom       | casadi::Func |
| e            | TOR          | )            | output       | tionInternal |
|              |              |              | scheme       |              |
+--------------+--------------+--------------+--------------+--------------+
| regularity_c | OT_BOOLEAN   | true         | Throw        | casadi::Func |
| heck         |              |              | exceptions   | tionInternal |
|              |              |              | when NaN or  |              |
|              |              |              | Inf appears  |              |
|              |              |              | during       |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| user_data    | OT_VOIDPTR   | GenericType( | A user-      | casadi::Func |
|              |              | )            | defined      | tionInternal |
|              |              |              | field that   |              |
|              |              |              | can be used  |              |
|              |              |              | to identify  |              |
|              |              |              | the function |              |
|              |              |              | or pass      |              |
|              |              |              | additional   |              |
|              |              |              | information  |              |
+--------------+--------------+--------------+--------------+--------------+
| verbose      | OT_BOOLEAN   | false        | Verbose      | casadi::Func |
|              |              |              | evaluation   | tionInternal |
|              |              |              | for          |              |
|              |              |              | debugging    |              |
+--------------+--------------+--------------+--------------+--------------+

Diagrams
--------



C++ includes: callback.hpp ";

%feature("docstring") casadi::Callback::addMonitor "

Add modules to be monitored.

";

%feature("docstring") casadi::Callback::n_out "

Get the number of function outputs.

";

%feature("docstring") casadi::Callback::removeMonitor "

Remove modules to be monitored.

";

%feature("docstring") casadi::Callback::getSanitizedName "

get function name with all non alphanumeric characters converted to '_'

";

%feature("docstring") casadi::Callback::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::Callback::free_mx "

Get all the free variables of the function.

";

%feature("docstring") casadi::Callback::index_in "

Find the index for a string describing a particular entry of an input
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::Callback::sparsity_out "

Get sparsity of a given output.

";

%feature("docstring") casadi::Callback::setJacSparsity "

Generate the sparsity of a Jacobian block

";

%feature("docstring") casadi::Callback::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::Callback::callReverse "

Create call to (cached) derivative function, reverse mode.

";

%feature("docstring") casadi::Callback::size_out "

Get output dimension.

";

%feature("docstring") casadi::Callback::fullJacobian "

Generate a Jacobian function of all the inputs elements with respect to all
the output elements).

";

%feature("docstring") casadi::Callback::derReverse "

Get a function that calculates nadj adjoint derivatives.

Returns a function with n_in + n_out +nadj*n_out inputs and nadj*n_in
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nadj*n_out inputs correspond to adjoint seeds, one direction at a time The
nadj*n_in outputs correspond to adjoint sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

(n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::Callback::getInput "

>  DMatrix  IOInterface< Function  >.getInput(int iind=0) const
------------------------------------------------------------------------

Get an input by index.

Parameters:
-----------

iind:  index within the range [0..n_in()-1]

>  DMatrix  IOInterface< Function  >.getInput(str iname) const
------------------------------------------------------------------------

Get an input by name.

Parameters:
-----------

iname:  input name. Only allowed when an input scheme is set.

>  void IOInterface< Function  >.getInput(T val, int iind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.getInput(T val, str iname)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::Callback::mx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::Callback::setDerForward "

Set a function that calculates nfwd forward derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";


// File: classcasadi_1_1casadi__limits.xml
%feature("docstring") casadi::casadi_limits "

casadi_limits class

The following class, which acts as a complements to the standard
numeric_limits class, allows specifying certain properties of scalar
objects. The template can be specialized for e.g. symbolic scalars Joel
Andersson

C++ includes: casadi_limits.hpp ";


// File: classcasadi_1_1CasadiException.xml
%feature("docstring") casadi::CasadiException::what "throw () Display
error.

";

%feature("docstring") casadi::CasadiException::CasadiException "

>  CasadiException()
------------------------------------------------------------------------

Default constructor.

>  CasadiException(str msg)
------------------------------------------------------------------------

Form message string.

";

%feature("docstring") casadi::CasadiException "

Casadi exception class.

Joel Andersson

C++ includes: casadi_exception.hpp ";

%feature("docstring") casadi::CasadiException::~CasadiException "throw ()
Destructor.

";


// File: classcasadi_1_1CasadiMeta.xml
%feature("docstring") casadi::CasadiMeta "

Collects global CasADi meta information.

Joris Gillis

C++ includes: casadi_meta.hpp ";


// File: classcasadi_1_1CasadiOptions.xml
%feature("docstring") casadi::CasadiOptions "

Collects global CasADi options.

Note to developers: use sparingly. Global options are - in general - a
rather bad idea

this class must never be instantiated. Access its static members directly
Joris Gillis

C++ includes: casadi_options.hpp ";


// File: classcasadi_1_1ClangCompiler.xml


// File: classcasadi_1_1CodeGenerator.xml
%feature("docstring") casadi::CodeGenerator::addInclude "

Add an include file optionally using a relative path \"...\" instead of an
absolute path <...>

";

%feature("docstring") casadi::CodeGenerator "C++ includes:
code_generator.hpp ";

%feature("docstring") casadi::CodeGenerator::compile "

Compile and load function.

";

%feature("docstring") casadi::CodeGenerator::add "

>  void CodeGenerator.add(Function f)
------------------------------------------------------------------------

Add a function (name generated)

>  void CodeGenerator.add(Function f, str fname)
------------------------------------------------------------------------

Add a function.

";

%feature("docstring") casadi::CodeGenerator::CodeGenerator "

Constructor.

";

%feature("docstring") casadi::CodeGenerator::generate "

>  void CodeGenerator.generate(str name) const 
------------------------------------------------------------------------

Generate a file.

>  str CodeGenerator.generate() const 
------------------------------------------------------------------------

Generate a file, return code as string.

";


// File: classcasadi_1_1CollocationIntegrator.xml


// File: classcasadi_1_1CommonExternal.xml


// File: classcasadi_1_1Compiler.xml


/*  Option Functionality  */ %feature("docstring")
casadi::Compiler::plugin_name "

Query plugin name.

";

%feature("docstring") casadi::Compiler::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::Compiler::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::Compiler::isNull "

Is a null pointer?

";

%feature("docstring") casadi::Compiler::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::Compiler::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::Compiler::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::Compiler::getOptionTypeName "

Get the type name of a certain option.

";

%feature("docstring") casadi::Compiler::setOptionByEnumValue "[INTERNAL]
Set a certain option by giving an enum value.

";

%feature("docstring") casadi::Compiler::print "

Print a description of the object.

";

%feature("docstring") casadi::Compiler::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::Compiler "

Compiler.

Just-in-time compilation of code

General information
===================



>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode     |
|              |              |              | according to |              |
|              |              |              | a given      |              |
|              |              |              | recipe (low- |              |
|              |              |              | level)       |              |
+--------------+--------------+--------------+--------------+--------------+

List of plugins
===============



- clang

- shell

Note: some of the plugins in this list might not be available on your
system. Also, there might be extra plugins available to you that are not
listed here. You can obtain their documentation with
Compiler.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

clang
-----



Interface to the JIT compiler CLANG

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| flags           | OT_STRINGVECTOR | GenericType()   | Compile flags   |
|                 |                 |                 | for the JIT     |
|                 |                 |                 | compiler.       |
|                 |                 |                 | Default: None   |
+-----------------+-----------------+-----------------+-----------------+
| include_path    | OT_STRING       | \"\"              | Include paths   |
|                 |                 |                 | for the JIT     |
|                 |                 |                 | compiler. The   |
|                 |                 |                 | include         |
|                 |                 |                 | directory       |
|                 |                 |                 | shipped with    |
|                 |                 |                 | CasADi will be  |
|                 |                 |                 | automatically   |
|                 |                 |                 | appended.       |
+-----------------+-----------------+-----------------+-----------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

shell
-----



Interface to the JIT compiler SHELL

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| compiler        | OT_STRING       | \"gcc\"           | Compiler        |
|                 |                 |                 | command         |
+-----------------+-----------------+-----------------+-----------------+
| compiler_setup  | OT_STRING       | \"-fPIC -shared\" | Compiler setup  |
|                 |                 |                 | command.        |
|                 |                 |                 | Intended to be  |
|                 |                 |                 | fixed. The      |
|                 |                 |                 | 'flag' option   |
|                 |                 |                 | is the prefered |
|                 |                 |                 | way to set      |
|                 |                 |                 | custom flags.   |
+-----------------+-----------------+-----------------+-----------------+
| flags           | OT_STRINGVECTOR | GenericType()   | Compile flags   |
|                 |                 |                 | for the JIT     |
|                 |                 |                 | compiler.       |
|                 |                 |                 | Default: None   |
+-----------------+-----------------+-----------------+-----------------+

--------------------------------------------------------------------------------



Joris Gillis
Diagrams
--------



C++ includes: compiler.hpp ";

%feature("docstring") casadi::Compiler::getOptionAllowedIndex "[INTERNAL]
Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::Compiler::setOptionByAllowedIndex "[INTERNAL]
Set a certain option by giving its index into the allowed values.

";

%feature("docstring") casadi::Compiler::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::Compiler::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::Compiler::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::Compiler::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::Compiler::getOptionEnumValue "[INTERNAL]  Get
the enum value corresponding to th certain option.

";

%feature("docstring") casadi::Compiler::repr "

Print a representation of the object.

";

%feature("docstring") casadi::Compiler::Compiler "

>  Compiler()
------------------------------------------------------------------------

Default constructor.

>  Compiler(str name, str compiler, Dict opts=Dict())
------------------------------------------------------------------------

Compiler factory (new syntax, includes initialization)

";

%feature("docstring") casadi::Compiler::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::Compiler::dictionary "

Get the dictionary.

";


// File: classcasadi_1_1Concat.xml


// File: classcasadi_1_1Constant.xml


// File: classcasadi_1_1ConstantDMatrix.xml


// File: classcasadi_1_1ConstantMX.xml


// File: classcasadi_1_1ConstantSX.xml


// File: classcasadi_1_1DaeBuilder.xml


/*  Variables and equations  */

/* Public data members

*/

/*  Symbolic modeling  */

/* Formulate an optimal control problem

*/

/*  Manipulation  */

/* Reformulate the dynamic optimization problem.

*/

/*  Import and export  */ %feature("docstring") casadi::DaeBuilder::setMin "

>  void DaeBuilder.setMin(str name, double val, bool normalized=false)
------------------------------------------------------------------------

Set the lower bound by name.

>  void DaeBuilder.setMin(MX var, [double ] val, bool normalized=false)
------------------------------------------------------------------------

Set the lower bound(s) by expression.

";

%feature("docstring") casadi::DaeBuilder::add_s "

Add a implicit state.

";

%feature("docstring") casadi::DaeBuilder::derivativeStart "

>  double DaeBuilder.derivativeStart(str name, bool normalized=false) const 
------------------------------------------------------------------------

Get the (optionally normalized) derivative value at time 0 by name.

>  [double] DaeBuilder.derivativeStart(MX var, bool normalized=false) const 
------------------------------------------------------------------------

Get the (optionally normalized) derivative value(s) at time 0 by expression.

";

%feature("docstring") casadi::DaeBuilder::sanity_check "

Check if dimensions match.

";

%feature("docstring") casadi::DaeBuilder::DaeBuilder "

Default constructor.

";

%feature("docstring") casadi::DaeBuilder::sort_dae "

Sort the DAE and implicitly defined states.

";

%feature("docstring") casadi::DaeBuilder::setStart "

>  void DaeBuilder.setStart(str name, double val, bool normalized=false)
------------------------------------------------------------------------

Set the (optionally normalized) value at time 0 by name.

>  void DaeBuilder.setStart(MX var, [double ] val, bool normalized=false)
------------------------------------------------------------------------

Set the (optionally normalized) value(s) at time 0 by expression.

";

%feature("docstring") casadi::DaeBuilder::min "

>  double DaeBuilder.min(str name, bool normalized=false) const 
------------------------------------------------------------------------

Get the lower bound by name.

>  [double] DaeBuilder.min(MX var, bool normalized=false) const 
------------------------------------------------------------------------

Get the lower bound(s) by expression.

";

%feature("docstring") casadi::DaeBuilder::add_quad "

Add a quadrature equation.

";

%feature("docstring") casadi::DaeBuilder::nominal "

>  double DaeBuilder.nominal(str name) const 
------------------------------------------------------------------------

Get the nominal value by name.

>  [double] DaeBuilder.nominal(MX var) const 
------------------------------------------------------------------------

Get the nominal value(s) by expression.

";

%feature("docstring") casadi::DaeBuilder::max "

>  double DaeBuilder.max(str name, bool normalized=false) const 
------------------------------------------------------------------------

Get the upper bound by name.

>  [double] DaeBuilder.max(MX var, bool normalized=false) const 
------------------------------------------------------------------------

Get the upper bound(s) by expression.

";

%feature("docstring") casadi::DaeBuilder::variable "

Access a variable by name

";

%feature("docstring") casadi::DaeBuilder::initialGuess "

>  double DaeBuilder.initialGuess(str name, bool normalized=false) const 
------------------------------------------------------------------------

Get the initial guess by name.

>  [double] DaeBuilder.initialGuess(MX var, bool normalized=false) const 
------------------------------------------------------------------------

Get the initial guess(es) by expression.

";

%feature("docstring") casadi::DaeBuilder::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::DaeBuilder::makeSemiExplicit "

Transform the implicit DAE to a semi-explicit DAE.

";

%feature("docstring") casadi::DaeBuilder::sort_d "

Sort dependent parameters.

";

%feature("docstring") casadi::DaeBuilder::print "

Print description.

";

%feature("docstring") casadi::DaeBuilder::sort_alg "

Sort the algebraic equations and algebraic states.

";

%feature("docstring") casadi::DaeBuilder::add_dae "

Add a differential-algebraic equation.

";

%feature("docstring") casadi::DaeBuilder::addVariable "

>  void DaeBuilder.addVariable(str name, Variable var)
------------------------------------------------------------------------

Add a variable.

>  MX DaeBuilder.addVariable(str name, int n=1)

>  MX DaeBuilder.addVariable(str name, Sparsity sp)
------------------------------------------------------------------------

Add a new variable: returns corresponding symbolic expression.

";

%feature("docstring") casadi::DaeBuilder::setMax "

>  void DaeBuilder.setMax(str name, double val, bool normalized=false)
------------------------------------------------------------------------

Set the upper bound by name.

>  void DaeBuilder.setMax(MX var, [double ] val, bool normalized=false)
------------------------------------------------------------------------

Set the upper bound(s) by expression.

";

%feature("docstring") casadi::DaeBuilder "

An initial-value problem in differential-algebraic equations.

Independent variables:
======================





::

  t:      time
  



Time-continuous variables:
==========================





::

  x:      states defined by ODE
  s:      implicitly defined states
  z:      algebraic variables
  u:      control signals
  q:      quadrature states
  y:      outputs
  



Time-constant variables:
========================





::

  p:      free parameters
  d:      dependent parameters
  



Dynamic constraints (imposed everywhere):
=========================================





::

  ODE                    \\\\dot{x} ==  ode(t, x, s, z, u, p, d)
  DAE or implicit ODE:         0 ==  dae(t, x, s, z, u, p, d, sdot)
  algebraic equations:         0 ==  alg(t, x, s, z, u, p, d)
  quadrature equations:  \\\\dot{q} == quad(t, x, s, z, u, p, d)
  deppendent parameters:       d == ddef(t, x, s, z, u, p, d)
  output equations:            y == ydef(t, x, s, z, u, p, d)
  



Point constraints (imposed pointwise):
======================================





::

  Initial equations:           0 == init(t, x, s, z, u, p, d, sdot)
  



Joel Andersson

C++ includes: dae_builder.hpp ";

%feature("docstring") casadi::DaeBuilder::parseFMI "

Import existing problem from FMI/XML

";

%feature("docstring") casadi::DaeBuilder::add_q "

Add a new quadrature state.

";

%feature("docstring") casadi::DaeBuilder::add_p "

Add a new parameter

";

%feature("docstring") casadi::DaeBuilder::add_u "

Add a new control.

";

%feature("docstring") casadi::DaeBuilder::setDerivativeStart "

>  void DaeBuilder.setDerivativeStart(str name, double val, bool normalized=false)
------------------------------------------------------------------------

Set the (optionally normalized) derivative value at time 0 by name.

>  void DaeBuilder.setDerivativeStart(MX var, [double ] val, bool normalized=false)
------------------------------------------------------------------------

Set the (optionally normalized) derivative value(s) at time 0 by expression.

";

%feature("docstring") casadi::DaeBuilder::add_z "

Add a new algebraic variable.

";

%feature("docstring") casadi::DaeBuilder::add_y "

Add a new output.

";

%feature("docstring") casadi::DaeBuilder::add_x "

Add a new differential state.

";

%feature("docstring") casadi::DaeBuilder::eliminate_quad "

Eliminate quadrature states and turn them into ODE states.

";

%feature("docstring") casadi::DaeBuilder::makeExplicit "

Transform the implicit DAE or semi-explicit DAE into an explicit ODE.

";

%feature("docstring") casadi::DaeBuilder::add_d "

Add a new dependent parameter.

";

%feature("docstring") casadi::DaeBuilder::split_dae "

Identify and separate the algebraic variables and equations in the DAE.

";

%feature("docstring") casadi::DaeBuilder::addLinearCombination "

Add a named linear combination of output expressions.

";

%feature("docstring") casadi::DaeBuilder::add_alg "

Add an algebraic equation.

";

%feature("docstring") casadi::DaeBuilder::scaleEquations "

Scale the implicit equations.

";

%feature("docstring") casadi::DaeBuilder::repr "

Print representation.

";

%feature("docstring") casadi::DaeBuilder::create "

Construct a function object.

";

%feature("docstring") casadi::DaeBuilder::scaleVariables "

Scale the variables.

";

%feature("docstring") casadi::DaeBuilder::setInitialGuess "

>  void DaeBuilder.setInitialGuess(str name, double val, bool normalized=false)
------------------------------------------------------------------------

Set the initial guess by name.

>  void DaeBuilder.setInitialGuess(MX var, [double ] val, bool normalized=false)
------------------------------------------------------------------------

Set the initial guess(es) by expression.

";

%feature("docstring") casadi::DaeBuilder::eliminate_alg "

Eliminate algebraic variables and equations transforming them into outputs.

";

%feature("docstring") casadi::DaeBuilder::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::DaeBuilder::unit "

>  str DaeBuilder.unit(str name) const 
------------------------------------------------------------------------

Get the unit for a component.

>  str DaeBuilder.unit(MX var) const 
------------------------------------------------------------------------

Get the unit given a vector of symbolic variables (all units must be
identical)

";

%feature("docstring") casadi::DaeBuilder::start "

>  double DaeBuilder.start(str name, bool normalized=false) const 
------------------------------------------------------------------------

Get the (optionally normalized) value at time 0 by name.

>  [double] DaeBuilder.start(MX var, bool normalized=false) const 
------------------------------------------------------------------------

Get the (optionally normalized) value(s) at time 0 by expression.

";

%feature("docstring") casadi::DaeBuilder::split_d "

Eliminate interdependencies amongst dependent parameters.

";

%feature("docstring") casadi::DaeBuilder::eliminate_d "

Eliminate dependent parameters.

";

%feature("docstring") casadi::DaeBuilder::der "

>  MX DaeBuilder.der(str name) const 
------------------------------------------------------------------------

Get a derivative expression by name.

>  MX DaeBuilder.der(MX var) const 
------------------------------------------------------------------------

Get a derivative expression by non-differentiated expression.

";

%feature("docstring") casadi::DaeBuilder::setNominal "

>  void DaeBuilder.setNominal(str name, double val)
------------------------------------------------------------------------

Set the nominal value by name.

>  void DaeBuilder.setNominal(MX var, [double ] val)
------------------------------------------------------------------------

Set the nominal value(s) by expression.

";

%feature("docstring") casadi::DaeBuilder::setUnit "

Set the unit for a component.

";

%feature("docstring") casadi::DaeBuilder::add_ode "

Add an ordinary differential equation.

";


// File: classcasadi_1_1DenseMultiplication.xml


// File: classcasadi_1_1DenseTranspose.xml


// File: classcasadi_1_1Determinant.xml


// File: classcasadi_1_1Diagcat.xml


// File: classcasadi_1_1Diagsplit.xml


// File: classcasadi_1_1Find.xml


// File: classcasadi_1_1FixedStepIntegrator.xml


// File: classcasadi_1_1Function.xml


/*  Simple Getters & Setters  */

/*  Advanced Getters  */

/*  Option Functionality  */ %feature("docstring") casadi::Function::repr "

Print a representation of the object.

";

%feature("docstring") casadi::Function::gradient "

Generate a gradient function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the output must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::Function::getOptionTypeName "

Get the type name of a certain option.

";

%feature("docstring") casadi::Function::nnz_in "

Get of number of input nonzeros For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::Function::sz_res "[INTERNAL]  Get required
length of res field.

";

%feature("docstring") casadi::Function::name "

Name of the function.

";

%feature("docstring") casadi::Function::checkInputs "[INTERNAL]  Check if
the numerical values of the supplied bounds make sense.

";

%feature("docstring") casadi::Function::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::Function::jacobian "

Generate a Jacobian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. If
compact is set to true, only the nonzeros of the input and output
expressions are considered. If symmetric is set to true, the Jacobian being
calculated is known to be symmetric (usually a Hessian), which can be
exploited by the algorithm.

The generated Jacobian has one more output than the calling function
corresponding to the Jacobian and the same number of inputs.

";

%feature("docstring") casadi::Function::size1_in "

Get input dimension.

";

%feature("docstring") casadi::Function::setFullJacobian "

Set the Jacobian of all the input nonzeros with respect to all output
nonzeros NOTE: Does not take ownership, only weak references to the Jacobian
are kept internally

";

%feature("docstring") casadi::Function::sparsity_out "

Get sparsity of a given output.

";

%feature("docstring") casadi::Function::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::Function::setOptionByAllowedIndex "[INTERNAL]
Set a certain option by giving its index into the allowed values.

";

%feature("docstring") casadi::Function::spEvaluate "[INTERNAL]  Propagate
the sparsity pattern through a set of directional.

derivatives forward or backward (for usage, see the example
propagating_sparsity.cpp)

";

%feature("docstring") casadi::Function::derForward "

Get a function that calculates nfwd forward derivatives.

Returns a function with n_in + n_out +nfwd*n_in inputs and nfwd*n_out
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nfwd*n_in inputs correspond to forward seeds, one direction at a time The
nfwd*n_out outputs correspond to forward sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::Function::getAtomicOperation "

Get an atomic operation operator index.

";

%feature("docstring") casadi::Function::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::Function::getAlgorithmSize "

Get the number of atomic operations.

";

%feature("docstring") casadi::Function::getOutput "

>  DMatrix  IOInterface< Function  >.getOutput(int oind=0) const
------------------------------------------------------------------------

Get an output by index.

Parameters:
-----------

oind:  index within the range [0..n_out()-1]

>  DMatrix  IOInterface< Function  >.getOutput(str oname) const
------------------------------------------------------------------------

Get an output by name.

Parameters:
-----------

oname:  output name. Only allowed when an output scheme is set.

>  void IOInterface< Function  >.getOutput(T val, int oind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.getOutput(T val, str oname)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::Function::callReverse "

Create call to (cached) derivative function, reverse mode.

";

%feature("docstring") casadi::Function::numel_out "

Get of number of output elements For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::Function::getSanitizedName "

get function name with all non alphanumeric characters converted to '_'

";

%feature("docstring") casadi::Function::derReverse "

Get a function that calculates nadj adjoint derivatives.

Returns a function with n_in + n_out +nadj*n_out inputs and nadj*n_in
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nadj*n_out inputs correspond to adjoint seeds, one direction at a time The
nadj*n_in outputs correspond to adjoint sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

(n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::Function::printDimensions "

Print dimensions of inputs and outputs.

";

%feature("docstring") casadi::Function::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::Function::name_out "

>  [str] Function.name_out() const 
------------------------------------------------------------------------

Get output scheme.

>  str Function.name_out(int ind) const 
------------------------------------------------------------------------

Get output scheme name by index.

";

%feature("docstring") casadi::Function::setDerReverse "

Set a function that calculates nadj adjoint derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::Function::n_in "

Get the number of function inputs.

";

%feature("docstring") casadi::Function::mx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::Function::countNodes "

Number of nodes in the algorithm.

";

%feature("docstring") casadi::Function::free_sx "

Get all the free variables of the function.

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
  



";

%feature("docstring") casadi::Function::derivative "

Get a function that calculates nfwd forward derivatives and nadj adjoint
derivatives Legacy function: Use derForward and derReverse instead.

Returns a function with (1+nfwd)*n_in+nadj*n_out inputs and (1+nfwd)*n_out +
nadj*n_in outputs. The first n_in inputs correspond to nondifferentiated
inputs. The next nfwd*n_in inputs correspond to forward seeds, one direction
at a time and the last nadj*n_out inputs correspond to adjoint seeds, one
direction at a time. The first n_out outputs correspond to nondifferentiated
outputs. The next nfwd*n_out outputs correspond to forward sensitivities,
one direction at a time and the last nadj*n_in outputs corresponds to
adjoint sensitivities, one direction at a time.

(n_in = n_in(), n_out = n_out())

";

%feature("docstring") casadi::Function::Function "

default constructor

";

%feature("docstring") casadi::Function::setDerForward "

Set a function that calculates nfwd forward derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::Function::type_name "

Get type name.

";

%feature("docstring") casadi::Function::numel_in "

Get of number of input elements For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::Function::n_out "

Get the number of function outputs.

";

%feature("docstring") casadi::Function::getAtomicInput "

Get the (integer) input arguments of an atomic operation.

";

%feature("docstring") casadi::Function::spInit "[INTERNAL]  Reset the
sparsity propagation.

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::Function::index_in "

Find the index for a string describing a particular entry of an input
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::Function::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::Function::~Function "

Destructor.

";

%feature("docstring") casadi::Function::index_out "

Find the index for a string describing a particular entry of an output
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::Function::free_mx "

Get all the free variables of the function.

";

%feature("docstring") casadi::Function::size_out "

Get output dimension.

";

%feature("docstring") casadi::Function::addMonitor "

Add modules to be monitored.

";

%feature("docstring") casadi::Function::description_out "

Get output scheme description by index.

";

%feature("docstring") casadi::Function::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

";

%feature("docstring") casadi::Function::setOutput "

>  void IOInterface< Function  >.setOutput(T val, int oind=0)
------------------------------------------------------------------------

Set an output by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.setOutput(T val, str oname)
------------------------------------------------------------------------

Set an output by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::Function::getStats "

Get all statistics obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::Function::print "

Print a description of the object.

";

%feature("docstring") casadi::Function::setJacSparsity "

Generate the sparsity of a Jacobian block

";

%feature("docstring") casadi::Function::hessian "

Generate a Hessian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The generated Hessian has two more outputs than the calling function
corresponding to the Hessian and the gradients.

";

%feature("docstring") casadi::Function::isNull "

Is a null pointer?

";

%feature("docstring") casadi::Function::getOptionEnumValue "[INTERNAL]  Get
the enum value corresponding to th certain option.

";

%feature("docstring") casadi::Function::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::Function::fullJacobian "

Generate a Jacobian function of all the inputs elements with respect to all
the output elements).

";

%feature("docstring") casadi::Function::map "

>  [[MX] ] Function.map([[MX ] ] arg, str parallelization=\"serial\")

>  [MX] Function.map([MX ] arg, str parallelization=\"serial\")
------------------------------------------------------------------------

Evaluate symbolically in parallel (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

>  Function Function.map(str name, int N, Dict options=Dict()) const 
------------------------------------------------------------------------

Create a mapped version of this function.

Suppose the function has a signature of:

::

     f: (a, p) -> ( s )
  



The the mapaccumulated version has the signature:

::

     F: (A, P) -> (S )
  
      with
          a: horzcat([a0, a1, ..., a_(N-1)])
          p: horzcat([p0, p1, ..., p_(N-1)])
          s: horzcat([s0, s1, ..., s_(N-1)])
      and
          s0 <- f(a0, p0)
          s1 <- f(a1, p1)
          ...
          s_(N-1) <- f(a_(N-1), p_(N-1))
  



>  Function Function.map(str name, int n, [bool ] repeat_in, [bool ] repeat_out, Dict opts=Dict()) const 
------------------------------------------------------------------------

Generic map.

";

%feature("docstring") casadi::Function::sx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::Function::evaluate "

Evaluate.

";

%feature("docstring") casadi::Function::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::Function::getStat "

Get a single statistic obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::Function::jacSparsity "

Get, if necessary generate, the sparsity of a Jacobian block

";

%feature("docstring") casadi::Function::description_in "

Get input scheme description by index.

";

%feature("docstring") casadi::Function::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::Function::getOptionAllowedIndex "[INTERNAL]
Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::Function::getInput "

>  DMatrix  IOInterface< Function  >.getInput(int iind=0) const
------------------------------------------------------------------------

Get an input by index.

Parameters:
-----------

iind:  index within the range [0..n_in()-1]

>  DMatrix  IOInterface< Function  >.getInput(str iname) const
------------------------------------------------------------------------

Get an input by name.

Parameters:
-----------

iname:  input name. Only allowed when an input scheme is set.

>  void IOInterface< Function  >.getInput(T val, int iind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.getInput(T val, str iname)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::Function::sz_iw "[INTERNAL]  Get required
length of iw field.

";

%feature("docstring") casadi::Function::is_a "

Check if the function is of a particular type Optionally check if name
matches one of the base classes (default true)

";

%feature("docstring") casadi::Function::name_in "

>  [str] Function.name_in() const 
------------------------------------------------------------------------

Get input scheme.

>  str Function.name_in(int ind) const 
------------------------------------------------------------------------

Get input scheme name by index.

";

%feature("docstring") casadi::Function::mx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::Function "

General function.

A general function $f$ in casadi can be multi-input, multi-output. Number of
inputs: nin n_in() Number of outputs: nout n_out()  We can view this
function as a being composed of a ( nin, nout) grid of single-input, single-
output primitive functions. Each such primitive function $f_ {i, j}
\\\\forall i \\\\in [0, nin-1], j \\\\in [0, nout-1]$ can map as $\\\\mathbf
{R}^{n, m}\\\\to\\\\mathbf{R}^{p, q}$, in which n, m, p, q can take
different values for every (i, j) pair.  When passing input, you specify
which partition $i$ is active. You pass the numbers vectorized, as a vector
of size $(n*m)$. When requesting output, you specify which partition $j$ is
active. You get the numbers vectorized, as a vector of size $(p*q)$.  To
calculate Jacobians, you need to have $(m=1, q=1)$.

Write the Jacobian as $J_ {i, j} = \\\\nabla f_{i, j} = \\\\frac
{\\\\partial f_{i, j}(\\\\vec{x})}{\\\\partial \\\\vec{x}}$.

We have the following relationships for function mapping from a row vector
to a row vector:

$ \\\\vec {s}_f = \\\\nabla f_{i, j} . \\\\vec{v}$ $ \\\\vec {s}_a =
(\\\\nabla f_{i, j})^T . \\\\vec{w}$

Some quantities in these formulas must be transposed: input col: transpose $
\\\\vec {v} $ and $\\\\vec{s}_a$ output col: transpose $ \\\\vec {w} $ and
$\\\\vec{s}_f$  NOTE: Functions are allowed to modify their input arguments
when evaluating: implicitFunction, IDAS solver Further releases may disallow
this.

Joel Andersson

>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| ad_weight    | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | derivative c |              |
|              |              |              | alculation.W |              |
|              |              |              | hen there is |              |
|              |              |              | an option of |              |
|              |              |              | either using |              |
|              |              |              | forward or   |              |
|              |              |              | reverse mode |              |
|              |              |              | directional  |              |
|              |              |              | derivatives, |              |
|              |              |              | the          |              |
|              |              |              | condition ad |              |
|              |              |              | _weight*nf<= |              |
|              |              |              | (1-ad_weight |              |
|              |              |              | )*na is used |              |
|              |              |              | where nf and |              |
|              |              |              | na are       |              |
|              |              |              | estimates of |              |
|              |              |              | the number   |              |
|              |              |              | of forward/r |              |
|              |              |              | everse mode  |              |
|              |              |              | directional  |              |
|              |              |              | derivatives  |              |
|              |              |              | needed. By   |              |
|              |              |              | default,     |              |
|              |              |              | ad_weight is |              |
|              |              |              | calculated a |              |
|              |              |              | utomatically |              |
|              |              |              | , but this   |              |
|              |              |              | can be       |              |
|              |              |              | overridden   |              |
|              |              |              | by setting   |              |
|              |              |              | this option. |              |
|              |              |              | In           |              |
|              |              |              | particular,  |              |
|              |              |              | 0 means      |              |
|              |              |              | forcing      |              |
|              |              |              | forward mode |              |
|              |              |              | and 1        |              |
|              |              |              | forcing      |              |
|              |              |              | reverse      |              |
|              |              |              | mode. Leave  |              |
|              |              |              | unset for    |              |
|              |              |              | (class       |              |
|              |              |              | specific)    |              |
|              |              |              | heuristics.  |              |
+--------------+--------------+--------------+--------------+--------------+
| ad_weight_sp | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | sparsity     |              |
|              |              |              | pattern      |              |
|              |              |              | calculation  |              |
|              |              |              | calculation. |              |
|              |              |              | Overrides    |              |
|              |              |              | default      |              |
|              |              |              | behavior.    |              |
|              |              |              | Set to 0 and |              |
|              |              |              | 1 to force   |              |
|              |              |              | forward and  |              |
|              |              |              | reverse mode |              |
|              |              |              | respectively |              |
|              |              |              | . Cf. option |              |
|              |              |              | \"ad_weight\". |              |
+--------------+--------------+--------------+--------------+--------------+
| compiler     | OT_STRING    | \"clang\"      | Just-in-time | casadi::Func |
|              |              |              | compiler     | tionInternal |
|              |              |              | plugin to be |              |
|              |              |              | used.        |              |
+--------------+--------------+--------------+--------------+--------------+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode     |
|              |              |              | according to |              |
|              |              |              | a given      |              |
|              |              |              | recipe (low- |              |
|              |              |              | level)       |              |
+--------------+--------------+--------------+--------------+--------------+
| gather_stats | OT_BOOLEAN   | false        | Flag to      | casadi::Func |
|              |              |              | indicate     | tionInternal |
|              |              |              | whether      |              |
|              |              |              | statistics   |              |
|              |              |              | must be      |              |
|              |              |              | gathered     |              |
+--------------+--------------+--------------+--------------+--------------+
| input_scheme | OT_STRINGVEC | GenericType( | Custom input | casadi::Func |
|              | TOR          | )            | scheme       | tionInternal |
+--------------+--------------+--------------+--------------+--------------+
| inputs_check | OT_BOOLEAN   | true         | Throw        | casadi::Func |
|              |              |              | exceptions   | tionInternal |
|              |              |              | when the     |              |
|              |              |              | numerical    |              |
|              |              |              | values of    |              |
|              |              |              | the inputs   |              |
|              |              |              | don't make   |              |
|              |              |              | sense        |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_penalty  | OT_REAL      | 2            | When         | casadi::Func |
|              |              |              | requested    | tionInternal |
|              |              |              | for a number |              |
|              |              |              | of forward/r |              |
|              |              |              | everse       |              |
|              |              |              | directions,  |              |
|              |              |              | it may be    |              |
|              |              |              | cheaper to   |              |
|              |              |              | compute      |              |
|              |              |              | first the    |              |
|              |              |              | full         |              |
|              |              |              | jacobian and |              |
|              |              |              | then         |              |
|              |              |              | multiply     |              |
|              |              |              | with seeds,  |              |
|              |              |              | rather than  |              |
|              |              |              | obtain the   |              |
|              |              |              | requested    |              |
|              |              |              | directions   |              |
|              |              |              | in a straigh |              |
|              |              |              | tforward     |              |
|              |              |              | manner.      |              |
|              |              |              | Casadi uses  |              |
|              |              |              | a heuristic  |              |
|              |              |              | to decide    |              |
|              |              |              | which is     |              |
|              |              |              | cheaper. A   |              |
|              |              |              | high value   |              |
|              |              |              | of 'jac_pena |              |
|              |              |              | lty' makes   |              |
|              |              |              | it less      |              |
|              |              |              | likely for   |              |
|              |              |              | the heurstic |              |
|              |              |              | to chose the |              |
|              |              |              | full         |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy.    |              |
|              |              |              | The special  |              |
|              |              |              | value -1     |              |
|              |              |              | indicates    |              |
|              |              |              | never to use |              |
|              |              |              | the full     |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy     |              |
+--------------+--------------+--------------+--------------+--------------+
| jit          | OT_BOOLEAN   | false        | Use just-in- | casadi::Func |
|              |              |              | time         | tionInternal |
|              |              |              | compiler to  |              |
|              |              |              | speed up the |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| jit_options  | OT_DICT      | GenericType( | Options to   | casadi::Func |
|              |              | )            | be passed to | tionInternal |
|              |              |              | the jit      |              |
|              |              |              | compiler.    |              |
+--------------+--------------+--------------+--------------+--------------+
| monitor      | OT_STRINGVEC | GenericType( | Monitors to  | casadi::Func |
|              | TOR          | )            | be activated | tionInternal |
|              |              |              | (inputs|outp |              |
|              |              |              | uts)         |              |
+--------------+--------------+--------------+--------------+--------------+
| output_schem | OT_STRINGVEC | GenericType( | Custom       | casadi::Func |
| e            | TOR          | )            | output       | tionInternal |
|              |              |              | scheme       |              |
+--------------+--------------+--------------+--------------+--------------+
| regularity_c | OT_BOOLEAN   | true         | Throw        | casadi::Func |
| heck         |              |              | exceptions   | tionInternal |
|              |              |              | when NaN or  |              |
|              |              |              | Inf appears  |              |
|              |              |              | during       |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| user_data    | OT_VOIDPTR   | GenericType( | A user-      | casadi::Func |
|              |              | )            | defined      | tionInternal |
|              |              |              | field that   |              |
|              |              |              | can be used  |              |
|              |              |              | to identify  |              |
|              |              |              | the function |              |
|              |              |              | or pass      |              |
|              |              |              | additional   |              |
|              |              |              | information  |              |
+--------------+--------------+--------------+--------------+--------------+
| verbose      | OT_BOOLEAN   | false        | Verbose      | casadi::Func |
|              |              |              | evaluation   | tionInternal |
|              |              |              | for          |              |
|              |              |              | debugging    |              |
+--------------+--------------+--------------+--------------+--------------+

>List of available monitors

+---------+--------------------------+
|   Id    |         Used in          |
+=========+==========================+
| inputs  | casadi::FunctionInternal |
+---------+--------------------------+
| outputs | casadi::FunctionInternal |
+---------+--------------------------+

Diagrams
--------



C++ includes: function.hpp ";

%feature("docstring") casadi::Function::setOptionByEnumValue "[INTERNAL]
Set a certain option by giving an enum value.

";

%feature("docstring") casadi::Function::tangent "

Generate a tangent function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the input must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::Function::getAtomicInputReal "

Get the floating point output argument of an atomic operation.

";

%feature("docstring") casadi::Function::setJacobian "

Set the Jacobian function of output oind with respect to input iind NOTE:
Does not take ownership, only weak references to the Jacobians are kept
internally

";

%feature("docstring") casadi::Function::generate "

Export / Generate C code for the function.

";

%feature("docstring") casadi::Function::getAtomicOutput "

Get the (integer) output argument of an atomic operation.

";

%feature("docstring") casadi::Function::nnz_out "

Get of number of output nonzeros For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::Function::dictionary "

Get the dictionary.

";

%feature("docstring") casadi::Function::sz_w "[INTERNAL]  Get required
length of w field.

";

%feature("docstring") casadi::Function::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::Function::size_in "

Get input dimension.

";

%feature("docstring") casadi::Function::callDerivative "[INTERNAL]
Evaluate the function symbolically or numerically with directional
derivatives The first two arguments are the nondifferentiated inputs and
results of the evaluation, the next two arguments are a set of forward
directional seeds and the resulting forward directional derivatives, the
length of the vector being the number of forward directions. The next two
arguments are a set of adjoint directional seeds and the resulting adjoint
directional derivatives, the length of the vector being the number of
adjoint directions.

";

%feature("docstring") casadi::Function::call "

Evaluate the function symbolically or numerically.

";

%feature("docstring") casadi::Function::size2_out "

Get output dimension.

";

%feature("docstring") casadi::Function::setInput "

>  void IOInterface< Function  >.setInput(T val, int iind=0)
------------------------------------------------------------------------

Set an input by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.setInput(T val, str iname)
------------------------------------------------------------------------

Set an input by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::Function::callForward "

Create call to (cached) derivative function, forward mode.

";

%feature("docstring") casadi::Function::size1_out "

Get output dimension.

";

%feature("docstring") casadi::Function::sparsity_in "

Get sparsity of a given input.

";

%feature("docstring") casadi::Function::generateLiftingFunctions "

Extract the functions needed for the Lifted Newton method.

";

%feature("docstring") casadi::Function::getWorkSize "

Get the length of the work vector.

";

%feature("docstring") casadi::Function::sz_arg "[INTERNAL]  Get required
length of arg field.

";

%feature("docstring") casadi::Function::size2_in "

Get input dimension.

";

%feature("docstring") casadi::Function::removeMonitor "

Remove modules to be monitored.

";

%feature("docstring") casadi::Function::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::Function::spCanEvaluate "[INTERNAL]  Is the
class able to propagate seeds through the algorithm?

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::Function::sx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::Function::default_in "

Get default input value (NOTE: constant reference)

";

%feature("docstring") casadi::Function::getDescription "

Return a string with a description (for SWIG)

";


// File: classcasadi_1_1GenericCall.xml


// File: classcasadi_1_1GenericExpression.xml
%feature("docstring") friendwrap_floor "

Round down to nearest integer.

";

%feature("docstring") friendwrap_acos "

Arc cosine.

";

%feature("docstring") friendwrap_copysign "

Copy sign.

";

%feature("docstring") friendwrap_exp "

Exponential function.

";

%feature("docstring") friendwrap_ceil "

Round up to nearest integer.

";

%feature("docstring") friendwrap_cos "

Cosine.

";

%feature("docstring") friendwrap_asinh "

Inverse hyperbolic sine.

";

%feature("docstring") friendwrap_atanh "

Inverse hyperbolic tangent.

";

%feature("docstring") friendwrap_iszero "

Addition.

";

%feature("docstring") friendwrap_tan "

Tangent.

";

%feature("docstring") friendwrap_acosh "

Inverse hyperbolic cosine.

";

%feature("docstring") friendwrap_erfinv "

Invers error function.

";

%feature("docstring") friendwrap_fmod "

Remainder after division.

";

%feature("docstring") friendwrap_log "

Natural logarithm.

";

%feature("docstring") friendwrap_log10 "

Base-10 logarithm.

";

%feature("docstring") friendwrap_constpow "

Elementwise power with const power.

";

%feature("docstring") friendwrap_abs "

Absolute value.

";

%feature("docstring") friendwrap_fmax "

Largest of two values.

";

%feature("docstring") friendwrap_sqrt "

Square root.

";

%feature("docstring") friendwrap_sign "

Sine function sign(x) := -1 for x<0 sign(x) := 1 for x>0, sign(0) := 0
sign(NaN) := NaN

";

%feature("docstring") friendwrap_logic_and "

Logical and, alternative syntax.

";

%feature("docstring") friendwrap_fmin "

Smallest of two values.

";

%feature("docstring") friendwrap_erf "

Error function.

";

%feature("docstring") friendwrap_pow "

Elementwise power.

";

%feature("docstring") friendwrap_atan2 "

Two argument arc tangent.

";

%feature("docstring") friendwrap_logic_or "

Logical or, alterntive syntax.

";

%feature("docstring") friendwrap_fabs "

Absolute value.

";

%feature("docstring") friendwrap_simplify "

Simplify an expression.

";

%feature("docstring") friendwrap_sinh "

Hyperbolic sine.

";

%feature("docstring") friendwrap_tanh "

Hyperbolic tangent.

";

%feature("docstring") friendwrap_cosh "

Hyperbolic cosine.

";

%feature("docstring") friendwrap_logic_not "

Logical not, alternative syntax.

";

%feature("docstring") friendwrap_atan "

Arc tangent.

";

%feature("docstring") casadi::GenericExpression "

Expression interface.

This is a common base class for SX, MX and Matrix<>, introducing a uniform
syntax and implementing common functionality using the curiously recurring
template pattern (CRTP) idiom. Joel Andersson

C++ includes: generic_expression.hpp ";

%feature("docstring") friendwrap_is_equal "

Check if two nodes are equivalent up to a given depth. Depth=0 checks if the
expressions are identical, i.e. points to the same node.

a = x*x b = x*x

a.is_equal(b, 0) will return false, but a.is_equal(b, 1) will return true

";

%feature("docstring") friendwrap_sin "

Sine.

";

%feature("docstring") friendwrap_asin "

Arc sine.

";


// File: classcasadi_1_1GenericExternal.xml


// File: classcasadi_1_1GenericMatrix.xml


/*  Construct symbolic primitives  */

/* The \"sym\" function is intended to work in a similar way as \"sym\" used
in the Symbolic Toolbox for Matlab but instead creating a CasADi symbolic
primitive.

*/ %feature("docstring") friendwrap_sumRows "

Return a row-wise summation of elements.

";

%feature("docstring") friendwrap_substitute "

>  MatType substitute(MatType ex, MatType v, MatType vdef)
------------------------------------------------------------------------

Substitute variable v with expression vdef in an expression ex.

>  [MatType] substitute([MatType ] ex, [MatType ] v, [MatType ] vdef)
------------------------------------------------------------------------

Substitute variable var with expression expr in multiple expressions.

";

%feature("docstring") casadi::GenericMatrix::numel "

>  int MatType .numel() const 
------------------------------------------------------------------------

Get the number of elements.

>  int MatType .numel(int i) const 
------------------------------------------------------------------------

Get the number of elements in slice (cf. MATLAB)

";

%feature("docstring") friendwrap_norm_inf "

Infinity-norm.

";

%feature("docstring") friendwrap_pinv "

>  MatType pinv(MatType A)
------------------------------------------------------------------------

Computes the Moore-Penrose pseudo-inverse.

If the matrix A is fat (size1<size2), mul(A, pinv(A)) is unity.

pinv(A)' = (AA')^(-1) A

If the matrix A is slender (size1>size2), mul(pinv(A), A) is unity.

pinv(A) = (A'A)^(-1) A'

>  MatType pinv(MatType A, str lsolver, Dict dict=Dict())
------------------------------------------------------------------------

Computes the Moore-Penrose pseudo-inverse.

If the matrix A is fat (size1>size2), mul(A, pinv(A)) is unity. If the
matrix A is slender (size2<size1), mul(pinv(A), A) is unity.

";

%feature("docstring") friendwrap_countNodes "

Count number of nodes

";

%feature("docstring") friendwrap_hessian "";

%feature("docstring") friendwrap_trace "

Matrix trace.

";

%feature("docstring") casadi::GenericMatrix::nnz_lower "

Get the number of non-zeros in the lower triangular half.

";

%feature("docstring") friendwrap_mldivide "

Matrix divide (cf. backslash '\\\\' in MATLAB)

";

%feature("docstring") friendwrap_inner_prod "

Inner product of two matrices with x and y matrices of the same dimension.

";

%feature("docstring") casadi::GenericMatrix::get_colind "

Get the sparsity pattern. See the Sparsity class for details.

";

%feature("docstring") friendwrap_cross "

Matlab's cross command.

";

%feature("docstring") friendwrap_tril2symm "

Convert a lower triangular matrix to a symmetric one.

";

%feature("docstring") friendwrap_substituteInPlace "

Inplace substitution with piggyback expressions Substitute variables v out
of the expressions vdef sequentially, as well as out of a number of other
expressions piggyback.

";

%feature("docstring") casadi::GenericMatrix::issquare "

Check if the matrix expression is square.

";

%feature("docstring") friendwrap_quad_form "

>  MatType quad_form(MatType X, MatType A)
------------------------------------------------------------------------

Calculate quadratic form X^T A X.

>  MatType quad_form(MatType X)
------------------------------------------------------------------------

Calculate quadratic form X^T X.

";

%feature("docstring") friendwrap_det "

Matrix determinant (experimental)

";

%feature("docstring") casadi::GenericMatrix::sparsity "

Get the sparsity pattern.

";

%feature("docstring") casadi::GenericMatrix::size2 "

Get the second dimension (i.e. number of columns)

";

%feature("docstring") friendwrap_dependsOn "

Check if expression depends on the argument The argument must be symbolic.

";

%feature("docstring") friendwrap_mpower "

Matrix power x^n.

";

%feature("docstring") casadi::GenericMatrix::nnz "

Get the number of (structural) non-zero elements.

";

%feature("docstring") friendwrap_norm_F "

Frobenius norm.

";

%feature("docstring") friendwrap_repsum "

Given a repeated matrix, computes the sum of repeated parts.

";

%feature("docstring") friendwrap_getOperatorRepresentation "

Get a string representation for a binary MatType, using custom arguments.

";

%feature("docstring") casadi::GenericMatrix::colind "

Get the sparsity pattern. See the Sparsity class for details.

";

%feature("docstring") casadi::GenericMatrix::nnz_diag "

Get get the number of non-zeros on the diagonal.

";

%feature("docstring") casadi::GenericMatrix::isempty "

Check if the sparsity is empty, i.e. if one of the dimensions is zero (or
optionally both dimensions)

";

%feature("docstring") casadi::GenericMatrix::size "

>  (int,int) MatType .size() const 
------------------------------------------------------------------------

Get the shape.

>  int MatType .size(int axis) const 
------------------------------------------------------------------------

Get the size along a particular dimensions.

";

%feature("docstring") casadi::GenericMatrix::iscolumn "

Check if the matrix is a column vector (i.e. size2()==1)

";

%feature("docstring") casadi::GenericMatrix "

Matrix base class.

This is a common base class for MX and Matrix<>, introducing a uniform
syntax and implementing common functionality using the curiously recurring
template pattern (CRTP) idiom.  The class is designed with the idea that
\"everything is a matrix\", that is, also scalars and vectors. This
philosophy makes it easy to use and to interface in particularly with Python
and Matlab/Octave.  The syntax tries to stay as close as possible to the
ublas syntax when it comes to vector/matrix operations.  Index starts with
0. Index vec happens as follows: (rr, cc) -> k = rr+cc*size1() Vectors are
column vectors.  The storage format is Compressed Column Storage (CCS),
similar to that used for sparse matrices in Matlab, but unlike this format,
we do allow for elements to be structurally non-zero but numerically zero.
The sparsity pattern, which is reference counted and cached, can be accessed
with Sparsity& sparsity() Joel Andersson

C++ includes: generic_matrix.hpp ";

%feature("docstring") friendwrap_tangent "

Matrix power x^n.

";

%feature("docstring") casadi::GenericMatrix::find "

Get the location of all non-zero elements as they would appear in a Dense
matrix A : DenseMatrix 4 x 3 B : SparseMatrix 4 x 3 , 5 structural non-
zeros.

k = A.find() A[k] will contain the elements of A that are non-zero in B

";

%feature("docstring") friendwrap_polyval "

Evaluate a polynomial with coefficients p in x.

";

%feature("docstring") friendwrap_triu2symm "

Convert a upper triangular matrix to a symmetric one.

";

%feature("docstring") friendwrap_outer_prod "

Take the outer product of two vectors Equals.

with x and y vectors

";

%feature("docstring") friendwrap_symvar "

Get all symbols contained in the supplied expression Get all symbols on
which the supplied expression depends.

See:  SXFunction::getFree(), MXFunction::getFree()

";

%feature("docstring") friendwrap_sum_square "

Calculate some of squares: sum_ij X_ij^2.

";

%feature("docstring") friendwrap_solve "

>  MatType solve(MatType A, MatType b)
------------------------------------------------------------------------

Solve a system of equations: A*x = b The solve routine works similar to
Matlab's backslash when A is square and nonsingular. The algorithm used is
the following:

A simple forward or backward substitution if A is upper or lower triangular

If the linear system is at most 3-by-3, form the inverse via minor expansion
and multiply

Permute the variables and equations as to get a (structurally) nonzero
diagonal, then perform a QR factorization without pivoting and solve the
factorized system.

Note 1: If there are entries of the linear system known to be zero, these
will be removed. Elements that are very small, or will evaluate to be zero,
can still cause numerical errors, due to the lack of pivoting (which is not
possible since cannot compare the size of entries)

Note 2: When permuting the linear system, a BLT (block lower triangular)
transformation is formed. Only the permutation part of this is however used.
An improvement would be to solve block-by-block if there are multiple BLT
blocks.

>  MatType solve(MatType A, MatType b, str lsolver, Dict dict=Dict())
------------------------------------------------------------------------

Solve a system of equations: A*x = b.

";

%feature("docstring") friendwrap_unite "

Unite two matrices no overlapping sparsity.

";

%feature("docstring") casadi::GenericMatrix::get_row "

Get the sparsity pattern. See the Sparsity class for details.

";

%feature("docstring") friendwrap_gradient "

Matrix power x^n.

";

%feature("docstring") friendwrap_sumCols "

Return a col-wise summation of elements.

";

%feature("docstring") casadi::GenericMatrix::nnz_upper "

Get the number of non-zeros in the upper triangular half.

";

%feature("docstring") casadi::GenericMatrix::isvector "

Check if the matrix is a row or column vector.

";

%feature("docstring") casadi::GenericMatrix::isrow "

Check if the matrix is a row vector (i.e. size1()==1)

";

%feature("docstring") friendwrap_inv "

Matrix inverse (experimental)

";

%feature("docstring") friendwrap_mrdivide "

Matrix divide (cf. slash '/' in MATLAB)

";

%feature("docstring") friendwrap_jacobian "

Calculate jacobian via source code transformation.

";

%feature("docstring") casadi::GenericMatrix::zeros "

Create a dense matrix or a matrix with specified sparsity with all entries
zero.

";

%feature("docstring") friendwrap_conditional "

Create a switch.

If the condition

Parameters:
-----------

ind:  evaluates to the integer k, where 0<=k<f.size(), then x[k] will be
returned, otherwise

x_default:  will be returned.

";

%feature("docstring") friendwrap_norm_1 "

1-norm

";

%feature("docstring") friendwrap_nullspace "

Computes the nullspace of a matrix A.

Finds Z m-by-(m-n) such that AZ = 0 with A n-by-m with m > n

Assumes A is full rank

Inspired by Numerical Methods in Scientific Computing by Ake Bjorck

";

%feature("docstring") friendwrap_norm_2 "

2-norm

";

%feature("docstring") friendwrap_if_else "

Branching on MX nodes Ternary operator, \"cond ? if_true : if_false\".

";

%feature("docstring") friendwrap_diag "

Get the diagonal of a matrix or construct a diagonal When the input is
square, the diagonal elements are returned. If the input is vector- like, a
diagonal matrix is constructed with it.

";

%feature("docstring") casadi::GenericMatrix::isdense "

Check if the matrix expression is dense.

";

%feature("docstring") casadi::GenericMatrix::isscalar "

Check if the matrix expression is scalar.

";

%feature("docstring") friendwrap_linspace "

Matlab's linspace command.

";

%feature("docstring") casadi::GenericMatrix::dim "

Get string representation of dimensions. The representation is (nrow x ncol
= numel | size)

";

%feature("docstring") friendwrap_densify "

>  MatType densify(MatType x)
------------------------------------------------------------------------

Make the matrix dense if not already.

>  MatType densify(MatType x, MatType val)
------------------------------------------------------------------------

Make the matrix dense and assign nonzeros to a value.

";

%feature("docstring") casadi::GenericMatrix::sym "

>  static MatType MatType .sym(str name, int nrow=1, int ncol=1)
------------------------------------------------------------------------

Create an nrow-by-ncol symbolic primitive.

>  static MatType MatType .sym(str name, (int,int) rc)
------------------------------------------------------------------------

Construct a symbolic primitive with given dimensions.

>  MatType MatType .sym(str name, Sparsity sp)
------------------------------------------------------------------------

Create symbolic primitive with a given sparsity pattern.

>  [MatType ] MatType .sym(str name, Sparsity sp, int p)
------------------------------------------------------------------------

Create a vector of length p with with matrices with symbolic primitives of
given sparsity.

>  static[MatType ] MatType .sym(str name, int nrow, int ncol, int p)
------------------------------------------------------------------------

Create a vector of length p with nrow-by-ncol symbolic primitives.

>  [[MatType ] ] MatType .sym(str name, Sparsity sp, int p, int r)
------------------------------------------------------------------------

Create a vector of length r of vectors of length p with symbolic primitives
with given sparsity.

>  static[[MatType] ] MatType .sym(str name, int nrow, int ncol, int p, int r)
------------------------------------------------------------------------

Create a vector of length r of vectors of length p with nrow-by-ncol
symbolic primitives.

>  SX SX .sym(str name, Sparsity sp)
------------------------------------------------------------------------
[INTERNAL] 
";

%feature("docstring") casadi::GenericMatrix::istril "

Check if the matrix is lower triangular.

";

%feature("docstring") friendwrap_project "

Create a new matrix with a given sparsity pattern but with the nonzeros
taken from an existing matrix.

";

%feature("docstring") casadi::GenericMatrix::row "

Get the sparsity pattern. See the Sparsity class for details.

";

%feature("docstring") friendwrap_extractShared "

Extract shared subexpressions from an set of expressions.

";

%feature("docstring") casadi::GenericMatrix::ones "

Create a dense matrix or a matrix with specified sparsity with all entries
one.

";

%feature("docstring") casadi::GenericMatrix::istriu "

Check if the matrix is upper triangular.

";

%feature("docstring") casadi::GenericMatrix::size1 "

Get the first dimension (i.e. number of rows)

";


// File: classcasadi_1_1GenericType.xml
%feature("docstring") casadi::GenericType::isDouble "

Check if a particular type.

";

%feature("docstring") casadi::GenericType::isBool "

Check if a particular type.

";

%feature("docstring") casadi::GenericType::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::GenericType::isemptyVector "

Check if a particular type.

";

%feature("docstring") casadi::GenericType::toIntVectorVector "

Convert to a type.

";

%feature("docstring") casadi::GenericType::print "

Print a description of the object.

";

%feature("docstring") casadi::GenericType::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::GenericType::toDict "

Convert to a type.

";

%feature("docstring") casadi::GenericType::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::GenericType::isString "

Check if a particular type.

";

%feature("docstring") casadi::GenericType::toFunction "

Convert to a type.

";

%feature("docstring") casadi::GenericType::isDoubleVector "

Check if a particular type.

";

%feature("docstring") casadi::GenericType "

Generic data type, can hold different types such as bool, int, string etc.

Joel Andersson

C++ includes: generic_type.hpp ";

%feature("docstring") casadi::GenericType::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::GenericType::getType "";

%feature("docstring") casadi::GenericType::isIntVectorVector "

Check if a particular type.

";

%feature("docstring") casadi::GenericType::isNull "

Is a null pointer?

";

%feature("docstring") casadi::GenericType::isInt "

Check if a particular type.

";

%feature("docstring") casadi::GenericType::isVoidPointer "

Check if a particular type.

";

%feature("docstring") casadi::GenericType::can_cast_to "";

%feature("docstring") casadi::GenericType::toStringVector "

Convert to a type.

";

%feature("docstring") casadi::GenericType::isDict "

Check if a particular type.

";

%feature("docstring") casadi::GenericType::toDoubleVector "

Convert to a type.

";

%feature("docstring") casadi::GenericType::toInt "

Convert to a type.

";

%feature("docstring") casadi::GenericType::repr "

Print a representation of the object.

";

%feature("docstring") casadi::GenericType::toBool "

Convert to a type.

";

%feature("docstring") casadi::GenericType::isFunction "

Check if a particular type.

";

%feature("docstring") casadi::GenericType::GenericType "

>  GenericType()
------------------------------------------------------------------------

Default constructor.

>  GenericType(bool b)
------------------------------------------------------------------------

Constructors (implicit type conversion)

";

%feature("docstring") casadi::GenericType::toString "

Convert to a type.

";

%feature("docstring") casadi::GenericType::isStringVector "

Check if a particular type.

";

%feature("docstring") casadi::GenericType::isIntVector "

Check if a particular type.

";

%feature("docstring") casadi::GenericType::get_description "

Get a description of the object's type.

";

%feature("docstring") casadi::GenericType::toDouble "

Convert to a type.

";

%feature("docstring") casadi::GenericType::toVoidPointer "

Convert to a type.

";

%feature("docstring") casadi::GenericType::toIntVector "

Convert to a type.

";


// File: classcasadi_1_1GenericTypeBase.xml


// File: classcasadi_1_1GetNonzeros.xml


// File: classcasadi_1_1GetNonzerosSlice.xml


// File: classcasadi_1_1GetNonzerosSlice2.xml


// File: classcasadi_1_1GetNonzerosVector.xml


// File: classcasadi_1_1Horzcat.xml


// File: classcasadi_1_1HorzRepmat.xml


// File: classcasadi_1_1HorzRepsum.xml


// File: classcasadi_1_1Horzsplit.xml


// File: classcasadi_1_1ImplicitFixedStepIntegrator.xml


// File: classcasadi_1_1ImplicitFunction.xml


/*  Simple Getters & Setters  */

/*  Advanced Getters  */

/*  Option Functionality  */ %feature("docstring")
casadi::ImplicitFunction::getInput "

>  DMatrix  IOInterface< Function  >.getInput(int iind=0) const
------------------------------------------------------------------------

Get an input by index.

Parameters:
-----------

iind:  index within the range [0..n_in()-1]

>  DMatrix  IOInterface< Function  >.getInput(str iname) const
------------------------------------------------------------------------

Get an input by name.

Parameters:
-----------

iname:  input name. Only allowed when an input scheme is set.

>  void IOInterface< Function  >.getInput(T val, int iind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.getInput(T val, str iname)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::ImplicitFunction::generateLiftingFunctions "

Extract the functions needed for the Lifted Newton method.

";

%feature("docstring") casadi::ImplicitFunction::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::ImplicitFunction::tangent "

Generate a tangent function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the input must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::ImplicitFunction::getOptionAllowedIndex "[INTERNAL]  Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::ImplicitFunction::spInit "[INTERNAL]  Reset
the sparsity propagation.

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::ImplicitFunction::description_out "

Get output scheme description by index.

";

%feature("docstring") casadi::ImplicitFunction::getStat "

Get a single statistic obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::ImplicitFunction::getAtomicInput "

Get the (integer) input arguments of an atomic operation.

";

%feature("docstring") casadi::ImplicitFunction::n_in "

Get the number of function inputs.

";

%feature("docstring") casadi::ImplicitFunction::setFullJacobian "

Set the Jacobian of all the input nonzeros with respect to all output
nonzeros NOTE: Does not take ownership, only weak references to the Jacobian
are kept internally

";

%feature("docstring") casadi::ImplicitFunction::spCanEvaluate "[INTERNAL]
Is the class able to propagate seeds through the algorithm?

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::ImplicitFunction::checkInputs "[INTERNAL]
Check if the numerical values of the supplied bounds make sense.

";

%feature("docstring") casadi::ImplicitFunction::name "

Name of the function.

";

%feature("docstring") casadi::ImplicitFunction::n_out "

Get the number of function outputs.

";

%feature("docstring") casadi::ImplicitFunction::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::ImplicitFunction::getOutput "

>  DMatrix  IOInterface< Function  >.getOutput(int oind=0) const
------------------------------------------------------------------------

Get an output by index.

Parameters:
-----------

oind:  index within the range [0..n_out()-1]

>  DMatrix  IOInterface< Function  >.getOutput(str oname) const
------------------------------------------------------------------------

Get an output by name.

Parameters:
-----------

oname:  output name. Only allowed when an output scheme is set.

>  void IOInterface< Function  >.getOutput(T val, int oind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.getOutput(T val, str oname)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::ImplicitFunction::mx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::ImplicitFunction::size1_out "

Get output dimension.

";

%feature("docstring") casadi::ImplicitFunction::type_name "

Get type name.

";

%feature("docstring") casadi::ImplicitFunction::hessian "

Generate a Hessian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The generated Hessian has two more outputs than the calling function
corresponding to the Hessian and the gradients.

";

%feature("docstring") casadi::ImplicitFunction::print "

Print a description of the object.

";

%feature("docstring") casadi::ImplicitFunction::countNodes "

Number of nodes in the algorithm.

";

%feature("docstring") casadi::ImplicitFunction::printDimensions "

Print dimensions of inputs and outputs.

";

%feature("docstring") casadi::ImplicitFunction::removeMonitor "

Remove modules to be monitored.

";

%feature("docstring") casadi::ImplicitFunction::jacSparsity "

Get, if necessary generate, the sparsity of a Jacobian block

";

%feature("docstring") casadi::ImplicitFunction::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

";

%feature("docstring") casadi::ImplicitFunction::derivative "

Get a function that calculates nfwd forward derivatives and nadj adjoint
derivatives Legacy function: Use derForward and derReverse instead.

Returns a function with (1+nfwd)*n_in+nadj*n_out inputs and (1+nfwd)*n_out +
nadj*n_in outputs. The first n_in inputs correspond to nondifferentiated
inputs. The next nfwd*n_in inputs correspond to forward seeds, one direction
at a time and the last nadj*n_out inputs correspond to adjoint seeds, one
direction at a time. The first n_out outputs correspond to nondifferentiated
outputs. The next nfwd*n_out outputs correspond to forward sensitivities,
one direction at a time and the last nadj*n_in outputs corresponds to
adjoint sensitivities, one direction at a time.

(n_in = n_in(), n_out = n_out())

";

%feature("docstring") casadi::ImplicitFunction::size2_in "

Get input dimension.

";

%feature("docstring") casadi::ImplicitFunction "

Abstract base class for the implicit function classes.

Mathematically, the equation:

F(z, x1, x2, ..., xn) == 0

where d_F/dz is invertible, implicitly defines the equation:

z := G(x1, x2, ..., xn)

In CasADi, F is a Function. The first input presents the variables that need
to be solved for. The first output is the residual that needs to attain
zero. Optional remaining outputs can be supplied; they are output
expressions.

In pseudo-code, we can write:

G* = ImplicitFunction('solver',F)

Here, G* is a Function with one extra input over the pure mathematical G:

z := G*(z0, x1, x2, ..., xn)

The first input to the ImplicitFunction is the intial guess for z.

General information
===================



>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| ad_weight    | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | derivative c |              |
|              |              |              | alculation.W |              |
|              |              |              | hen there is |              |
|              |              |              | an option of |              |
|              |              |              | either using |              |
|              |              |              | forward or   |              |
|              |              |              | reverse mode |              |
|              |              |              | directional  |              |
|              |              |              | derivatives, |              |
|              |              |              | the          |              |
|              |              |              | condition ad |              |
|              |              |              | _weight*nf<= |              |
|              |              |              | (1-ad_weight |              |
|              |              |              | )*na is used |              |
|              |              |              | where nf and |              |
|              |              |              | na are       |              |
|              |              |              | estimates of |              |
|              |              |              | the number   |              |
|              |              |              | of forward/r |              |
|              |              |              | everse mode  |              |
|              |              |              | directional  |              |
|              |              |              | derivatives  |              |
|              |              |              | needed. By   |              |
|              |              |              | default,     |              |
|              |              |              | ad_weight is |              |
|              |              |              | calculated a |              |
|              |              |              | utomatically |              |
|              |              |              | , but this   |              |
|              |              |              | can be       |              |
|              |              |              | overridden   |              |
|              |              |              | by setting   |              |
|              |              |              | this option. |              |
|              |              |              | In           |              |
|              |              |              | particular,  |              |
|              |              |              | 0 means      |              |
|              |              |              | forcing      |              |
|              |              |              | forward mode |              |
|              |              |              | and 1        |              |
|              |              |              | forcing      |              |
|              |              |              | reverse      |              |
|              |              |              | mode. Leave  |              |
|              |              |              | unset for    |              |
|              |              |              | (class       |              |
|              |              |              | specific)    |              |
|              |              |              | heuristics.  |              |
+--------------+--------------+--------------+--------------+--------------+
| ad_weight_sp | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | sparsity     |              |
|              |              |              | pattern      |              |
|              |              |              | calculation  |              |
|              |              |              | calculation. |              |
|              |              |              | Overrides    |              |
|              |              |              | default      |              |
|              |              |              | behavior.    |              |
|              |              |              | Set to 0 and |              |
|              |              |              | 1 to force   |              |
|              |              |              | forward and  |              |
|              |              |              | reverse mode |              |
|              |              |              | respectively |              |
|              |              |              | . Cf. option |              |
|              |              |              | \"ad_weight\". |              |
+--------------+--------------+--------------+--------------+--------------+
| compiler     | OT_STRING    | \"clang\"      | Just-in-time | casadi::Func |
|              |              |              | compiler     | tionInternal |
|              |              |              | plugin to be |              |
|              |              |              | used.        |              |
+--------------+--------------+--------------+--------------+--------------+
| constraints  | OT_INTEGERVE | GenericType( | Constrain    | casadi::Impl |
|              | CTOR         | )            | the          | icitFunction |
|              |              |              | unknowns. 0  | Internal     |
|              |              |              | (default):   |              |
|              |              |              | no           |              |
|              |              |              | constraint   |              |
|              |              |              | on ui, 1: ui |              |
|              |              |              | >= 0.0, -1:  |              |
|              |              |              | ui <= 0.0,   |              |
|              |              |              | 2: ui > 0.0, |              |
|              |              |              | -2: ui <     |              |
|              |              |              | 0.0.         |              |
+--------------+--------------+--------------+--------------+--------------+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode     |
|              |              |              | according to |              |
|              |              |              | a given      |              |
|              |              |              | recipe (low- |              |
|              |              |              | level)       |              |
+--------------+--------------+--------------+--------------+--------------+
| gather_stats | OT_BOOLEAN   | false        | Flag to      | casadi::Func |
|              |              |              | indicate     | tionInternal |
|              |              |              | whether      |              |
|              |              |              | statistics   |              |
|              |              |              | must be      |              |
|              |              |              | gathered     |              |
+--------------+--------------+--------------+--------------+--------------+
| implicit_inp | OT_INTEGER   | 0            | Index of the | casadi::Impl |
| ut           |              |              | input that   | icitFunction |
|              |              |              | corresponds  | Internal     |
|              |              |              | to the       |              |
|              |              |              | actual root- |              |
|              |              |              | finding      |              |
+--------------+--------------+--------------+--------------+--------------+
| implicit_out | OT_INTEGER   | 0            | Index of the | casadi::Impl |
| put          |              |              | output that  | icitFunction |
|              |              |              | corresponds  | Internal     |
|              |              |              | to the       |              |
|              |              |              | actual root- |              |
|              |              |              | finding      |              |
+--------------+--------------+--------------+--------------+--------------+
| input_scheme | OT_STRINGVEC | GenericType( | Custom input | casadi::Func |
|              | TOR          | )            | scheme       | tionInternal |
+--------------+--------------+--------------+--------------+--------------+
| inputs_check | OT_BOOLEAN   | true         | Throw        | casadi::Func |
|              |              |              | exceptions   | tionInternal |
|              |              |              | when the     |              |
|              |              |              | numerical    |              |
|              |              |              | values of    |              |
|              |              |              | the inputs   |              |
|              |              |              | don't make   |              |
|              |              |              | sense        |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_penalty  | OT_REAL      | 2            | When         | casadi::Func |
|              |              |              | requested    | tionInternal |
|              |              |              | for a number |              |
|              |              |              | of forward/r |              |
|              |              |              | everse       |              |
|              |              |              | directions,  |              |
|              |              |              | it may be    |              |
|              |              |              | cheaper to   |              |
|              |              |              | compute      |              |
|              |              |              | first the    |              |
|              |              |              | full         |              |
|              |              |              | jacobian and |              |
|              |              |              | then         |              |
|              |              |              | multiply     |              |
|              |              |              | with seeds,  |              |
|              |              |              | rather than  |              |
|              |              |              | obtain the   |              |
|              |              |              | requested    |              |
|              |              |              | directions   |              |
|              |              |              | in a straigh |              |
|              |              |              | tforward     |              |
|              |              |              | manner.      |              |
|              |              |              | Casadi uses  |              |
|              |              |              | a heuristic  |              |
|              |              |              | to decide    |              |
|              |              |              | which is     |              |
|              |              |              | cheaper. A   |              |
|              |              |              | high value   |              |
|              |              |              | of 'jac_pena |              |
|              |              |              | lty' makes   |              |
|              |              |              | it less      |              |
|              |              |              | likely for   |              |
|              |              |              | the heurstic |              |
|              |              |              | to chose the |              |
|              |              |              | full         |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy.    |              |
|              |              |              | The special  |              |
|              |              |              | value -1     |              |
|              |              |              | indicates    |              |
|              |              |              | never to use |              |
|              |              |              | the full     |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy     |              |
+--------------+--------------+--------------+--------------+--------------+
| jacobian_fun | OT_FUNCTION  | GenericType( | Function     | casadi::Impl |
| ction        |              | )            | object for   | icitFunction |
|              |              |              | calculating  | Internal     |
|              |              |              | the Jacobian |              |
|              |              |              | (autogenerat |              |
|              |              |              | ed by        |              |
|              |              |              | default)     |              |
+--------------+--------------+--------------+--------------+--------------+
| jit          | OT_BOOLEAN   | false        | Use just-in- | casadi::Func |
|              |              |              | time         | tionInternal |
|              |              |              | compiler to  |              |
|              |              |              | speed up the |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| jit_options  | OT_DICT      | GenericType( | Options to   | casadi::Func |
|              |              | )            | be passed to | tionInternal |
|              |              |              | the jit      |              |
|              |              |              | compiler.    |              |
+--------------+--------------+--------------+--------------+--------------+
| linear_solve | OT_STRING    | \"csparse\"    | User-defined | casadi::Impl |
| r            |              |              | linear       | icitFunction |
|              |              |              | solver       | Internal     |
|              |              |              | class.       |              |
|              |              |              | Needed for s |              |
|              |              |              | ensitivities |              |
|              |              |              | .            |              |
+--------------+--------------+--------------+--------------+--------------+
| linear_solve | OT_FUNCTION  | GenericType( | Function     | casadi::Impl |
| r_function   |              | )            | object for   | icitFunction |
|              |              |              | solving the  | Internal     |
|              |              |              | linearized   |              |
|              |              |              | problem (aut |              |
|              |              |              | ogenerated   |              |
|              |              |              | by default)  |              |
+--------------+--------------+--------------+--------------+--------------+
| linear_solve | OT_DICT      | GenericType( | Options to   | casadi::Impl |
| r_options    |              | )            | be passed to | icitFunction |
|              |              |              | the linear   | Internal     |
|              |              |              | solver.      |              |
+--------------+--------------+--------------+--------------+--------------+
| monitor      | OT_STRINGVEC | GenericType( | Monitors to  | casadi::Func |
|              | TOR          | )            | be activated | tionInternal |
|              |              |              | (inputs|outp |              |
|              |              |              | uts)         |              |
+--------------+--------------+--------------+--------------+--------------+
| output_schem | OT_STRINGVEC | GenericType( | Custom       | casadi::Func |
| e            | TOR          | )            | output       | tionInternal |
|              |              |              | scheme       |              |
+--------------+--------------+--------------+--------------+--------------+
| regularity_c | OT_BOOLEAN   | true         | Throw        | casadi::Func |
| heck         |              |              | exceptions   | tionInternal |
|              |              |              | when NaN or  |              |
|              |              |              | Inf appears  |              |
|              |              |              | during       |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| user_data    | OT_VOIDPTR   | GenericType( | A user-      | casadi::Func |
|              |              | )            | defined      | tionInternal |
|              |              |              | field that   |              |
|              |              |              | can be used  |              |
|              |              |              | to identify  |              |
|              |              |              | the function |              |
|              |              |              | or pass      |              |
|              |              |              | additional   |              |
|              |              |              | information  |              |
+--------------+--------------+--------------+--------------+--------------+
| verbose      | OT_BOOLEAN   | false        | Verbose      | casadi::Func |
|              |              |              | evaluation   | tionInternal |
|              |              |              | for          |              |
|              |              |              | debugging    |              |
+--------------+--------------+--------------+--------------+--------------+

List of plugins
===============



- kinsol

- nlp

- newton

Note: some of the plugins in this list might not be available on your
system. Also, there might be extra plugins available to you that are not
listed here. You can obtain their documentation with
ImplicitFunction.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

kinsol
------



KINSOL interface from the Sundials suite

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| abstol          | OT_REAL         | 0.000           | Stopping        |
|                 |                 |                 | criterion       |
|                 |                 |                 | tolerance       |
+-----------------+-----------------+-----------------+-----------------+
| disable_interna | OT_BOOLEAN      | false           | Disable KINSOL  |
| l_warnings      |                 |                 | internal        |
|                 |                 |                 | warning         |
|                 |                 |                 | messages        |
+-----------------+-----------------+-----------------+-----------------+
| exact_jacobian  | OT_BOOLEAN      | true            |                 |
+-----------------+-----------------+-----------------+-----------------+
| f_scale         | OT_REALVECTOR   |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| iterative_solve | OT_STRING       | \"gmres\"         | gmres|bcgstab|t |
| r               |                 |                 | fqmr            |
+-----------------+-----------------+-----------------+-----------------+
| linear_solver_t | OT_STRING       | \"dense\"         | dense|banded|it |
| ype             |                 |                 | erative|user_de |
|                 |                 |                 | fined           |
+-----------------+-----------------+-----------------+-----------------+
| lower_bandwidth | OT_INTEGER      |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| max_iter        | OT_INTEGER      | 0               | Maximum number  |
|                 |                 |                 | of Newton       |
|                 |                 |                 | iterations.     |
|                 |                 |                 | Putting 0 sets  |
|                 |                 |                 | the default     |
|                 |                 |                 | value of        |
|                 |                 |                 | KinSol.         |
+-----------------+-----------------+-----------------+-----------------+
| max_krylov      | OT_INTEGER      | 0               |                 |
+-----------------+-----------------+-----------------+-----------------+
| pretype         | OT_STRING       | \"none\"          | (none|left|righ |
|                 |                 |                 | t|both)         |
+-----------------+-----------------+-----------------+-----------------+
| strategy        | OT_STRING       | \"none\"          | Globalization   |
|                 |                 |                 | strategy (none| |
|                 |                 |                 | linesearch)     |
+-----------------+-----------------+-----------------+-----------------+
| u_scale         | OT_REALVECTOR   |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| upper_bandwidth | OT_INTEGER      |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| use_preconditio | OT_BOOLEAN      | false           | precondition an |
| ner             |                 |                 | iterative       |
|                 |                 |                 | solver          |
+-----------------+-----------------+-----------------+-----------------+

>List of available monitors

+-----------+
|    Id     |
+===========+
| eval_djac |
+-----------+
| eval_f    |
+-----------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

nlp
---



Use an NlpSolver as ImplicitFunction solver

>List of available options

+----+------+---------+-------------+
| Id | Type | Default | Description |
+====+======+=========+=============+
+----+------+---------+-------------+

>List of available stats

+--------------+
|      Id      |
+==============+
| solver_stats |
+--------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

newton
------



Implements simple newton iterations to solve an implicit function.

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| abstol          | OT_REAL         | 0.000           | Stopping        |
|                 |                 |                 | criterion       |
|                 |                 |                 | tolerance on    |
|                 |                 |                 | max(|F|)        |
+-----------------+-----------------+-----------------+-----------------+
| abstolStep      | OT_REAL         | 0.000           | Stopping        |
|                 |                 |                 | criterion       |
|                 |                 |                 | tolerance on    |
|                 |                 |                 | step size       |
+-----------------+-----------------+-----------------+-----------------+
| max_iter        | OT_INTEGER      | 1000            | Maximum number  |
|                 |                 |                 | of Newton       |
|                 |                 |                 | iterations to   |
|                 |                 |                 | perform before  |
|                 |                 |                 | returning.      |
+-----------------+-----------------+-----------------+-----------------+
| print_iteration | OT_BOOLEAN      | false           | Print           |
|                 |                 |                 | information     |
|                 |                 |                 | about each      |
|                 |                 |                 | iteration       |
+-----------------+-----------------+-----------------+-----------------+

>List of available monitors

+----------+
|    Id    |
+==========+
| F        |
+----------+
| J        |
+----------+
| normF    |
+----------+
| step     |
+----------+
| stepsize |
+----------+

>List of available stats

+---------------+
|      Id       |
+===============+
| iter          |
+---------------+
| return_status |
+---------------+

--------------------------------------------------------------------------------



Joel Andersson
Diagrams
--------



C++ includes: implicit_function.hpp ";

%feature("docstring") casadi::ImplicitFunction::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::ImplicitFunction::callForward "

Create call to (cached) derivative function, forward mode.

";

%feature("docstring") casadi::ImplicitFunction::call "

Evaluate the function symbolically or numerically.

";

%feature("docstring") casadi::ImplicitFunction::spEvaluate "[INTERNAL]
Propagate the sparsity pattern through a set of directional.

derivatives forward or backward (for usage, see the example
propagating_sparsity.cpp)

";

%feature("docstring") casadi::ImplicitFunction::free_mx "

Get all the free variables of the function.

";

%feature("docstring") casadi::ImplicitFunction::getAtomicOperation "

Get an atomic operation operator index.

";

%feature("docstring") casadi::ImplicitFunction::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::ImplicitFunction::size_in "

Get input dimension.

";

%feature("docstring") casadi::ImplicitFunction::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::ImplicitFunction::ImplicitFunction "

>  ImplicitFunction()
------------------------------------------------------------------------

Default constructor.

>  ImplicitFunction(str name, str solver, Function f, Dict opts=Dict())
------------------------------------------------------------------------

Create an implicit function solver (new syntax, includes initialization)

Parameters:
-----------

solver:

Name of a solver. It might be one of:

- kinsol

- nlp

- newton

Note: some of the plugins in this list might not be available on your
system. Also, there might be extra plugins available to you that are not
listed here. You can obtain their documentation with
ImplicitFunction.doc(\"myextraplugin\")

f:   Function where one of the inputs (by default the first) is an unknown
and one of the outputs (by default the first) is a residual.

";

%feature("docstring") casadi::ImplicitFunction::setDerReverse "

Set a function that calculates nadj adjoint derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::ImplicitFunction::nnz_in "

Get of number of input nonzeros For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::ImplicitFunction::sparsity_out "

Get sparsity of a given output.

";

%feature("docstring") casadi::ImplicitFunction::size_out "

Get output dimension.

";

%feature("docstring") casadi::ImplicitFunction::setInput "

>  void IOInterface< Function  >.setInput(T val, int iind=0)
------------------------------------------------------------------------

Set an input by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.setInput(T val, str iname)
------------------------------------------------------------------------

Set an input by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::ImplicitFunction::repr "

Print a representation of the object.

";

%feature("docstring") casadi::ImplicitFunction::sz_arg "[INTERNAL]  Get
required length of arg field.

";

%feature("docstring") casadi::ImplicitFunction::name_in "

>  [str] Function.name_in() const 
------------------------------------------------------------------------

Get input scheme.

>  str Function.name_in(int ind) const 
------------------------------------------------------------------------

Get input scheme name by index.

";

%feature("docstring") casadi::ImplicitFunction::free_sx "

Get all the free variables of the function.

";

%feature("docstring") casadi::ImplicitFunction::is_a "

Check if the function is of a particular type Optionally check if name
matches one of the base classes (default true)

";

%feature("docstring") casadi::ImplicitFunction::setOutput "

>  void IOInterface< Function  >.setOutput(T val, int oind=0)
------------------------------------------------------------------------

Set an output by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.setOutput(T val, str oname)
------------------------------------------------------------------------

Set an output by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::ImplicitFunction::generate "

Export / Generate C code for the function.

";

%feature("docstring") casadi::ImplicitFunction::sx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::ImplicitFunction::setJacobian "

Set the Jacobian function of output oind with respect to input iind NOTE:
Does not take ownership, only weak references to the Jacobians are kept
internally

";

%feature("docstring") casadi::ImplicitFunction::getLinsol "

Access linear solver.

";

%feature("docstring") casadi::ImplicitFunction::map "

>  [[MX] ] Function.map([[MX ] ] arg, str parallelization=\"serial\")

>  [MX] Function.map([MX ] arg, str parallelization=\"serial\")
------------------------------------------------------------------------

Evaluate symbolically in parallel (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

>  Function Function.map(str name, int N, Dict options=Dict()) const 
------------------------------------------------------------------------

Create a mapped version of this function.

Suppose the function has a signature of:

::

     f: (a, p) -> ( s )
  



The the mapaccumulated version has the signature:

::

     F: (A, P) -> (S )
  
      with
          a: horzcat([a0, a1, ..., a_(N-1)])
          p: horzcat([p0, p1, ..., p_(N-1)])
          s: horzcat([s0, s1, ..., s_(N-1)])
      and
          s0 <- f(a0, p0)
          s1 <- f(a1, p1)
          ...
          s_(N-1) <- f(a_(N-1), p_(N-1))
  



>  Function Function.map(str name, int n, [bool ] repeat_in, [bool ] repeat_out, Dict opts=Dict()) const 
------------------------------------------------------------------------

Generic map.

";

%feature("docstring") casadi::ImplicitFunction::sparsity_in "

Get sparsity of a given input.

";

%feature("docstring") casadi::ImplicitFunction::derReverse "

Get a function that calculates nadj adjoint derivatives.

Returns a function with n_in + n_out +nadj*n_out inputs and nadj*n_in
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nadj*n_out inputs correspond to adjoint seeds, one direction at a time The
nadj*n_in outputs correspond to adjoint sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

(n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::ImplicitFunction::nnz_out "

Get of number of output nonzeros For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::ImplicitFunction::setDerForward "

Set a function that calculates nfwd forward derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::ImplicitFunction::sz_res "[INTERNAL]  Get
required length of res field.

";

%feature("docstring") casadi::ImplicitFunction::getAtomicOutput "

Get the (integer) output argument of an atomic operation.

";

%feature("docstring") casadi::ImplicitFunction::index_out "

Find the index for a string describing a particular entry of an output
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::ImplicitFunction::default_in "

Get default input value (NOTE: constant reference)

";

%feature("docstring") casadi::ImplicitFunction::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::ImplicitFunction::name_out "

>  [str] Function.name_out() const 
------------------------------------------------------------------------

Get output scheme.

>  str Function.name_out(int ind) const 
------------------------------------------------------------------------

Get output scheme name by index.

";

%feature("docstring") casadi::ImplicitFunction::jacobian "

Generate a Jacobian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. If
compact is set to true, only the nonzeros of the input and output
expressions are considered. If symmetric is set to true, the Jacobian being
calculated is known to be symmetric (usually a Hessian), which can be
exploited by the algorithm.

The generated Jacobian has one more output than the calling function
corresponding to the Jacobian and the same number of inputs.

";

%feature("docstring") casadi::ImplicitFunction::getOptionTypeName "

Get the type name of a certain option.

";

%feature("docstring") casadi::ImplicitFunction::setOptionByEnumValue "[INTERNAL]  Set a certain option by giving an enum value.

";

%feature("docstring") casadi::ImplicitFunction::getAlgorithmSize "

Get the number of atomic operations.

";

%feature("docstring") casadi::ImplicitFunction::setJacSparsity "

Generate the sparsity of a Jacobian block

";

%feature("docstring") casadi::ImplicitFunction::getF "

Access F.

";

%feature("docstring") casadi::ImplicitFunction::derForward "

Get a function that calculates nfwd forward derivatives.

Returns a function with n_in + n_out +nfwd*n_in inputs and nfwd*n_out
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nfwd*n_in inputs correspond to forward seeds, one direction at a time The
nfwd*n_out outputs correspond to forward sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::ImplicitFunction::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::ImplicitFunction::evaluate "

Evaluate.

";

%feature("docstring") casadi::ImplicitFunction::dictionary "

Get the dictionary.

";

%feature("docstring") casadi::ImplicitFunction::index_in "

Find the index for a string describing a particular entry of an input
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::ImplicitFunction::isNull "

Is a null pointer?

";

%feature("docstring") casadi::ImplicitFunction::addMonitor "

Add modules to be monitored.

";

%feature("docstring") casadi::ImplicitFunction::getWorkSize "

Get the length of the work vector.

";

%feature("docstring") casadi::ImplicitFunction::getJac "

Access Jacobian.

";

%feature("docstring") casadi::ImplicitFunction::mx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::ImplicitFunction::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::ImplicitFunction::printPtr "[INTERNAL]  Print
the pointer to the internal class

";

%feature("docstring") casadi::ImplicitFunction::sz_iw "[INTERNAL]  Get
required length of iw field.

";

%feature("docstring") casadi::ImplicitFunction::sz_w "[INTERNAL]  Get
required length of w field.

";

%feature("docstring") casadi::ImplicitFunction::setOptionByAllowedIndex "[INTERNAL]  Set a certain option by giving its index into the allowed
values.

";

%feature("docstring") casadi::ImplicitFunction::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::ImplicitFunction::gradient "

Generate a gradient function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the output must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::ImplicitFunction::mapaccum "

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
  



";

%feature("docstring") casadi::ImplicitFunction::getAtomicInputReal "

Get the floating point output argument of an atomic operation.

";

%feature("docstring") casadi::ImplicitFunction::getSanitizedName "

get function name with all non alphanumeric characters converted to '_'

";

%feature("docstring") casadi::ImplicitFunction::size2_out "

Get output dimension.

";

%feature("docstring") casadi::ImplicitFunction::numel_in "

Get of number of input elements For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::ImplicitFunction::callReverse "

Create call to (cached) derivative function, reverse mode.

";

%feature("docstring") casadi::ImplicitFunction::fullJacobian "

Generate a Jacobian function of all the inputs elements with respect to all
the output elements).

";

%feature("docstring") casadi::ImplicitFunction::numel_out "

Get of number of output elements For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::ImplicitFunction::getOptionEnumValue "[INTERNAL]  Get the enum value corresponding to th certain option.

";

%feature("docstring") casadi::ImplicitFunction::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::ImplicitFunction::callDerivative "[INTERNAL]
Evaluate the function symbolically or numerically with directional
derivatives The first two arguments are the nondifferentiated inputs and
results of the evaluation, the next two arguments are a set of forward
directional seeds and the resulting forward directional derivatives, the
length of the vector being the number of forward directions. The next two
arguments are a set of adjoint directional seeds and the resulting adjoint
directional derivatives, the length of the vector being the number of
adjoint directions.

";

%feature("docstring") casadi::ImplicitFunction::getStats "

Get all statistics obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::ImplicitFunction::size1_in "

Get input dimension.

";

%feature("docstring") casadi::ImplicitFunction::description_in "

Get input scheme description by index.

";

%feature("docstring") casadi::ImplicitFunction::sx_in "

Get symbolic primitives equivalent to the input expressions.

";


// File: classcasadi_1_1InfSX.xml


// File: classcasadi_1_1InnerProd.xml


// File: classcasadi_1_1IntegerSX.xml


// File: classcasadi_1_1Integrator.xml


/*  Simple Getters & Setters  */

/*  Advanced Getters  */

/*  Option Functionality  */ %feature("docstring")
casadi::Integrator::name_out "

>  [str] Function.name_out() const 
------------------------------------------------------------------------

Get output scheme.

>  str Function.name_out(int ind) const 
------------------------------------------------------------------------

Get output scheme name by index.

";

%feature("docstring") casadi::Integrator::size_in "

Get input dimension.

";

%feature("docstring") casadi::Integrator::description_in "

Get input scheme description by index.

";

%feature("docstring") casadi::Integrator::sz_w "[INTERNAL]  Get required
length of w field.

";

%feature("docstring") casadi::Integrator::spCanEvaluate "[INTERNAL]  Is the
class able to propagate seeds through the algorithm?

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::Integrator::sz_iw "[INTERNAL]  Get required
length of iw field.

";

%feature("docstring") casadi::Integrator::hessian "

Generate a Hessian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The generated Hessian has two more outputs than the calling function
corresponding to the Hessian and the gradients.

";

%feature("docstring") casadi::Integrator::jacobian "

Generate a Jacobian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. If
compact is set to true, only the nonzeros of the input and output
expressions are considered. If symmetric is set to true, the Jacobian being
calculated is known to be symmetric (usually a Hessian), which can be
exploited by the algorithm.

The generated Jacobian has one more output than the calling function
corresponding to the Jacobian and the same number of inputs.

";

%feature("docstring") casadi::Integrator::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::Integrator::size1_in "

Get input dimension.

";

%feature("docstring") casadi::Integrator::setOptionByAllowedIndex "[INTERNAL]  Set a certain option by giving its index into the allowed
values.

";

%feature("docstring") casadi::Integrator::n_out "

Get the number of function outputs.

";

%feature("docstring") casadi::Integrator::generate "

Export / Generate C code for the function.

";

%feature("docstring") casadi::Integrator::setDerReverse "

Set a function that calculates nadj adjoint derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::Integrator::printDimensions "

Print dimensions of inputs and outputs.

";

%feature("docstring") casadi::Integrator::sx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::Integrator::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::Integrator::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::Integrator::is_a "

Check if the function is of a particular type Optionally check if name
matches one of the base classes (default true)

";

%feature("docstring") casadi::Integrator::setJacSparsity "

Generate the sparsity of a Jacobian block

";

%feature("docstring") casadi::Integrator::addMonitor "

Add modules to be monitored.

";

%feature("docstring") casadi::Integrator::tangent "

Generate a tangent function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the input must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::Integrator::derForward "

Get a function that calculates nfwd forward derivatives.

Returns a function with n_in + n_out +nfwd*n_in inputs and nfwd*n_out
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nfwd*n_in inputs correspond to forward seeds, one direction at a time The
nfwd*n_out outputs correspond to forward sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::Integrator::call "

Evaluate the function symbolically or numerically.

";

%feature("docstring") casadi::Integrator::getStats "

Get all statistics obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::Integrator::getDAE "

Get the DAE.

";

%feature("docstring") casadi::Integrator::getOptionAllowedIndex "[INTERNAL]
Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::Integrator::fullJacobian "

Generate a Jacobian function of all the inputs elements with respect to all
the output elements).

";

%feature("docstring") casadi::Integrator::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::Integrator::print "

Print a description of the object.

";

%feature("docstring") casadi::Integrator::getStat "

Get a single statistic obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::Integrator::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::Integrator::Integrator "

>  Integrator()
------------------------------------------------------------------------

Default constructor.

>  Integrator(str name, str solver, Function f, Dict opts=Dict())
------------------------------------------------------------------------

Integrator factory (new syntax, includes initialization)

Parameters:
-----------

solver:

Name of a solver. It might be one of:

- cvodes

- idas

- collocation

- oldcollocation

- rk

Note: some of the plugins in this list might not be available on your
system. Also, there might be extra plugins available to you that are not
listed here. You can obtain their documentation with
Integrator.doc(\"myextraplugin\")

f:  dynamical system

>Input scheme: casadi::DAEInput (DAE_NUM_IN = 4) [daeIn]

+-----------+-------+----------------------------+
| Full name | Short |        Description         |
+===========+=======+============================+
| DAE_X     | x     | Differential state .       |
+-----------+-------+----------------------------+
| DAE_Z     | z     | Algebraic state .          |
+-----------+-------+----------------------------+
| DAE_P     | p     | Parameter .                |
+-----------+-------+----------------------------+
| DAE_T     | t     | Explicit time dependence . |
+-----------+-------+----------------------------+

>Output scheme: casadi::DAEOutput (DAE_NUM_OUT = 3) [daeOut]

+-----------+-------+--------------------------------------------+
| Full name | Short |                Description                 |
+===========+=======+============================================+
| DAE_ODE   | ode   | Right hand side of the implicit ODE .      |
+-----------+-------+--------------------------------------------+
| DAE_ALG   | alg   | Right hand side of algebraic equations .   |
+-----------+-------+--------------------------------------------+
| DAE_QUAD  | quad  | Right hand side of quadratures equations . |
+-----------+-------+--------------------------------------------+

";

%feature("docstring") casadi::Integrator::setOutput "

>  void IOInterface< Function  >.setOutput(T val, int oind=0)
------------------------------------------------------------------------

Set an output by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.setOutput(T val, str oname)
------------------------------------------------------------------------

Set an output by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::Integrator::jacSparsity "

Get, if necessary generate, the sparsity of a Jacobian block

";

%feature("docstring") casadi::Integrator::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::Integrator::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::Integrator::dictionary "

Get the dictionary.

";

%feature("docstring") casadi::Integrator::map "

>  [[MX] ] Function.map([[MX ] ] arg, str parallelization=\"serial\")

>  [MX] Function.map([MX ] arg, str parallelization=\"serial\")
------------------------------------------------------------------------

Evaluate symbolically in parallel (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

>  Function Function.map(str name, int N, Dict options=Dict()) const 
------------------------------------------------------------------------

Create a mapped version of this function.

Suppose the function has a signature of:

::

     f: (a, p) -> ( s )
  



The the mapaccumulated version has the signature:

::

     F: (A, P) -> (S )
  
      with
          a: horzcat([a0, a1, ..., a_(N-1)])
          p: horzcat([p0, p1, ..., p_(N-1)])
          s: horzcat([s0, s1, ..., s_(N-1)])
      and
          s0 <- f(a0, p0)
          s1 <- f(a1, p1)
          ...
          s_(N-1) <- f(a_(N-1), p_(N-1))
  



>  Function Function.map(str name, int n, [bool ] repeat_in, [bool ] repeat_out, Dict opts=Dict()) const 
------------------------------------------------------------------------

Generic map.

";

%feature("docstring") casadi::Integrator "

Base class for integrators.

Integrator abstract base class

Solves an initial value problem (IVP) coupled to a terminal value problem
with differential equation given as an implicit ODE coupled to an algebraic
equation and a set of quadratures:

::

   Initial conditions at t=t0
   x(t0)  = x0
   q(t0)  = 0
  
   Forward integration from t=t0 to t=tf
   der(x) = function(x, z, p, t)                  Forward ODE
   0 = fz(x, z, p, t)                  Forward algebraic equations
   der(q) = fq(x, z, p, t)                  Forward quadratures
  
   Terminal conditions at t=tf
   rx(tf)  = rx0
   rq(tf)  = 0
  
   Backward integration from t=tf to t=t0
   der(rx) = gx(rx, rz, rp, x, z, p, t)        Backward ODE
   0 = gz(rx, rz, rp, x, z, p, t)        Backward algebraic equations
   der(rq) = gq(rx, rz, rp, x, z, p, t)        Backward quadratures
  
   where we assume that both the forward and backwards integrations are index-1
   (i.e. dfz/dz, dgz/drz are invertible) and furthermore that
   gx, gz and gq have a linear dependency on rx, rz and rp.



The Integrator class provides some additional functionality, such as getting
the value of the state and/or sensitivities at certain time points.

General information
===================



>Input scheme: casadi::IntegratorInput (INTEGRATOR_NUM_IN = 6) [integratorIn]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| INTEGRATOR_X0          | x0                     | Differential state at  |
|                        |                        | the initial time .     |
+------------------------+------------------------+------------------------+
| INTEGRATOR_P           | p                      | Parameters .           |
+------------------------+------------------------+------------------------+
| INTEGRATOR_Z0          | z0                     | Initial guess for the  |
|                        |                        | algebraic variable .   |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RX0         | rx0                    | Backward differential  |
|                        |                        | state at the final     |
|                        |                        | time .                 |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RP          | rp                     | Backward parameter     |
|                        |                        | vector .               |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RZ0         | rz0                    | Initial guess for the  |
|                        |                        | backwards algebraic    |
|                        |                        | variable .             |
+------------------------+------------------------+------------------------+

>Output scheme: casadi::IntegratorOutput (INTEGRATOR_NUM_OUT = 6) [integratorOut]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| INTEGRATOR_XF          | xf                     | Differential state at  |
|                        |                        | the final time .       |
+------------------------+------------------------+------------------------+
| INTEGRATOR_QF          | qf                     | Quadrature state at    |
|                        |                        | the final time .       |
+------------------------+------------------------+------------------------+
| INTEGRATOR_ZF          | zf                     | Algebraic variable at  |
|                        |                        | the final time .       |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RXF         | rxf                    | Backward differential  |
|                        |                        | state at the initial   |
|                        |                        | time .                 |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RQF         | rqf                    | Backward quadrature    |
|                        |                        | state at the initial   |
|                        |                        | time .                 |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RZF         | rzf                    | Backward algebraic     |
|                        |                        | variable at the        |
|                        |                        | initial time .         |
+------------------------+------------------------+------------------------+

>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| ad_weight    | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | derivative c |              |
|              |              |              | alculation.W |              |
|              |              |              | hen there is |              |
|              |              |              | an option of |              |
|              |              |              | either using |              |
|              |              |              | forward or   |              |
|              |              |              | reverse mode |              |
|              |              |              | directional  |              |
|              |              |              | derivatives, |              |
|              |              |              | the          |              |
|              |              |              | condition ad |              |
|              |              |              | _weight*nf<= |              |
|              |              |              | (1-ad_weight |              |
|              |              |              | )*na is used |              |
|              |              |              | where nf and |              |
|              |              |              | na are       |              |
|              |              |              | estimates of |              |
|              |              |              | the number   |              |
|              |              |              | of forward/r |              |
|              |              |              | everse mode  |              |
|              |              |              | directional  |              |
|              |              |              | derivatives  |              |
|              |              |              | needed. By   |              |
|              |              |              | default,     |              |
|              |              |              | ad_weight is |              |
|              |              |              | calculated a |              |
|              |              |              | utomatically |              |
|              |              |              | , but this   |              |
|              |              |              | can be       |              |
|              |              |              | overridden   |              |
|              |              |              | by setting   |              |
|              |              |              | this option. |              |
|              |              |              | In           |              |
|              |              |              | particular,  |              |
|              |              |              | 0 means      |              |
|              |              |              | forcing      |              |
|              |              |              | forward mode |              |
|              |              |              | and 1        |              |
|              |              |              | forcing      |              |
|              |              |              | reverse      |              |
|              |              |              | mode. Leave  |              |
|              |              |              | unset for    |              |
|              |              |              | (class       |              |
|              |              |              | specific)    |              |
|              |              |              | heuristics.  |              |
+--------------+--------------+--------------+--------------+--------------+
| ad_weight_sp | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | sparsity     |              |
|              |              |              | pattern      |              |
|              |              |              | calculation  |              |
|              |              |              | calculation. |              |
|              |              |              | Overrides    |              |
|              |              |              | default      |              |
|              |              |              | behavior.    |              |
|              |              |              | Set to 0 and |              |
|              |              |              | 1 to force   |              |
|              |              |              | forward and  |              |
|              |              |              | reverse mode |              |
|              |              |              | respectively |              |
|              |              |              | . Cf. option |              |
|              |              |              | \"ad_weight\". |              |
+--------------+--------------+--------------+--------------+--------------+
| augmented_op | OT_DICT      | GenericType( | Options to   | casadi::Inte |
| tions        |              | )            | be passed    | gratorIntern |
|              |              |              | down to the  | al           |
|              |              |              | augmented    |              |
|              |              |              | integrator,  |              |
|              |              |              | if one is    |              |
|              |              |              | constructed. |              |
+--------------+--------------+--------------+--------------+--------------+
| compiler     | OT_STRING    | \"clang\"      | Just-in-time | casadi::Func |
|              |              |              | compiler     | tionInternal |
|              |              |              | plugin to be |              |
|              |              |              | used.        |              |
+--------------+--------------+--------------+--------------+--------------+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode     |
|              |              |              | according to |              |
|              |              |              | a given      |              |
|              |              |              | recipe (low- |              |
|              |              |              | level)       |              |
+--------------+--------------+--------------+--------------+--------------+
| expand_augme | OT_BOOLEAN   | true         | If DAE       | casadi::Inte |
| nted         |              |              | callback     | gratorIntern |
|              |              |              | functions    | al           |
|              |              |              | are          |              |
|              |              |              | SXFunction,  |              |
|              |              |              | have         |              |
|              |              |              | augmented    |              |
|              |              |              | DAE callback |              |
|              |              |              | function     |              |
|              |              |              | also be      |              |
|              |              |              | SXFunction.  |              |
+--------------+--------------+--------------+--------------+--------------+
| gather_stats | OT_BOOLEAN   | false        | Flag to      | casadi::Func |
|              |              |              | indicate     | tionInternal |
|              |              |              | whether      |              |
|              |              |              | statistics   |              |
|              |              |              | must be      |              |
|              |              |              | gathered     |              |
+--------------+--------------+--------------+--------------+--------------+
| input_scheme | OT_STRINGVEC | GenericType( | Custom input | casadi::Func |
|              | TOR          | )            | scheme       | tionInternal |
+--------------+--------------+--------------+--------------+--------------+
| inputs_check | OT_BOOLEAN   | true         | Throw        | casadi::Func |
|              |              |              | exceptions   | tionInternal |
|              |              |              | when the     |              |
|              |              |              | numerical    |              |
|              |              |              | values of    |              |
|              |              |              | the inputs   |              |
|              |              |              | don't make   |              |
|              |              |              | sense        |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_penalty  | OT_REAL      | 2            | When         | casadi::Func |
|              |              |              | requested    | tionInternal |
|              |              |              | for a number |              |
|              |              |              | of forward/r |              |
|              |              |              | everse       |              |
|              |              |              | directions,  |              |
|              |              |              | it may be    |              |
|              |              |              | cheaper to   |              |
|              |              |              | compute      |              |
|              |              |              | first the    |              |
|              |              |              | full         |              |
|              |              |              | jacobian and |              |
|              |              |              | then         |              |
|              |              |              | multiply     |              |
|              |              |              | with seeds,  |              |
|              |              |              | rather than  |              |
|              |              |              | obtain the   |              |
|              |              |              | requested    |              |
|              |              |              | directions   |              |
|              |              |              | in a straigh |              |
|              |              |              | tforward     |              |
|              |              |              | manner.      |              |
|              |              |              | Casadi uses  |              |
|              |              |              | a heuristic  |              |
|              |              |              | to decide    |              |
|              |              |              | which is     |              |
|              |              |              | cheaper. A   |              |
|              |              |              | high value   |              |
|              |              |              | of 'jac_pena |              |
|              |              |              | lty' makes   |              |
|              |              |              | it less      |              |
|              |              |              | likely for   |              |
|              |              |              | the heurstic |              |
|              |              |              | to chose the |              |
|              |              |              | full         |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy.    |              |
|              |              |              | The special  |              |
|              |              |              | value -1     |              |
|              |              |              | indicates    |              |
|              |              |              | never to use |              |
|              |              |              | the full     |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy     |              |
+--------------+--------------+--------------+--------------+--------------+
| jit          | OT_BOOLEAN   | false        | Use just-in- | casadi::Func |
|              |              |              | time         | tionInternal |
|              |              |              | compiler to  |              |
|              |              |              | speed up the |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| jit_options  | OT_DICT      | GenericType( | Options to   | casadi::Func |
|              |              | )            | be passed to | tionInternal |
|              |              |              | the jit      |              |
|              |              |              | compiler.    |              |
+--------------+--------------+--------------+--------------+--------------+
| monitor      | OT_STRINGVEC | GenericType( | Monitors to  | casadi::Func |
|              | TOR          | )            | be activated | tionInternal |
|              |              |              | (inputs|outp |              |
|              |              |              | uts)         |              |
+--------------+--------------+--------------+--------------+--------------+
| output_schem | OT_STRINGVEC | GenericType( | Custom       | casadi::Func |
| e            | TOR          | )            | output       | tionInternal |
|              |              |              | scheme       |              |
+--------------+--------------+--------------+--------------+--------------+
| print_stats  | OT_BOOLEAN   | false        | Print out    | casadi::Inte |
|              |              |              | statistics   | gratorIntern |
|              |              |              | after        | al           |
|              |              |              | integration  |              |
+--------------+--------------+--------------+--------------+--------------+
| regularity_c | OT_BOOLEAN   | true         | Throw        | casadi::Func |
| heck         |              |              | exceptions   | tionInternal |
|              |              |              | when NaN or  |              |
|              |              |              | Inf appears  |              |
|              |              |              | during       |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| t0           | OT_REAL      | 0            | Beginning of | casadi::Inte |
|              |              |              | the time     | gratorIntern |
|              |              |              | horizon      | al           |
+--------------+--------------+--------------+--------------+--------------+
| tf           | OT_REAL      | 1            | End of the   | casadi::Inte |
|              |              |              | time horizon | gratorIntern |
|              |              |              |              | al           |
+--------------+--------------+--------------+--------------+--------------+
| user_data    | OT_VOIDPTR   | GenericType( | A user-      | casadi::Func |
|              |              | )            | defined      | tionInternal |
|              |              |              | field that   |              |
|              |              |              | can be used  |              |
|              |              |              | to identify  |              |
|              |              |              | the function |              |
|              |              |              | or pass      |              |
|              |              |              | additional   |              |
|              |              |              | information  |              |
+--------------+--------------+--------------+--------------+--------------+
| verbose      | OT_BOOLEAN   | false        | Verbose      | casadi::Func |
|              |              |              | evaluation   | tionInternal |
|              |              |              | for          |              |
|              |              |              | debugging    |              |
+--------------+--------------+--------------+--------------+--------------+

List of plugins
===============



- cvodes

- idas

- collocation

- oldcollocation

- rk

Note: some of the plugins in this list might not be available on your
system. Also, there might be extra plugins available to you that are not
listed here. You can obtain their documentation with
Integrator.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

cvodes
------



Interface to CVodes from the Sundials suite.

A call to evaluate will integrate to the end.

You can retrieve the entire state trajectory as follows, after the evaluate
call: Call reset. Then call integrate(t_i) and getOuput for a series of
times t_i.

Note: depending on the dimension and structure of your problem, you may
experience a dramatic speed-up by using a sparse linear solver:



::

     intg.setOption(\"linear_solver\",\"csparse\")
     intg.setOption(\"linear_solver_type\",\"user_defined\")



>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| abstol          | OT_REAL         | 0.000           | Absolute        |
|                 |                 |                 | tolerence for   |
|                 |                 |                 | the IVP         |
|                 |                 |                 | solution        |
+-----------------+-----------------+-----------------+-----------------+
| abstolB         | OT_REAL         | GenericType()   | Absolute        |
|                 |                 |                 | tolerence for   |
|                 |                 |                 | the adjoint     |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | solution        |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to abstol]      |
+-----------------+-----------------+-----------------+-----------------+
| disable_interna | OT_BOOLEAN      | false           | Disable CVodes  |
| l_warnings      |                 |                 | internal        |
|                 |                 |                 | warning         |
|                 |                 |                 | messages        |
+-----------------+-----------------+-----------------+-----------------+
| exact_jacobian  | OT_BOOLEAN      | true            | Use exact       |
|                 |                 |                 | Jacobian        |
|                 |                 |                 | information for |
|                 |                 |                 | the forward     |
|                 |                 |                 | integration     |
+-----------------+-----------------+-----------------+-----------------+
| exact_jacobianB | OT_BOOLEAN      | GenericType()   | Use exact       |
|                 |                 |                 | Jacobian        |
|                 |                 |                 | information for |
|                 |                 |                 | the backward    |
|                 |                 |                 | integration     |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to              |
|                 |                 |                 | exact_jacobian] |
+-----------------+-----------------+-----------------+-----------------+
| finite_differen | OT_BOOLEAN      | false           | Use finite      |
| ce_fsens        |                 |                 | differences to  |
|                 |                 |                 | approximate the |
|                 |                 |                 | forward         |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | equations (if   |
|                 |                 |                 | AD is not       |
|                 |                 |                 | available)      |
+-----------------+-----------------+-----------------+-----------------+
| fsens_abstol    | OT_REAL         | GenericType()   | Absolute        |
|                 |                 |                 | tolerence for   |
|                 |                 |                 | the forward     |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | solution        |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to abstol]      |
+-----------------+-----------------+-----------------+-----------------+
| fsens_all_at_on | OT_BOOLEAN      | true            | Calculate all   |
| ce              |                 |                 | right hand      |
|                 |                 |                 | sides of the    |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | equations at    |
|                 |                 |                 | once            |
+-----------------+-----------------+-----------------+-----------------+
| fsens_err_con   | OT_BOOLEAN      | true            | include the     |
|                 |                 |                 | forward         |
|                 |                 |                 | sensitivities   |
|                 |                 |                 | in all error    |
|                 |                 |                 | controls        |
+-----------------+-----------------+-----------------+-----------------+
| fsens_reltol    | OT_REAL         | GenericType()   | Relative        |
|                 |                 |                 | tolerence for   |
|                 |                 |                 | the forward     |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | solution        |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to reltol]      |
+-----------------+-----------------+-----------------+-----------------+
| fsens_scaling_f | OT_REALVECTOR   | GenericType()   | Scaling factor  |
| actors          |                 |                 | for the         |
|                 |                 |                 | components if   |
|                 |                 |                 | finite          |
|                 |                 |                 | differences is  |
|                 |                 |                 | used            |
+-----------------+-----------------+-----------------+-----------------+
| fsens_sensitivi | OT_INTEGERVECTO | GenericType()   | Specifies which |
| y_parameters    | R               |                 | components will |
|                 |                 |                 | be used when    |
|                 |                 |                 | estimating the  |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | equations       |
+-----------------+-----------------+-----------------+-----------------+
| interpolation_t | OT_STRING       | \"hermite\"       | Type of         |
| ype             |                 |                 | interpolation   |
|                 |                 |                 | for the adjoint |
|                 |                 |                 | sensitivities ( |
|                 |                 |                 | hermite|polynom |
|                 |                 |                 | ial)            |
+-----------------+-----------------+-----------------+-----------------+
| iterative_solve | OT_STRING       | \"gmres\"         | (gmres|bcgstab| |
| r               |                 |                 | tfqmr)          |
+-----------------+-----------------+-----------------+-----------------+
| iterative_solve | OT_STRING       | GenericType()   | (gmres|bcgstab| |
| rB              |                 |                 | tfqmr)          |
+-----------------+-----------------+-----------------+-----------------+
| linear_multiste | OT_STRING       | \"bdf\"           | Integrator      |
| p_method        |                 |                 | scheme          |
|                 |                 |                 | (bdf|adams)     |
+-----------------+-----------------+-----------------+-----------------+
| linear_solver   | OT_STRING       | GenericType()   | A custom linear |
|                 |                 |                 | solver creator  |
|                 |                 |                 | function        |
+-----------------+-----------------+-----------------+-----------------+
| linear_solverB  | OT_STRING       | GenericType()   | A custom linear |
|                 |                 |                 | solver creator  |
|                 |                 |                 | function for    |
|                 |                 |                 | backwards       |
|                 |                 |                 | integration     |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to              |
|                 |                 |                 | linear_solver]  |
+-----------------+-----------------+-----------------+-----------------+
| linear_solver_o | OT_DICT         | GenericType()   | Options to be   |
| ptions          |                 |                 | passed to the   |
|                 |                 |                 | linear solver   |
+-----------------+-----------------+-----------------+-----------------+
| linear_solver_o | OT_DICT         | GenericType()   | Options to be   |
| ptionsB         |                 |                 | passed to the   |
|                 |                 |                 | linear solver   |
|                 |                 |                 | for backwards   |
|                 |                 |                 | integration     |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to linear_solve |
|                 |                 |                 | r_options]      |
+-----------------+-----------------+-----------------+-----------------+
| linear_solver_t | OT_STRING       | \"dense\"         | (user_defined|d |
| ype             |                 |                 | ense|banded|ite |
|                 |                 |                 | rative)         |
+-----------------+-----------------+-----------------+-----------------+
| linear_solver_t | OT_STRING       | GenericType()   | (user_defined|d |
| ypeB            |                 |                 | ense|banded|ite |
|                 |                 |                 | rative)         |
+-----------------+-----------------+-----------------+-----------------+
| lower_bandwidth | OT_INTEGER      | GenericType()   | Lower band-     |
|                 |                 |                 | width of banded |
|                 |                 |                 | Jacobian        |
|                 |                 |                 | (estimations)   |
+-----------------+-----------------+-----------------+-----------------+
| lower_bandwidth | OT_INTEGER      | GenericType()   | lower band-     |
| B               |                 |                 | width of banded |
|                 |                 |                 | jacobians for   |
|                 |                 |                 | backward        |
|                 |                 |                 | integration     |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to lower_bandwi |
|                 |                 |                 | dth]            |
+-----------------+-----------------+-----------------+-----------------+
| max_krylov      | OT_INTEGER      | 10              | Maximum Krylov  |
|                 |                 |                 | subspace size   |
+-----------------+-----------------+-----------------+-----------------+
| max_krylovB     | OT_INTEGER      | GenericType()   | Maximum krylov  |
|                 |                 |                 | subspace size   |
+-----------------+-----------------+-----------------+-----------------+
| max_multistep_o | OT_INTEGER      | 5               |                 |
| rder            |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| max_num_steps   | OT_INTEGER      | 10000           | Maximum number  |
|                 |                 |                 | of integrator   |
|                 |                 |                 | steps           |
+-----------------+-----------------+-----------------+-----------------+
| nonlinear_solve | OT_STRING       | \"newton\"        | (newton|functio |
| r_iteration     |                 |                 | nal)            |
+-----------------+-----------------+-----------------+-----------------+
| pretype         | OT_STRING       | \"none\"          | (none|left|righ |
|                 |                 |                 | t|both)         |
+-----------------+-----------------+-----------------+-----------------+
| pretypeB        | OT_STRING       | GenericType()   | (none|left|righ |
|                 |                 |                 | t|both)         |
+-----------------+-----------------+-----------------+-----------------+
| quad_err_con    | OT_BOOLEAN      | false           | Should the      |
|                 |                 |                 | quadratures     |
|                 |                 |                 | affect the step |
|                 |                 |                 | size control    |
+-----------------+-----------------+-----------------+-----------------+
| reltol          | OT_REAL         | 0.000           | Relative        |
|                 |                 |                 | tolerence for   |
|                 |                 |                 | the IVP         |
|                 |                 |                 | solution        |
+-----------------+-----------------+-----------------+-----------------+
| reltolB         | OT_REAL         | GenericType()   | Relative        |
|                 |                 |                 | tolerence for   |
|                 |                 |                 | the adjoint     |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | solution        |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to reltol]      |
+-----------------+-----------------+-----------------+-----------------+
| sensitivity_met | OT_STRING       | \"simultaneous\"  | (simultaneous|s |
| hod             |                 |                 | taggered)       |
+-----------------+-----------------+-----------------+-----------------+
| steps_per_check | OT_INTEGER      | 20              | Number of steps |
| point           |                 |                 | between two     |
|                 |                 |                 | consecutive     |
|                 |                 |                 | checkpoints     |
+-----------------+-----------------+-----------------+-----------------+
| stop_at_end     | OT_BOOLEAN      | true            | Stop the        |
|                 |                 |                 | integrator at   |
|                 |                 |                 | the end of the  |
|                 |                 |                 | interval        |
+-----------------+-----------------+-----------------+-----------------+
| upper_bandwidth | OT_INTEGER      | GenericType()   | Upper band-     |
|                 |                 |                 | width of banded |
|                 |                 |                 | Jacobian        |
|                 |                 |                 | (estimations)   |
+-----------------+-----------------+-----------------+-----------------+
| upper_bandwidth | OT_INTEGER      | GenericType()   | Upper band-     |
| B               |                 |                 | width of banded |
|                 |                 |                 | jacobians for   |
|                 |                 |                 | backward        |
|                 |                 |                 | integration     |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to upper_bandwi |
|                 |                 |                 | dth]            |
+-----------------+-----------------+-----------------+-----------------+
| use_preconditio | OT_BOOLEAN      | false           | Precondition an |
| ner             |                 |                 | iterative       |
|                 |                 |                 | solver          |
+-----------------+-----------------+-----------------+-----------------+
| use_preconditio | OT_BOOLEAN      | GenericType()   | Precondition an |
| nerB            |                 |                 | iterative       |
|                 |                 |                 | solver for the  |
|                 |                 |                 | backwards       |
|                 |                 |                 | problem         |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to use_precondi |
|                 |                 |                 | tioner]         |
+-----------------+-----------------+-----------------+-----------------+

>List of available monitors

+---------+
|   Id    |
+=========+
| djacB   |
+---------+
| psetupB |
+---------+
| res     |
+---------+
| resB    |
+---------+
| resQB   |
+---------+
| reset   |
+---------+

>List of available stats

+-------------+
|     Id      |
+=============+
| nlinsetups  |
+-------------+
| nlinsetupsB |
+-------------+
| nsteps      |
+-------------+
| nstepsB     |
+-------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

idas
----



Interface to IDAS from the Sundials suite.

Note: depending on the dimension and structure of your problem, you may
experience a dramatic speed-up by using a sparse linear solver:



::

     intg.setOption(\"linear_solver\",\"csparse\")
     intg.setOption(\"linear_solver_type\",\"user_defined\")



>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| abstol          | OT_REAL         | 0.000           | Absolute        |
|                 |                 |                 | tolerence for   |
|                 |                 |                 | the IVP         |
|                 |                 |                 | solution        |
+-----------------+-----------------+-----------------+-----------------+
| abstolB         | OT_REAL         | GenericType()   | Absolute        |
|                 |                 |                 | tolerence for   |
|                 |                 |                 | the adjoint     |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | solution        |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to abstol]      |
+-----------------+-----------------+-----------------+-----------------+
| abstolv         | OT_REALVECTOR   |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| calc_ic         | OT_BOOLEAN      | true            | Use IDACalcIC   |
|                 |                 |                 | to get          |
|                 |                 |                 | consistent      |
|                 |                 |                 | initial         |
|                 |                 |                 | conditions.     |
+-----------------+-----------------+-----------------+-----------------+
| calc_icB        | OT_BOOLEAN      | GenericType()   | Use IDACalcIC   |
|                 |                 |                 | to get          |
|                 |                 |                 | consistent      |
|                 |                 |                 | initial         |
|                 |                 |                 | conditions for  |
|                 |                 |                 | backwards       |
|                 |                 |                 | system          |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to calc_ic].    |
+-----------------+-----------------+-----------------+-----------------+
| cj_scaling      | OT_BOOLEAN      | false           | IDAS scaling on |
|                 |                 |                 | cj for the      |
|                 |                 |                 | user-defined    |
|                 |                 |                 | linear solver   |
|                 |                 |                 | module          |
+-----------------+-----------------+-----------------+-----------------+
| disable_interna | OT_BOOLEAN      | false           | Disable IDAS    |
| l_warnings      |                 |                 | internal        |
|                 |                 |                 | warning         |
|                 |                 |                 | messages        |
+-----------------+-----------------+-----------------+-----------------+
| exact_jacobian  | OT_BOOLEAN      | true            | Use exact       |
|                 |                 |                 | Jacobian        |
|                 |                 |                 | information for |
|                 |                 |                 | the forward     |
|                 |                 |                 | integration     |
+-----------------+-----------------+-----------------+-----------------+
| exact_jacobianB | OT_BOOLEAN      | GenericType()   | Use exact       |
|                 |                 |                 | Jacobian        |
|                 |                 |                 | information for |
|                 |                 |                 | the backward    |
|                 |                 |                 | integration     |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to              |
|                 |                 |                 | exact_jacobian] |
+-----------------+-----------------+-----------------+-----------------+
| extra_fsens_cal | OT_BOOLEAN      | false           | Call calc ic an |
| c_ic            |                 |                 | extra time,     |
|                 |                 |                 | with fsens=0    |
+-----------------+-----------------+-----------------+-----------------+
| finite_differen | OT_BOOLEAN      | false           | Use finite      |
| ce_fsens        |                 |                 | differences to  |
|                 |                 |                 | approximate the |
|                 |                 |                 | forward         |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | equations (if   |
|                 |                 |                 | AD is not       |
|                 |                 |                 | available)      |
+-----------------+-----------------+-----------------+-----------------+
| first_time      | OT_REAL         | GenericType()   | First requested |
|                 |                 |                 | time as a       |
|                 |                 |                 | fraction of the |
|                 |                 |                 | time interval   |
+-----------------+-----------------+-----------------+-----------------+
| fsens_abstol    | OT_REAL         | GenericType()   | Absolute        |
|                 |                 |                 | tolerence for   |
|                 |                 |                 | the forward     |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | solution        |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to abstol]      |
+-----------------+-----------------+-----------------+-----------------+
| fsens_abstolv   | OT_REALVECTOR   |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| fsens_err_con   | OT_BOOLEAN      | true            | include the     |
|                 |                 |                 | forward         |
|                 |                 |                 | sensitivities   |
|                 |                 |                 | in all error    |
|                 |                 |                 | controls        |
+-----------------+-----------------+-----------------+-----------------+
| fsens_reltol    | OT_REAL         | GenericType()   | Relative        |
|                 |                 |                 | tolerence for   |
|                 |                 |                 | the forward     |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | solution        |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to reltol]      |
+-----------------+-----------------+-----------------+-----------------+
| fsens_scaling_f | OT_REALVECTOR   | GenericType()   | Scaling factor  |
| actors          |                 |                 | for the         |
|                 |                 |                 | components if   |
|                 |                 |                 | finite          |
|                 |                 |                 | differences is  |
|                 |                 |                 | used            |
+-----------------+-----------------+-----------------+-----------------+
| fsens_sensitivi | OT_INTEGERVECTO | GenericType()   | Specifies which |
| y_parameters    | R               |                 | components will |
|                 |                 |                 | be used when    |
|                 |                 |                 | estimating the  |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | equations       |
+-----------------+-----------------+-----------------+-----------------+
| init_xdot       | OT_REALVECTOR   | GenericType()   | Initial values  |
|                 |                 |                 | for the state   |
|                 |                 |                 | derivatives     |
+-----------------+-----------------+-----------------+-----------------+
| interpolation_t | OT_STRING       | \"hermite\"       | Type of         |
| ype             |                 |                 | interpolation   |
|                 |                 |                 | for the adjoint |
|                 |                 |                 | sensitivities ( |
|                 |                 |                 | hermite|polynom |
|                 |                 |                 | ial)            |
+-----------------+-----------------+-----------------+-----------------+
| iterative_solve | OT_STRING       | \"gmres\"         | (gmres|bcgstab| |
| r               |                 |                 | tfqmr)          |
+-----------------+-----------------+-----------------+-----------------+
| iterative_solve | OT_STRING       | GenericType()   | (gmres|bcgstab| |
| rB              |                 |                 | tfqmr)          |
+-----------------+-----------------+-----------------+-----------------+
| linear_solver   | OT_STRING       | GenericType()   | A custom linear |
|                 |                 |                 | solver creator  |
|                 |                 |                 | function        |
+-----------------+-----------------+-----------------+-----------------+
| linear_solverB  | OT_STRING       | GenericType()   | A custom linear |
|                 |                 |                 | solver creator  |
|                 |                 |                 | function for    |
|                 |                 |                 | backwards       |
|                 |                 |                 | integration     |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to              |
|                 |                 |                 | linear_solver]  |
+-----------------+-----------------+-----------------+-----------------+
| linear_solver_o | OT_DICT         | GenericType()   | Options to be   |
| ptions          |                 |                 | passed to the   |
|                 |                 |                 | linear solver   |
+-----------------+-----------------+-----------------+-----------------+
| linear_solver_o | OT_DICT         | GenericType()   | Options to be   |
| ptionsB         |                 |                 | passed to the   |
|                 |                 |                 | linear solver   |
|                 |                 |                 | for backwards   |
|                 |                 |                 | integration     |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to linear_solve |
|                 |                 |                 | r_options]      |
+-----------------+-----------------+-----------------+-----------------+
| linear_solver_t | OT_STRING       | \"dense\"         | (user_defined|d |
| ype             |                 |                 | ense|banded|ite |
|                 |                 |                 | rative)         |
+-----------------+-----------------+-----------------+-----------------+
| linear_solver_t | OT_STRING       | GenericType()   | (user_defined|d |
| ypeB            |                 |                 | ense|banded|ite |
|                 |                 |                 | rative)         |
+-----------------+-----------------+-----------------+-----------------+
| lower_bandwidth | OT_INTEGER      | GenericType()   | Lower band-     |
|                 |                 |                 | width of banded |
|                 |                 |                 | Jacobian        |
|                 |                 |                 | (estimations)   |
+-----------------+-----------------+-----------------+-----------------+
| lower_bandwidth | OT_INTEGER      | GenericType()   | lower band-     |
| B               |                 |                 | width of banded |
|                 |                 |                 | jacobians for   |
|                 |                 |                 | backward        |
|                 |                 |                 | integration     |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to lower_bandwi |
|                 |                 |                 | dth]            |
+-----------------+-----------------+-----------------+-----------------+
| max_krylov      | OT_INTEGER      | 10              | Maximum Krylov  |
|                 |                 |                 | subspace size   |
+-----------------+-----------------+-----------------+-----------------+
| max_krylovB     | OT_INTEGER      | GenericType()   | Maximum krylov  |
|                 |                 |                 | subspace size   |
+-----------------+-----------------+-----------------+-----------------+
| max_multistep_o | OT_INTEGER      | 5               |                 |
| rder            |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| max_num_steps   | OT_INTEGER      | 10000           | Maximum number  |
|                 |                 |                 | of integrator   |
|                 |                 |                 | steps           |
+-----------------+-----------------+-----------------+-----------------+
| max_step_size   | OT_REAL         | 0               | Maximim step    |
|                 |                 |                 | size            |
+-----------------+-----------------+-----------------+-----------------+
| pretype         | OT_STRING       | \"none\"          | (none|left|righ |
|                 |                 |                 | t|both)         |
+-----------------+-----------------+-----------------+-----------------+
| pretypeB        | OT_STRING       | GenericType()   | (none|left|righ |
|                 |                 |                 | t|both)         |
+-----------------+-----------------+-----------------+-----------------+
| quad_err_con    | OT_BOOLEAN      | false           | Should the      |
|                 |                 |                 | quadratures     |
|                 |                 |                 | affect the step |
|                 |                 |                 | size control    |
+-----------------+-----------------+-----------------+-----------------+
| reltol          | OT_REAL         | 0.000           | Relative        |
|                 |                 |                 | tolerence for   |
|                 |                 |                 | the IVP         |
|                 |                 |                 | solution        |
+-----------------+-----------------+-----------------+-----------------+
| reltolB         | OT_REAL         | GenericType()   | Relative        |
|                 |                 |                 | tolerence for   |
|                 |                 |                 | the adjoint     |
|                 |                 |                 | sensitivity     |
|                 |                 |                 | solution        |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to reltol]      |
+-----------------+-----------------+-----------------+-----------------+
| sensitivity_met | OT_STRING       | \"simultaneous\"  | (simultaneous|s |
| hod             |                 |                 | taggered)       |
+-----------------+-----------------+-----------------+-----------------+
| steps_per_check | OT_INTEGER      | 20              | Number of steps |
| point           |                 |                 | between two     |
|                 |                 |                 | consecutive     |
|                 |                 |                 | checkpoints     |
+-----------------+-----------------+-----------------+-----------------+
| stop_at_end     | OT_BOOLEAN      | true            | Stop the        |
|                 |                 |                 | integrator at   |
|                 |                 |                 | the end of the  |
|                 |                 |                 | interval        |
+-----------------+-----------------+-----------------+-----------------+
| suppress_algebr | OT_BOOLEAN      | false           | Suppress        |
| aic             |                 |                 | algebraic       |
|                 |                 |                 | variables in    |
|                 |                 |                 | the error       |
|                 |                 |                 | testing         |
+-----------------+-----------------+-----------------+-----------------+
| upper_bandwidth | OT_INTEGER      | GenericType()   | Upper band-     |
|                 |                 |                 | width of banded |
|                 |                 |                 | Jacobian        |
|                 |                 |                 | (estimations)   |
+-----------------+-----------------+-----------------+-----------------+
| upper_bandwidth | OT_INTEGER      | GenericType()   | Upper band-     |
| B               |                 |                 | width of banded |
|                 |                 |                 | jacobians for   |
|                 |                 |                 | backward        |
|                 |                 |                 | integration     |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to upper_bandwi |
|                 |                 |                 | dth]            |
+-----------------+-----------------+-----------------+-----------------+
| use_preconditio | OT_BOOLEAN      | false           | Precondition an |
| ner             |                 |                 | iterative       |
|                 |                 |                 | solver          |
+-----------------+-----------------+-----------------+-----------------+
| use_preconditio | OT_BOOLEAN      | GenericType()   | Precondition an |
| nerB            |                 |                 | iterative       |
|                 |                 |                 | solver for the  |
|                 |                 |                 | backwards       |
|                 |                 |                 | problem         |
|                 |                 |                 | [default: equal |
|                 |                 |                 | to use_precondi |
|                 |                 |                 | tioner]         |
+-----------------+-----------------+-----------------+-----------------+

>List of available monitors

+--------------------------+
|            Id            |
+==========================+
| bjacB                    |
+--------------------------+
| correctInitialConditions |
+--------------------------+
| jtimesB                  |
+--------------------------+
| psetup                   |
+--------------------------+
| psetupB                  |
+--------------------------+
| psolveB                  |
+--------------------------+
| res                      |
+--------------------------+
| resB                     |
+--------------------------+
| resS                     |
+--------------------------+
| rhsQB                    |
+--------------------------+

>List of available stats

+-------------+
|     Id      |
+=============+
| nlinsetups  |
+-------------+
| nlinsetupsB |
+-------------+
| nsteps      |
+-------------+
| nstepsB     |
+-------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

collocation
-----------



Fixed-step implicit Runge-Kutta integrator ODE/DAE integrator based on
collocation schemes

The method is still under development

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| collocation_sch | OT_STRING       | \"radau\"         | Collocation     |
| eme             |                 |                 | scheme (radau|l |
|                 |                 |                 | egendre)        |
+-----------------+-----------------+-----------------+-----------------+
| implicit_solver | OT_STRING       | GenericType()   | An implicit     |
|                 |                 |                 | function solver |
+-----------------+-----------------+-----------------+-----------------+
| implicit_solver | OT_DICT         | GenericType()   | Options to be   |
| _options        |                 |                 | passed to the   |
|                 |                 |                 | NLP Solver      |
+-----------------+-----------------+-----------------+-----------------+
| interpolation_o | OT_INTEGER      | 3               | Order of the    |
| rder            |                 |                 | interpolating   |
|                 |                 |                 | polynomials     |
+-----------------+-----------------+-----------------+-----------------+
| number_of_finit | OT_INTEGER      | 20              | Number of       |
| e_elements      |                 |                 | finite elements |
+-----------------+-----------------+-----------------+-----------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

oldcollocation
--------------



Collocation integrator ODE/DAE integrator based on collocation

The method is still under development

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| collocation_sch | OT_STRING       | \"radau\"         | Collocation     |
| eme             |                 |                 | scheme (radau|l |
|                 |                 |                 | egendre)        |
+-----------------+-----------------+-----------------+-----------------+
| expand_f        | OT_BOOLEAN      | false           | Expand the      |
|                 |                 |                 | ODE/DAE         |
|                 |                 |                 | residual        |
|                 |                 |                 | function in an  |
|                 |                 |                 | SX graph        |
+-----------------+-----------------+-----------------+-----------------+
| expand_q        | OT_BOOLEAN      | false           | Expand the      |
|                 |                 |                 | quadrature      |
|                 |                 |                 | function in an  |
|                 |                 |                 | SX graph        |
+-----------------+-----------------+-----------------+-----------------+
| hotstart        | OT_BOOLEAN      | true            | Initialize the  |
|                 |                 |                 | trajectory at   |
|                 |                 |                 | the previous    |
|                 |                 |                 | solution        |
+-----------------+-----------------+-----------------+-----------------+
| implicit_solver | OT_STRING       | GenericType()   | An implicit     |
|                 |                 |                 | function solver |
+-----------------+-----------------+-----------------+-----------------+
| implicit_solver | OT_DICT         | GenericType()   | Options to be   |
| _options        |                 |                 | passed to the   |
|                 |                 |                 | implicit solver |
+-----------------+-----------------+-----------------+-----------------+
| interpolation_o | OT_INTEGER      | 3               | Order of the    |
| rder            |                 |                 | interpolating   |
|                 |                 |                 | polynomials     |
+-----------------+-----------------+-----------------+-----------------+
| number_of_finit | OT_INTEGER      | 20              | Number of       |
| e_elements      |                 |                 | finite elements |
+-----------------+-----------------+-----------------+-----------------+
| startup_integra | OT_STRING       | GenericType()   | An ODE/DAE      |
| tor             |                 |                 | integrator that |
|                 |                 |                 | can be used to  |
|                 |                 |                 | generate a      |
|                 |                 |                 | startup         |
|                 |                 |                 | trajectory      |
+-----------------+-----------------+-----------------+-----------------+
| startup_integra | OT_DICT         | GenericType()   | Options to be   |
| tor_options     |                 |                 | passed to the   |
|                 |                 |                 | startup         |
|                 |                 |                 | integrator      |
+-----------------+-----------------+-----------------+-----------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

rk --



Fixed-step explicit Runge-Kutta integrator for ODEs Currently implements
RK4.

The method is still under development

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| number_of_finit | OT_INTEGER      | 20              | Number of       |
| e_elements      |                 |                 | finite elements |
+-----------------+-----------------+-----------------+-----------------+

--------------------------------------------------------------------------------



Joel Andersson
Diagrams
--------



C++ includes: integrator.hpp ";

%feature("docstring") casadi::Integrator::free_sx "

Get all the free variables of the function.

";

%feature("docstring") casadi::Integrator::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::Integrator::setStopTime "

Set a stop time for the forward integration.

";

%feature("docstring") casadi::Integrator::type_name "

Get type name.

";

%feature("docstring") casadi::Integrator::getAtomicOperation "

Get an atomic operation operator index.

";

%feature("docstring") casadi::Integrator::getAtomicOutput "

Get the (integer) output argument of an atomic operation.

";

%feature("docstring") casadi::Integrator::removeMonitor "

Remove modules to be monitored.

";

%feature("docstring") casadi::Integrator::setDerForward "

Set a function that calculates nfwd forward derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::Integrator::sparsity_in "

Get sparsity of a given input.

";

%feature("docstring") casadi::Integrator::size_out "

Get output dimension.

";

%feature("docstring") casadi::Integrator::spInit "[INTERNAL]  Reset the
sparsity propagation.

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::Integrator::nnz_out "

Get of number of output nonzeros For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::Integrator::callReverse "

Create call to (cached) derivative function, reverse mode.

";

%feature("docstring") casadi::Integrator::numel_out "

Get of number of output elements For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::Integrator::repr "

Print a representation of the object.

";

%feature("docstring") casadi::Integrator::getAtomicInput "

Get the (integer) input arguments of an atomic operation.

";

%feature("docstring") casadi::Integrator::evaluate "

Evaluate.

";

%feature("docstring") casadi::Integrator::setOptionByEnumValue "[INTERNAL]
Set a certain option by giving an enum value.

";

%feature("docstring") casadi::Integrator::isNull "

Is a null pointer?

";

%feature("docstring") casadi::Integrator::getAlgorithmSize "

Get the number of atomic operations.

";

%feature("docstring") casadi::Integrator::generateLiftingFunctions "

Extract the functions needed for the Lifted Newton method.

";

%feature("docstring") casadi::Integrator::getAugmented "

Generate a augmented DAE system with nfwd forward sensitivities and nadj
adjoint sensitivities.

";

%feature("docstring") casadi::Integrator::getOutput "

>  DMatrix  IOInterface< Function  >.getOutput(int oind=0) const
------------------------------------------------------------------------

Get an output by index.

Parameters:
-----------

oind:  index within the range [0..n_out()-1]

>  DMatrix  IOInterface< Function  >.getOutput(str oname) const
------------------------------------------------------------------------

Get an output by name.

Parameters:
-----------

oname:  output name. Only allowed when an output scheme is set.

>  void IOInterface< Function  >.getOutput(T val, int oind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.getOutput(T val, str oname)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::Integrator::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

";

%feature("docstring") casadi::Integrator::getAtomicInputReal "

Get the floating point output argument of an atomic operation.

";

%feature("docstring") casadi::Integrator::free_mx "

Get all the free variables of the function.

";

%feature("docstring") casadi::Integrator::setJacobian "

Set the Jacobian function of output oind with respect to input iind NOTE:
Does not take ownership, only weak references to the Jacobians are kept
internally

";

%feature("docstring") casadi::Integrator::derReverse "

Get a function that calculates nadj adjoint derivatives.

Returns a function with n_in + n_out +nadj*n_out inputs and nadj*n_in
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nadj*n_out inputs correspond to adjoint seeds, one direction at a time The
nadj*n_in outputs correspond to adjoint sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

(n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::Integrator::getOptionEnumValue "[INTERNAL]
Get the enum value corresponding to th certain option.

";

%feature("docstring") casadi::Integrator::callForward "

Create call to (cached) derivative function, forward mode.

";

%feature("docstring") casadi::Integrator::getWorkSize "

Get the length of the work vector.

";

%feature("docstring") casadi::Integrator::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::Integrator::sz_arg "[INTERNAL]  Get required
length of arg field.

";

%feature("docstring") casadi::Integrator::mx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::Integrator::nnz_in "

Get of number of input nonzeros For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::Integrator::size2_out "

Get output dimension.

";

%feature("docstring") casadi::Integrator::mapaccum "

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
  



";

%feature("docstring") casadi::Integrator::gradient "

Generate a gradient function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the output must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::Integrator::sparsity_out "

Get sparsity of a given output.

";

%feature("docstring") casadi::Integrator::callDerivative "[INTERNAL]
Evaluate the function symbolically or numerically with directional
derivatives The first two arguments are the nondifferentiated inputs and
results of the evaluation, the next two arguments are a set of forward
directional seeds and the resulting forward directional derivatives, the
length of the vector being the number of forward directions. The next two
arguments are a set of adjoint directional seeds and the resulting adjoint
directional derivatives, the length of the vector being the number of
adjoint directions.

";

%feature("docstring") casadi::Integrator::setFullJacobian "

Set the Jacobian of all the input nonzeros with respect to all output
nonzeros NOTE: Does not take ownership, only weak references to the Jacobian
are kept internally

";

%feature("docstring") casadi::Integrator::name "

Name of the function.

";

%feature("docstring") casadi::Integrator::size2_in "

Get input dimension.

";

%feature("docstring") casadi::Integrator::getInput "

>  DMatrix  IOInterface< Function  >.getInput(int iind=0) const
------------------------------------------------------------------------

Get an input by index.

Parameters:
-----------

iind:  index within the range [0..n_in()-1]

>  DMatrix  IOInterface< Function  >.getInput(str iname) const
------------------------------------------------------------------------

Get an input by name.

Parameters:
-----------

iname:  input name. Only allowed when an input scheme is set.

>  void IOInterface< Function  >.getInput(T val, int iind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.getInput(T val, str iname)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::Integrator::index_in "

Find the index for a string describing a particular entry of an input
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::Integrator::sx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::Integrator::setInput "

>  void IOInterface< Function  >.setInput(T val, int iind=0)
------------------------------------------------------------------------

Set an input by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.setInput(T val, str iname)
------------------------------------------------------------------------

Set an input by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::Integrator::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::Integrator::getSanitizedName "

get function name with all non alphanumeric characters converted to '_'

";

%feature("docstring") casadi::Integrator::index_out "

Find the index for a string describing a particular entry of an output
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::Integrator::name_in "

>  [str] Function.name_in() const 
------------------------------------------------------------------------

Get input scheme.

>  str Function.name_in(int ind) const 
------------------------------------------------------------------------

Get input scheme name by index.

";

%feature("docstring") casadi::Integrator::numel_in "

Get of number of input elements For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::Integrator::n_in "

Get the number of function inputs.

";

%feature("docstring") casadi::Integrator::mx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::Integrator::countNodes "

Number of nodes in the algorithm.

";

%feature("docstring") casadi::Integrator::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::Integrator::derivative "

Get a function that calculates nfwd forward derivatives and nadj adjoint
derivatives Legacy function: Use derForward and derReverse instead.

Returns a function with (1+nfwd)*n_in+nadj*n_out inputs and (1+nfwd)*n_out +
nadj*n_in outputs. The first n_in inputs correspond to nondifferentiated
inputs. The next nfwd*n_in inputs correspond to forward seeds, one direction
at a time and the last nadj*n_out inputs correspond to adjoint seeds, one
direction at a time. The first n_out outputs correspond to nondifferentiated
outputs. The next nfwd*n_out outputs correspond to forward sensitivities,
one direction at a time and the last nadj*n_in outputs corresponds to
adjoint sensitivities, one direction at a time.

(n_in = n_in(), n_out = n_out())

";

%feature("docstring") casadi::Integrator::printStats "

Print solver statistics.

";

%feature("docstring") casadi::Integrator::default_in "

Get default input value (NOTE: constant reference)

";

%feature("docstring") casadi::Integrator::sz_res "[INTERNAL]  Get required
length of res field.

";

%feature("docstring") casadi::Integrator::size1_out "

Get output dimension.

";

%feature("docstring") casadi::Integrator::getOptionTypeName "

Get the type name of a certain option.

";

%feature("docstring") casadi::Integrator::spEvaluate "[INTERNAL]  Propagate
the sparsity pattern through a set of directional.

derivatives forward or backward (for usage, see the example
propagating_sparsity.cpp)

";

%feature("docstring") casadi::Integrator::description_out "

Get output scheme description by index.

";

%feature("docstring") casadi::Integrator::checkInputs "[INTERNAL]  Check if
the numerical values of the supplied bounds make sense.

";


// File: classcasadi_1_1InterruptHandler.xml
%feature("docstring") casadi::InterruptHandler "[INTERNAL]  Takes care of
user interrupts (Ctrl+C)

This is an internal class.

Joris Gillis

C++ includes: casadi_interrupt.hpp ";


// File: classcasadi_1_1Inverse.xml


// File: classcasadi_1_1IpoptUserClass.xml
%feature("docstring") casadi::IpoptUserClass::get_starting_point "[INTERNAL]  Method to return the starting point for the algorithm

";

%feature("docstring") casadi::IpoptUserClass::finalize_solution "[INTERNAL]
This method is called when the algorithm is complete so the TNLP can
store/write the solution

";

%feature("docstring") casadi::IpoptUserClass "[INTERNAL] C++ includes:
ipopt_nlp.hpp ";

%feature("docstring")
casadi::IpoptUserClass::get_list_of_nonlinear_variables "[INTERNAL]
Specify which variables that appear in the Hessian

";

%feature("docstring") casadi::IpoptUserClass::eval_grad_f "[INTERNAL]
Method to return the gradient of the objective

";

%feature("docstring") casadi::IpoptUserClass::get_var_con_metadata "[INTERNAL]  Allows setting information about variables and constraints

";

%feature("docstring") casadi::IpoptUserClass::~IpoptUserClass "[INTERNAL]
";

%feature("docstring") casadi::IpoptUserClass::eval_g "[INTERNAL]  Method to
return the constraint residuals

";

%feature("docstring") casadi::IpoptUserClass::get_nlp_info "[INTERNAL]
Method to return some info about the nlp

";

%feature("docstring") casadi::IpoptUserClass::eval_f "[INTERNAL]  Method to
return the objective value

";

%feature("docstring")
casadi::IpoptUserClass::get_number_of_nonlinear_variables "[INTERNAL]
Specify the number of variables that appear in the Hessian

";

%feature("docstring") casadi::IpoptUserClass::eval_jac_g "[INTERNAL]
Method to return: 1) The structure of the Jacobian (if \"values\" is NULL)
2) The values of the Jacobian (if \"values\" is not NULL)

";

%feature("docstring") casadi::IpoptUserClass::finalize_metadata "[INTERNAL]
Retrieve information about variables and constraints

";

%feature("docstring") casadi::IpoptUserClass::get_bounds_info "[INTERNAL]
Method to return the bounds for my problem

";

%feature("docstring") casadi::IpoptUserClass::eval_h "[INTERNAL]  Method to
return: 1) The structure of the hessian of the Lagrangian (if \"values\" is
NULL) 2) The values of the hessian of the Lagrangian (if \"values\" is not
NULL)

";

%feature("docstring") casadi::IpoptUserClass::intermediate_callback "[INTERNAL]  This method is called at every iteration

";

%feature("docstring") casadi::IpoptUserClass::IpoptUserClass "[INTERNAL] ";


// File: classcasadi_1_1KernelSum2D.xml


/*  Simple Getters & Setters  */

/*  Advanced Getters  */

/*  Option Functionality  */ %feature("docstring")
casadi::KernelSum2D::tangent "

Generate a tangent function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the input must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::KernelSum2D::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::KernelSum2D::name_in "

>  [str] Function.name_in() const 
------------------------------------------------------------------------

Get input scheme.

>  str Function.name_in(int ind) const 
------------------------------------------------------------------------

Get input scheme name by index.

";

%feature("docstring") casadi::KernelSum2D::nnz_in "

Get of number of input nonzeros For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::KernelSum2D::printDimensions "

Print dimensions of inputs and outputs.

";

%feature("docstring") casadi::KernelSum2D::generateLiftingFunctions "

Extract the functions needed for the Lifted Newton method.

";

%feature("docstring") casadi::KernelSum2D::sparsity_in "

Get sparsity of a given input.

";

%feature("docstring") casadi::KernelSum2D::index_out "

Find the index for a string describing a particular entry of an output
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::KernelSum2D::map "

>  [[MX] ] Function.map([[MX ] ] arg, str parallelization=\"serial\")

>  [MX] Function.map([MX ] arg, str parallelization=\"serial\")
------------------------------------------------------------------------

Evaluate symbolically in parallel (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

>  Function Function.map(str name, int N, Dict options=Dict()) const 
------------------------------------------------------------------------

Create a mapped version of this function.

Suppose the function has a signature of:

::

     f: (a, p) -> ( s )
  



The the mapaccumulated version has the signature:

::

     F: (A, P) -> (S )
  
      with
          a: horzcat([a0, a1, ..., a_(N-1)])
          p: horzcat([p0, p1, ..., p_(N-1)])
          s: horzcat([s0, s1, ..., s_(N-1)])
      and
          s0 <- f(a0, p0)
          s1 <- f(a1, p1)
          ...
          s_(N-1) <- f(a_(N-1), p_(N-1))
  



>  Function Function.map(str name, int n, [bool ] repeat_in, [bool ] repeat_out, Dict opts=Dict()) const 
------------------------------------------------------------------------

Generic map.

";

%feature("docstring") casadi::KernelSum2D::print "

Print a description of the object.

";

%feature("docstring") casadi::KernelSum2D::description_in "

Get input scheme description by index.

";

%feature("docstring") casadi::KernelSum2D::sz_res "[INTERNAL]  Get required
length of res field.

";

%feature("docstring") casadi::KernelSum2D::default_in "

Get default input value (NOTE: constant reference)

";

%feature("docstring") casadi::KernelSum2D::is_a "

Check if the function is of a particular type Optionally check if name
matches one of the base classes (default true)

";

%feature("docstring") casadi::KernelSum2D::spInit "[INTERNAL]  Reset the
sparsity propagation.

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::KernelSum2D::callForward "

Create call to (cached) derivative function, forward mode.

";

%feature("docstring") casadi::KernelSum2D::spCanEvaluate "[INTERNAL]  Is
the class able to propagate seeds through the algorithm?

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::KernelSum2D::countNodes "

Number of nodes in the algorithm.

";

%feature("docstring") casadi::KernelSum2D::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::KernelSum2D::getStats "

Get all statistics obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::KernelSum2D::setJacobian "

Set the Jacobian function of output oind with respect to input iind NOTE:
Does not take ownership, only weak references to the Jacobians are kept
internally

";

%feature("docstring") casadi::KernelSum2D::callDerivative "[INTERNAL]
Evaluate the function symbolically or numerically with directional
derivatives The first two arguments are the nondifferentiated inputs and
results of the evaluation, the next two arguments are a set of forward
directional seeds and the resulting forward directional derivatives, the
length of the vector being the number of forward directions. The next two
arguments are a set of adjoint directional seeds and the resulting adjoint
directional derivatives, the length of the vector being the number of
adjoint directions.

";

%feature("docstring") casadi::KernelSum2D::getAtomicInput "

Get the (integer) input arguments of an atomic operation.

";

%feature("docstring") casadi::KernelSum2D::getAtomicOperation "

Get an atomic operation operator index.

";

%feature("docstring") casadi::KernelSum2D::getInput "

>  DMatrix  IOInterface< Function  >.getInput(int iind=0) const
------------------------------------------------------------------------

Get an input by index.

Parameters:
-----------

iind:  index within the range [0..n_in()-1]

>  DMatrix  IOInterface< Function  >.getInput(str iname) const
------------------------------------------------------------------------

Get an input by name.

Parameters:
-----------

iname:  input name. Only allowed when an input scheme is set.

>  void IOInterface< Function  >.getInput(T val, int iind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.getInput(T val, str iname)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::KernelSum2D::size2_out "

Get output dimension.

";

%feature("docstring") casadi::KernelSum2D::n_in "

Get the number of function inputs.

";

%feature("docstring") casadi::KernelSum2D::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::KernelSum2D::KernelSum2D "

>  KernelSum2D()
------------------------------------------------------------------------

Default constructor.

>  KernelSum2D(str name, Function f, (int,int) size, double r, int n, Dict opts=Dict())
------------------------------------------------------------------------

Constructor (generic kernel_sum_2d)

";

%feature("docstring") casadi::KernelSum2D::setFullJacobian "

Set the Jacobian of all the input nonzeros with respect to all output
nonzeros NOTE: Does not take ownership, only weak references to the Jacobian
are kept internally

";

%feature("docstring") casadi::KernelSum2D::size1_out "

Get output dimension.

";

%feature("docstring") casadi::KernelSum2D::setOptionByAllowedIndex "[INTERNAL]  Set a certain option by giving its index into the allowed
values.

";

%feature("docstring") casadi::KernelSum2D::checkInputs "[INTERNAL]  Check
if the numerical values of the supplied bounds make sense.

";

%feature("docstring") casadi::KernelSum2D::jacobian "

Generate a Jacobian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. If
compact is set to true, only the nonzeros of the input and output
expressions are considered. If symmetric is set to true, the Jacobian being
calculated is known to be symmetric (usually a Hessian), which can be
exploited by the algorithm.

The generated Jacobian has one more output than the calling function
corresponding to the Jacobian and the same number of inputs.

";

%feature("docstring") casadi::KernelSum2D::call "

Evaluate the function symbolically or numerically.

";

%feature("docstring") casadi::KernelSum2D::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::KernelSum2D::mx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::KernelSum2D::name "

Name of the function.

";

%feature("docstring") casadi::KernelSum2D::mx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::KernelSum2D::getAtomicOutput "

Get the (integer) output argument of an atomic operation.

";

%feature("docstring") casadi::KernelSum2D::derReverse "

Get a function that calculates nadj adjoint derivatives.

Returns a function with n_in + n_out +nadj*n_out inputs and nadj*n_in
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nadj*n_out inputs correspond to adjoint seeds, one direction at a time The
nadj*n_in outputs correspond to adjoint sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

(n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::KernelSum2D::derivative "

Get a function that calculates nfwd forward derivatives and nadj adjoint
derivatives Legacy function: Use derForward and derReverse instead.

Returns a function with (1+nfwd)*n_in+nadj*n_out inputs and (1+nfwd)*n_out +
nadj*n_in outputs. The first n_in inputs correspond to nondifferentiated
inputs. The next nfwd*n_in inputs correspond to forward seeds, one direction
at a time and the last nadj*n_out inputs correspond to adjoint seeds, one
direction at a time. The first n_out outputs correspond to nondifferentiated
outputs. The next nfwd*n_out outputs correspond to forward sensitivities,
one direction at a time and the last nadj*n_in outputs corresponds to
adjoint sensitivities, one direction at a time.

(n_in = n_in(), n_out = n_out())

";

%feature("docstring") casadi::KernelSum2D::setDerReverse "

Set a function that calculates nadj adjoint derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::KernelSum2D "

KernelSum2D

Consider a dense matrix V.

KernelSum computes

F(V,X) = sum_i sum_j f ( [i;j], V(i,j), X)

with X: [x;y]

where the summation is taken for all entries (i,j) that are a distance r
away from X.

This function assumes that V is fixed: sensitivities with respect to it are
not computed.

This allows for improved speed of evaluation.

Having V fixed is a common use case: V may be a large bitmap (observation),
onto which a kernel is fitted.

Joris Gillis

>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| ad_weight    | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | derivative c |              |
|              |              |              | alculation.W |              |
|              |              |              | hen there is |              |
|              |              |              | an option of |              |
|              |              |              | either using |              |
|              |              |              | forward or   |              |
|              |              |              | reverse mode |              |
|              |              |              | directional  |              |
|              |              |              | derivatives, |              |
|              |              |              | the          |              |
|              |              |              | condition ad |              |
|              |              |              | _weight*nf<= |              |
|              |              |              | (1-ad_weight |              |
|              |              |              | )*na is used |              |
|              |              |              | where nf and |              |
|              |              |              | na are       |              |
|              |              |              | estimates of |              |
|              |              |              | the number   |              |
|              |              |              | of forward/r |              |
|              |              |              | everse mode  |              |
|              |              |              | directional  |              |
|              |              |              | derivatives  |              |
|              |              |              | needed. By   |              |
|              |              |              | default,     |              |
|              |              |              | ad_weight is |              |
|              |              |              | calculated a |              |
|              |              |              | utomatically |              |
|              |              |              | , but this   |              |
|              |              |              | can be       |              |
|              |              |              | overridden   |              |
|              |              |              | by setting   |              |
|              |              |              | this option. |              |
|              |              |              | In           |              |
|              |              |              | particular,  |              |
|              |              |              | 0 means      |              |
|              |              |              | forcing      |              |
|              |              |              | forward mode |              |
|              |              |              | and 1        |              |
|              |              |              | forcing      |              |
|              |              |              | reverse      |              |
|              |              |              | mode. Leave  |              |
|              |              |              | unset for    |              |
|              |              |              | (class       |              |
|              |              |              | specific)    |              |
|              |              |              | heuristics.  |              |
+--------------+--------------+--------------+--------------+--------------+
| ad_weight_sp | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | sparsity     |              |
|              |              |              | pattern      |              |
|              |              |              | calculation  |              |
|              |              |              | calculation. |              |
|              |              |              | Overrides    |              |
|              |              |              | default      |              |
|              |              |              | behavior.    |              |
|              |              |              | Set to 0 and |              |
|              |              |              | 1 to force   |              |
|              |              |              | forward and  |              |
|              |              |              | reverse mode |              |
|              |              |              | respectively |              |
|              |              |              | . Cf. option |              |
|              |              |              | \"ad_weight\". |              |
+--------------+--------------+--------------+--------------+--------------+
| compiler     | OT_STRING    | \"clang\"      | Just-in-time | casadi::Func |
|              |              |              | compiler     | tionInternal |
|              |              |              | plugin to be |              |
|              |              |              | used.        |              |
+--------------+--------------+--------------+--------------+--------------+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode     |
|              |              |              | according to |              |
|              |              |              | a given      |              |
|              |              |              | recipe (low- |              |
|              |              |              | level)       |              |
+--------------+--------------+--------------+--------------+--------------+
| gather_stats | OT_BOOLEAN   | false        | Flag to      | casadi::Func |
|              |              |              | indicate     | tionInternal |
|              |              |              | whether      |              |
|              |              |              | statistics   |              |
|              |              |              | must be      |              |
|              |              |              | gathered     |              |
+--------------+--------------+--------------+--------------+--------------+
| input_scheme | OT_STRINGVEC | GenericType( | Custom input | casadi::Func |
|              | TOR          | )            | scheme       | tionInternal |
+--------------+--------------+--------------+--------------+--------------+
| inputs_check | OT_BOOLEAN   | true         | Throw        | casadi::Func |
|              |              |              | exceptions   | tionInternal |
|              |              |              | when the     |              |
|              |              |              | numerical    |              |
|              |              |              | values of    |              |
|              |              |              | the inputs   |              |
|              |              |              | don't make   |              |
|              |              |              | sense        |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_penalty  | OT_REAL      | 2            | When         | casadi::Func |
|              |              |              | requested    | tionInternal |
|              |              |              | for a number |              |
|              |              |              | of forward/r |              |
|              |              |              | everse       |              |
|              |              |              | directions,  |              |
|              |              |              | it may be    |              |
|              |              |              | cheaper to   |              |
|              |              |              | compute      |              |
|              |              |              | first the    |              |
|              |              |              | full         |              |
|              |              |              | jacobian and |              |
|              |              |              | then         |              |
|              |              |              | multiply     |              |
|              |              |              | with seeds,  |              |
|              |              |              | rather than  |              |
|              |              |              | obtain the   |              |
|              |              |              | requested    |              |
|              |              |              | directions   |              |
|              |              |              | in a straigh |              |
|              |              |              | tforward     |              |
|              |              |              | manner.      |              |
|              |              |              | Casadi uses  |              |
|              |              |              | a heuristic  |              |
|              |              |              | to decide    |              |
|              |              |              | which is     |              |
|              |              |              | cheaper. A   |              |
|              |              |              | high value   |              |
|              |              |              | of 'jac_pena |              |
|              |              |              | lty' makes   |              |
|              |              |              | it less      |              |
|              |              |              | likely for   |              |
|              |              |              | the heurstic |              |
|              |              |              | to chose the |              |
|              |              |              | full         |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy.    |              |
|              |              |              | The special  |              |
|              |              |              | value -1     |              |
|              |              |              | indicates    |              |
|              |              |              | never to use |              |
|              |              |              | the full     |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy     |              |
+--------------+--------------+--------------+--------------+--------------+
| jit          | OT_BOOLEAN   | false        | Use just-in- | casadi::Func |
|              |              |              | time         | tionInternal |
|              |              |              | compiler to  |              |
|              |              |              | speed up the |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| jit_options  | OT_DICT      | GenericType( | Options to   | casadi::Func |
|              |              | )            | be passed to | tionInternal |
|              |              |              | the jit      |              |
|              |              |              | compiler.    |              |
+--------------+--------------+--------------+--------------+--------------+
| monitor      | OT_STRINGVEC | GenericType( | Monitors to  | casadi::Func |
|              | TOR          | )            | be activated | tionInternal |
|              |              |              | (inputs|outp |              |
|              |              |              | uts)         |              |
+--------------+--------------+--------------+--------------+--------------+
| output_schem | OT_STRINGVEC | GenericType( | Custom       | casadi::Func |
| e            | TOR          | )            | output       | tionInternal |
|              |              |              | scheme       |              |
+--------------+--------------+--------------+--------------+--------------+
| regularity_c | OT_BOOLEAN   | true         | Throw        | casadi::Func |
| heck         |              |              | exceptions   | tionInternal |
|              |              |              | when NaN or  |              |
|              |              |              | Inf appears  |              |
|              |              |              | during       |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| user_data    | OT_VOIDPTR   | GenericType( | A user-      | casadi::Func |
|              |              | )            | defined      | tionInternal |
|              |              |              | field that   |              |
|              |              |              | can be used  |              |
|              |              |              | to identify  |              |
|              |              |              | the function |              |
|              |              |              | or pass      |              |
|              |              |              | additional   |              |
|              |              |              | information  |              |
+--------------+--------------+--------------+--------------+--------------+
| verbose      | OT_BOOLEAN   | false        | Verbose      | casadi::Func |
|              |              |              | evaluation   | tionInternal |
|              |              |              | for          |              |
|              |              |              | debugging    |              |
+--------------+--------------+--------------+--------------+--------------+

Diagrams
--------



C++ includes: kernel_sum_2d.hpp ";

%feature("docstring") casadi::KernelSum2D::size_out "

Get output dimension.

";

%feature("docstring") casadi::KernelSum2D::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

";

%feature("docstring") casadi::KernelSum2D::numel_out "

Get of number of output elements For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::KernelSum2D::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::KernelSum2D::setJacSparsity "

Generate the sparsity of a Jacobian block

";

%feature("docstring") casadi::KernelSum2D::getAtomicInputReal "

Get the floating point output argument of an atomic operation.

";

%feature("docstring") casadi::KernelSum2D::getOutput "

>  DMatrix  IOInterface< Function  >.getOutput(int oind=0) const
------------------------------------------------------------------------

Get an output by index.

Parameters:
-----------

oind:  index within the range [0..n_out()-1]

>  DMatrix  IOInterface< Function  >.getOutput(str oname) const
------------------------------------------------------------------------

Get an output by name.

Parameters:
-----------

oname:  output name. Only allowed when an output scheme is set.

>  void IOInterface< Function  >.getOutput(T val, int oind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.getOutput(T val, str oname)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::KernelSum2D::index_in "

Find the index for a string describing a particular entry of an input
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::KernelSum2D::fullJacobian "

Generate a Jacobian function of all the inputs elements with respect to all
the output elements).

";

%feature("docstring") casadi::KernelSum2D::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::KernelSum2D::dictionary "

Get the dictionary.

";

%feature("docstring") casadi::KernelSum2D::getOptionTypeName "

Get the type name of a certain option.

";

%feature("docstring") casadi::KernelSum2D::type_name "

Get type name.

";

%feature("docstring") casadi::KernelSum2D::nnz_out "

Get of number of output nonzeros For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::KernelSum2D::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::KernelSum2D::gradient "

Generate a gradient function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the output must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::KernelSum2D::evaluate "

Evaluate.

";

%feature("docstring") casadi::KernelSum2D::free_mx "

Get all the free variables of the function.

";

%feature("docstring") casadi::KernelSum2D::callReverse "

Create call to (cached) derivative function, reverse mode.

";

%feature("docstring") casadi::KernelSum2D::name_out "

>  [str] Function.name_out() const 
------------------------------------------------------------------------

Get output scheme.

>  str Function.name_out(int ind) const 
------------------------------------------------------------------------

Get output scheme name by index.

";

%feature("docstring") casadi::KernelSum2D::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::KernelSum2D::sz_w "[INTERNAL]  Get required
length of w field.

";

%feature("docstring") casadi::KernelSum2D::setOutput "

>  void IOInterface< Function  >.setOutput(T val, int oind=0)
------------------------------------------------------------------------

Set an output by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.setOutput(T val, str oname)
------------------------------------------------------------------------

Set an output by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::KernelSum2D::numel_in "

Get of number of input elements For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::KernelSum2D::n_out "

Get the number of function outputs.

";

%feature("docstring") casadi::KernelSum2D::derForward "

Get a function that calculates nfwd forward derivatives.

Returns a function with n_in + n_out +nfwd*n_in inputs and nfwd*n_out
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nfwd*n_in inputs correspond to forward seeds, one direction at a time The
nfwd*n_out outputs correspond to forward sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::KernelSum2D::jacSparsity "

Get, if necessary generate, the sparsity of a Jacobian block

";

%feature("docstring") casadi::KernelSum2D::getAlgorithmSize "

Get the number of atomic operations.

";

%feature("docstring") casadi::KernelSum2D::mapaccum "

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
  



";

%feature("docstring") casadi::KernelSum2D::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::KernelSum2D::getSanitizedName "

get function name with all non alphanumeric characters converted to '_'

";

%feature("docstring") casadi::KernelSum2D::addMonitor "

Add modules to be monitored.

";

%feature("docstring") casadi::KernelSum2D::setDerForward "

Set a function that calculates nfwd forward derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::KernelSum2D::hessian "

Generate a Hessian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The generated Hessian has two more outputs than the calling function
corresponding to the Hessian and the gradients.

";

%feature("docstring") casadi::KernelSum2D::size_in "

Get input dimension.

";

%feature("docstring") casadi::KernelSum2D::isNull "

Is a null pointer?

";

%feature("docstring") casadi::KernelSum2D::repr "

Print a representation of the object.

";

%feature("docstring") casadi::KernelSum2D::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::KernelSum2D::setOptionByEnumValue "[INTERNAL]
Set a certain option by giving an enum value.

";

%feature("docstring") casadi::KernelSum2D::generate "

Export / Generate C code for the function.

";

%feature("docstring") casadi::KernelSum2D::getStat "

Get a single statistic obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::KernelSum2D::sx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::KernelSum2D::removeMonitor "

Remove modules to be monitored.

";

%feature("docstring") casadi::KernelSum2D::sz_arg "[INTERNAL]  Get required
length of arg field.

";

%feature("docstring") casadi::KernelSum2D::getOptionAllowedIndex "[INTERNAL]  Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::KernelSum2D::spEvaluate "[INTERNAL]
Propagate the sparsity pattern through a set of directional.

derivatives forward or backward (for usage, see the example
propagating_sparsity.cpp)

";

%feature("docstring") casadi::KernelSum2D::setInput "

>  void IOInterface< Function  >.setInput(T val, int iind=0)
------------------------------------------------------------------------

Set an input by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.setInput(T val, str iname)
------------------------------------------------------------------------

Set an input by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::KernelSum2D::getOptionEnumValue "[INTERNAL]
Get the enum value corresponding to th certain option.

";

%feature("docstring") casadi::KernelSum2D::free_sx "

Get all the free variables of the function.

";

%feature("docstring") casadi::KernelSum2D::sx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::KernelSum2D::size1_in "

Get input dimension.

";

%feature("docstring") casadi::KernelSum2D::description_out "

Get output scheme description by index.

";

%feature("docstring") casadi::KernelSum2D::sz_iw "[INTERNAL]  Get required
length of iw field.

";

%feature("docstring") casadi::KernelSum2D::size2_in "

Get input dimension.

";

%feature("docstring") casadi::KernelSum2D::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::KernelSum2D::sparsity_out "

Get sparsity of a given output.

";

%feature("docstring") casadi::KernelSum2D::getWorkSize "

Get the length of the work vector.

";


// File: classcasadi_1_1LapackLuDense.xml


// File: classcasadi_1_1LapackQrDense.xml


// File: classcasadi_1_1LibInfo.xml
%feature("docstring") casadi::LibInfo "[INTERNAL]  Structure with
information about the library.

C++ includes: external_function_internal.hpp ";


// File: classcasadi_1_1LibInfo_3_01Compiler_01_4.xml
%feature("docstring") casadi::LibInfo< Compiler >::LibInfo " [INTERNAL] ";

%feature("docstring") casadi::LibInfo< Compiler > " [INTERNAL]  Library that
has been just-in-time compiled.

C++ includes: external_function_internal.hpp ";

%feature("docstring") casadi::LibInfo< Compiler >::get " [INTERNAL] ";

%feature("docstring") casadi::LibInfo< Compiler >::clear " [INTERNAL] ";


// File: classcasadi_1_1LibInfo_3_01std_1_1string_01_4.xml
%feature("docstring") casadi::LibInfo< std::string > " [INTERNAL]  Library
given as a dynamically linked library.

C++ includes: external_function_internal.hpp ";

%feature("docstring") casadi::LibInfo< std::string >::clear " [INTERNAL] ";

%feature("docstring") casadi::LibInfo< std::string >::LibInfo " [INTERNAL]
";

%feature("docstring") casadi::LibInfo< std::string >::get " [INTERNAL] ";


// File: classcasadi_1_1LinearSolver.xml


/*  Simple Getters & Setters  */

/*  Advanced Getters  */

/*  Option Functionality  */ %feature("docstring")
casadi::LinearSolver::mx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::LinearSolver::getAlgorithmSize "

Get the number of atomic operations.

";

%feature("docstring") casadi::LinearSolver::getInput "

>  DMatrix  IOInterface< Function  >.getInput(int iind=0) const
------------------------------------------------------------------------

Get an input by index.

Parameters:
-----------

iind:  index within the range [0..n_in()-1]

>  DMatrix  IOInterface< Function  >.getInput(str iname) const
------------------------------------------------------------------------

Get an input by name.

Parameters:
-----------

iname:  input name. Only allowed when an input scheme is set.

>  void IOInterface< Function  >.getInput(T val, int iind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.getInput(T val, str iname)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::LinearSolver::getAtomicInputReal "

Get the floating point output argument of an atomic operation.

";

%feature("docstring") casadi::LinearSolver::hessian "

Generate a Hessian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The generated Hessian has two more outputs than the calling function
corresponding to the Hessian and the gradients.

";

%feature("docstring") casadi::LinearSolver::spEvaluate "[INTERNAL]
Propagate the sparsity pattern through a set of directional.

derivatives forward or backward (for usage, see the example
propagating_sparsity.cpp)

";

%feature("docstring") casadi::LinearSolver::default_in "

Get default input value (NOTE: constant reference)

";

%feature("docstring") casadi::LinearSolver::sparsity_out "

Get sparsity of a given output.

";

%feature("docstring") casadi::LinearSolver::jacobian "

Generate a Jacobian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. If
compact is set to true, only the nonzeros of the input and output
expressions are considered. If symmetric is set to true, the Jacobian being
calculated is known to be symmetric (usually a Hessian), which can be
exploited by the algorithm.

The generated Jacobian has one more output than the calling function
corresponding to the Jacobian and the same number of inputs.

";

%feature("docstring") casadi::LinearSolver::size1_in "

Get input dimension.

";

%feature("docstring") casadi::LinearSolver::numel_out "

Get of number of output elements For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::LinearSolver::setOutput "

>  void IOInterface< Function  >.setOutput(T val, int oind=0)
------------------------------------------------------------------------

Set an output by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.setOutput(T val, str oname)
------------------------------------------------------------------------

Set an output by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::LinearSolver::sz_iw "[INTERNAL]  Get required
length of iw field.

";

%feature("docstring") casadi::LinearSolver::callForward "

Create call to (cached) derivative function, forward mode.

";

%feature("docstring") casadi::LinearSolver::free_sx "

Get all the free variables of the function.

";

%feature("docstring") casadi::LinearSolver::sz_arg "[INTERNAL]  Get
required length of arg field.

";

%feature("docstring") casadi::LinearSolver::sz_w "[INTERNAL]  Get required
length of w field.

";

%feature("docstring") casadi::LinearSolver::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::LinearSolver::sz_res "[INTERNAL]  Get
required length of res field.

";

%feature("docstring") casadi::LinearSolver::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::LinearSolver::is_a "

Check if the function is of a particular type Optionally check if name
matches one of the base classes (default true)

";

%feature("docstring") casadi::LinearSolver::evaluate "

Evaluate.

";

%feature("docstring") casadi::LinearSolver::dictionary "

Get the dictionary.

";

%feature("docstring") casadi::LinearSolver::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::LinearSolver::setOptionByAllowedIndex "[INTERNAL]  Set a certain option by giving its index into the allowed
values.

";

%feature("docstring") casadi::LinearSolver::n_out "

Get the number of function outputs.

";

%feature("docstring") casadi::LinearSolver::free_mx "

Get all the free variables of the function.

";

%feature("docstring") casadi::LinearSolver::prepare "

Factorize the matrix.

";

%feature("docstring") casadi::LinearSolver::size1_out "

Get output dimension.

";

%feature("docstring") casadi::LinearSolver::fullJacobian "

Generate a Jacobian function of all the inputs elements with respect to all
the output elements).

";

%feature("docstring") casadi::LinearSolver::map "

>  [[MX] ] Function.map([[MX ] ] arg, str parallelization=\"serial\")

>  [MX] Function.map([MX ] arg, str parallelization=\"serial\")
------------------------------------------------------------------------

Evaluate symbolically in parallel (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

>  Function Function.map(str name, int N, Dict options=Dict()) const 
------------------------------------------------------------------------

Create a mapped version of this function.

Suppose the function has a signature of:

::

     f: (a, p) -> ( s )
  



The the mapaccumulated version has the signature:

::

     F: (A, P) -> (S )
  
      with
          a: horzcat([a0, a1, ..., a_(N-1)])
          p: horzcat([p0, p1, ..., p_(N-1)])
          s: horzcat([s0, s1, ..., s_(N-1)])
      and
          s0 <- f(a0, p0)
          s1 <- f(a1, p1)
          ...
          s_(N-1) <- f(a_(N-1), p_(N-1))
  



>  Function Function.map(str name, int n, [bool ] repeat_in, [bool ] repeat_out, Dict opts=Dict()) const 
------------------------------------------------------------------------

Generic map.

";

%feature("docstring") casadi::LinearSolver::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::LinearSolver::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

";

%feature("docstring") casadi::LinearSolver::print "

Print a description of the object.

";

%feature("docstring") casadi::LinearSolver::gradient "

Generate a gradient function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the output must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::LinearSolver::callReverse "

Create call to (cached) derivative function, reverse mode.

";

%feature("docstring") casadi::LinearSolver::getSanitizedName "

get function name with all non alphanumeric characters converted to '_'

";

%feature("docstring") casadi::LinearSolver::name_out "

>  [str] Function.name_out() const 
------------------------------------------------------------------------

Get output scheme.

>  str Function.name_out(int ind) const 
------------------------------------------------------------------------

Get output scheme name by index.

";

%feature("docstring") casadi::LinearSolver::derReverse "

Get a function that calculates nadj adjoint derivatives.

Returns a function with n_in + n_out +nadj*n_out inputs and nadj*n_in
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nadj*n_out inputs correspond to adjoint seeds, one direction at a time The
nadj*n_in outputs correspond to adjoint sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

(n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::LinearSolver::removeMonitor "

Remove modules to be monitored.

";

%feature("docstring") casadi::LinearSolver::setDerReverse "

Set a function that calculates nadj adjoint derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::LinearSolver::printDimensions "

Print dimensions of inputs and outputs.

";

%feature("docstring") casadi::LinearSolver::getAtomicOperation "

Get an atomic operation operator index.

";

%feature("docstring") casadi::LinearSolver::prepared "

Check if prepared.

";

%feature("docstring") casadi::LinearSolver::type_name "

Get type name.

";

%feature("docstring") casadi::LinearSolver::nnz_out "

Get of number of output nonzeros For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::LinearSolver::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::LinearSolver::index_out "

Find the index for a string describing a particular entry of an output
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::LinearSolver::tangent "

Generate a tangent function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the input must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::LinearSolver::countNodes "

Number of nodes in the algorithm.

";

%feature("docstring") casadi::LinearSolver::index_in "

Find the index for a string describing a particular entry of an input
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::LinearSolver::repr "

Print a representation of the object.

";

%feature("docstring") casadi::LinearSolver::description_in "

Get input scheme description by index.

";

%feature("docstring") casadi::LinearSolver::isNull "

Is a null pointer?

";

%feature("docstring") casadi::LinearSolver "

Base class for the linear solver classes.

Solves the linear system A*X = B or A^T*X = B for X with A square and non-
singular

If A is structurally singular, an error will be thrown during init. If A is
numerically singular, the prepare step will fail.

The usual procedure to use LinearSolver is: init()

set the first input (A)

prepare()

set the second input (b)

solve()

Repeat steps 4 and 5 to work with other b vectors.

The method evaluate() combines the prepare() and solve() step and is
therefore more expensive if A is invariant.

General information
===================



>Input scheme: casadi::LinsolInput (LINSOL_NUM_IN = 2) [linsolIn]

+-----------+-------+------------------------------------------------+
| Full name | Short |                  Description                   |
+===========+=======+================================================+
| LINSOL_A  | A     | The square matrix A: sparse, (n x n). .        |
+-----------+-------+------------------------------------------------+
| LINSOL_B  | B     | The right-hand-side matrix b: dense, (n x m) . |
+-----------+-------+------------------------------------------------+

>Output scheme: casadi::LinsolOutput (LINSOL_NUM_OUT = 1) [linsolOut]

+-----------+-------+----------------------------------------------+
| Full name | Short |                 Description                  |
+===========+=======+==============================================+
| LINSOL_X  | X     | Solution to the linear system of equations . |
+-----------+-------+----------------------------------------------+

>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| ad_weight    | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | derivative c |              |
|              |              |              | alculation.W |              |
|              |              |              | hen there is |              |
|              |              |              | an option of |              |
|              |              |              | either using |              |
|              |              |              | forward or   |              |
|              |              |              | reverse mode |              |
|              |              |              | directional  |              |
|              |              |              | derivatives, |              |
|              |              |              | the          |              |
|              |              |              | condition ad |              |
|              |              |              | _weight*nf<= |              |
|              |              |              | (1-ad_weight |              |
|              |              |              | )*na is used |              |
|              |              |              | where nf and |              |
|              |              |              | na are       |              |
|              |              |              | estimates of |              |
|              |              |              | the number   |              |
|              |              |              | of forward/r |              |
|              |              |              | everse mode  |              |
|              |              |              | directional  |              |
|              |              |              | derivatives  |              |
|              |              |              | needed. By   |              |
|              |              |              | default,     |              |
|              |              |              | ad_weight is |              |
|              |              |              | calculated a |              |
|              |              |              | utomatically |              |
|              |              |              | , but this   |              |
|              |              |              | can be       |              |
|              |              |              | overridden   |              |
|              |              |              | by setting   |              |
|              |              |              | this option. |              |
|              |              |              | In           |              |
|              |              |              | particular,  |              |
|              |              |              | 0 means      |              |
|              |              |              | forcing      |              |
|              |              |              | forward mode |              |
|              |              |              | and 1        |              |
|              |              |              | forcing      |              |
|              |              |              | reverse      |              |
|              |              |              | mode. Leave  |              |
|              |              |              | unset for    |              |
|              |              |              | (class       |              |
|              |              |              | specific)    |              |
|              |              |              | heuristics.  |              |
+--------------+--------------+--------------+--------------+--------------+
| ad_weight_sp | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | sparsity     |              |
|              |              |              | pattern      |              |
|              |              |              | calculation  |              |
|              |              |              | calculation. |              |
|              |              |              | Overrides    |              |
|              |              |              | default      |              |
|              |              |              | behavior.    |              |
|              |              |              | Set to 0 and |              |
|              |              |              | 1 to force   |              |
|              |              |              | forward and  |              |
|              |              |              | reverse mode |              |
|              |              |              | respectively |              |
|              |              |              | . Cf. option |              |
|              |              |              | \"ad_weight\". |              |
+--------------+--------------+--------------+--------------+--------------+
| compiler     | OT_STRING    | \"clang\"      | Just-in-time | casadi::Func |
|              |              |              | compiler     | tionInternal |
|              |              |              | plugin to be |              |
|              |              |              | used.        |              |
+--------------+--------------+--------------+--------------+--------------+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode     |
|              |              |              | according to |              |
|              |              |              | a given      |              |
|              |              |              | recipe (low- |              |
|              |              |              | level)       |              |
+--------------+--------------+--------------+--------------+--------------+
| gather_stats | OT_BOOLEAN   | false        | Flag to      | casadi::Func |
|              |              |              | indicate     | tionInternal |
|              |              |              | whether      |              |
|              |              |              | statistics   |              |
|              |              |              | must be      |              |
|              |              |              | gathered     |              |
+--------------+--------------+--------------+--------------+--------------+
| input_scheme | OT_STRINGVEC | GenericType( | Custom input | casadi::Func |
|              | TOR          | )            | scheme       | tionInternal |
+--------------+--------------+--------------+--------------+--------------+
| inputs_check | OT_BOOLEAN   | true         | Throw        | casadi::Func |
|              |              |              | exceptions   | tionInternal |
|              |              |              | when the     |              |
|              |              |              | numerical    |              |
|              |              |              | values of    |              |
|              |              |              | the inputs   |              |
|              |              |              | don't make   |              |
|              |              |              | sense        |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_penalty  | OT_REAL      | 2            | When         | casadi::Func |
|              |              |              | requested    | tionInternal |
|              |              |              | for a number |              |
|              |              |              | of forward/r |              |
|              |              |              | everse       |              |
|              |              |              | directions,  |              |
|              |              |              | it may be    |              |
|              |              |              | cheaper to   |              |
|              |              |              | compute      |              |
|              |              |              | first the    |              |
|              |              |              | full         |              |
|              |              |              | jacobian and |              |
|              |              |              | then         |              |
|              |              |              | multiply     |              |
|              |              |              | with seeds,  |              |
|              |              |              | rather than  |              |
|              |              |              | obtain the   |              |
|              |              |              | requested    |              |
|              |              |              | directions   |              |
|              |              |              | in a straigh |              |
|              |              |              | tforward     |              |
|              |              |              | manner.      |              |
|              |              |              | Casadi uses  |              |
|              |              |              | a heuristic  |              |
|              |              |              | to decide    |              |
|              |              |              | which is     |              |
|              |              |              | cheaper. A   |              |
|              |              |              | high value   |              |
|              |              |              | of 'jac_pena |              |
|              |              |              | lty' makes   |              |
|              |              |              | it less      |              |
|              |              |              | likely for   |              |
|              |              |              | the heurstic |              |
|              |              |              | to chose the |              |
|              |              |              | full         |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy.    |              |
|              |              |              | The special  |              |
|              |              |              | value -1     |              |
|              |              |              | indicates    |              |
|              |              |              | never to use |              |
|              |              |              | the full     |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy     |              |
+--------------+--------------+--------------+--------------+--------------+
| jit          | OT_BOOLEAN   | false        | Use just-in- | casadi::Func |
|              |              |              | time         | tionInternal |
|              |              |              | compiler to  |              |
|              |              |              | speed up the |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| jit_options  | OT_DICT      | GenericType( | Options to   | casadi::Func |
|              |              | )            | be passed to | tionInternal |
|              |              |              | the jit      |              |
|              |              |              | compiler.    |              |
+--------------+--------------+--------------+--------------+--------------+
| monitor      | OT_STRINGVEC | GenericType( | Monitors to  | casadi::Func |
|              | TOR          | )            | be activated | tionInternal |
|              |              |              | (inputs|outp |              |
|              |              |              | uts)         |              |
+--------------+--------------+--------------+--------------+--------------+
| output_schem | OT_STRINGVEC | GenericType( | Custom       | casadi::Func |
| e            | TOR          | )            | output       | tionInternal |
|              |              |              | scheme       |              |
+--------------+--------------+--------------+--------------+--------------+
| regularity_c | OT_BOOLEAN   | true         | Throw        | casadi::Func |
| heck         |              |              | exceptions   | tionInternal |
|              |              |              | when NaN or  |              |
|              |              |              | Inf appears  |              |
|              |              |              | during       |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| user_data    | OT_VOIDPTR   | GenericType( | A user-      | casadi::Func |
|              |              | )            | defined      | tionInternal |
|              |              |              | field that   |              |
|              |              |              | can be used  |              |
|              |              |              | to identify  |              |
|              |              |              | the function |              |
|              |              |              | or pass      |              |
|              |              |              | additional   |              |
|              |              |              | information  |              |
+--------------+--------------+--------------+--------------+--------------+
| verbose      | OT_BOOLEAN   | false        | Verbose      | casadi::Func |
|              |              |              | evaluation   | tionInternal |
|              |              |              | for          |              |
|              |              |              | debugging    |              |
+--------------+--------------+--------------+--------------+--------------+

List of plugins
===============



- csparsecholesky

- csparse

- lapacklu

- lapackqr

- symbolicqr

Note: some of the plugins in this list might not be available on your
system. Also, there might be extra plugins available to you that are not
listed here. You can obtain their documentation with
LinearSolver.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

csparsecholesky
---------------



LinearSolver with CSparseCholesky Interface

>List of available options

+----+------+---------+-------------+
| Id | Type | Default | Description |
+====+======+=========+=============+
+----+------+---------+-------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

csparse
-------



LinearSolver with CSparse Interface

>List of available options

+----+------+---------+-------------+
| Id | Type | Default | Description |
+====+======+=========+=============+
+----+------+---------+-------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

lapacklu
--------



This class solves the linear system A.x=b by making an LU factorization of
A: A = L.U, with L lower and U upper triangular

>List of available options

+-----------------------------+------------+---------+-------------+
|             Id              |    Type    | Default | Description |
+=============================+============+=========+=============+
| allow_equilibration_failure | OT_BOOLEAN | false   |             |
+-----------------------------+------------+---------+-------------+
| equilibration               | OT_BOOLEAN | true    |             |
+-----------------------------+------------+---------+-------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

lapackqr
--------



This class solves the linear system A.x=b by making an QR factorization of
A: A = Q.R, with Q orthogonal and R upper triangular

>List of available options

+----+------+---------+-------------+
| Id | Type | Default | Description |
+====+======+=========+=============+
+----+------+---------+-------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

symbolicqr
----------



LinearSolver based on QR factorization with sparsity pattern based
reordering without partial pivoting

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| codegen         | OT_BOOLEAN      | false           | C-code          |
|                 |                 |                 | generation      |
+-----------------+-----------------+-----------------+-----------------+
| compiler        | OT_STRING       | \"gcc -fPIC -O2\" | Compiler        |
|                 |                 |                 | command to be   |
|                 |                 |                 | used for        |
|                 |                 |                 | compiling       |
|                 |                 |                 | generated code  |
+-----------------+-----------------+-----------------+-----------------+

--------------------------------------------------------------------------------



Joel Andersson
Diagrams
--------



C++ includes: linear_solver.hpp ";

%feature("docstring") casadi::LinearSolver::mapaccum "

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
  



";

%feature("docstring") casadi::LinearSolver::nnz_in "

Get of number of input nonzeros For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::LinearSolver::spInit "[INTERNAL]  Reset the
sparsity propagation.

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::LinearSolver::description_out "

Get output scheme description by index.

";

%feature("docstring") casadi::LinearSolver::getOptionAllowedIndex "[INTERNAL]  Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::LinearSolver::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::LinearSolver::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::LinearSolver::callDerivative "[INTERNAL]
Evaluate the function symbolically or numerically with directional
derivatives The first two arguments are the nondifferentiated inputs and
results of the evaluation, the next two arguments are a set of forward
directional seeds and the resulting forward directional derivatives, the
length of the vector being the number of forward directions. The next two
arguments are a set of adjoint directional seeds and the resulting adjoint
directional derivatives, the length of the vector being the number of
adjoint directions.

";

%feature("docstring") casadi::LinearSolver::sx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::LinearSolver::size_in "

Get input dimension.

";

%feature("docstring") casadi::LinearSolver::derivative "

Get a function that calculates nfwd forward derivatives and nadj adjoint
derivatives Legacy function: Use derForward and derReverse instead.

Returns a function with (1+nfwd)*n_in+nadj*n_out inputs and (1+nfwd)*n_out +
nadj*n_in outputs. The first n_in inputs correspond to nondifferentiated
inputs. The next nfwd*n_in inputs correspond to forward seeds, one direction
at a time and the last nadj*n_out inputs correspond to adjoint seeds, one
direction at a time. The first n_out outputs correspond to nondifferentiated
outputs. The next nfwd*n_out outputs correspond to forward sensitivities,
one direction at a time and the last nadj*n_in outputs corresponds to
adjoint sensitivities, one direction at a time.

(n_in = n_in(), n_out = n_out())

";

%feature("docstring") casadi::LinearSolver::derForward "

Get a function that calculates nfwd forward derivatives.

Returns a function with n_in + n_out +nfwd*n_in inputs and nfwd*n_out
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nfwd*n_in inputs correspond to forward seeds, one direction at a time The
nfwd*n_out outputs correspond to forward sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::LinearSolver::getFactorizationSparsity "

Obtain a symbolic Cholesky factorization Only for Cholesky solvers.

";

%feature("docstring") casadi::LinearSolver::getOptionEnumValue "[INTERNAL]
Get the enum value corresponding to th certain option.

";

%feature("docstring") casadi::LinearSolver::LinearSolver "

>  LinearSolver(str name, str solver, Sparsity sp, Dict opts=Dict())

>  LinearSolver(str name, str solver, Sparsity sp, int nrhs, Dict opts=Dict())
------------------------------------------------------------------------

Create a linear solver given a sparsity pattern (new syntax, includes
initialization)

Parameters:
-----------

solver:

Name of a solver. It might be one of:

- csparsecholesky

- csparse

- lapacklu

- lapackqr

- symbolicqr

Note: some of the plugins in this list might not be available on your
system. Also, there might be extra plugins available to you that are not
listed here. You can obtain their documentation with
LinearSolver.doc(\"myextraplugin\")

>  LinearSolver()
------------------------------------------------------------------------
[INTERNAL] 
Default (empty) constructor

";

%feature("docstring") casadi::LinearSolver::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::LinearSolver::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::LinearSolver::getStats "

Get all statistics obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::LinearSolver::getAtomicOutput "

Get the (integer) output argument of an atomic operation.

";

%feature("docstring") casadi::LinearSolver::setOptionByEnumValue "[INTERNAL]  Set a certain option by giving an enum value.

";

%feature("docstring") casadi::LinearSolver::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::LinearSolver::sparsity_in "

Get sparsity of a given input.

";

%feature("docstring") casadi::LinearSolver::name "

Name of the function.

";

%feature("docstring") casadi::LinearSolver::setDerForward "

Set a function that calculates nfwd forward derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::LinearSolver::size_out "

Get output dimension.

";

%feature("docstring") casadi::LinearSolver::size2_in "

Get input dimension.

";

%feature("docstring") casadi::LinearSolver::spCanEvaluate "[INTERNAL]  Is
the class able to propagate seeds through the algorithm?

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::LinearSolver::jacSparsity "

Get, if necessary generate, the sparsity of a Jacobian block

";

%feature("docstring") casadi::LinearSolver::generate "

Export / Generate C code for the function.

";

%feature("docstring") casadi::LinearSolver::call "

Evaluate the function symbolically or numerically.

";

%feature("docstring") casadi::LinearSolver::solve "

>  void LinearSolver.solve(bool transpose=false)
------------------------------------------------------------------------

Solve the system of equations, internal vector.

>  MX LinearSolver.solve(MX A, MX B, bool transpose=false)
------------------------------------------------------------------------

Create a solve node.

";

%feature("docstring") casadi::LinearSolver::setJacSparsity "

Generate the sparsity of a Jacobian block

";

%feature("docstring") casadi::LinearSolver::n_in "

Get the number of function inputs.

";

%feature("docstring") casadi::LinearSolver::size2_out "

Get output dimension.

";

%feature("docstring") casadi::LinearSolver::addMonitor "

Add modules to be monitored.

";

%feature("docstring") casadi::LinearSolver::setJacobian "

Set the Jacobian function of output oind with respect to input iind NOTE:
Does not take ownership, only weak references to the Jacobians are kept
internally

";

%feature("docstring") casadi::LinearSolver::generateLiftingFunctions "

Extract the functions needed for the Lifted Newton method.

";

%feature("docstring") casadi::LinearSolver::setFullJacobian "

Set the Jacobian of all the input nonzeros with respect to all output
nonzeros NOTE: Does not take ownership, only weak references to the Jacobian
are kept internally

";

%feature("docstring") casadi::LinearSolver::getOutput "

>  DMatrix  IOInterface< Function  >.getOutput(int oind=0) const
------------------------------------------------------------------------

Get an output by index.

Parameters:
-----------

oind:  index within the range [0..n_out()-1]

>  DMatrix  IOInterface< Function  >.getOutput(str oname) const
------------------------------------------------------------------------

Get an output by name.

Parameters:
-----------

oname:  output name. Only allowed when an output scheme is set.

>  void IOInterface< Function  >.getOutput(T val, int oind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.getOutput(T val, str oname)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::LinearSolver::numel_in "

Get of number of input elements For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::LinearSolver::name_in "

>  [str] Function.name_in() const 
------------------------------------------------------------------------

Get input scheme.

>  str Function.name_in(int ind) const 
------------------------------------------------------------------------

Get input scheme name by index.

";

%feature("docstring") casadi::LinearSolver::checkInputs "[INTERNAL]  Check
if the numerical values of the supplied bounds make sense.

";

%feature("docstring") casadi::LinearSolver::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::LinearSolver::getStat "

Get a single statistic obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::LinearSolver::getWorkSize "

Get the length of the work vector.

";

%feature("docstring") casadi::LinearSolver::getOptionTypeName "

Get the type name of a certain option.

";

%feature("docstring") casadi::LinearSolver::setInput "

>  void IOInterface< Function  >.setInput(T val, int iind=0)
------------------------------------------------------------------------

Set an input by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.setInput(T val, str iname)
------------------------------------------------------------------------

Set an input by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::LinearSolver::sx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::LinearSolver::getFactorization "

Obtain a numeric Cholesky factorization Only for Cholesky solvers.

";

%feature("docstring") casadi::LinearSolver::mx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::LinearSolver::getAtomicInput "

Get the (integer) input arguments of an atomic operation.

";


// File: classcasadi_1_1Logger.xml
%feature("docstring") casadi::Logger "

Keeps track of logging output to screen and/or files. All printout from
CasADi routines should go through this files.

Joel Andersson

C++ includes: casadi_logger.hpp ";


// File: classcasadi_1_1MapBase.xml


// File: classcasadi_1_1MapReduce.xml


// File: classcasadi_1_1MapSerial.xml


// File: singletoncasadi_1_1Matrix.xml


/*  Construct symbolic primitives  */

/* The \"sym\" function is intended to work in a similar way as \"sym\" used
in the Symbolic Toolbox for Matlab but instead creating a CasADi symbolic
primitive.

*/ %feature("docstring") casadi::Matrix::getIntValue "

Get double value (only if integer constant)

";

%feature("docstring") casadi::Matrix::nnz_upper "

Get the number of non-zeros in the upper triangular half.

";

%feature("docstring") casadi::Matrix::repr "

Print a representation of the object.

";

%feature("docstring") casadi::Matrix::getNZ "

>  void array(DataType) .getNZ([double ] output_m) const 
------------------------------------------------------------------------

Get the elements numerically.

>  void array(DataType) .getNZ(array(DataType) output_m, bool ind1, Slice k) const

>  void array(DataType) .getNZ(array(DataType) output_m, bool ind1, IMatrix k) const 
------------------------------------------------------------------------

Get a set of nonzeros

";

%feature("docstring") friendwrap_expand "

Expand the expression as a weighted sum (with constant weights)

";

%feature("docstring") friendwrap_mtaylor "

>  array(DataType)  mtaylor(array(DataType) ex, array(DataType) x, array(DataType) a, int order=1)
------------------------------------------------------------------------

multivariate Taylor series expansion

Do Taylor expansions until the aggregated order of a term is equal to
'order'. The aggregated order of $x^n y^m$ equals $n+m$.

>  array(DataType)  mtaylor(array(DataType) ex, array(DataType) x, array(DataType) a, int order, [int ] order_contributions)
------------------------------------------------------------------------

multivariate Taylor series expansion

Do Taylor expansions until the aggregated order of a term is equal to
'order'. The aggregated order of $x^n y^m$ equals $n+m$.

The argument order_contributions can denote how match each variable
contributes to the aggregated order. If x=[x, y] and order_contributions=[1,
2], then the aggregated order of $x^n y^m$ equals $1n+2m$.

Example usage

$ \\\\sin(b+a)+\\\\cos(b+a)(x-a)+\\\\cos(b+a)(y-b) $ $ y+x-(x^3+3y x^2+3 y^2
x+y^3)/6 $ $ (-3 x^2 y-x^3)/6+y+x $

";

%feature("docstring") casadi::Matrix::grad "

Gradient expression.

";

%feature("docstring") casadi::Matrix::hasNZ "

Returns true if the matrix has a non-zero at location rr, cc.

";

%feature("docstring") casadi::Matrix::set "

>  void array(DataType) .set(double val)

>  void array(DataType) .set(const double *val, bool tr=false)

>  void array(DataType) .set([double ] val, bool tr=false)
------------------------------------------------------------------------

Get the elements numerically.

>  void array(DataType) .set(array(DataType) m, bool ind1, Slice rr)

>  void array(DataType) .set(array(DataType) m, bool ind1, IMatrix rr)

>  void array(DataType) .set(array(DataType) m, bool ind1, Sparsity sp)
------------------------------------------------------------------------

Set a submatrix, single argument

>  void array(DataType) .set(array(DataType) m, bool ind1, Slice rr, Slice cc)

>  void array(DataType) .set(array(DataType) m, bool ind1, Slice rr, IMatrix cc)

>  void array(DataType) .set(array(DataType) m, bool ind1, IMatrix rr, Slice cc)

>  void array(DataType) .set(array(DataType) m, bool ind1, IMatrix rr, IMatrix cc)
------------------------------------------------------------------------

Set a submatrix, two arguments

>  void array(DataType) .set(array(DataType) val)
------------------------------------------------------------------------

Set all the entries without changing sparsity pattern.

";

%feature("docstring") casadi::Matrix::hasNonStructuralZeros "

Check if the matrix has any zero entries which are not structural zeros.

";

%feature("docstring") casadi::Matrix::nnz "

Get the number of (structural) non-zero elements.

";

%feature("docstring") casadi::Matrix::remove "

Remove columns and rows Remove/delete rows and/or columns of a matrix.

";

%feature("docstring") casadi::Matrix::getDep "

Get expressions of the children of the expression Only defined if symbolic
scalar. Wraps SXElement SXElement::getDep(int ch=0) const.

";

%feature("docstring") casadi::Matrix::get "

>  void array(DataType) .get([double ] output_m) const 
------------------------------------------------------------------------

Get the elements numerically.

>  void array(DataType) .get(array(DataType) output_m, bool ind1, Slice rr) const

>  void array(DataType) .get(array(DataType) output_m, bool ind1, IMatrix rr) const

>  void array(DataType) .get(array(DataType) output_m, bool ind1, Sparsity sp) const 
------------------------------------------------------------------------

Get a submatrix, single argument

>  void array(DataType) .get(array(DataType) output_m, bool ind1, Slice rr, Slice cc) const

>  void array(DataType) .get(array(DataType) output_m, bool ind1, Slice rr, IMatrix cc) const

>  void array(DataType) .get(array(DataType) output_m, bool ind1, IMatrix rr, Slice cc) const

>  void array(DataType) .get(array(DataType) output_m, bool ind1, IMatrix rr, IMatrix cc) const 
------------------------------------------------------------------------

Get a submatrix, two arguments

";

%feature("docstring") friendwrap_triangle "

triangle function

\\\\[ \\\\begin {cases} \\\\Lambda(x) = 0 & |x| >= 1 \\\\\\\\ \\\\Lambda(x)
= 1-|x| & |x| < 1 \\\\end {cases} \\\\]

";

%feature("docstring") friendwrap_adj "

Matrix adjoint.

";

%feature("docstring") casadi::Matrix::triplet "";

%feature("docstring") casadi::Matrix::printDense "

Print dense matrix-stype.

";

%feature("docstring") casadi::Matrix::getSparsity "

Get an owning reference to the sparsity pattern.

";

%feature("docstring") casadi::Matrix::isIdentity "

check if the matrix is an identity matrix (note that false negative answers
are possible)

";

%feature("docstring") casadi::Matrix::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::Matrix::setNZ "

>  void array(DataType) .setNZ(double val)

>  void array(DataType) .setNZ(const double *val)

>  void array(DataType) .setNZ([double ] val)
------------------------------------------------------------------------

Set the elements numerically.

>  void array(DataType) .setNZ(array(DataType) m, bool ind1, Slice k)

>  void array(DataType) .setNZ(array(DataType) m, bool ind1, IMatrix k)
------------------------------------------------------------------------

Set a set of nonzeros

";

%feature("docstring") casadi::Matrix::sanity_check "

Check if the dimensions and colind, row vectors are compatible.

Parameters:
-----------

complete:  set to true to also check elementwise throws an error as possible
result

";

%feature("docstring") casadi::Matrix::getSym "

Get upper triangular elements.

";

%feature("docstring") casadi::Matrix::unary "[INTERNAL]  Create nodes by
their ID.

";

%feature("docstring") friendwrap_norm_inf_mul "

Inf-norm of a Matrix-Matrix product.

";

%feature("docstring") casadi::Matrix::numel "

>  int array(DataType) .numel() const
------------------------------------------------------------------------

Get the number of elements.

>  int array(DataType) .numel(int i) const
------------------------------------------------------------------------

Get the number of elements in slice (cf. MATLAB)

";

%feature("docstring") friendwrap_all "

Returns true only if every element in the matrix is true.

";

%feature("docstring") casadi::Matrix::find "

Get the location of all non-zero elements as they would appear in a Dense
matrix A : DenseMatrix 4 x 3 B : SparseMatrix 4 x 3 , 5 structural non-
zeros.

k = A.find() A[k] will contain the elements of A that are non-zero in B

";

%feature("docstring") casadi::Matrix::toSlice "

>  Slice array(DataType) .toSlice(bool ind1=false) const 
------------------------------------------------------------------------

Convert to Slice (only for IMatrix)

";

%feature("docstring") casadi::Matrix::nnz_diag "

Get get the number of non-zeros on the diagonal.

";

%feature("docstring") casadi::Matrix::sparsity "

Get the sparsity pattern.

";

%feature("docstring") casadi::Matrix::isCommutative "

Check whether a binary SX is commutative.

Only defined if symbolic scalar.

";

%feature("docstring") casadi::Matrix::isvector "

Check if the matrix is a row or column vector.

";

%feature("docstring") casadi::Matrix::isSlice "

>  bool array(DataType) .isSlice(bool ind1=false) const 
------------------------------------------------------------------------

Is the Matrix a Slice (only for IMatrix)

";

%feature("docstring") casadi::Matrix::hasDuplicates "[INTERNAL]  Detect
duplicate symbolic expressions If there are symbolic primitives appearing
more than once, the function will return true and the names of the duplicate
expressions will be printed to userOut<true, PL_WARN>(). Note: Will mark the
node using SXElement::setTemp. Make sure to call resetInput() after usage.

";

%feature("docstring") casadi::Matrix::isempty "

Check if the sparsity is empty, i.e. if one of the dimensions is zero (or
optionally both dimensions)

";

%feature("docstring") casadi::Matrix::dim "

Get string representation of dimensions. The representation is (nrow x ncol
= numel | size)

";

%feature("docstring") casadi::Matrix::printSparse "

Print sparse matrix style.

";

%feature("docstring") casadi::Matrix::T "

Transpose the matrix.

";

%feature("docstring") friendwrap_any "

Returns true only if any element in the matrix is true.

";

%feature("docstring") casadi::Matrix::isscalar "

Check if the matrix expression is scalar.

";

%feature("docstring") casadi::Matrix::clear "";

%feature("docstring") casadi::Matrix::nonzeros_int "

Get all nonzeros.

";

%feature("docstring") friendwrap_poly_roots "

Attempts to find the roots of a polynomial.

This will only work for polynomials up to order 3 It is assumed that the
roots are real.

";

%feature("docstring") casadi::Matrix::setScientific "

Set the 'precision, width & scientific' used in printing and serializing to
streams.

";

%feature("docstring") casadi::Matrix::isrow "

Check if the matrix is a row vector (i.e. size1()==1)

";

%feature("docstring") casadi::Matrix::hess "

Hessian expression

";

%feature("docstring") casadi::Matrix::getValue "

>  double array(DataType) .getValue() const 
------------------------------------------------------------------------

Get double value (only if constant)

>  double array(DataType) .getValue(int k) const 
------------------------------------------------------------------------

Get double value (particular nonzero)

";

%feature("docstring") casadi::Matrix::isOne "

check if the matrix is 1 (note that false negative answers are possible)

";

%feature("docstring") casadi::Matrix::nnz_lower "

Get the number of non-zeros in the lower triangular half.

";

%feature("docstring") casadi::Matrix::reserve "";

%feature("docstring") casadi::Matrix::printScalar "

Print scalar.

";

%feature("docstring") casadi::Matrix::tang "

Tangent expression.

";

%feature("docstring") casadi::Matrix::isSmooth "

Check if smooth.

";

%feature("docstring") casadi::Matrix::erase "

>  void array(DataType) .erase([int ] rr, [int ] cc, bool ind1=false)
------------------------------------------------------------------------

Erase a submatrix (leaving structural zeros in its place) Erase rows and/or
columns of a matrix.

>  void array(DataType) .erase([int ] rr, bool ind1=false)
------------------------------------------------------------------------

Erase a submatrix (leaving structural zeros in its place) Erase elements of
a matrix.

";

%feature("docstring") casadi::Matrix::inf "

create a matrix with all inf

";

%feature("docstring") casadi::Matrix::scalar_matrix "[INTERNAL]  Create
nodes by their ID.

";

%feature("docstring") friendwrap_ramp "

ramp function

\\\\[ \\\\begin {cases} R(x) = 0 & x <= 1 \\\\\\\\ R(x) = x & x > 1 \\\\\\\\
\\\\end {cases} \\\\]

Also called: slope function

";

%feature("docstring") casadi::Matrix::ones "

Create a dense matrix or a matrix with specified sparsity with all entries
one.

";

%feature("docstring") friendwrap_jacobianTimesVector "

Calculate the Jacobian and multiply by a vector from the right This is
equivalent to mul(jacobian(ex, arg), v) or mul(jacobian(ex, arg).T, v) for
transpose_jacobian set to false and true respectively. If contrast to these
expressions, it will use directional derivatives which is typically (but not
necessarily) more efficient if the complete Jacobian is not needed and v has
few rows.

";

%feature("docstring") casadi::Matrix::addSub "

Add a submatrix to an existing matrix (TODO: remove memory allocation)

";

%feature("docstring") casadi::Matrix::resize "";

%feature("docstring") friendwrap_qr "

QR factorization using the modified Gram-Schmidt algorithm More stable than
the classical Gram-Schmidt, but may break down if the rows of A are nearly
linearly dependent See J. Demmel: Applied Numerical Linear Algebra
(algorithm 3.1.). Note that in SWIG, Q and R are returned by value.

";

%feature("docstring") casadi::Matrix::row "

Get the sparsity pattern. See the Sparsity class for details.

";

%feature("docstring") casadi::Matrix::zeros "

Create a dense matrix or a matrix with specified sparsity with all entries
zero.

";

%feature("docstring") casadi::Matrix::isConstant "

Check if the matrix is constant (note that false negative answers are
possible)

";

%feature("docstring") casadi::Matrix::getNdeps "

Get the number of dependencies of a binary SXElement Only defined if
symbolic scalar.

";

%feature("docstring") casadi::Matrix::printSplit "

Get strings corresponding to the nonzeros and the interdependencies.

";

%feature("docstring") casadi::Matrix "

Sparse matrix class. SX and DMatrix are specializations.

General sparse matrix class that is designed with the idea that \"everything
is a matrix\", that is, also scalars and vectors. This philosophy makes it
easy to use and to interface in particularly with Python and Matlab/Octave.
Index starts with 0. Index vec happens as follows: (rr, cc) -> k =
rr+cc*size1() Vectors are column vectors.  The storage format is Compressed
Column Storage (CCS), similar to that used for sparse matrices in Matlab,
but unlike this format, we do allow for elements to be structurally non-zero
but numerically zero.  Matrix<DataType> is polymorphic with a
std::vector<DataType> that contain all non-identical-zero elements. The
sparsity can be accessed with Sparsity& sparsity() Joel Andersson

C++ includes: casadi_types.hpp ";

%feature("docstring") casadi::Matrix::setValue "

>  void array(DataType) .setValue(double m)
------------------------------------------------------------------------

Set double value (only if constant)

>  void array(DataType) .setValue(double m, int k)
------------------------------------------------------------------------

Set double value (particular nonzero)

";

%feature("docstring") casadi::Matrix::isMinusOne "

check if the matrix is -1 (note that false negative answers are possible)

";

%feature("docstring") casadi::Matrix::printVector "

Print vector-style.

";

%feature("docstring") casadi::Matrix::isRegular "

Checks if expression does not contain NaN or Inf.

";

%feature("docstring") casadi::Matrix::istril "

Check if the matrix is lower triangular.

";

%feature("docstring") friendwrap_eig_symbolic "

Attempts to find the eigenvalues of a symbolic matrix This will only work
for up to 3x3 matrices.

";

%feature("docstring") casadi::Matrix::jac "

Jacobian expression.

";

%feature("docstring") casadi::Matrix::istriu "

Check if the matrix is upper triangular.

";

%feature("docstring") friendwrap_sparsify "

Make a matrix sparse by removing numerical zeros.

";

%feature("docstring") casadi::Matrix::nonzeros "

Get all nonzeros.

";

%feature("docstring") friendwrap_cofactor "

Get the (i,j) cofactor matrix.

";

%feature("docstring") casadi::Matrix::colind "

Get the sparsity pattern. See the Sparsity class for details.

";

%feature("docstring") casadi::Matrix::issquare "

Check if the matrix expression is square.

";

%feature("docstring") casadi::Matrix::isLeaf "

Check if SX is a leaf of the SX graph.

Only defined if symbolic scalar.

";

%feature("docstring") casadi::Matrix::matrix_scalar "[INTERNAL]  Create
nodes by their ID.

";

%feature("docstring") casadi::Matrix::nan "

create a matrix with all nan

";

%feature("docstring") casadi::Matrix::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") friendwrap_poly_coeff "

extracts polynomial coefficients from an expression

Parameters:
-----------

ex:  Scalar expression that represents a polynomial

x:  Scalar symbol that the polynomial is build up with

";

%feature("docstring") friendwrap_getMinor "

Get the (i,j) minor matrix.

";

%feature("docstring") friendwrap_heaviside "

Heaviside function.

\\\\[ \\\\begin {cases} H(x) = 0 & x<0 \\\\\\\\ H(x) = 1/2 & x=0 \\\\\\\\
H(x) = 1 & x>0 \\\\\\\\ \\\\end {cases} \\\\]

";

%feature("docstring") casadi::Matrix::setPrecision "

Set the 'precision, width & scientific' used in printing and serializing to
streams.

";

%feature("docstring") casadi::Matrix::isValidInput "

Check if matrix can be used to define function inputs. Sparse matrices can
return true if all non-zero elements are symbolic.

";

%feature("docstring") casadi::Matrix::isZero "

check if the matrix is 0 (note that false negative answers are possible)

";

%feature("docstring") casadi::Matrix::setSym "

Set upper triangular elements.

";

%feature("docstring") casadi::Matrix::size "

>  (int,int) array(DataType) .size() const
------------------------------------------------------------------------

Get the shape.

>  int array(DataType) .size(int axis) const
------------------------------------------------------------------------

Get the size along a particular dimensions.

";

%feature("docstring") casadi::Matrix::isdense "

Check if the matrix expression is dense.

";

%feature("docstring") casadi::Matrix::__nonzero__ "

>  bool array(DataType) .__nonzero__() const 
------------------------------------------------------------------------

Returns the truth value of a Matrix.

>  bool SX.__nonzero__() const
------------------------------------------------------------------------
[INTERNAL] 
";

%feature("docstring") friendwrap_taylor "

univariate Taylor series expansion

Calculate the Taylor expansion of expression 'ex' up to order 'order' with
respect to variable 'x' around the point 'a'

$(x)=f(a)+f'(a)(x-a)+f''(a)\\\\frac
{(x-a)^2}{2!}+f'''(a)\\\\frac{(x-a)^3}{3!}+\\\\ldots$

Example usage:

::

>>   x



";

%feature("docstring") casadi::Matrix::isSymbolic "

Check if symbolic (Dense) Sparse matrices invariable return false.

";

%feature("docstring") casadi::Matrix::size1 "

Get the first dimension (i.e. number of rows)

";

%feature("docstring") friendwrap_chol "

Obtain a Cholesky factorisation of a matrix Returns an upper triangular R
such that R'R = A. Matrix A must be positive definite.

At the moment, the algorithm is dense (Cholesky-Banachiewicz). There is an
open ticket #1212 to make it sparse.

";

%feature("docstring") casadi::Matrix::getName "

Get name (only if symbolic scalar)

";

%feature("docstring") casadi::Matrix::printme "";

%feature("docstring") friendwrap_gauss_quadrature "

>  array(DataType)  gauss_quadrature(array(DataType) f, array(DataType) x, array(DataType) a, array(DataType) b, int order=5)
------------------------------------------------------------------------

Integrate f from a to b using Gaussian quadrature with n points.

>  array(DataType)  gauss_quadrature(array(DataType) f, array(DataType) x, array(DataType) a, array(DataType) b, int order, array(DataType) w)
------------------------------------------------------------------------

Matrix adjoint.

";

%feature("docstring") casadi::Matrix::get_row "

Get the sparsity pattern. See the Sparsity class for details.

";

%feature("docstring") friendwrap_pw_const "

Create a piecewise constant function Create a piecewise constant function
with n=val.size() intervals.

Inputs:

Parameters:
-----------

t:  a scalar variable (e.g. time)

tval:  vector with the discrete values of t at the interval transitions
(length n-1)

val:  vector with the value of the function for each interval (length n)

";

%feature("docstring") casadi::Matrix::enlarge "

Enlarge matrix Make the matrix larger by inserting empty rows and columns,
keeping the existing non-zeros.

";

%feature("docstring") casadi::Matrix::iscolumn "

Check if the matrix is a column vector (i.e. size2()==1)

";

%feature("docstring") casadi::Matrix::size2 "

Get the second dimension (i.e. number of columns)

";

%feature("docstring") friendwrap_pw_lin "

t a scalar variable (e.g. time)

Create a piecewise linear function Create a piecewise linear function:

Inputs: tval vector with the the discrete values of t (monotonically
increasing) val vector with the corresponding function values (same length
as tval)

";

%feature("docstring") casadi::Matrix::setWidth "

Set the 'precision, width & scientific' used in printing and serializing to
streams.

";

%feature("docstring") casadi::Matrix::Matrix "

>  array(DataType) ()
------------------------------------------------------------------------

constructors

empty 0-by-0 matrix constructor

>  array(DataType) (array(DataType) m)
------------------------------------------------------------------------

Copy constructor.

>  array(DataType) (int nrow, int ncol)
------------------------------------------------------------------------

Create a sparse matrix with all structural zeros.

>  array(DataType) (Sparsity sp)
------------------------------------------------------------------------

Create a sparse matrix from a sparsity pattern. Same as
Matrix::ones(sparsity)

>  array(DataType) (Sparsity sp, array(DataType) d)
------------------------------------------------------------------------

Construct matrix with a given sparsity and nonzeros.

>  array(DataType) (double val)
------------------------------------------------------------------------

This constructor enables implicit type conversion from a numeric type.

>  array(DataType) ([[double ] ] m)
------------------------------------------------------------------------

Dense matrix constructor with data given as vector of vectors.

>  array(DataType) (array(A) x)
------------------------------------------------------------------------

Create a matrix from another matrix with a different entry type Assumes that
the scalar conversion is valid.

>  array(DataType) ([A ] x)
------------------------------------------------------------------------

Create an expression from a vector.

>  array(DataType) ([DataType ] x)

>  array(DataType) ((int,int) rc)

>  array(DataType) (Sparsity sp, DataType val, bool dummy)

>  array(DataType) (Sparsity sp, [DataType ] d, bool dummy)
------------------------------------------------------------------------
[INTERNAL] 
";

%feature("docstring") casadi::Matrix::getElementHash "

Returns a number that is unique for a given symbolic scalar.

Only defined if symbolic scalar.

";

%feature("docstring") friendwrap_rectangle "

rectangle function

\\\\[ \\\\begin {cases} \\\\Pi(x) = 1 & |x| < 1/2 \\\\\\\\ \\\\Pi(x) = 1/2 &
|x| = 1/2 \\\\\\\\ \\\\Pi(x) = 0 & |x| > 1/2 \\\\\\\\ \\\\end {cases} \\\\]

Also called: gate function, block function, band function, pulse function,
window function

";

%feature("docstring") casadi::Matrix::sym "

>  static array(DataType)   array(DataType) .sym(str name, int nrow=1, int ncol=1)
------------------------------------------------------------------------

Create an nrow-by-ncol symbolic primitive.

>  static array(DataType)   array(DataType) .sym(str name, (int,int) rc)
------------------------------------------------------------------------

Construct a symbolic primitive with given dimensions.

>  static array(DataType)   array(DataType) .sym(str name, Sparsity sp)
------------------------------------------------------------------------

Create symbolic primitive with a given sparsity pattern.

>  static[array(DataType)   ] array(DataType) .sym(str name, Sparsity sp, int p)
------------------------------------------------------------------------

Create a vector of length p with with matrices with symbolic primitives of
given sparsity.

>  static[array(DataType)   ] array(DataType) .sym(str name, int nrow, int ncol, int p)
------------------------------------------------------------------------

Create a vector of length p with nrow-by-ncol symbolic primitives.

>  static[[array(DataType)  ] ] array(DataType) .sym(str name, Sparsity sp, int p, int r)
------------------------------------------------------------------------

Create a vector of length r of vectors of length p with symbolic primitives
with given sparsity.

>  static[[array(DataType)  ] ] array(DataType) .sym(str name, int nrow, int ncol, int p, int r)
------------------------------------------------------------------------

Create a vector of length r of vectors of length p with nrow-by-ncol
symbolic primitives.

";

%feature("docstring") casadi::Matrix::matrix_matrix "[INTERNAL]  Create
nodes by their ID.

";

%feature("docstring") casadi::Matrix::binary "[INTERNAL]  Create nodes by
their ID.

";

%feature("docstring") casadi::Matrix::isInteger "

Check if the matrix is integer-valued (note that false negative answers are
possible)

";

%feature("docstring") casadi::Matrix::get_colind "

Get the sparsity pattern. See the Sparsity class for details.

";

%feature("docstring") casadi::Matrix::print "

Print a description of the object.

";

%feature("docstring") casadi::Matrix::resetInput "[INTERNAL]  Reset the
marker for an input expression.

";


// File: classcasadi_1_1MinusInfSX.xml


// File: classcasadi_1_1MinusOneSX.xml


// File: classcasadi_1_1Monitor.xml


// File: classcasadi_1_1MultipleOutput.xml


// File: classcasadi_1_1Multiplication.xml


// File: classcasadi_1_1MX.xml


/*  Construct symbolic primitives  */

/* The \"sym\" function is intended to work in a similar way as \"sym\" used
in the Symbolic Toolbox for Matlab but instead creating a CasADi symbolic
primitive.

*/ %feature("docstring") casadi::MX::attachAssert "

returns itself, but with an assertion attached

If y does not evaluate to 1, a runtime error is raised

";

%feature("docstring") casadi::MX "

MX - Matrix expression.

The MX class is used to build up trees made up from MXNodes. It is a more
general graph representation than the scalar expression, SX, and much less
efficient for small objects. On the other hand, the class allows much more
general operations than does SX, in particular matrix valued operations and
calls to arbitrary differentiable functions.

The MX class is designed to have identical syntax with the Matrix<> template
class, and uses Matrix<double> as its internal representation of the values
at a node. By keeping the syntaxes identical, it is possible to switch from
one class to the other, as well as inlining MX functions to SXElement
functions.

Note that an operation is always \"lazy\", making a matrix multiplication
will create a matrix multiplication node, not perform the actual
multiplication.

Joel Andersson

C++ includes: mx.hpp ";

%feature("docstring") casadi::MX::grad "

Gradient expression.

";

%feature("docstring") casadi::MX::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::MX::erase "

>  void MX.erase([int ] rr, [int ] cc, bool ind1=false)
------------------------------------------------------------------------

Erase a submatrix (leaving structural zeros in its place) Erase rows and/or
columns of a matrix.

>  void MX.erase([int ] rr, bool ind1=false)
------------------------------------------------------------------------

Erase a submatrix (leaving structural zeros in its place) Erase elements of
a matrix.

";

%feature("docstring") casadi::MX::monitor "

Monitor an expression Returns itself, but with the side effect of printing
the nonzeros along with a comment.

";

%feature("docstring") casadi::MX::isValidInput "

Check if matrix can be used to define function inputs. Valid inputs for
MXFunctions are combinations of Reshape, concatenations and SymbolicMX.

";

%feature("docstring") casadi::MX::set "

>  void MX.set(MX m, bool ind1, Slice rr)

>  void MX.set(MX m, bool ind1, IMatrix rr)

>  void MX.set(MX m, bool ind1, Sparsity sp)
------------------------------------------------------------------------

Set a submatrix, single argument

";

%feature("docstring") casadi::MX::getTemp "[INTERNAL]  Get the temporary
variable

";

%feature("docstring") casadi::MX::hasDuplicates "[INTERNAL]  Detect
duplicate symbolic expressions If there are symbolic primitives appearing
more than once, the function will return true and the names of the duplicate
expressions will be printed to userOut<true, PL_WARN>(). Note: Will mark the
node using MX::setTemp. Make sure to call resetInput() after usage.

";

%feature("docstring") casadi::MX::binary "

Create nodes by their ID.

";

%feature("docstring") casadi::MX::sparsity "

Get the sparsity pattern.

";

%feature("docstring") friendwrap_lift "

Lift the expression Experimental feature.

";

%feature("docstring") casadi::MX::size2 "

Get the second dimension (i.e. number of columns)

";

%feature("docstring") casadi::MX::size1 "

Get the first dimension (i.e. number of rows)

";

%feature("docstring") casadi::MX::iscolumn "

Check if the matrix is a column vector (i.e. size2()==1)

";

%feature("docstring") casadi::MX::isMultiplication "

Check if multiplication.

";

%feature("docstring") casadi::MX::~MX "[INTERNAL]  Destructor.

";

%feature("docstring") casadi::MX::isscalar "

Check if the matrix expression is scalar.

";

%feature("docstring") casadi::MX::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::MX::n_out "

Number of outputs.

";

%feature("docstring") casadi::MX::isrow "

Check if the matrix is a row vector (i.e. size1()==1)

";

%feature("docstring") casadi::MX::isEvaluationOutput "

Check if evaluation output.

";

%feature("docstring") casadi::MX::print "

Print a description of the object.

";

%feature("docstring") casadi::MX::splitPrimitives "

Split up an expression along symbolic primitives.

";

%feature("docstring") casadi::MX::ones "

Create a dense matrix or a matrix with specified sparsity with all entries
one.

";

%feature("docstring") casadi::MX::getSparsity "

Get an owning reference to the sparsity pattern.

";

%feature("docstring") casadi::MX::isIdentity "

check if identity

";

%feature("docstring") casadi::MX::isdense "

Check if the matrix expression is dense.

";

%feature("docstring") casadi::MX::isMinusOne "

check if zero (note that false negative answers are possible)

";

%feature("docstring") casadi::MX::inf "

create a matrix with all inf

";

%feature("docstring") friendwrap_matrix_expand "

Expand MX graph to SXFunction call.

Expand the given expression e, optionally supplying expressions contained in
it at which expansion should stop.

";

%feature("docstring") casadi::MX::colind "

Get the sparsity pattern. See the Sparsity class for details.

";

%feature("docstring") casadi::MX::isRegular "

Checks if expression does not contain NaN or Inf.

";

%feature("docstring") casadi::MX::getName "

Get the name.

";

%feature("docstring") casadi::MX::get_colind "

Get the sparsity pattern. See the Sparsity class for details.

";

%feature("docstring") casadi::MX::getDep "

Get the nth dependency as MX.

";

%feature("docstring") casadi::MX::isvector "

Check if the matrix is a row or column vector.

";

%feature("docstring") casadi::MX::printPtr "[INTERNAL]  Print the pointer
to the internal class

";

%feature("docstring") casadi::MX::nnz "

Get the number of (structural) non-zero elements.

";

%feature("docstring") casadi::MX::find "

Get the location of all non-zero elements as they would appear in a Dense
matrix A : DenseMatrix 4 x 3 B : SparseMatrix 4 x 3 , 5 structural non-
zeros.

k = A.find() A[k] will contain the elements of A that are non-zero in B

";

%feature("docstring") casadi::MX::dim "

Get string representation of dimensions. The representation is (nrow x ncol
= numel | size)

";

%feature("docstring") friendwrap_find "

Find first nonzero If failed, returns the number of rows.

";

%feature("docstring") casadi::MX::isZero "

check if zero (note that false negative answers are possible)

";

%feature("docstring") friendwrap_graph_substitute "

>  MX graph_substitute(MX ex, [MX ] v, [MX ] vdef)
------------------------------------------------------------------------

Substitute single expression in graph Substitute variable v with expression
vdef in an expression ex, preserving nodes.

>  [MX] graph_substitute([MX ] ex, [MX ] v, [MX ] vdef)
------------------------------------------------------------------------

Substitute multiple expressions in graph Substitute variable var with
expression expr in multiple expressions, preserving nodes.

";

%feature("docstring") casadi::MX::getValue "

Get the value (only for scalar constant nodes)

";

%feature("docstring") casadi::MX::getEvaluationOutput "

Get the index of evaluation output - only valid when isEvaluationoutput() is
true.

";

%feature("docstring") casadi::MX::get_row "

Get the sparsity pattern. See the Sparsity class for details.

";

%feature("docstring") casadi::MX::nnz_upper "

Get the number of non-zeros in the upper triangular half.

";

%feature("docstring") casadi::MX::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::MX::jac "

Jacobian expression.

";

%feature("docstring") casadi::MX::getOutput "

Get an output.

";

%feature("docstring") casadi::MX::zeros "

Create a dense matrix or a matrix with specified sparsity with all entries
zero.

";

%feature("docstring") casadi::MX::size "

>  (int,int) MX .size() const
------------------------------------------------------------------------

Get the shape.

>  int MX .size(int axis) const
------------------------------------------------------------------------

Get the size along a particular dimensions.

";

%feature("docstring") casadi::MX::T "

Transpose the matrix.

";

%feature("docstring") casadi::MX::printme "";

%feature("docstring") casadi::MX::joinPrimitives "

Join an expression along symbolic primitives.

";

%feature("docstring") casadi::MX::nan "

create a matrix with all nan

";

%feature("docstring") casadi::MX::isTranspose "

Is the expression a transpose?

";

%feature("docstring") casadi::MX::row "

Get the sparsity pattern. See the Sparsity class for details.

";

%feature("docstring") casadi::MX::isCommutative "

Check if commutative operation.

";

%feature("docstring") casadi::MX::numel "

>  int MX .numel() const
------------------------------------------------------------------------

Get the number of elements.

>  int MX .numel(int i) const
------------------------------------------------------------------------

Get the number of elements in slice (cf. MATLAB)

";

%feature("docstring") casadi::MX::enlarge "

Enlarge matrix Make the matrix larger by inserting empty rows and columns,
keeping the existing non-zeros.

";

%feature("docstring") casadi::MX::sym "

>  static MX  MX .sym(str name, int nrow=1, int ncol=1)
------------------------------------------------------------------------

Create an nrow-by-ncol symbolic primitive.

>  static MX  MX .sym(str name, (int,int) rc)
------------------------------------------------------------------------

Construct a symbolic primitive with given dimensions.

>  static MX  MX .sym(str name, Sparsity sp)
------------------------------------------------------------------------

Create symbolic primitive with a given sparsity pattern.

>  static[MX  ] MX .sym(str name, Sparsity sp, int p)
------------------------------------------------------------------------

Create a vector of length p with with matrices with symbolic primitives of
given sparsity.

>  static[MX  ] MX .sym(str name, int nrow, int ncol, int p)
------------------------------------------------------------------------

Create a vector of length p with nrow-by-ncol symbolic primitives.

>  static[[MX ] ] MX .sym(str name, Sparsity sp, int p, int r)
------------------------------------------------------------------------

Create a vector of length r of vectors of length p with symbolic primitives
with given sparsity.

>  static[[MX ] ] MX .sym(str name, int nrow, int ncol, int p, int r)
------------------------------------------------------------------------

Create a vector of length r of vectors of length p with nrow-by-ncol
symbolic primitives.

";

%feature("docstring") casadi::MX::setNZ "

Set a set of nonzeros

";

%feature("docstring") casadi::MX::isOne "

check if zero (note that false negative answers are possible)

";

%feature("docstring") casadi::MX::get "

>  void MX.get(MX &output_m, bool ind1, Slice rr) const

>  void MX.get(MX &output_m, bool ind1, IMatrix rr) const

>  void MX.get(MX &output_m, bool ind1, Sparsity sp) const 
------------------------------------------------------------------------

Get a submatrix, single argument

>  void MX.get(MX &output_m, bool ind1, Slice rr, Slice cc) const

>  void MX.get(MX &output_m, bool ind1, Slice rr, IMatrix cc) const

>  void MX.get(MX &output_m, bool ind1, IMatrix rr, Slice cc) const

>  void MX.get(MX &output_m, bool ind1, IMatrix rr, IMatrix cc) const 
------------------------------------------------------------------------

Get a submatrix, two arguments

";

%feature("docstring") casadi::MX::repr "

Print a representation of the object.

";

%feature("docstring") casadi::MX::getFunction "

Get function.

";

%feature("docstring") casadi::MX::numFunctions "

Number of functions.

";

%feature("docstring") casadi::MX::__nonzero__ "

Returns the truth value of an MX expression.

";

%feature("docstring") casadi::MX::resetInput "[INTERNAL]  Reset the marker
for an input expression.

";

%feature("docstring") casadi::MX::tang "

Tangent expression.

";

%feature("docstring") casadi::MX::istril "

Check if the matrix is lower triangular.

";

%feature("docstring") casadi::MX::isempty "

Check if the sparsity is empty, i.e. if one of the dimensions is zero (or
optionally both dimensions)

";

%feature("docstring") casadi::MX::zz_project "

Set sparse.

";

%feature("docstring") casadi::MX::istriu "

Check if the matrix is upper triangular.

";

%feature("docstring") casadi::MX::MX "

>  MX()
------------------------------------------------------------------------

Default constructor.

>  MX(int nrow, int ncol)
------------------------------------------------------------------------

Create a sparse matrix with all structural zeros.

>  MX(Sparsity sp)
------------------------------------------------------------------------

Create a sparse matrix from a sparsity pattern. Same as MX::ones(sparsity)

>  MX(Sparsity sp, MX val)
------------------------------------------------------------------------

Construct matrix with a given sparsity and nonzeros.

>  MX(double x)
------------------------------------------------------------------------

Create scalar constant (also implicit type conversion)

>  MX(MX x)
------------------------------------------------------------------------

Copy constructor.

>  MX([double ] x)
------------------------------------------------------------------------

Create vector constant (also implicit type conversion)

>  MX(DMatrix x)
------------------------------------------------------------------------

Create sparse matrix constant (also implicit type conversion)

";

%feature("docstring") casadi::MX::numPrimitives "

Get the number of symbolic primitive Assumes isValidInput() returns true.

";

%feature("docstring") casadi::MX::isBinary "

Is binary operation.

";

%feature("docstring") casadi::MX::isUnary "

Is unary operation.

";

%feature("docstring") casadi::MX::getOp "

Get operation type.

";

%feature("docstring") casadi::MX::nnz_lower "

Get the number of non-zeros in the lower triangular half.

";

%feature("docstring") casadi::MX::isNorm "

Check if norm.

";

%feature("docstring") casadi::MX::setTemp "[INTERNAL]  Set the temporary
variable.

";

%feature("docstring") casadi::MX::getPrimitives "

Get symbolic primitives.

";

%feature("docstring") casadi::MX::getMatrixValue "

Get the value (only for constant nodes)

";

%feature("docstring") casadi::MX::isConstant "

Check if constant.

";

%feature("docstring") casadi::MX::mapping "

Get an IMatrix representation of a GetNonzeros or SetNonzeros node.

";

%feature("docstring") casadi::MX::getNZ "

Get a set of nonzeros

";

%feature("docstring") casadi::MX::isEvaluation "

Check if evaluation.

";

%feature("docstring") casadi::MX::isOperation "

Is it a certain operation.

";

%feature("docstring") casadi::MX::getNdeps "

Get the number of dependencies of a binary SXElement.

";

%feature("docstring") casadi::MX::isSymbolic "

Check if symbolic.

";

%feature("docstring") casadi::MX::issquare "

Check if the matrix expression is square.

";

%feature("docstring") casadi::MX::isNull "

Is a null pointer?

";

%feature("docstring") casadi::MX::nnz_diag "

Get get the number of non-zeros on the diagonal.

";

%feature("docstring") casadi::MX::unary "

Create nodes by their ID.

";


// File: classcasadi_1_1NanSX.xml


// File: classcasadi_1_1Newton.xml


// File: classcasadi_1_1NlpBuilder.xml


/*  Symbolic representation of the NLP  */

/* Data members

*/ %feature("docstring") casadi::NlpBuilder "

A symbolic NLP representation.

Joel Andersson

C++ includes: nlp_builder.hpp ";

%feature("docstring") casadi::NlpBuilder::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::NlpBuilder::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::NlpBuilder::parseNL "

Parse an AMPL och PyOmo NL-file.

";

%feature("docstring") casadi::NlpBuilder::repr "

Print a representation of the object.

";

%feature("docstring") casadi::NlpBuilder::print "

Print a description of the object.

";


// File: classcasadi_1_1NlpSolver.xml


/*  Simple Getters & Setters  */

/*  Advanced Getters  */

/*  Option Functionality  */ %feature("docstring")
casadi::NlpSolver::getAtomicInput "

Get the (integer) input arguments of an atomic operation.

";

%feature("docstring") casadi::NlpSolver::size_out "

Get output dimension.

";

%feature("docstring") casadi::NlpSolver::mx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::NlpSolver::getInput "

>  DMatrix  IOInterface< Function  >.getInput(int iind=0) const
------------------------------------------------------------------------

Get an input by index.

Parameters:
-----------

iind:  index within the range [0..n_in()-1]

>  DMatrix  IOInterface< Function  >.getInput(str iname) const
------------------------------------------------------------------------

Get an input by name.

Parameters:
-----------

iname:  input name. Only allowed when an input scheme is set.

>  void IOInterface< Function  >.getInput(T val, int iind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.getInput(T val, str iname)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::NlpSolver::sparsity_out "

Get sparsity of a given output.

";

%feature("docstring") casadi::NlpSolver::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::NlpSolver::callForward "

Create call to (cached) derivative function, forward mode.

";

%feature("docstring") casadi::NlpSolver::getOutput "

>  DMatrix  IOInterface< Function  >.getOutput(int oind=0) const
------------------------------------------------------------------------

Get an output by index.

Parameters:
-----------

oind:  index within the range [0..n_out()-1]

>  DMatrix  IOInterface< Function  >.getOutput(str oname) const
------------------------------------------------------------------------

Get an output by name.

Parameters:
-----------

oname:  output name. Only allowed when an output scheme is set.

>  void IOInterface< Function  >.getOutput(T val, int oind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.getOutput(T val, str oname)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::NlpSolver::setDerReverse "

Set a function that calculates nadj adjoint derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::NlpSolver::setInput "

>  void IOInterface< Function  >.setInput(T val, int iind=0)
------------------------------------------------------------------------

Set an input by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.setInput(T val, str iname)
------------------------------------------------------------------------

Set an input by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::NlpSolver::tangent "

Generate a tangent function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the input must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::NlpSolver::NlpSolver "

>  NlpSolver()
------------------------------------------------------------------------

Default constructor.

>  NlpSolver(str name, str solver, Function nlp, Dict opts=Dict())
------------------------------------------------------------------------

NLP solver factory (new syntax, includes initialization)

>  NlpSolver(str name, str solver, const str:SX &nlp, Dict opts=Dict())
------------------------------------------------------------------------

Create an NLP solver from a dictionary with SX expressions.

>  NlpSolver(str name, str solver, const str:MX &nlp, Dict opts=Dict())
------------------------------------------------------------------------

Create an NLP solver from a dictionary with MX expressions.

";

%feature("docstring") casadi::NlpSolver::numel_in "

Get of number of input elements For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::NlpSolver::getAtomicOperation "

Get an atomic operation operator index.

";

%feature("docstring") casadi::NlpSolver::getAtomicInputReal "

Get the floating point output argument of an atomic operation.

";

%feature("docstring") casadi::NlpSolver::removeMonitor "

Remove modules to be monitored.

";

%feature("docstring") casadi::NlpSolver::size2_out "

Get output dimension.

";

%feature("docstring") casadi::NlpSolver::jacG "

Access the Hessian of the Lagrangian function.

>Input scheme: casadi::JacGInput (JACG_NUM_IN = 2) [jacGIn]

+-----------+-------+---------------------+
| Full name | Short |     Description     |
+===========+=======+=====================+
| JACG_X    | x     | Decision variable . |
+-----------+-------+---------------------+
| JACG_P    | p     | Fixed parameter .   |
+-----------+-------+---------------------+

>Output scheme: casadi::JacGOutput (JACG_NUM_OUT = 3) [jacGOut]

+-----------+-------+-------------------------------+
| Full name | Short |          Description          |
+===========+=======+===============================+
| JACG_JAC  | jac   | Jacobian of the constraints . |
+-----------+-------+-------------------------------+
| JACG_F    | f     | Objective function .          |
+-----------+-------+-------------------------------+
| JACG_G    | g     | Constraint function .         |
+-----------+-------+-------------------------------+

";

%feature("docstring") casadi::NlpSolver::setOutput "

>  void IOInterface< Function  >.setOutput(T val, int oind=0)
------------------------------------------------------------------------

Set an output by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.setOutput(T val, str oname)
------------------------------------------------------------------------

Set an output by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::NlpSolver::size1_in "

Get input dimension.

";

%feature("docstring") casadi::NlpSolver::setDerForward "

Set a function that calculates nfwd forward derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::NlpSolver::addMonitor "

Add modules to be monitored.

";

%feature("docstring") casadi::NlpSolver::getOptionTypeName "

Get the type name of a certain option.

";

%feature("docstring") casadi::NlpSolver::nlp "

Access the NLP.

>Input scheme: casadi::NlpSolverInput (NLP_SOLVER_NUM_IN = 8) [nlpSolverIn]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| NLP_SOLVER_X0          | x0                     | Decision variables,    |
|                        |                        | initial guess (nx x 1) |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_P           | p                      | Value of fixed         |
|                        |                        | parameters (np x 1) .  |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LBX         | lbx                    | Decision variables     |
|                        |                        | lower bound (nx x 1),  |
|                        |                        | default -inf .         |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_UBX         | ubx                    | Decision variables     |
|                        |                        | upper bound (nx x 1),  |
|                        |                        | default +inf .         |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LBG         | lbg                    | Constraints lower      |
|                        |                        | bound (ng x 1),        |
|                        |                        | default -inf .         |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_UBG         | ubg                    | Constraints upper      |
|                        |                        | bound (ng x 1),        |
|                        |                        | default +inf .         |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_X0      | lam_x0                 | Lagrange multipliers   |
|                        |                        | for bounds on X,       |
|                        |                        | initial guess (nx x 1) |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_G0      | lam_g0                 | Lagrange multipliers   |
|                        |                        | for bounds on G,       |
|                        |                        | initial guess (ng x 1) |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+

>Output scheme: casadi::NlpSolverOutput (NLP_SOLVER_NUM_OUT = 6) [nlpSolverOut]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| NLP_SOLVER_X           | x                      | Decision variables at  |
|                        |                        | the optimal solution   |
|                        |                        | (nx x 1) .             |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_F           | f                      | Cost function value at |
|                        |                        | the optimal solution   |
|                        |                        | (1 x 1) .              |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_G           | g                      | Constraints function   |
|                        |                        | at the optimal         |
|                        |                        | solution (ng x 1) .    |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_X       | lam_x                  | Lagrange multipliers   |
|                        |                        | for bounds on X at the |
|                        |                        | solution (nx x 1) .    |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_G       | lam_g                  | Lagrange multipliers   |
|                        |                        | for bounds on G at the |
|                        |                        | solution (ng x 1) .    |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_P       | lam_p                  | Lagrange multipliers   |
|                        |                        | for bounds on P at the |
|                        |                        | solution (np x 1) .    |
+------------------------+------------------------+------------------------+

";

%feature("docstring") casadi::NlpSolver::default_in "

Get default input value (NOTE: constant reference)

";

%feature("docstring") casadi::NlpSolver::printDimensions "

Print dimensions of inputs and outputs.

";

%feature("docstring") casadi::NlpSolver::getSanitizedName "

get function name with all non alphanumeric characters converted to '_'

";

%feature("docstring") casadi::NlpSolver::sz_iw "[INTERNAL]  Get required
length of iw field.

";

%feature("docstring") casadi::NlpSolver::map "

>  [[MX] ] Function.map([[MX ] ] arg, str parallelization=\"serial\")

>  [MX] Function.map([MX ] arg, str parallelization=\"serial\")
------------------------------------------------------------------------

Evaluate symbolically in parallel (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

>  Function Function.map(str name, int N, Dict options=Dict()) const 
------------------------------------------------------------------------

Create a mapped version of this function.

Suppose the function has a signature of:

::

     f: (a, p) -> ( s )
  



The the mapaccumulated version has the signature:

::

     F: (A, P) -> (S )
  
      with
          a: horzcat([a0, a1, ..., a_(N-1)])
          p: horzcat([p0, p1, ..., p_(N-1)])
          s: horzcat([s0, s1, ..., s_(N-1)])
      and
          s0 <- f(a0, p0)
          s1 <- f(a1, p1)
          ...
          s_(N-1) <- f(a_(N-1), p_(N-1))
  



>  Function Function.map(str name, int n, [bool ] repeat_in, [bool ] repeat_out, Dict opts=Dict()) const 
------------------------------------------------------------------------

Generic map.

";

%feature("docstring") casadi::NlpSolver::getReducedHessian "

Get the reduced Hessian. Requires a patched sIPOPT installation, see CasADi
documentation.

";

%feature("docstring") casadi::NlpSolver::setFullJacobian "

Set the Jacobian of all the input nonzeros with respect to all output
nonzeros NOTE: Does not take ownership, only weak references to the Jacobian
are kept internally

";

%feature("docstring") casadi::NlpSolver::derForward "

Get a function that calculates nfwd forward derivatives.

Returns a function with n_in + n_out +nfwd*n_in inputs and nfwd*n_out
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nfwd*n_in inputs correspond to forward seeds, one direction at a time The
nfwd*n_out outputs correspond to forward sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::NlpSolver::description_out "

Get output scheme description by index.

";

%feature("docstring") casadi::NlpSolver::hessian "

Generate a Hessian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The generated Hessian has two more outputs than the calling function
corresponding to the Hessian and the gradients.

";

%feature("docstring") casadi::NlpSolver::type_name "

Get type name.

";

%feature("docstring") casadi::NlpSolver::mx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::NlpSolver::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

";

%feature("docstring") casadi::NlpSolver::countNodes "

Number of nodes in the algorithm.

";

%feature("docstring") casadi::NlpSolver::numel_out "

Get of number of output elements For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::NlpSolver::generate "

Export / Generate C code for the function.

";

%feature("docstring") casadi::NlpSolver::spCanEvaluate "[INTERNAL]  Is the
class able to propagate seeds through the algorithm?

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::NlpSolver::setJacobian "

Set the Jacobian function of output oind with respect to input iind NOTE:
Does not take ownership, only weak references to the Jacobians are kept
internally

";

%feature("docstring") casadi::NlpSolver::derivative "

Get a function that calculates nfwd forward derivatives and nadj adjoint
derivatives Legacy function: Use derForward and derReverse instead.

Returns a function with (1+nfwd)*n_in+nadj*n_out inputs and (1+nfwd)*n_out +
nadj*n_in outputs. The first n_in inputs correspond to nondifferentiated
inputs. The next nfwd*n_in inputs correspond to forward seeds, one direction
at a time and the last nadj*n_out inputs correspond to adjoint seeds, one
direction at a time. The first n_out outputs correspond to nondifferentiated
outputs. The next nfwd*n_out outputs correspond to forward sensitivities,
one direction at a time and the last nadj*n_in outputs corresponds to
adjoint sensitivities, one direction at a time.

(n_in = n_in(), n_out = n_out())

";

%feature("docstring") casadi::NlpSolver::sz_res "[INTERNAL]  Get required
length of res field.

";

%feature("docstring") casadi::NlpSolver::gradient "

Generate a gradient function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the output must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::NlpSolver::getAtomicOutput "

Get the (integer) output argument of an atomic operation.

";

%feature("docstring") casadi::NlpSolver::generateLiftingFunctions "

Extract the functions needed for the Lifted Newton method.

";

%feature("docstring") casadi::NlpSolver::reportConstraints "

Prints out a human readable report about possible constraint violations,
after solving.

";

%feature("docstring") casadi::NlpSolver::getOptionAllowedIndex "[INTERNAL]
Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::NlpSolver::jacobian "

Generate a Jacobian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. If
compact is set to true, only the nonzeros of the input and output
expressions are considered. If symmetric is set to true, the Jacobian being
calculated is known to be symmetric (usually a Hessian), which can be
exploited by the algorithm.

The generated Jacobian has one more output than the calling function
corresponding to the Jacobian and the same number of inputs.

";

%feature("docstring") casadi::NlpSolver::name_out "

>  [str] Function.name_out() const 
------------------------------------------------------------------------

Get output scheme.

>  str Function.name_out(int ind) const 
------------------------------------------------------------------------

Get output scheme name by index.

";

%feature("docstring") casadi::NlpSolver::jacSparsity "

Get, if necessary generate, the sparsity of a Jacobian block

";

%feature("docstring") casadi::NlpSolver::setOptionByAllowedIndex "[INTERNAL]  Set a certain option by giving its index into the allowed
values.

";

%feature("docstring") casadi::NlpSolver::free_mx "

Get all the free variables of the function.

";

%feature("docstring") casadi::NlpSolver::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::NlpSolver::getWorkSize "

Get the length of the work vector.

";

%feature("docstring") casadi::NlpSolver::sx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::NlpSolver::sz_w "[INTERNAL]  Get required
length of w field.

";

%feature("docstring") casadi::NlpSolver::fullJacobian "

Generate a Jacobian function of all the inputs elements with respect to all
the output elements).

";

%feature("docstring") casadi::NlpSolver::getStats "

Get all statistics obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::NlpSolver::derReverse "

Get a function that calculates nadj adjoint derivatives.

Returns a function with n_in + n_out +nadj*n_out inputs and nadj*n_in
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nadj*n_out inputs correspond to adjoint seeds, one direction at a time The
nadj*n_in outputs correspond to adjoint sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

(n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::NlpSolver::callReverse "

Create call to (cached) derivative function, reverse mode.

";

%feature("docstring") casadi::NlpSolver::print "

Print a description of the object.

";

%feature("docstring") casadi::NlpSolver::getAlgorithmSize "

Get the number of atomic operations.

";

%feature("docstring") casadi::NlpSolver::description_in "

Get input scheme description by index.

";

%feature("docstring") casadi::NlpSolver::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::NlpSolver "

NlpSolver.

Solves the following parametric nonlinear program (NLP):

::

  min          F(x, p)
   x
  
  subject to
              LBX <=   x    <= UBX
              LBG <= G(x, p) <= UBG
                         p  == P
  
      nx: number of decision variables
      ng: number of constraints
      np: number of parameters
  



General information
===================



>Input scheme: casadi::NlpSolverInput (NLP_SOLVER_NUM_IN = 8) [nlpSolverIn]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| NLP_SOLVER_X0          | x0                     | Decision variables,    |
|                        |                        | initial guess (nx x 1) |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_P           | p                      | Value of fixed         |
|                        |                        | parameters (np x 1) .  |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LBX         | lbx                    | Decision variables     |
|                        |                        | lower bound (nx x 1),  |
|                        |                        | default -inf .         |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_UBX         | ubx                    | Decision variables     |
|                        |                        | upper bound (nx x 1),  |
|                        |                        | default +inf .         |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LBG         | lbg                    | Constraints lower      |
|                        |                        | bound (ng x 1),        |
|                        |                        | default -inf .         |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_UBG         | ubg                    | Constraints upper      |
|                        |                        | bound (ng x 1),        |
|                        |                        | default +inf .         |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_X0      | lam_x0                 | Lagrange multipliers   |
|                        |                        | for bounds on X,       |
|                        |                        | initial guess (nx x 1) |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_G0      | lam_g0                 | Lagrange multipliers   |
|                        |                        | for bounds on G,       |
|                        |                        | initial guess (ng x 1) |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+

>Output scheme: casadi::NlpSolverOutput (NLP_SOLVER_NUM_OUT = 6) [nlpSolverOut]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| NLP_SOLVER_X           | x                      | Decision variables at  |
|                        |                        | the optimal solution   |
|                        |                        | (nx x 1) .             |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_F           | f                      | Cost function value at |
|                        |                        | the optimal solution   |
|                        |                        | (1 x 1) .              |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_G           | g                      | Constraints function   |
|                        |                        | at the optimal         |
|                        |                        | solution (ng x 1) .    |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_X       | lam_x                  | Lagrange multipliers   |
|                        |                        | for bounds on X at the |
|                        |                        | solution (nx x 1) .    |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_G       | lam_g                  | Lagrange multipliers   |
|                        |                        | for bounds on G at the |
|                        |                        | solution (ng x 1) .    |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_P       | lam_p                  | Lagrange multipliers   |
|                        |                        | for bounds on P at the |
|                        |                        | solution (np x 1) .    |
+------------------------+------------------------+------------------------+

>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| ad_weight    | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | derivative c |              |
|              |              |              | alculation.W |              |
|              |              |              | hen there is |              |
|              |              |              | an option of |              |
|              |              |              | either using |              |
|              |              |              | forward or   |              |
|              |              |              | reverse mode |              |
|              |              |              | directional  |              |
|              |              |              | derivatives, |              |
|              |              |              | the          |              |
|              |              |              | condition ad |              |
|              |              |              | _weight*nf<= |              |
|              |              |              | (1-ad_weight |              |
|              |              |              | )*na is used |              |
|              |              |              | where nf and |              |
|              |              |              | na are       |              |
|              |              |              | estimates of |              |
|              |              |              | the number   |              |
|              |              |              | of forward/r |              |
|              |              |              | everse mode  |              |
|              |              |              | directional  |              |
|              |              |              | derivatives  |              |
|              |              |              | needed. By   |              |
|              |              |              | default,     |              |
|              |              |              | ad_weight is |              |
|              |              |              | calculated a |              |
|              |              |              | utomatically |              |
|              |              |              | , but this   |              |
|              |              |              | can be       |              |
|              |              |              | overridden   |              |
|              |              |              | by setting   |              |
|              |              |              | this option. |              |
|              |              |              | In           |              |
|              |              |              | particular,  |              |
|              |              |              | 0 means      |              |
|              |              |              | forcing      |              |
|              |              |              | forward mode |              |
|              |              |              | and 1        |              |
|              |              |              | forcing      |              |
|              |              |              | reverse      |              |
|              |              |              | mode. Leave  |              |
|              |              |              | unset for    |              |
|              |              |              | (class       |              |
|              |              |              | specific)    |              |
|              |              |              | heuristics.  |              |
+--------------+--------------+--------------+--------------+--------------+
| ad_weight_sp | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | sparsity     |              |
|              |              |              | pattern      |              |
|              |              |              | calculation  |              |
|              |              |              | calculation. |              |
|              |              |              | Overrides    |              |
|              |              |              | default      |              |
|              |              |              | behavior.    |              |
|              |              |              | Set to 0 and |              |
|              |              |              | 1 to force   |              |
|              |              |              | forward and  |              |
|              |              |              | reverse mode |              |
|              |              |              | respectively |              |
|              |              |              | . Cf. option |              |
|              |              |              | \"ad_weight\". |              |
+--------------+--------------+--------------+--------------+--------------+
| compiler     | OT_STRING    | \"clang\"      | Just-in-time | casadi::Func |
|              |              |              | compiler     | tionInternal |
|              |              |              | plugin to be |              |
|              |              |              | used.        |              |
+--------------+--------------+--------------+--------------+--------------+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode   c |
|              |              |              | according to | asadi::NlpSo |
|              |              |              | a given      | lverInternal |
|              |              |              | recipe (low- |              |
|              |              |              | level)  (qp) |              |
+--------------+--------------+--------------+--------------+--------------+
| eval_errors_ | OT_BOOLEAN   | false        | When errors  | casadi::NlpS |
| fatal        |              |              | occur during | olverInterna |
|              |              |              | evaluation   | l            |
|              |              |              | of           |              |
|              |              |              | f,g,...,stop |              |
|              |              |              | the          |              |
|              |              |              | iterations   |              |
+--------------+--------------+--------------+--------------+--------------+
| expand       | OT_BOOLEAN   | false        | Expand the   | casadi::NlpS |
|              |              |              | NLP function | olverInterna |
|              |              |              | in terms of  | l            |
|              |              |              | scalar       |              |
|              |              |              | operations,  |              |
|              |              |              | i.e. MX->SX  |              |
+--------------+--------------+--------------+--------------+--------------+
| gather_stats | OT_BOOLEAN   | false        | Flag to      | casadi::Func |
|              |              |              | indicate     | tionInternal |
|              |              |              | whether      |              |
|              |              |              | statistics   |              |
|              |              |              | must be      |              |
|              |              |              | gathered     |              |
+--------------+--------------+--------------+--------------+--------------+
| grad_f       | OT_FUNCTION  | GenericType( | Function for | casadi::NlpS |
|              |              | )            | calculating  | olverInterna |
|              |              |              | the gradient | l            |
|              |              |              | of the       |              |
|              |              |              | objective    |              |
|              |              |              | (column, aut |              |
|              |              |              | ogenerated   |              |
|              |              |              | by default)  |              |
+--------------+--------------+--------------+--------------+--------------+
| grad_f_optio | OT_DICT      | GenericType( | Options for  | casadi::NlpS |
| ns           |              | )            | the autogene | olverInterna |
|              |              |              | rated        | l            |
|              |              |              | gradient of  |              |
|              |              |              | the          |              |
|              |              |              | objective.   |              |
+--------------+--------------+--------------+--------------+--------------+
| grad_lag     | OT_FUNCTION  | GenericType( | Function for | casadi::NlpS |
|              |              | )            | calculating  | olverInterna |
|              |              |              | the gradient | l            |
|              |              |              | of the       |              |
|              |              |              | Lagrangian ( |              |
|              |              |              | autogenerate |              |
|              |              |              | d by         |              |
|              |              |              | default)     |              |
+--------------+--------------+--------------+--------------+--------------+
| grad_lag_opt | OT_DICT      | GenericType( | Options for  | casadi::NlpS |
| ions         |              | )            | the autogene | olverInterna |
|              |              |              | rated        | l            |
|              |              |              | gradient of  |              |
|              |              |              | the          |              |
|              |              |              | Lagrangian.  |              |
+--------------+--------------+--------------+--------------+--------------+
| hess_lag     | OT_FUNCTION  | GenericType( | Function for | casadi::NlpS |
|              |              | )            | calculating  | olverInterna |
|              |              |              | the Hessian  | l            |
|              |              |              | of the       |              |
|              |              |              | Lagrangian ( |              |
|              |              |              | autogenerate |              |
|              |              |              | d by         |              |
|              |              |              | default)     |              |
+--------------+--------------+--------------+--------------+--------------+
| hess_lag_opt | OT_DICT      | GenericType( | Options for  | casadi::NlpS |
| ions         |              | )            | the autogene | olverInterna |
|              |              |              | rated        | l            |
|              |              |              | Hessian of   |              |
|              |              |              | the          |              |
|              |              |              | Lagrangian.  |              |
+--------------+--------------+--------------+--------------+--------------+
| ignore_check | OT_BOOLEAN   | false        | If set to    | casadi::NlpS |
| _vec         |              |              | true, the    | olverInterna |
|              |              |              | input shape  | l            |
|              |              |              | of F will    |              |
|              |              |              | not be       |              |
|              |              |              | checked.     |              |
+--------------+--------------+--------------+--------------+--------------+
| input_scheme | OT_STRINGVEC | GenericType( | Custom input | casadi::Func |
|              | TOR          | )            | scheme       | tionInternal |
+--------------+--------------+--------------+--------------+--------------+
| inputs_check | OT_BOOLEAN   | true         | Throw        | casadi::Func |
|              |              |              | exceptions   | tionInternal |
|              |              |              | when the     |              |
|              |              |              | numerical    |              |
|              |              |              | values of    |              |
|              |              |              | the inputs   |              |
|              |              |              | don't make   |              |
|              |              |              | sense        |              |
+--------------+--------------+--------------+--------------+--------------+
| iteration_ca | OT_FUNCTION  | GenericType( | A function   | casadi::NlpS |
| llback       |              | )            | that will be | olverInterna |
|              |              |              | called at    | l            |
|              |              |              | each         |              |
|              |              |              | iteration    |              |
|              |              |              | with the     |              |
|              |              |              | solver as    |              |
|              |              |              | input. Check |              |
|              |              |              | documentatio |              |
|              |              |              | n of         |              |
|              |              |              | Callback .   |              |
+--------------+--------------+--------------+--------------+--------------+
| iteration_ca | OT_BOOLEAN   | false        | If set to    | casadi::NlpS |
| llback_ignor |              |              | true, errors | olverInterna |
| e_errors     |              |              | thrown by it | l            |
|              |              |              | eration_call |              |
|              |              |              | back will be |              |
|              |              |              | ignored.     |              |
+--------------+--------------+--------------+--------------+--------------+
| iteration_ca | OT_INTEGER   | 1            | Only call    | casadi::NlpS |
| llback_step  |              |              | the callback | olverInterna |
|              |              |              | function     | l            |
|              |              |              | every few    |              |
|              |              |              | iterations.  |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_f        | OT_FUNCTION  | GenericType( | Function for | casadi::NlpS |
|              |              | )            | calculating  | olverInterna |
|              |              |              | the Jacobian | l            |
|              |              |              | of the       |              |
|              |              |              | objective    |              |
|              |              |              | (sparse row, |              |
|              |              |              | autogenerate |              |
|              |              |              | d by         |              |
|              |              |              | default)     |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_f_option | OT_DICT      | GenericType( | Options for  | casadi::NlpS |
| s            |              | )            | the autogene | olverInterna |
|              |              |              | rated        | l            |
|              |              |              | Jacobian of  |              |
|              |              |              | the          |              |
|              |              |              | objective.   |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_g        | OT_FUNCTION  | GenericType( | Function for | casadi::NlpS |
|              |              | )            | calculating  | olverInterna |
|              |              |              | the Jacobian | l            |
|              |              |              | of the       |              |
|              |              |              | constraints  |              |
|              |              |              | (autogenerat |              |
|              |              |              | ed by        |              |
|              |              |              | default)     |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_g_option | OT_DICT      | GenericType( | Options for  | casadi::NlpS |
| s            |              | )            | the autogene | olverInterna |
|              |              |              | rated        | l            |
|              |              |              | Jacobian of  |              |
|              |              |              | the          |              |
|              |              |              | constraints. |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_penalty  | OT_REAL      | 2            | When         | casadi::Func |
|              |              |              | requested    | tionInternal |
|              |              |              | for a number |              |
|              |              |              | of forward/r |              |
|              |              |              | everse       |              |
|              |              |              | directions,  |              |
|              |              |              | it may be    |              |
|              |              |              | cheaper to   |              |
|              |              |              | compute      |              |
|              |              |              | first the    |              |
|              |              |              | full         |              |
|              |              |              | jacobian and |              |
|              |              |              | then         |              |
|              |              |              | multiply     |              |
|              |              |              | with seeds,  |              |
|              |              |              | rather than  |              |
|              |              |              | obtain the   |              |
|              |              |              | requested    |              |
|              |              |              | directions   |              |
|              |              |              | in a straigh |              |
|              |              |              | tforward     |              |
|              |              |              | manner.      |              |
|              |              |              | Casadi uses  |              |
|              |              |              | a heuristic  |              |
|              |              |              | to decide    |              |
|              |              |              | which is     |              |
|              |              |              | cheaper. A   |              |
|              |              |              | high value   |              |
|              |              |              | of 'jac_pena |              |
|              |              |              | lty' makes   |              |
|              |              |              | it less      |              |
|              |              |              | likely for   |              |
|              |              |              | the heurstic |              |
|              |              |              | to chose the |              |
|              |              |              | full         |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy.    |              |
|              |              |              | The special  |              |
|              |              |              | value -1     |              |
|              |              |              | indicates    |              |
|              |              |              | never to use |              |
|              |              |              | the full     |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy     |              |
+--------------+--------------+--------------+--------------+--------------+
| jit          | OT_BOOLEAN   | false        | Use just-in- | casadi::Func |
|              |              |              | time         | tionInternal |
|              |              |              | compiler to  |              |
|              |              |              | speed up the |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| jit_options  | OT_DICT      | GenericType( | Options to   | casadi::Func |
|              |              | )            | be passed to | tionInternal |
|              |              |              | the jit      |              |
|              |              |              | compiler.    |              |
+--------------+--------------+--------------+--------------+--------------+
| monitor      | OT_STRINGVEC | GenericType( | Monitors to  | casadi::Func |
|              | TOR          | )            | be activated | tionInternal |
|              |              |              | (inputs|outp |              |
|              |              |              | uts)         |              |
+--------------+--------------+--------------+--------------+--------------+
| output_schem | OT_STRINGVEC | GenericType( | Custom       | casadi::Func |
| e            | TOR          | )            | output       | tionInternal |
|              |              |              | scheme       |              |
+--------------+--------------+--------------+--------------+--------------+
| regularity_c | OT_BOOLEAN   | true         | Throw        | casadi::Func |
| heck         |              |              | exceptions   | tionInternal |
|              |              |              | when NaN or  |              |
|              |              |              | Inf appears  |              |
|              |              |              | during       |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| user_data    | OT_VOIDPTR   | GenericType( | A user-      | casadi::Func |
|              |              | )            | defined      | tionInternal |
|              |              |              | field that   |              |
|              |              |              | can be used  |              |
|              |              |              | to identify  |              |
|              |              |              | the function |              |
|              |              |              | or pass      |              |
|              |              |              | additional   |              |
|              |              |              | information  |              |
+--------------+--------------+--------------+--------------+--------------+
| verbose      | OT_BOOLEAN   | false        | Verbose      | casadi::Func |
|              |              |              | evaluation   | tionInternal |
|              |              |              | for          |              |
|              |              |              | debugging    |              |
+--------------+--------------+--------------+--------------+--------------+
| verbose_init | OT_BOOLEAN   | false        | Print out    | casadi::NlpS |
|              |              |              | timing       | olverInterna |
|              |              |              | information  | l            |
|              |              |              | about the    |              |
|              |              |              | different    |              |
|              |              |              | stages of in |              |
|              |              |              | itialization |              |
+--------------+--------------+--------------+--------------+--------------+
| warn_initial | OT_BOOLEAN   | false        | Warn if the  | casadi::NlpS |
| _bounds      |              |              | initial      | olverInterna |
|              |              |              | guess does   | l            |
|              |              |              | not satisfy  |              |
|              |              |              | LBX and UBX  |              |
+--------------+--------------+--------------+--------------+--------------+

>List of available stats

+------------------------------+---------------------------+
|              Id              |          Used in          |
+==============================+===========================+
| base class init time         | casadi::NlpSolverInternal |
+------------------------------+---------------------------+
| constraint jacobian gen time | casadi::NlpSolverInternal |
+------------------------------+---------------------------+
| grad lag gen time            | casadi::NlpSolverInternal |
+------------------------------+---------------------------+
| hess lag gen time            | casadi::NlpSolverInternal |
+------------------------------+---------------------------+
| hess lag sparsity time       | casadi::NlpSolverInternal |
+------------------------------+---------------------------+
| objective gradient gen time  | casadi::NlpSolverInternal |
+------------------------------+---------------------------+
| objective jacobian gen time  | casadi::NlpSolverInternal |
+------------------------------+---------------------------+

List of plugins
===============



- ipopt

- knitro

- snopt

- worhp

- scpgen

- sqpmethod

- stabilizedsqp

Note: some of the plugins in this list might not be available on your
system. Also, there might be extra plugins available to you that are not
listed here. You can obtain their documentation with
NlpSolver.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

ipopt
-----



When in warmstart mode, output NLP_SOLVER_LAM_X may be used as input

NOTE: Even when max_iter == 0, it is not guaranteed that
input(NLP_SOLVER_X0) == output(NLP_SOLVER_X). Indeed if bounds on X or
constraints are unmet, they will differ.

For a good tutorial on IPOPT,
seehttp://drops.dagstuhl.de/volltexte/2009/2089/pdf/09061.WaechterAndreas.Paper.2089.pdf

A good resource about the algorithms in IPOPT is: Wachter and L. T. Biegler,
On the Implementation of an Interior-Point Filter Line-Search Algorithm for
Large-Scale Nonlinear Programming, Mathematical Programming 106(1), pp.
25-57, 2006 (As Research Report RC 23149, IBM T. J. Watson Research Center,
Yorktown, USA

Caveats: with default options, multipliers for the decision variables are
wrong for equality constraints. Change the 'fixed_variable_treatment' to
'make_constraint' or 'relax_bounds' to obtain correct results.

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| accept_after_ma | OT_INTEGER      | -1              | Accept a trial  |
| x_steps         |                 |                 | point after     |
|                 |                 |                 | maximal this    |
|                 |                 |                 | number of       |
|                 |                 |                 | steps. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| accept_every_tr | OT_STRING       | no              | Always accept   |
| ial_step        |                 |                 | the first trial |
|                 |                 |                 | step. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| acceptable_comp | OT_REAL         | 0.010           | \"Acceptance\"    |
| l_inf_tol       |                 |                 | threshold for   |
|                 |                 |                 | the             |
|                 |                 |                 | complementarity |
|                 |                 |                 | conditions.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| acceptable_cons | OT_REAL         | 0.010           | \"Acceptance\"    |
| tr_viol_tol     |                 |                 | threshold for   |
|                 |                 |                 | the constraint  |
|                 |                 |                 | violation. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| acceptable_dual | OT_REAL         | 1.000e+10       | \"Acceptance\"    |
| _inf_tol        |                 |                 | threshold for   |
|                 |                 |                 | the dual        |
|                 |                 |                 | infeasibility.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| acceptable_iter | OT_INTEGER      | 15              | Number of       |
|                 |                 |                 | \"acceptable\"    |
|                 |                 |                 | iterates before |
|                 |                 |                 | triggering      |
|                 |                 |                 | termination.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| acceptable_obj_ | OT_REAL         | 1.000e+20       | \"Acceptance\"    |
| change_tol      |                 |                 | stopping        |
|                 |                 |                 | criterion based |
|                 |                 |                 | on objective    |
|                 |                 |                 | function        |
|                 |                 |                 | change. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| acceptable_tol  | OT_REAL         | 0.000           | \"Acceptable\"    |
|                 |                 |                 | convergence     |
|                 |                 |                 | tolerance       |
|                 |                 |                 | (relative).     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| adaptive_mu_glo | OT_STRING       | obj-constr-     | Globalization   |
| balization      |                 | filter          | strategy for    |
|                 |                 |                 | the adaptive mu |
|                 |                 |                 | selection mode. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| adaptive_mu_kkt | OT_STRING       | 2-norm-squared  | Norm used for   |
| _norm_type      |                 |                 | the KKT error   |
|                 |                 |                 | in the adaptive |
|                 |                 |                 | mu              |
|                 |                 |                 | globalization   |
|                 |                 |                 | strategies.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| adaptive_mu_kkt | OT_REAL         | 1.000           | Sufficient      |
| error_red_fact  |                 |                 | decrease factor |
|                 |                 |                 | for \"kkt-error\" |
|                 |                 |                 | globalization   |
|                 |                 |                 | strategy. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| adaptive_mu_kkt | OT_INTEGER      | 4               | Maximum number  |
| error_red_iters |                 |                 | of iterations   |
|                 |                 |                 | requiring       |
|                 |                 |                 | sufficient      |
|                 |                 |                 | progress. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| adaptive_mu_mon | OT_REAL         | 0.800           | Determines the  |
| otone_init_fact |                 |                 | initial value   |
| or              |                 |                 | of the barrier  |
|                 |                 |                 | parameter when  |
|                 |                 |                 | switching to    |
|                 |                 |                 | the monotone    |
|                 |                 |                 | mode. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| adaptive_mu_res | OT_STRING       | no              | Indicates if    |
| tore_previous_i |                 |                 | the previous    |
| terate          |                 |                 | iterate should  |
|                 |                 |                 | be restored if  |
|                 |                 |                 | the monotone    |
|                 |                 |                 | mode is         |
|                 |                 |                 | entered. (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| adaptive_mu_saf | OT_REAL         | 0               | (see IPOPT      |
| eguard_factor   |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| alpha_for_y     | OT_STRING       | primal          | Method to       |
|                 |                 |                 | determine the   |
|                 |                 |                 | step size for   |
|                 |                 |                 | constraint      |
|                 |                 |                 | multipliers.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| alpha_for_y_tol | OT_REAL         | 10              | Tolerance for   |
|                 |                 |                 | switching to    |
|                 |                 |                 | full equality   |
|                 |                 |                 | multiplier      |
|                 |                 |                 | steps. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| alpha_min_frac  | OT_REAL         | 0.050           | Safety factor   |
|                 |                 |                 | for the minimal |
|                 |                 |                 | step size       |
|                 |                 |                 | (before         |
|                 |                 |                 | switching to    |
|                 |                 |                 | restoration     |
|                 |                 |                 | phase). (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| alpha_red_facto | OT_REAL         | 0.500           | Fractional      |
| r               |                 |                 | reduction of    |
|                 |                 |                 | the trial step  |
|                 |                 |                 | size in the     |
|                 |                 |                 | backtracking    |
|                 |                 |                 | line search.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| barrier_tol_fac | OT_REAL         | 10              | Factor for mu   |
| tor             |                 |                 | in barrier stop |
|                 |                 |                 | test. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| bound_frac      | OT_REAL         | 0.010           | Desired minimum |
|                 |                 |                 | relative        |
|                 |                 |                 | distance from   |
|                 |                 |                 | the initial     |
|                 |                 |                 | point to bound. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| bound_mult_init | OT_STRING       | constant        | Initialization  |
| _method         |                 |                 | method for      |
|                 |                 |                 | bound           |
|                 |                 |                 | multipliers     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| bound_mult_init | OT_REAL         | 1               | Initial value   |
| _val            |                 |                 | for the bound   |
|                 |                 |                 | multipliers.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| bound_mult_rese | OT_REAL         | 1000            | Threshold for   |
| t_threshold     |                 |                 | resetting bound |
|                 |                 |                 | multipliers     |
|                 |                 |                 | after the       |
|                 |                 |                 | restoration     |
|                 |                 |                 | phase. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| bound_push      | OT_REAL         | 0.010           | Desired minimum |
|                 |                 |                 | absolute        |
|                 |                 |                 | distance from   |
|                 |                 |                 | the initial     |
|                 |                 |                 | point to bound. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| bound_relax_fac | OT_REAL         | 0.000           | Factor for      |
| tor             |                 |                 | initial         |
|                 |                 |                 | relaxation of   |
|                 |                 |                 | the bounds.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| check_derivativ | OT_STRING       | no              | Indicates       |
| es_for_naninf   |                 |                 | whether it is   |
|                 |                 |                 | desired to      |
|                 |                 |                 | check for       |
|                 |                 |                 | Nan/Inf in      |
|                 |                 |                 | derivative      |
|                 |                 |                 | matrices (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| chi_cup         | OT_REAL         | 1.500           | LIFENG WRITES   |
|                 |                 |                 | THIS. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| chi_hat         | OT_REAL         | 2               | LIFENG WRITES   |
|                 |                 |                 | THIS. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| chi_tilde       | OT_REAL         | 5               | LIFENG WRITES   |
|                 |                 |                 | THIS. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| compl_inf_tol   | OT_REAL         | 0.000           | Desired         |
|                 |                 |                 | threshold for   |
|                 |                 |                 | the             |
|                 |                 |                 | complementarity |
|                 |                 |                 | conditions.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| con_integer_md  | OT_DICT         | None            | Integer         |
|                 |                 |                 | metadata (a     |
|                 |                 |                 | dictionary with |
|                 |                 |                 | lists of        |
|                 |                 |                 | integers) about |
|                 |                 |                 | constraints to  |
|                 |                 |                 | be passed to    |
|                 |                 |                 | IPOPT           |
+-----------------+-----------------+-----------------+-----------------+
| con_numeric_md  | OT_DICT         | None            | Numeric         |
|                 |                 |                 | metadata (a     |
|                 |                 |                 | dictionary with |
|                 |                 |                 | lists of reals) |
|                 |                 |                 | about           |
|                 |                 |                 | constraints to  |
|                 |                 |                 | be passed to    |
|                 |                 |                 | IPOPT           |
+-----------------+-----------------+-----------------+-----------------+
| con_string_md   | OT_DICT         | None            | String metadata |
|                 |                 |                 | (a dictionary   |
|                 |                 |                 | with lists of   |
|                 |                 |                 | strings) about  |
|                 |                 |                 | constraints to  |
|                 |                 |                 | be passed to    |
|                 |                 |                 | IPOPT           |
+-----------------+-----------------+-----------------+-----------------+
| constr_mult_ini | OT_REAL         | 1000            | Maximum allowed |
| t_max           |                 |                 | least-square    |
|                 |                 |                 | guess of        |
|                 |                 |                 | constraint      |
|                 |                 |                 | multipliers.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| constr_mult_res | OT_REAL         | 0               | Threshold for   |
| et_threshold    |                 |                 | resetting       |
|                 |                 |                 | equality and    |
|                 |                 |                 | inequality      |
|                 |                 |                 | multipliers     |
|                 |                 |                 | after           |
|                 |                 |                 | restoration     |
|                 |                 |                 | phase. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| constr_viol_tol | OT_REAL         | 0.000           | Desired         |
|                 |                 |                 | threshold for   |
|                 |                 |                 | the constraint  |
|                 |                 |                 | violation. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| constraint_viol | OT_STRING       | 1-norm          | Norm to be used |
| ation_norm_type |                 |                 | for the         |
|                 |                 |                 | constraint      |
|                 |                 |                 | violation in    |
|                 |                 |                 | the line        |
|                 |                 |                 | search. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| corrector_compl | OT_REAL         | 1               | Complementarity |
| _avrg_red_fact  |                 |                 | tolerance       |
|                 |                 |                 | factor for      |
|                 |                 |                 | accepting       |
|                 |                 |                 | corrector step  |
|                 |                 |                 | (unsupported!). |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| corrector_type  | OT_STRING       | none            | The type of     |
|                 |                 |                 | corrector steps |
|                 |                 |                 | that should be  |
|                 |                 |                 | taken           |
|                 |                 |                 | (unsupported!). |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| delta           | OT_REAL         | 1               | Multiplier for  |
|                 |                 |                 | constraint      |
|                 |                 |                 | violation in    |
|                 |                 |                 | the switching   |
|                 |                 |                 | rule. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| delta_y_max     | OT_REAL         | 1.000e+12       | a parameter     |
|                 |                 |                 | used to check   |
|                 |                 |                 | if the fast     |
|                 |                 |                 | direction can   |
|                 |                 |                 | be used asthe   |
|                 |                 |                 | line search     |
|                 |                 |                 | direction (for  |
|                 |                 |                 | Chen-Goldfarb   |
|                 |                 |                 | line search).   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| dependency_dete | OT_STRING       | no              | Indicates if    |
| ction_with_rhs  |                 |                 | the right hand  |
|                 |                 |                 | sides of the    |
|                 |                 |                 | constraints     |
|                 |                 |                 | should be       |
|                 |                 |                 | considered      |
|                 |                 |                 | during          |
|                 |                 |                 | dependency      |
|                 |                 |                 | detection (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| dependency_dete | OT_STRING       | none            | Indicates which |
| ctor            |                 |                 | linear solver   |
|                 |                 |                 | should be used  |
|                 |                 |                 | to detect       |
|                 |                 |                 | linearly        |
|                 |                 |                 | dependent       |
|                 |                 |                 | equality        |
|                 |                 |                 | constraints.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| derivative_test | OT_STRING       | none            | Enable          |
|                 |                 |                 | derivative      |
|                 |                 |                 | checker (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| derivative_test | OT_INTEGER      | -2              | Index of first  |
| _first_index    |                 |                 | quantity to be  |
|                 |                 |                 | checked by      |
|                 |                 |                 | derivative      |
|                 |                 |                 | checker (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| derivative_test | OT_REAL         | 0.000           | Size of the     |
| _perturbation   |                 |                 | finite          |
|                 |                 |                 | difference      |
|                 |                 |                 | perturbation in |
|                 |                 |                 | derivative      |
|                 |                 |                 | test. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| derivative_test | OT_STRING       | no              | Indicates       |
| _print_all      |                 |                 | whether         |
|                 |                 |                 | information for |
|                 |                 |                 | all estimated   |
|                 |                 |                 | derivatives     |
|                 |                 |                 | should be       |
|                 |                 |                 | printed. (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| derivative_test | OT_REAL         | 0.000           | Threshold for   |
| _tol            |                 |                 | indicating      |
|                 |                 |                 | wrong           |
|                 |                 |                 | derivative.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| diverging_itera | OT_REAL         | 1.000e+20       | Threshold for   |
| tes_tol         |                 |                 | maximal value   |
|                 |                 |                 | of primal       |
|                 |                 |                 | iterates. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| dual_inf_tol    | OT_REAL         | 1               | Desired         |
|                 |                 |                 | threshold for   |
|                 |                 |                 | the dual        |
|                 |                 |                 | infeasibility.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| epsilon_c       | OT_REAL         | 0.010           | LIFENG WRITES   |
|                 |                 |                 | THIS. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| eta_min         | OT_REAL         | 10              | LIFENG WRITES   |
|                 |                 |                 | THIS. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| eta_penalty     | OT_REAL         | 0.000           | Relaxation      |
|                 |                 |                 | factor in the   |
|                 |                 |                 | Armijo          |
|                 |                 |                 | condition for   |
|                 |                 |                 | the penalty     |
|                 |                 |                 | function. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| eta_phi         | OT_REAL         | 0.000           | Relaxation      |
|                 |                 |                 | factor in the   |
|                 |                 |                 | Armijo          |
|                 |                 |                 | condition. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| evaluate_orig_o | OT_STRING       | yes             | Determines if   |
| bj_at_resto_tri |                 |                 | the original    |
| al              |                 |                 | objective       |
|                 |                 |                 | function should |
|                 |                 |                 | be evaluated at |
|                 |                 |                 | restoration     |
|                 |                 |                 | phase trial     |
|                 |                 |                 | points. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| expect_infeasib | OT_STRING       | no              | Enable          |
| le_problem      |                 |                 | heuristics to   |
|                 |                 |                 | quickly detect  |
|                 |                 |                 | an infeasible   |
|                 |                 |                 | problem. (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| expect_infeasib | OT_REAL         | 0.001           | Threshold for   |
| le_problem_ctol |                 |                 | disabling \"expe |
|                 |                 |                 | ct_infeasible_p |
|                 |                 |                 | roblem\" option. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| expect_infeasib | OT_REAL         | 100000000       | Multiplier      |
| le_problem_ytol |                 |                 | threshold for   |
|                 |                 |                 | activating \"exp |
|                 |                 |                 | ect_infeasible_ |
|                 |                 |                 | problem\"        |
|                 |                 |                 | option. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| fast_des_fact   | OT_REAL         | 0.100           | a parameter     |
|                 |                 |                 | used to check   |
|                 |                 |                 | if the fast     |
|                 |                 |                 | direction can   |
|                 |                 |                 | be used asthe   |
|                 |                 |                 | line search     |
|                 |                 |                 | direction (for  |
|                 |                 |                 | Chen-Goldfarb   |
|                 |                 |                 | line search).   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| fast_step_compu | OT_STRING       | no              | Indicates if    |
| tation          |                 |                 | the linear      |
|                 |                 |                 | system should   |
|                 |                 |                 | be solved       |
|                 |                 |                 | quickly. (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| file_print_leve | OT_INTEGER      | 5               | Verbosity level |
| l               |                 |                 | for output      |
|                 |                 |                 | file. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| filter_margin_f | OT_REAL         | 0.000           | Factor          |
| act             |                 |                 | determining     |
|                 |                 |                 | width of margin |
|                 |                 |                 | for obj-constr- |
|                 |                 |                 | filter adaptive |
|                 |                 |                 | globalization   |
|                 |                 |                 | strategy. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| filter_max_marg | OT_REAL         | 1               | Maximum width   |
| in              |                 |                 | of margin in    |
|                 |                 |                 | obj-constr-     |
|                 |                 |                 | filter adaptive |
|                 |                 |                 | globalization   |
|                 |                 |                 | strategy. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| filter_reset_tr | OT_INTEGER      | 5               | Number of       |
| igger           |                 |                 | iterations that |
|                 |                 |                 | trigger the     |
|                 |                 |                 | filter reset.   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| findiff_perturb | OT_REAL         | 0.000           | Size of the     |
| ation           |                 |                 | finite          |
|                 |                 |                 | difference      |
|                 |                 |                 | perturbation    |
|                 |                 |                 | for derivative  |
|                 |                 |                 | approximation.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| first_hessian_p | OT_REAL         | 0.000           | Size of first   |
| erturbation     |                 |                 | x-s             |
|                 |                 |                 | perturbation    |
|                 |                 |                 | tried. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| fixed_mu_oracle | OT_STRING       | average_compl   | Oracle for the  |
|                 |                 |                 | barrier         |
|                 |                 |                 | parameter when  |
|                 |                 |                 | switching to    |
|                 |                 |                 | fixed mode.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| fixed_variable_ | OT_STRING       | make_parameter  | Determines how  |
| treatment       |                 |                 | fixed variables |
|                 |                 |                 | should be       |
|                 |                 |                 | handled. (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| gamma_hat       | OT_REAL         | 0.040           | LIFENG WRITES   |
|                 |                 |                 | THIS. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| gamma_phi       | OT_REAL         | 0.000           | Relaxation      |
|                 |                 |                 | factor in the   |
|                 |                 |                 | filter margin   |
|                 |                 |                 | for the barrier |
|                 |                 |                 | function. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| gamma_theta     | OT_REAL         | 0.000           | Relaxation      |
|                 |                 |                 | factor in the   |
|                 |                 |                 | filter margin   |
|                 |                 |                 | for the         |
|                 |                 |                 | constraint      |
|                 |                 |                 | violation. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| gamma_tilde     | OT_REAL         | 4               | LIFENG WRITES   |
|                 |                 |                 | THIS. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| hessian_approxi | OT_STRING       | exact           | Indicates what  |
| mation          |                 |                 | Hessian         |
|                 |                 |                 | information is  |
|                 |                 |                 | to be used.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| hessian_approxi | OT_STRING       | nonlinear-      | Indicates in    |
| mation_space    |                 | variables       | which subspace  |
|                 |                 |                 | the Hessian     |
|                 |                 |                 | information is  |
|                 |                 |                 | to be           |
|                 |                 |                 | approximated.   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| hessian_constan | OT_STRING       | no              | Indicates       |
| t               |                 |                 | whether the     |
|                 |                 |                 | problem is a    |
|                 |                 |                 | quadratic       |
|                 |                 |                 | problem (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| honor_original_ | OT_STRING       | yes             | Indicates       |
| bounds          |                 |                 | whether final   |
|                 |                 |                 | points should   |
|                 |                 |                 | be projected    |
|                 |                 |                 | into original   |
|                 |                 |                 | bounds. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| inf_pr_output   | OT_STRING       | original        | Determines what |
|                 |                 |                 | value is        |
|                 |                 |                 | printed in the  |
|                 |                 |                 | \"inf_pr\" output |
|                 |                 |                 | column. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| jac_c_constant  | OT_STRING       | no              | Indicates       |
|                 |                 |                 | whether all     |
|                 |                 |                 | equality        |
|                 |                 |                 | constraints are |
|                 |                 |                 | linear (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| jac_d_constant  | OT_STRING       | no              | Indicates       |
|                 |                 |                 | whether all     |
|                 |                 |                 | inequality      |
|                 |                 |                 | constraints are |
|                 |                 |                 | linear (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| jacobian_approx | OT_STRING       | exact           | Specifies       |
| imation         |                 |                 | technique to    |
|                 |                 |                 | compute         |
|                 |                 |                 | constraint      |
|                 |                 |                 | Jacobian (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| jacobian_regula | OT_REAL         | 0.250           | Exponent for mu |
| rization_expone |                 |                 | in the          |
| nt              |                 |                 | regularization  |
|                 |                 |                 | for rank-       |
|                 |                 |                 | deficient       |
|                 |                 |                 | constraint      |
|                 |                 |                 | Jacobians. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| jacobian_regula | OT_REAL         | 0.000           | Size of the     |
| rization_value  |                 |                 | regularization  |
|                 |                 |                 | for rank-       |
|                 |                 |                 | deficient       |
|                 |                 |                 | constraint      |
|                 |                 |                 | Jacobians. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| kappa_d         | OT_REAL         | 0.000           | Weight for      |
|                 |                 |                 | linear damping  |
|                 |                 |                 | term (to handle |
|                 |                 |                 | one-sided       |
|                 |                 |                 | bounds). (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| kappa_sigma     | OT_REAL         | 1.000e+10       | Factor limiting |
|                 |                 |                 | the deviation   |
|                 |                 |                 | of dual         |
|                 |                 |                 | variables from  |
|                 |                 |                 | primal          |
|                 |                 |                 | estimates. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| kappa_soc       | OT_REAL         | 0.990           | Factor in the   |
|                 |                 |                 | sufficient      |
|                 |                 |                 | reduction rule  |
|                 |                 |                 | for second      |
|                 |                 |                 | order           |
|                 |                 |                 | correction.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| kappa_x_dis     | OT_REAL         | 100             | a parameter     |
|                 |                 |                 | used to check   |
|                 |                 |                 | if the fast     |
|                 |                 |                 | direction can   |
|                 |                 |                 | be used asthe   |
|                 |                 |                 | line search     |
|                 |                 |                 | direction (for  |
|                 |                 |                 | Chen-Goldfarb   |
|                 |                 |                 | line search).   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| kappa_y_dis     | OT_REAL         | 10000           | a parameter     |
|                 |                 |                 | used to check   |
|                 |                 |                 | if the fast     |
|                 |                 |                 | direction can   |
|                 |                 |                 | be used asthe   |
|                 |                 |                 | line search     |
|                 |                 |                 | direction (for  |
|                 |                 |                 | Chen-Goldfarb   |
|                 |                 |                 | line search).   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| least_square_in | OT_STRING       | no              | Least square    |
| it_duals        |                 |                 | initialization  |
|                 |                 |                 | of all dual     |
|                 |                 |                 | variables (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| least_square_in | OT_STRING       | no              | Least square    |
| it_primal       |                 |                 | initialization  |
|                 |                 |                 | of the primal   |
|                 |                 |                 | variables (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| limited_memory_ | OT_STRING       | sherman-        | Strategy for    |
| aug_solver      |                 | morrison        | solving the     |
|                 |                 |                 | augmented       |
|                 |                 |                 | system for low- |
|                 |                 |                 | rank Hessian.   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| limited_memory_ | OT_REAL         | 1               | Value for B0 in |
| init_val        |                 |                 | low-rank        |
|                 |                 |                 | update. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| limited_memory_ | OT_REAL         | 100000000       | Upper bound on  |
| init_val_max    |                 |                 | value for B0 in |
|                 |                 |                 | low-rank        |
|                 |                 |                 | update. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| limited_memory_ | OT_REAL         | 0.000           | Lower bound on  |
| init_val_min    |                 |                 | value for B0 in |
|                 |                 |                 | low-rank        |
|                 |                 |                 | update. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| limited_memory_ | OT_STRING       | scalar1         | Initialization  |
| initialization  |                 |                 | strategy for    |
|                 |                 |                 | the limited     |
|                 |                 |                 | memory quasi-   |
|                 |                 |                 | Newton          |
|                 |                 |                 | approximation.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| limited_memory_ | OT_INTEGER      | 6               | Maximum size of |
| max_history     |                 |                 | the history for |
|                 |                 |                 | the limited     |
|                 |                 |                 | quasi-Newton    |
|                 |                 |                 | Hessian         |
|                 |                 |                 | approximation.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| limited_memory_ | OT_INTEGER      | 2               | Threshold for   |
| max_skipping    |                 |                 | successive      |
|                 |                 |                 | iterations      |
|                 |                 |                 | where update is |
|                 |                 |                 | skipped. (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| limited_memory_ | OT_STRING       | no              | Determines if   |
| special_for_res |                 |                 | the quasi-      |
| to              |                 |                 | Newton updates  |
|                 |                 |                 | should be       |
|                 |                 |                 | special during  |
|                 |                 |                 | the restoration |
|                 |                 |                 | phase. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| limited_memory_ | OT_STRING       | bfgs            | Quasi-Newton    |
| update_type     |                 |                 | update formula  |
|                 |                 |                 | for the limited |
|                 |                 |                 | memory          |
|                 |                 |                 | approximation.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| line_search_met | OT_STRING       | filter          | Globalization   |
| hod             |                 |                 | method used in  |
|                 |                 |                 | backtracking    |
|                 |                 |                 | line search     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| linear_scaling_ | OT_STRING       | yes             | Flag indicating |
| on_demand       |                 |                 | that linear     |
|                 |                 |                 | scaling is only |
|                 |                 |                 | done if it      |
|                 |                 |                 | seems required. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| linear_solver   | OT_STRING       | mumps           | Linear solver   |
|                 |                 |                 | used for step   |
|                 |                 |                 | computations.   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| linear_system_s | OT_STRING       | none            | Method for      |
| caling          |                 |                 | scaling the     |
|                 |                 |                 | linear system.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma27_ignore_sin | OT_STRING       | no              | Enables MA27's  |
| gularity        |                 |                 | ability to      |
|                 |                 |                 | solve a linear  |
|                 |                 |                 | system even if  |
|                 |                 |                 | the matrix is   |
|                 |                 |                 | singular. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma27_la_init_fa | OT_REAL         | 5               | Real workspace  |
| ctor            |                 |                 | memory for      |
|                 |                 |                 | MA27. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma27_liw_init_f | OT_REAL         | 5               | Integer         |
| actor           |                 |                 | workspace       |
|                 |                 |                 | memory for      |
|                 |                 |                 | MA27. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma27_meminc_fac | OT_REAL         | 10              | Increment       |
| tor             |                 |                 | factor for      |
|                 |                 |                 | workspace size  |
|                 |                 |                 | for MA27. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma27_pivtol     | OT_REAL         | 0.000           | Pivot tolerance |
|                 |                 |                 | for the linear  |
|                 |                 |                 | solver MA27.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma27_pivtolmax  | OT_REAL         | 0.000           | Maximum pivot   |
|                 |                 |                 | tolerance for   |
|                 |                 |                 | the linear      |
|                 |                 |                 | solver MA27.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma27_skip_inert | OT_STRING       | no              | Always pretend  |
| ia_check        |                 |                 | inertia is      |
|                 |                 |                 | correct. (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma28_pivtol     | OT_REAL         | 0.010           | Pivot tolerance |
|                 |                 |                 | for linear      |
|                 |                 |                 | solver MA28.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma57_automatic_ | OT_STRING       | yes             | Controls MA57   |
| scaling         |                 |                 | automatic       |
|                 |                 |                 | scaling (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma57_block_size | OT_INTEGER      | 16              | Controls block  |
|                 |                 |                 | size used by    |
|                 |                 |                 | Level 3 BLAS in |
|                 |                 |                 | MA57BD (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma57_node_amalg | OT_INTEGER      | 16              | Node            |
| amation         |                 |                 | amalgamation    |
|                 |                 |                 | parameter (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma57_pivot_orde | OT_INTEGER      | 5               | Controls pivot  |
| r               |                 |                 | order in MA57   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma57_pivtol     | OT_REAL         | 0.000           | Pivot tolerance |
|                 |                 |                 | for the linear  |
|                 |                 |                 | solver MA57.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma57_pivtolmax  | OT_REAL         | 0.000           | Maximum pivot   |
|                 |                 |                 | tolerance for   |
|                 |                 |                 | the linear      |
|                 |                 |                 | solver MA57.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma57_pre_alloc  | OT_REAL         | 1.050           | Safety factor   |
|                 |                 |                 | for work space  |
|                 |                 |                 | memory          |
|                 |                 |                 | allocation for  |
|                 |                 |                 | the linear      |
|                 |                 |                 | solver MA57.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma57_small_pivo | OT_INTEGER      | 0               | If set to 1,    |
| t_flag          |                 |                 | then when small |
|                 |                 |                 | entries defined |
|                 |                 |                 | by CNTL(2) are  |
|                 |                 |                 | detected they   |
|                 |                 |                 | are removed and |
|                 |                 |                 | the             |
|                 |                 |                 | corresponding   |
|                 |                 |                 | pivots placed   |
|                 |                 |                 | at the end of   |
|                 |                 |                 | the             |
|                 |                 |                 | factorization.  |
|                 |                 |                 | This can be     |
|                 |                 |                 | particularly    |
|                 |                 |                 | efficient if    |
|                 |                 |                 | the matrix is   |
|                 |                 |                 | highly rank     |
|                 |                 |                 | deficient. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma86_nemin      | OT_INTEGER      | 32              | Node            |
|                 |                 |                 | Amalgamation    |
|                 |                 |                 | parameter (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma86_print_leve | OT_INTEGER      | 0               | Debug printing  |
| l               |                 |                 | level for the   |
|                 |                 |                 | linear solver   |
|                 |                 |                 | MA86 (see IPOPT |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma86_small      | OT_REAL         | 0.000           | Zero Pivot      |
|                 |                 |                 | Threshold (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma86_static     | OT_REAL         | 0               | Static Pivoting |
|                 |                 |                 | Threshold (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma86_u          | OT_REAL         | 0.000           | Pivoting        |
|                 |                 |                 | Threshold (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| ma86_umax       | OT_REAL         | 0.000           | Maximum         |
|                 |                 |                 | Pivoting        |
|                 |                 |                 | Threshold (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| magic_steps     | OT_STRING       | no              | Enables magic   |
|                 |                 |                 | steps. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| max_cpu_time    | OT_REAL         | 1000000         | Maximum number  |
|                 |                 |                 | of CPU seconds. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| max_filter_rese | OT_INTEGER      | 5               | Maximal allowed |
| ts              |                 |                 | number of       |
|                 |                 |                 | filter resets   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| max_hessian_per | OT_REAL         | 1.000e+20       | Maximum value   |
| turbation       |                 |                 | of              |
|                 |                 |                 | regularization  |
|                 |                 |                 | parameter for   |
|                 |                 |                 | handling        |
|                 |                 |                 | negative        |
|                 |                 |                 | curvature. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| max_iter        | OT_INTEGER      | 3000            | Maximum number  |
|                 |                 |                 | of iterations.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| max_refinement_ | OT_INTEGER      | 10              | Maximum number  |
| steps           |                 |                 | of iterative    |
|                 |                 |                 | refinement      |
|                 |                 |                 | steps per       |
|                 |                 |                 | linear system   |
|                 |                 |                 | solve. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| max_resto_iter  | OT_INTEGER      | 3000000         | Maximum number  |
|                 |                 |                 | of successive   |
|                 |                 |                 | iterations in   |
|                 |                 |                 | restoration     |
|                 |                 |                 | phase. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| max_soc         | OT_INTEGER      | 4               | Maximum number  |
|                 |                 |                 | of second order |
|                 |                 |                 | correction      |
|                 |                 |                 | trial steps at  |
|                 |                 |                 | each iteration. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| max_soft_resto_ | OT_INTEGER      | 10              | Maximum number  |
| iters           |                 |                 | of iterations   |
|                 |                 |                 | performed       |
|                 |                 |                 | successively in |
|                 |                 |                 | soft            |
|                 |                 |                 | restoration     |
|                 |                 |                 | phase. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mehrotra_algori | OT_STRING       | no              | Indicates if we |
| thm             |                 |                 | want to do      |
|                 |                 |                 | Mehrotra's      |
|                 |                 |                 | algorithm. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| min_alpha_prima | OT_REAL         | 0.000           | LIFENG WRITES   |
| l               |                 |                 | THIS. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| min_hessian_per | OT_REAL         | 0.000           | Smallest        |
| turbation       |                 |                 | perturbation of |
|                 |                 |                 | the Hessian     |
|                 |                 |                 | block. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| min_refinement_ | OT_INTEGER      | 1               | Minimum number  |
| steps           |                 |                 | of iterative    |
|                 |                 |                 | refinement      |
|                 |                 |                 | steps per       |
|                 |                 |                 | linear system   |
|                 |                 |                 | solve. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mu_allow_fast_m | OT_STRING       | yes             | Allow skipping  |
| onotone_decreas |                 |                 | of barrier      |
| e               |                 |                 | problem if      |
|                 |                 |                 | barrier test is |
|                 |                 |                 | already met.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mu_init         | OT_REAL         | 0.100           | Initial value   |
|                 |                 |                 | for the barrier |
|                 |                 |                 | parameter. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mu_linear_decre | OT_REAL         | 0.200           | Determines      |
| ase_factor      |                 |                 | linear decrease |
|                 |                 |                 | rate of barrier |
|                 |                 |                 | parameter. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mu_max          | OT_REAL         | 100000          | Maximum value   |
|                 |                 |                 | for barrier     |
|                 |                 |                 | parameter. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mu_max_fact     | OT_REAL         | 1000            | Factor for      |
|                 |                 |                 | initialization  |
|                 |                 |                 | of maximum      |
|                 |                 |                 | value for       |
|                 |                 |                 | barrier         |
|                 |                 |                 | parameter. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mu_min          | OT_REAL         | 0.000           | Minimum value   |
|                 |                 |                 | for barrier     |
|                 |                 |                 | parameter. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mu_oracle       | OT_STRING       | quality-        | Oracle for a    |
|                 |                 | function        | new barrier     |
|                 |                 |                 | parameter in    |
|                 |                 |                 | the adaptive    |
|                 |                 |                 | strategy. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mu_strategy     | OT_STRING       | monotone        | Update strategy |
|                 |                 |                 | for barrier     |
|                 |                 |                 | parameter. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mu_superlinear_ | OT_REAL         | 1.500           | Determines      |
| decrease_power  |                 |                 | superlinear     |
|                 |                 |                 | decrease rate   |
|                 |                 |                 | of barrier      |
|                 |                 |                 | parameter. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mu_target       | OT_REAL         | 0               | Desired value   |
|                 |                 |                 | of complementar |
|                 |                 |                 | ity. (see IPOPT |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mult_diverg_fea | OT_REAL         | 0.000           | tolerance for   |
| sibility_tol    |                 |                 | deciding if the |
|                 |                 |                 | multipliers are |
|                 |                 |                 | diverging (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mult_diverg_y_t | OT_REAL         | 100000000       | tolerance for   |
| ol              |                 |                 | deciding if the |
|                 |                 |                 | multipliers are |
|                 |                 |                 | diverging (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mumps_dep_tol   | OT_REAL         | -1              | Pivot threshold |
|                 |                 |                 | for detection   |
|                 |                 |                 | of linearly     |
|                 |                 |                 | dependent       |
|                 |                 |                 | constraints in  |
|                 |                 |                 | MUMPS. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mumps_mem_perce | OT_INTEGER      | 1000            | Percentage      |
| nt              |                 |                 | increase in the |
|                 |                 |                 | estimated       |
|                 |                 |                 | working space   |
|                 |                 |                 | for MUMPS. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mumps_permuting | OT_INTEGER      | 7               | Controls        |
| _scaling        |                 |                 | permuting and   |
|                 |                 |                 | scaling in      |
|                 |                 |                 | MUMPS (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mumps_pivot_ord | OT_INTEGER      | 7               | Controls pivot  |
| er              |                 |                 | order in MUMPS  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mumps_pivtol    | OT_REAL         | 0.000           | Pivot tolerance |
|                 |                 |                 | for the linear  |
|                 |                 |                 | solver MUMPS.   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mumps_pivtolmax | OT_REAL         | 0.100           | Maximum pivot   |
|                 |                 |                 | tolerance for   |
|                 |                 |                 | the linear      |
|                 |                 |                 | solver MUMPS.   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| mumps_scaling   | OT_INTEGER      | 77              | Controls        |
|                 |                 |                 | scaling in      |
|                 |                 |                 | MUMPS (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| neg_curv_test_t | OT_REAL         | 0               | Tolerance for   |
| ol              |                 |                 | heuristic to    |
|                 |                 |                 | ignore wrong    |
|                 |                 |                 | inertia. (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| never_use_fact_ | OT_STRING       | no              | Toggle to       |
| cgpen_direction |                 |                 | switch off the  |
|                 |                 |                 | fast Chen-      |
|                 |                 |                 | Goldfarb        |
|                 |                 |                 | direction (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| never_use_piece | OT_STRING       | no              | Toggle to       |
| wise_penalty_ls |                 |                 | switch off the  |
|                 |                 |                 | piecewise       |
|                 |                 |                 | penalty method  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| nlp_lower_bound | OT_REAL         | -1.000e+19      | any bound less  |
| _inf            |                 |                 | or equal this   |
|                 |                 |                 | value will be   |
|                 |                 |                 | considered -inf |
|                 |                 |                 | (i.e. not lower |
|                 |                 |                 | bounded). (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| nlp_scaling_con | OT_REAL         | 0               | Target value    |
| str_target_grad |                 |                 | for constraint  |
| ient            |                 |                 | function        |
|                 |                 |                 | gradient size.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| nlp_scaling_max | OT_REAL         | 100             | Maximum         |
| _gradient       |                 |                 | gradient after  |
|                 |                 |                 | NLP scaling.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| nlp_scaling_met | OT_STRING       | gradient-based  | Select the      |
| hod             |                 |                 | technique used  |
|                 |                 |                 | for scaling the |
|                 |                 |                 | NLP. (see IPOPT |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| nlp_scaling_min | OT_REAL         | 0.000           | Minimum value   |
| _value          |                 |                 | of gradient-    |
|                 |                 |                 | based scaling   |
|                 |                 |                 | values. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| nlp_scaling_obj | OT_REAL         | 0               | Target value    |
| _target_gradien |                 |                 | for objective   |
| t               |                 |                 | function        |
|                 |                 |                 | gradient size.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| nlp_upper_bound | OT_REAL         | 1.000e+19       | any bound       |
| _inf            |                 |                 | greater or this |
|                 |                 |                 | value will be   |
|                 |                 |                 | considered +inf |
|                 |                 |                 | (i.e. not upper |
|                 |                 |                 | bounded). (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| nu_inc          | OT_REAL         | 0.000           | Increment of    |
|                 |                 |                 | the penalty     |
|                 |                 |                 | parameter. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| nu_init         | OT_REAL         | 0.000           | Initial value   |
|                 |                 |                 | of the penalty  |
|                 |                 |                 | parameter. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| num_linear_vari | OT_INTEGER      | 0               | Number of       |
| ables           |                 |                 | linear          |
|                 |                 |                 | variables (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| obj_max_inc     | OT_REAL         | 5               | Determines the  |
|                 |                 |                 | upper bound on  |
|                 |                 |                 | the acceptable  |
|                 |                 |                 | increase of     |
|                 |                 |                 | barrier         |
|                 |                 |                 | objective       |
|                 |                 |                 | function. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| obj_scaling_fac | OT_REAL         | 1               | Scaling factor  |
| tor             |                 |                 | for the         |
|                 |                 |                 | objective       |
|                 |                 |                 | function. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| option_file_nam | OT_STRING       |                 | File name of    |
| e               |                 |                 | options file    |
|                 |                 |                 | (to overwrite   |
|                 |                 |                 | default). (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| output_file     | OT_STRING       |                 | File name of    |
|                 |                 |                 | desired output  |
|                 |                 |                 | file (leave     |
|                 |                 |                 | unset for no    |
|                 |                 |                 | file output).   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_iter_co | OT_INTEGER      | 5000            | Maximum Size of |
| arse_size       |                 |                 | Coarse Grid     |
|                 |                 |                 | Matrix (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_iter_dr | OT_REAL         | 0.500           | dropping value  |
| opping_factor   |                 |                 | for incomplete  |
|                 |                 |                 | factor (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_iter_dr | OT_REAL         | 0.100           | dropping value  |
| opping_schur    |                 |                 | for sparsify    |
|                 |                 |                 | schur           |
|                 |                 |                 | complement      |
|                 |                 |                 | factor (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_iter_in | OT_REAL         | 5000000         | (see IPOPT      |
| verse_norm_fact |                 |                 | documentation)  |
| or              |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_iter_ma | OT_INTEGER      | 10              | Maximum Size of |
| x_levels        |                 |                 | Grid Levels     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_iter_ma | OT_INTEGER      | 10000000        | max fill for    |
| x_row_fill      |                 |                 | each row (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_iter_re | OT_REAL         | 0.000           | Relative        |
| lative_tol      |                 |                 | Residual        |
|                 |                 |                 | Convergence     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_iterati | OT_STRING       | no              | Switch on       |
| ve              |                 |                 | iterative       |
|                 |                 |                 | solver in       |
|                 |                 |                 | Pardiso library |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_matchin | OT_STRING       | complete+2x2    | Matching        |
| g_strategy      |                 |                 | strategy to be  |
|                 |                 |                 | used by Pardiso |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_max_dro | OT_INTEGER      | 4               | Maximal number  |
| ptol_correction |                 |                 | of decreases of |
| s               |                 |                 | drop tolerance  |
|                 |                 |                 | during one      |
|                 |                 |                 | solve. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_max_ite | OT_INTEGER      | 500             | Maximum number  |
| r               |                 |                 | of Krylov-      |
|                 |                 |                 | Subspace        |
|                 |                 |                 | Iteration (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_msglvl  | OT_INTEGER      | 0               | Pardiso message |
|                 |                 |                 | level (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_out_of_ | OT_INTEGER      | 0               | Enables out-of- |
| core_power      |                 |                 | core variant of |
|                 |                 |                 | Pardiso (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_redo_sy | OT_STRING       | no              | Toggle for      |
| mbolic_fact_onl |                 |                 | handling case   |
| y_if_inertia_wr |                 |                 | when elements   |
| ong             |                 |                 | were perturbed  |
|                 |                 |                 | by Pardiso.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_repeate | OT_STRING       | no              | Interpretation  |
| d_perturbation_ |                 |                 | of perturbed    |
| means_singular  |                 |                 | elements. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pardiso_skip_in | OT_STRING       | no              | Always pretend  |
| ertia_check     |                 |                 | inertia is      |
|                 |                 |                 | correct. (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pass_nonlinear_ | OT_BOOLEAN      | False           | n/a             |
| variables       |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| pen_des_fact    | OT_REAL         | 0.200           | a parameter     |
|                 |                 |                 | used in penalty |
|                 |                 |                 | parameter       |
|                 |                 |                 | computation     |
|                 |                 |                 | (for Chen-      |
|                 |                 |                 | Goldfarb line   |
|                 |                 |                 | search). (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pen_init_fac    | OT_REAL         | 50              | a parameter     |
|                 |                 |                 | used to choose  |
|                 |                 |                 | initial penalty |
|                 |                 |                 | parameterswhen  |
|                 |                 |                 | the regularized |
|                 |                 |                 | Newton method   |
|                 |                 |                 | is used. (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| pen_theta_max_f | OT_REAL         | 10000           | Determines      |
| act             |                 |                 | upper bound for |
|                 |                 |                 | constraint      |
|                 |                 |                 | violation in    |
|                 |                 |                 | the filter.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| penalty_init_ma | OT_REAL         | 100000          | Maximal value   |
| x               |                 |                 | for the intial  |
|                 |                 |                 | penalty         |
|                 |                 |                 | parameter (for  |
|                 |                 |                 | Chen-Goldfarb   |
|                 |                 |                 | line search).   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| penalty_init_mi | OT_REAL         | 1               | Minimal value   |
| n               |                 |                 | for the intial  |
|                 |                 |                 | penalty         |
|                 |                 |                 | parameter for   |
|                 |                 |                 | line search(for |
|                 |                 |                 | Chen-Goldfarb   |
|                 |                 |                 | line search).   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| penalty_max     | OT_REAL         | 1.000e+30       | Maximal value   |
|                 |                 |                 | for the penalty |
|                 |                 |                 | parameter (for  |
|                 |                 |                 | Chen-Goldfarb   |
|                 |                 |                 | line search).   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| penalty_update_ | OT_REAL         | 10              | LIFENG WRITES   |
| compl_tol       |                 |                 | THIS. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| penalty_update_ | OT_REAL         | 0.000           | Threshold for   |
| infeasibility_t |                 |                 | infeasibility   |
| ol              |                 |                 | in penalty      |
|                 |                 |                 | parameter       |
|                 |                 |                 | update test.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| perturb_always_ | OT_STRING       | no              | Active          |
| cd              |                 |                 | permanent       |
|                 |                 |                 | perturbation of |
|                 |                 |                 | constraint      |
|                 |                 |                 | linearization.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| perturb_dec_fac | OT_REAL         | 0.333           | Decrease factor |
| t               |                 |                 | for x-s         |
|                 |                 |                 | perturbation.   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| perturb_inc_fac | OT_REAL         | 8               | Increase factor |
| t               |                 |                 | for x-s         |
|                 |                 |                 | perturbation.   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| perturb_inc_fac | OT_REAL         | 100             | Increase factor |
| t_first         |                 |                 | for x-s         |
|                 |                 |                 | perturbation    |
|                 |                 |                 | for very first  |
|                 |                 |                 | perturbation.   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| piecewisepenalt | OT_REAL         | 0.000           | LIFENG WRITES   |
| y_gamma_infeasi |                 |                 | THIS. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| piecewisepenalt | OT_REAL         | 0.000           | LIFENG WRITES   |
| y_gamma_obj     |                 |                 | THIS. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| point_perturbat | OT_REAL         | 10              | Maximal         |
| ion_radius      |                 |                 | perturbation of |
|                 |                 |                 | an evaluation   |
|                 |                 |                 | point. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| print_info_stri | OT_STRING       | no              | Enables         |
| ng              |                 |                 | printing of     |
|                 |                 |                 | additional info |
|                 |                 |                 | string at end   |
|                 |                 |                 | of iteration    |
|                 |                 |                 | output. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| print_level     | OT_INTEGER      | 5               | Output          |
|                 |                 |                 | verbosity       |
|                 |                 |                 | level. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| print_options_d | OT_STRING       | no              | Switch to print |
| ocumentation    |                 |                 | all algorithmic |
|                 |                 |                 | options. (see   |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| print_options_l | OT_STRING       | no              | Undocumented    |
| atex_mode       |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| print_time      | OT_BOOLEAN      | True            | print           |
|                 |                 |                 | information     |
|                 |                 |                 | about execution |
|                 |                 |                 | time            |
+-----------------+-----------------+-----------------+-----------------+
| print_timing_st | OT_STRING       | no              | Switch to print |
| atistics        |                 |                 | timing          |
|                 |                 |                 | statistics.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| print_user_opti | OT_STRING       | no              | Print all       |
| ons             |                 |                 | options set by  |
|                 |                 |                 | the user. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| quality_functio | OT_STRING       | none            | The balancing   |
| n_balancing_ter |                 |                 | term included   |
| m               |                 |                 | in the quality  |
|                 |                 |                 | function for    |
|                 |                 |                 | centrality.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| quality_functio | OT_STRING       | none            | The penalty     |
| n_centrality    |                 |                 | term for        |
|                 |                 |                 | centrality that |
|                 |                 |                 | is included in  |
|                 |                 |                 | quality         |
|                 |                 |                 | function. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| quality_functio | OT_INTEGER      | 8               | Maximum number  |
| n_max_section_s |                 |                 | of search steps |
| teps            |                 |                 | during direct   |
|                 |                 |                 | search          |
|                 |                 |                 | procedure       |
|                 |                 |                 | determining the |
|                 |                 |                 | optimal         |
|                 |                 |                 | centering       |
|                 |                 |                 | parameter. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| quality_functio | OT_STRING       | 2-norm-squared  | Norm used for   |
| n_norm_type     |                 |                 | components of   |
|                 |                 |                 | the quality     |
|                 |                 |                 | function. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| quality_functio | OT_REAL         | 0               | Tolerance for   |
| n_section_qf_to |                 |                 | the golden      |
| l               |                 |                 | section search  |
|                 |                 |                 | procedure       |
|                 |                 |                 | determining the |
|                 |                 |                 | optimal         |
|                 |                 |                 | centering       |
|                 |                 |                 | parameter (in   |
|                 |                 |                 | the function    |
|                 |                 |                 | value space).   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| quality_functio | OT_REAL         | 0.010           | Tolerance for   |
| n_section_sigma |                 |                 | the section     |
| _tol            |                 |                 | search          |
|                 |                 |                 | procedure       |
|                 |                 |                 | determining the |
|                 |                 |                 | optimal         |
|                 |                 |                 | centering       |
|                 |                 |                 | parameter (in   |
|                 |                 |                 | sigma space).   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| recalc_y        | OT_STRING       | no              | Tells the       |
|                 |                 |                 | algorithm to    |
|                 |                 |                 | recalculate the |
|                 |                 |                 | equality and    |
|                 |                 |                 | inequality      |
|                 |                 |                 | multipliers as  |
|                 |                 |                 | least square    |
|                 |                 |                 | estimates. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| recalc_y_feas_t | OT_REAL         | 0.000           | Feasibility     |
| ol              |                 |                 | threshold for   |
|                 |                 |                 | recomputation   |
|                 |                 |                 | of multipliers. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| replace_bounds  | OT_STRING       | no              | Indicates if    |
|                 |                 |                 | all variable    |
|                 |                 |                 | bounds should   |
|                 |                 |                 | be replaced by  |
|                 |                 |                 | inequality      |
|                 |                 |                 | constraints     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| required_infeas | OT_REAL         | 0.900           | Required        |
| ibility_reducti |                 |                 | reduction of    |
| on              |                 |                 | infeasibility   |
|                 |                 |                 | before leaving  |
|                 |                 |                 | restoration     |
|                 |                 |                 | phase. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| residual_improv | OT_REAL         | 1.000           | Minimal         |
| ement_factor    |                 |                 | required        |
|                 |                 |                 | reduction of    |
|                 |                 |                 | residual test   |
|                 |                 |                 | ratio in        |
|                 |                 |                 | iterative       |
|                 |                 |                 | refinement.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| residual_ratio_ | OT_REAL         | 0.000           | Iterative       |
| max             |                 |                 | refinement      |
|                 |                 |                 | tolerance (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| residual_ratio_ | OT_REAL         | 0.000           | Threshold for   |
| singular        |                 |                 | declaring       |
|                 |                 |                 | linear system   |
|                 |                 |                 | singular after  |
|                 |                 |                 | failed          |
|                 |                 |                 | iterative       |
|                 |                 |                 | refinement.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| resto_failure_f | OT_REAL         | 0               | Threshold for   |
| easibility_thre |                 |                 | primal          |
| shold           |                 |                 | infeasibility   |
|                 |                 |                 | to declare      |
|                 |                 |                 | failure of      |
|                 |                 |                 | restoration     |
|                 |                 |                 | phase. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| resto_penalty_p | OT_REAL         | 1000            | Penalty         |
| arameter        |                 |                 | parameter in    |
|                 |                 |                 | the restoration |
|                 |                 |                 | phase objective |
|                 |                 |                 | function. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| resto_proximity | OT_REAL         | 1               | Weighting       |
| _weight         |                 |                 | factor for the  |
|                 |                 |                 | proximity term  |
|                 |                 |                 | in restoration  |
|                 |                 |                 | phase           |
|                 |                 |                 | objective. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| rho             | OT_REAL         | 0.100           | Value in        |
|                 |                 |                 | penalty         |
|                 |                 |                 | parameter       |
|                 |                 |                 | update formula. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| s_max           | OT_REAL         | 100             | Scaling         |
|                 |                 |                 | threshold for   |
|                 |                 |                 | the NLP error.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| s_phi           | OT_REAL         | 2.300           | Exponent for    |
|                 |                 |                 | linear barrier  |
|                 |                 |                 | function model  |
|                 |                 |                 | in the          |
|                 |                 |                 | switching rule. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| s_theta         | OT_REAL         | 1.100           | Exponent for    |
|                 |                 |                 | current         |
|                 |                 |                 | constraint      |
|                 |                 |                 | violation in    |
|                 |                 |                 | the switching   |
|                 |                 |                 | rule. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| sb              | OT_STRING       | no              | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| sigma_max       | OT_REAL         | 100             | Maximum value   |
|                 |                 |                 | of the          |
|                 |                 |                 | centering       |
|                 |                 |                 | parameter. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| sigma_min       | OT_REAL         | 0.000           | Minimum value   |
|                 |                 |                 | of the          |
|                 |                 |                 | centering       |
|                 |                 |                 | parameter. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| skip_corr_if_ne | OT_STRING       | yes             | Skip the        |
| g_curv          |                 |                 | corrector step  |
|                 |                 |                 | in negative     |
|                 |                 |                 | curvature       |
|                 |                 |                 | iteration       |
|                 |                 |                 | (unsupported!). |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| skip_corr_in_mo | OT_STRING       | yes             | Skip the        |
| notone_mode     |                 |                 | corrector step  |
|                 |                 |                 | during monotone |
|                 |                 |                 | barrier         |
|                 |                 |                 | parameter mode  |
|                 |                 |                 | (unsupported!). |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| skip_finalize_s | OT_STRING       | no              | Indicates if    |
| olution_call    |                 |                 | call to NLP::Fi |
|                 |                 |                 | nalizeSolution  |
|                 |                 |                 | after           |
|                 |                 |                 | optimization    |
|                 |                 |                 | should be       |
|                 |                 |                 | suppressed (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| slack_bound_fra | OT_REAL         | 0.010           | Desired minimum |
| c               |                 |                 | relative        |
|                 |                 |                 | distance from   |
|                 |                 |                 | the initial     |
|                 |                 |                 | slack to bound. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| slack_bound_pus | OT_REAL         | 0.010           | Desired minimum |
| h               |                 |                 | absolute        |
|                 |                 |                 | distance from   |
|                 |                 |                 | the initial     |
|                 |                 |                 | slack to bound. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| slack_move      | OT_REAL         | 0.000           | Correction size |
|                 |                 |                 | for very small  |
|                 |                 |                 | slacks. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| soft_resto_pder | OT_REAL         | 1.000           | Required        |
| ror_reduction_f |                 |                 | reduction in    |
| actor           |                 |                 | primal-dual     |
|                 |                 |                 | error in the    |
|                 |                 |                 | soft            |
|                 |                 |                 | restoration     |
|                 |                 |                 | phase. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| start_with_rest | OT_STRING       | no              | Tells algorithm |
| o               |                 |                 | to switch to    |
|                 |                 |                 | restoration     |
|                 |                 |                 | phase in first  |
|                 |                 |                 | iteration. (see |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| suppress_all_ou | OT_STRING       | no              | Undocumented    |
| tput            |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| tau_min         | OT_REAL         | 0.990           | Lower bound on  |
|                 |                 |                 | fraction-to-    |
|                 |                 |                 | the-boundary    |
|                 |                 |                 | parameter tau.  |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| theta_max_fact  | OT_REAL         | 10000           | Determines      |
|                 |                 |                 | upper bound for |
|                 |                 |                 | constraint      |
|                 |                 |                 | violation in    |
|                 |                 |                 | the filter.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| theta_min       | OT_REAL         | 0.000           | LIFENG WRITES   |
|                 |                 |                 | THIS. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| theta_min_fact  | OT_REAL         | 0.000           | Determines      |
|                 |                 |                 | constraint      |
|                 |                 |                 | violation       |
|                 |                 |                 | threshold in    |
|                 |                 |                 | the switching   |
|                 |                 |                 | rule. (see      |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| tiny_step_tol   | OT_REAL         | 0.000           | Tolerance for   |
|                 |                 |                 | detecting       |
|                 |                 |                 | numerically     |
|                 |                 |                 | insignificant   |
|                 |                 |                 | steps. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| tiny_step_y_tol | OT_REAL         | 0.010           | Tolerance for   |
|                 |                 |                 | quitting        |
|                 |                 |                 | because of      |
|                 |                 |                 | numerically     |
|                 |                 |                 | insignificant   |
|                 |                 |                 | steps. (see     |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| tol             | OT_REAL         | 0.000           | Desired         |
|                 |                 |                 | convergence     |
|                 |                 |                 | tolerance       |
|                 |                 |                 | (relative).     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| var_integer_md  | OT_DICT         | None            | Integer         |
|                 |                 |                 | metadata (a     |
|                 |                 |                 | dictionary with |
|                 |                 |                 | lists of        |
|                 |                 |                 | integers) about |
|                 |                 |                 | variables to be |
|                 |                 |                 | passed to IPOPT |
+-----------------+-----------------+-----------------+-----------------+
| var_numeric_md  | OT_DICT         | None            | Numeric         |
|                 |                 |                 | metadata (a     |
|                 |                 |                 | dictionary with |
|                 |                 |                 | lists of reals) |
|                 |                 |                 | about variables |
|                 |                 |                 | to be passed to |
|                 |                 |                 | IPOPT           |
+-----------------+-----------------+-----------------+-----------------+
| var_string_md   | OT_DICT         | None            | String metadata |
|                 |                 |                 | (a dictionary   |
|                 |                 |                 | with lists of   |
|                 |                 |                 | strings) about  |
|                 |                 |                 | variables to be |
|                 |                 |                 | passed to IPOPT |
+-----------------+-----------------+-----------------+-----------------+
| vartheta        | OT_REAL         | 0.500           | a parameter     |
|                 |                 |                 | used to check   |
|                 |                 |                 | if the fast     |
|                 |                 |                 | direction can   |
|                 |                 |                 | be used asthe   |
|                 |                 |                 | line search     |
|                 |                 |                 | direction (for  |
|                 |                 |                 | Chen-Goldfarb   |
|                 |                 |                 | line search).   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| warm_start_boun | OT_REAL         | 0.001           | same as         |
| d_frac          |                 |                 | bound_frac for  |
|                 |                 |                 | the regular     |
|                 |                 |                 | initializer.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| warm_start_boun | OT_REAL         | 0.001           | same as         |
| d_push          |                 |                 | bound_push for  |
|                 |                 |                 | the regular     |
|                 |                 |                 | initializer.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| warm_start_enti | OT_STRING       | no              | Tells algorithm |
| re_iterate      |                 |                 | whether to use  |
|                 |                 |                 | the GetWarmStar |
|                 |                 |                 | tIterate method |
|                 |                 |                 | in the NLP.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| warm_start_init | OT_STRING       | no              | Warm-start for  |
| _point          |                 |                 | initial point   |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| warm_start_mult | OT_REAL         | 0.001           | same as         |
| _bound_push     |                 |                 | mult_bound_push |
|                 |                 |                 | for the regular |
|                 |                 |                 | initializer.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| warm_start_mult | OT_REAL         | 1000000         | Maximum initial |
| _init_max       |                 |                 | value for the   |
|                 |                 |                 | equality        |
|                 |                 |                 | multipliers.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| warm_start_same | OT_STRING       | no              | Indicates       |
| _structure      |                 |                 | whether a       |
|                 |                 |                 | problem with a  |
|                 |                 |                 | structure       |
|                 |                 |                 | identical to    |
|                 |                 |                 | the previous    |
|                 |                 |                 | one is to be    |
|                 |                 |                 | solved. (see    |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| warm_start_slac | OT_REAL         | 0.001           | same as slack_b |
| k_bound_frac    |                 |                 | ound_frac for   |
|                 |                 |                 | the regular     |
|                 |                 |                 | initializer.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| warm_start_slac | OT_REAL         | 0.001           | same as slack_b |
| k_bound_push    |                 |                 | ound_push for   |
|                 |                 |                 | the regular     |
|                 |                 |                 | initializer.    |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| warm_start_targ | OT_REAL         | 0               | Unsupported!    |
| et_mu           |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| watchdog_shorte | OT_INTEGER      | 10              | Number of       |
| ned_iter_trigge |                 |                 | shortened       |
| r               |                 |                 | iterations that |
|                 |                 |                 | trigger the     |
|                 |                 |                 | watchdog. (see  |
|                 |                 |                 | IPOPT           |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| watchdog_trial_ | OT_INTEGER      | 3               | Maximum number  |
| iter_max        |                 |                 | of watchdog     |
|                 |                 |                 | iterations.     |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+
| wsmp_iterative  | OT_STRING       | no              | Switches to     |
|                 |                 |                 | iterative       |
|                 |                 |                 | solver in WSMP. |
|                 |                 |                 | (see IPOPT      |
|                 |                 |                 | documentation)  |
+-----------------+-----------------+-----------------+-----------------+

>List of available monitors

+-------------+
|     Id      |
+=============+
| eval_f      |
+-------------+
| eval_g      |
+-------------+
| eval_grad_f |
+-------------+
| eval_h      |
+-------------+
| eval_jac_g  |
+-------------+

>List of available stats

+--------------------+
|         Id         |
+====================+
| con_integer_md     |
+--------------------+
| con_numeric_md     |
+--------------------+
| con_string_md      |
+--------------------+
| iter_count         |
+--------------------+
| iteration          |
+--------------------+
| iterations         |
+--------------------+
| n_eval_callback    |
+--------------------+
| n_eval_f           |
+--------------------+
| n_eval_g           |
+--------------------+
| n_eval_grad_f      |
+--------------------+
| n_eval_h           |
+--------------------+
| n_eval_jac_g       |
+--------------------+
| return_status      |
+--------------------+
| t_callback_fun     |
+--------------------+
| t_callback_prepare |
+--------------------+
| t_eval_f           |
+--------------------+
| t_eval_g           |
+--------------------+
| t_eval_grad_f      |
+--------------------+
| t_eval_h           |
+--------------------+
| t_eval_jac_g       |
+--------------------+
| t_mainloop         |
+--------------------+
| var_integer_md     |
+--------------------+
| var_numeric_md     |
+--------------------+
| var_string_md      |
+--------------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

knitro
------



KNITRO interface

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| BarRule         | OT_INTEGER      | 0               | Barrier Rule    |
+-----------------+-----------------+-----------------+-----------------+
| Debug           | OT_INTEGER      | 0               | Debug level     |
+-----------------+-----------------+-----------------+-----------------+
| Delta           | OT_REAL         | 1               | Initial region  |
|                 |                 |                 | scaling factor  |
+-----------------+-----------------+-----------------+-----------------+
| FeasModeTol     | OT_REAL         | 0.000           | Feasible mode   |
|                 |                 |                 | tolerance       |
+-----------------+-----------------+-----------------+-----------------+
| FeasTol         | OT_REAL         | 0.000           | Feasible        |
|                 |                 |                 | tolerance       |
+-----------------+-----------------+-----------------+-----------------+
| FeasTolAbs      | OT_REAL         | 0               | Absolute        |
|                 |                 |                 | feasible        |
|                 |                 |                 | tolerance       |
+-----------------+-----------------+-----------------+-----------------+
| Feasible        | OT_BOOLEAN      | 1               | Allow           |
|                 |                 |                 | infeasible      |
|                 |                 |                 | iterations      |
+-----------------+-----------------+-----------------+-----------------+
| GradOpt         | OT_INTEGER      | 1               | Gradient        |
|                 |                 |                 | calculation     |
|                 |                 |                 | method          |
+-----------------+-----------------+-----------------+-----------------+
| HessOpt         | OT_INTEGER      | 1               | Hessian         |
|                 |                 |                 | calculation     |
|                 |                 |                 | method          |
+-----------------+-----------------+-----------------+-----------------+
| HonorBnds       | OT_BOOLEAN      | 0               | Enforce bounds  |
+-----------------+-----------------+-----------------+-----------------+
| InitPt          | OT_BOOLEAN      | 0               | Use initial     |
|                 |                 |                 | point strategy  |
+-----------------+-----------------+-----------------+-----------------+
| LmSize          | OT_INTEGER      | 10              | Memory pairsize |
|                 |                 |                 | limit           |
+-----------------+-----------------+-----------------+-----------------+
| LpSolver        | OT_BOOLEAN      | 0               | Use LpSolver    |
+-----------------+-----------------+-----------------+-----------------+
| MaxCgIt         | OT_INTEGER      | 0               | Maximum         |
|                 |                 |                 | conjugate       |
|                 |                 |                 | gradient        |
|                 |                 |                 | iterations      |
+-----------------+-----------------+-----------------+-----------------+
| MaxIt           | OT_INTEGER      | 10000           | Iteration limit |
+-----------------+-----------------+-----------------+-----------------+
| Mu              | OT_REAL         | 0.100           | Initial barrier |
|                 |                 |                 | parameter       |
+-----------------+-----------------+-----------------+-----------------+
| Multistart      | OT_BOOLEAN      | 0               | Use multistart  |
+-----------------+-----------------+-----------------+-----------------+
| NewPoint        | OT_BOOLEAN      | 0               | Select new-     |
|                 |                 |                 | point feature   |
+-----------------+-----------------+-----------------+-----------------+
| ObjRange        | OT_REAL         | 1.000e+20       | Maximum         |
|                 |                 |                 | objective value |
+-----------------+-----------------+-----------------+-----------------+
| OptTol          | OT_REAL         | 0.000           | Relative        |
|                 |                 |                 | optimality      |
|                 |                 |                 | tolerance       |
+-----------------+-----------------+-----------------+-----------------+
| OptTolAbs       | OT_REAL         | 0               | Absolute        |
|                 |                 |                 | optimality      |
|                 |                 |                 | tolerance       |
+-----------------+-----------------+-----------------+-----------------+
| OutLev          | OT_INTEGER      | 2               | Log output      |
|                 |                 |                 | level           |
+-----------------+-----------------+-----------------+-----------------+
| Pivot           | OT_REAL         | 0.000           | Initial pivot   |
|                 |                 |                 | threshold       |
+-----------------+-----------------+-----------------+-----------------+
| Scale           | OT_BOOLEAN      | 1               | Perform scaling |
+-----------------+-----------------+-----------------+-----------------+
| ShiftInit       | OT_BOOLEAN      | 1               | Interior-point  |
|                 |                 |                 | shifting        |
|                 |                 |                 | initial point   |
+-----------------+-----------------+-----------------+-----------------+
| Soc             | OT_INTEGER      | 1               | Second order    |
|                 |                 |                 | correction      |
+-----------------+-----------------+-----------------+-----------------+
| XTol            | OT_REAL         | 0.000           | Relative        |
|                 |                 |                 | solution change |
|                 |                 |                 | tolerance       |
+-----------------+-----------------+-----------------+-----------------+
| contype         | OT_INTEGERVECTO |                 |                 |
|                 | R               |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

>List of available monitors

+-------------+
|     Id      |
+=============+
| eval_f      |
+-------------+
| eval_g      |
+-------------+
| eval_grad_f |
+-------------+
| eval_h      |
+-------------+
| eval_jac_g  |
+-------------+

>List of available stats

+---------------+
|      Id       |
+===============+
| return_status |
+---------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

snopt
-----



SNOPT interface

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| Backup basis    | OT_INTEGER      | None            | 0 * output      |
| file            |                 |                 | extra basis map |
+-----------------+-----------------+-----------------+-----------------+
| Central         | OT_REAL         | None            | 6.7e-5 * (      |
| difference      |                 |                 | Function        |
| interval        |                 |                 | precision)^1/3  |
+-----------------+-----------------+-----------------+-----------------+
| Check frequency | OT_INTEGER      | None            | 60 * test row   |
|                 |                 |                 | residuals kAx - |
|                 |                 |                 | sk              |
+-----------------+-----------------+-----------------+-----------------+
| Crash option    | OT_INTEGER      | None            | 3 * first basis |
|                 |                 |                 | is essentially  |
|                 |                 |                 | triangular      |
+-----------------+-----------------+-----------------+-----------------+
| Crash tolerance | OT_REAL         | None            | 0.100           |
+-----------------+-----------------+-----------------+-----------------+
| Debug level     | OT_INTEGER      | None            | 0 * for         |
|                 |                 |                 | developers      |
+-----------------+-----------------+-----------------+-----------------+
| Derivative      | OT_INTEGER      | None            | 3               |
| level           |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| Difference      | OT_REAL         | None            | 5.5e-7 * (      |
| interval        |                 |                 | Function        |
|                 |                 |                 | precision)^1/2  |
+-----------------+-----------------+-----------------+-----------------+
| Dump file       | OT_INTEGER      | None            | 0 * output Load |
|                 |                 |                 | data            |
+-----------------+-----------------+-----------------+-----------------+
| Elastic weight  | OT_REAL         | None            | 1.0e+4 * used   |
|                 |                 |                 | only during     |
|                 |                 |                 | elastic mode    |
+-----------------+-----------------+-----------------+-----------------+
| Expand          | OT_INTEGER      | None            | 10000 * for     |
| frequency       |                 |                 | anti-cycling    |
|                 |                 |                 | procedure       |
+-----------------+-----------------+-----------------+-----------------+
| Factorization   | OT_INTEGER      | None            | 50 * 100 for    |
| frequency       |                 |                 | LPs             |
+-----------------+-----------------+-----------------+-----------------+
| Function        | OT_REAL         | None            | 3.0e-13 * e^0.8 |
| precision       |                 |                 | (almost full    |
|                 |                 |                 | accuracy)       |
+-----------------+-----------------+-----------------+-----------------+
| Hessian         | OT_STRING       | None            | full memory *   |
|                 |                 |                 | default if n1   |
|                 |                 |                 | 75  limited     |
|                 |                 |                 | memory *        |
|                 |                 |                 | default if n1 > |
|                 |                 |                 | 75              |
+-----------------+-----------------+-----------------+-----------------+
| Hessian flush   | OT_INTEGER      | None            | 999999 * no     |
|                 |                 |                 | flushing        |
+-----------------+-----------------+-----------------+-----------------+
| Hessian         | OT_INTEGER      | None            | 999999 * for    |
| frequency       |                 |                 | full Hessian    |
|                 |                 |                 | (never reset)   |
+-----------------+-----------------+-----------------+-----------------+
| Hessian updates | OT_INTEGER      | None            | 10 * for        |
|                 |                 |                 | limited memory  |
|                 |                 |                 | Hessian         |
+-----------------+-----------------+-----------------+-----------------+
| Insert file     | OT_INTEGER      | None            | 0 * input in    |
|                 |                 |                 | industry format |
+-----------------+-----------------+-----------------+-----------------+
| Iterations      | OT_INTEGER      | None            | 10000 * or 20m  |
| limit           |                 |                 | if that is more |
+-----------------+-----------------+-----------------+-----------------+
| LU              | OT_STRING       | None            | LU partial      |
|                 |                 |                 | pivoting *      |
|                 |                 |                 | default         |
|                 |                 |                 | threshold       |
|                 |                 |                 | pivoting        |
|                 |                 |                 | strategy  LU    |
|                 |                 |                 | rook pivoting * |
|                 |                 |                 | threshold rook  |
|                 |                 |                 | pivoting  LU    |
|                 |                 |                 | complete        |
|                 |                 |                 | pivoting *      |
|                 |                 |                 | threshold       |
|                 |                 |                 | complete        |
|                 |                 |                 | pivoting        |
+-----------------+-----------------+-----------------+-----------------+
| LU factor       | OT_REAL         | None            | 3.99 * for NP   |
| tolerance       |                 |                 | (100.0 for LP)  |
+-----------------+-----------------+-----------------+-----------------+
| LU singularity  | OT_REAL         | None            | 0.000           |
| tolerance       |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| LU update       | OT_REAL         | None            | 3.99 * for NP ( |
| tolerance       |                 |                 | 10.0 for LP)    |
+-----------------+-----------------+-----------------+-----------------+
| Linesearch      | OT_REAL         | None            | 0.9 * smaller   |
| tolerance       |                 |                 | for more        |
|                 |                 |                 | accurate search |
+-----------------+-----------------+-----------------+-----------------+
| Load file       | OT_INTEGER      | None            | 0 * input names |
|                 |                 |                 | and values      |
+-----------------+-----------------+-----------------+-----------------+
| Major           | OT_REAL         | None            | 1.0e-6 * target |
| feasibility     |                 |                 | nonlinear       |
| tolerance       |                 |                 | constraint      |
|                 |                 |                 | violation       |
+-----------------+-----------------+-----------------+-----------------+
| Major           | OT_INTEGER      | None            | 1000 * or m if  |
| iterations      |                 |                 | that is more    |
| limit           |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| Major           | OT_REAL         | None            | 1.0e-6 * target |
| optimality      |                 |                 | complementarity |
| tolerance       |                 |                 | gap             |
+-----------------+-----------------+-----------------+-----------------+
| Major print     | OT_INTEGER      | None            | 1 * 1-line      |
| level           |                 |                 | major iteration |
|                 |                 |                 | log             |
+-----------------+-----------------+-----------------+-----------------+
| Major step      | OT_REAL         | None            | 2               |
| limit           |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| Minor           | OT_REAL         | None            | 1.0e-6 * for    |
| feasibility     |                 |                 | satisfying the  |
| tolerance       |                 |                 | QP bounds       |
+-----------------+-----------------+-----------------+-----------------+
| Minor           | OT_INTEGER      | None            | 500 * or 3m if  |
| iterations      |                 |                 | that is more    |
| limit           |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| Minor print     | OT_INTEGER      | None            | 1 * 1-line      |
| level           |                 |                 | minor iteration |
|                 |                 |                 | log             |
+-----------------+-----------------+-----------------+-----------------+
| New basis file  | OT_INTEGER      | None            | 0 * output      |
|                 |                 |                 | basis map       |
+-----------------+-----------------+-----------------+-----------------+
| New superbasics | OT_INTEGER      | None            | 99 * controls   |
| limit           |                 |                 | early           |
|                 |                 |                 | termination of  |
|                 |                 |                 | QPs             |
+-----------------+-----------------+-----------------+-----------------+
| Old basis file  | OT_INTEGER      | None            | 0 * input basis |
|                 |                 |                 | map             |
+-----------------+-----------------+-----------------+-----------------+
| Partial price   | OT_INTEGER      | None            | 1 * 10 for      |
|                 |                 |                 | large LPs       |
+-----------------+-----------------+-----------------+-----------------+
| Penalty         | OT_REAL         | None            | 0.0 * initial   |
| parameter       |                 |                 | penalty         |
|                 |                 |                 | parameter       |
+-----------------+-----------------+-----------------+-----------------+
| Pivot tolerance | OT_REAL         | None            | 3.7e-11 * e^2/3 |
+-----------------+-----------------+-----------------+-----------------+
| Print frequency | OT_INTEGER      | None            | 100 * minor     |
|                 |                 |                 | iterations log  |
|                 |                 |                 | on Print file   |
+-----------------+-----------------+-----------------+-----------------+
| Proximal point  | OT_INTEGER      | None            | 1 * satisfies   |
| method          |                 |                 | linear          |
|                 |                 |                 | constraints     |
|                 |                 |                 | near x0         |
+-----------------+-----------------+-----------------+-----------------+
| Punch file      | OT_INTEGER      | None            | 0 * output      |
|                 |                 |                 | Insert data     |
+-----------------+-----------------+-----------------+-----------------+
| QPSolver        | OT_STRING       | None            | Cholesky *      |
|                 |                 |                 | default         |
+-----------------+-----------------+-----------------+-----------------+
| Reduced Hessian | OT_INTEGER      | None            | 2000 * or       |
| dimension       |                 |                 | Superbasics     |
|                 |                 |                 | limit if that   |
|                 |                 |                 | is less         |
+-----------------+-----------------+-----------------+-----------------+
| Save frequency  | OT_INTEGER      | None            | 100 * save      |
|                 |                 |                 | basis map       |
+-----------------+-----------------+-----------------+-----------------+
| Scale option    | OT_INTEGER      | None            | 1 * linear      |
|                 |                 |                 | constraints and |
|                 |                 |                 | variables       |
+-----------------+-----------------+-----------------+-----------------+
| Scale tolerance | OT_REAL         | None            | 0.900           |
+-----------------+-----------------+-----------------+-----------------+
| Solution        | OT_STRING       | None            | Yes * on the    |
|                 |                 |                 | Print file      |
+-----------------+-----------------+-----------------+-----------------+
| Solution file   | OT_INTEGER      | None            | 0 * different   |
|                 |                 |                 | from printed    |
|                 |                 |                 | solution        |
+-----------------+-----------------+-----------------+-----------------+
| Sticky          | OT_STRING       | None            | No * Yes makes  |
| parameters      |                 |                 | parameter       |
|                 |                 |                 | values persist  |
+-----------------+-----------------+-----------------+-----------------+
| Summary         | OT_INTEGER      | None            | 100 * minor     |
| frequency       |                 |                 | iterations log  |
|                 |                 |                 | on Summary file |
+-----------------+-----------------+-----------------+-----------------+
| Superbasics     | OT_INTEGER      | None            | n1 + 1 * n1 =   |
| limit           |                 |                 | number of       |
|                 |                 |                 | nonlinear       |
|                 |                 |                 | variables       |
+-----------------+-----------------+-----------------+-----------------+
| System          | OT_STRING       | None            | No * Yes prints |
| information     |                 |                 | more system     |
|                 |                 |                 | information     |
+-----------------+-----------------+-----------------+-----------------+
| Timing level    | OT_INTEGER      | None            | 3 * print cpu   |
|                 |                 |                 | times           |
+-----------------+-----------------+-----------------+-----------------+
| Unbounded       | OT_REAL         | None            | 1.000e+15       |
| objective       |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| Unbounded step  | OT_REAL         | None            | 1.000e+18       |
| size            |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| Verify level    | OT_INTEGER      | None            | 0 * cheap check |
|                 |                 |                 | on gradients    |
+-----------------+-----------------+-----------------+-----------------+
| Violation limit | OT_REAL         | None            | 10.0 * unscaled |
|                 |                 |                 | constraint      |
|                 |                 |                 | violation limit |
+-----------------+-----------------+-----------------+-----------------+
| detect_linear   | OT_BOOLEAN      | True            | Make an effort  |
|                 |                 |                 | to treat linear |
|                 |                 |                 | constraints and |
|                 |                 |                 | linear          |
|                 |                 |                 | variables       |
|                 |                 |                 | specially.      |
+-----------------+-----------------+-----------------+-----------------+
| print file      | OT_STRING       | None            | n/a             |
+-----------------+-----------------+-----------------+-----------------+
| print_time      | OT_BOOLEAN      | True            | print           |
|                 |                 |                 | information     |
|                 |                 |                 | about execution |
|                 |                 |                 | time            |
+-----------------+-----------------+-----------------+-----------------+
| specs file      | OT_STRING       | None            | n/a             |
+-----------------+-----------------+-----------------+-----------------+
| start           | OT_STRING       | Cold            |                 |
+-----------------+-----------------+-----------------+-----------------+
| summary         | OT_BOOLEAN      | True            | n/a             |
+-----------------+-----------------+-----------------+-----------------+

>List of available monitors

+-----------+
|    Id     |
+===========+
| eval_nlp  |
+-----------+
| setup_nlp |
+-----------+

>List of available stats

+----------------+
|       Id       |
+================+
| iter_count     |
+----------------+
| iterations     |
+----------------+
| n_callback_fun |
+----------------+
| n_eval_grad_f  |
+----------------+
| n_eval_jac_g   |
+----------------+
| return_status  |
+----------------+
| t_callback_fun |
+----------------+
| t_eval_grad_f  |
+----------------+
| t_eval_jac_g   |
+----------------+
| t_mainloop     |
+----------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

worhp
-----



WORHP interface

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| AcceptTolFeas   | OT_REAL         | 0.001           | Tolerance for   |
|                 |                 |                 | acceptable      |
|                 |                 |                 | feasibility     |
+-----------------+-----------------+-----------------+-----------------+
| AcceptTolOpti   | OT_REAL         | 0.001           | Tolerance for   |
|                 |                 |                 | acceptable      |
|                 |                 |                 | optimality      |
+-----------------+-----------------+-----------------+-----------------+
| AlphaMinConst   | OT_BOOLEAN      | False           | Use a constant  |
|                 |                 |                 | lower bound on  |
|                 |                 |                 | Armijo stepsize |
|                 |                 |                 | in Filter       |
+-----------------+-----------------+-----------------+-----------------+
| Ares            | OT_INTEGERVECTO | [42, 41, 42,    | Armijo recovery |
|                 | R               | 43, 44, 41, 50] | strategies.     |
|                 |                 |                 | Vector of size  |
|                 |                 |                 | 7               |
+-----------------+-----------------+-----------------+-----------------+
| ArmijoBeta      | OT_REAL         | 0.712           | Trial stepsize  |
|                 |                 |                 | decrease factor |
|                 |                 |                 | for Armijo rule |
+-----------------+-----------------+-----------------+-----------------+
| ArmijoMaxAlpha  | OT_REAL         | 1               | Initial alpha   |
|                 |                 |                 | for Armijo rule |
+-----------------+-----------------+-----------------+-----------------+
| ArmijoMinAlpha  | OT_REAL         | 0.000           | Lower bound on  |
|                 |                 |                 | alpha for       |
|                 |                 |                 | Armijo rule     |
+-----------------+-----------------+-----------------+-----------------+
| ArmijoMinAlphaR | OT_REAL         | 0.000           | Lower bound on  |
| ec              |                 |                 | alpha for       |
|                 |                 |                 | Armijo rule     |
|                 |                 |                 | during recovery |
+-----------------+-----------------+-----------------+-----------------+
| ArmijoSigma     | OT_REAL         | 0.005           | Scale factor    |
|                 |                 |                 | for linearised  |
|                 |                 |                 | descent check   |
|                 |                 |                 | in Armijo rule  |
+-----------------+-----------------+-----------------+-----------------+
| AutoQPRecovery  | OT_BOOLEAN      | True            | Enable          |
|                 |                 |                 | automatic QP    |
|                 |                 |                 | recovery        |
+-----------------+-----------------+-----------------+-----------------+
| BFGSmaxblockSiz | OT_INTEGER      | 300             | Block size      |
| e               |                 |                 | parameter used  |
|                 |                 |                 | by certain BFGS |
|                 |                 |                 | methods         |
+-----------------+-----------------+-----------------+-----------------+
| BFGSmethod      | OT_INTEGER      | 0               | Choose BFGS     |
|                 |                 |                 | method (0:      |
|                 |                 |                 | dense, 1-3:     |
|                 |                 |                 | block, 100+:    |
|                 |                 |                 | sparse)         |
+-----------------+-----------------+-----------------+-----------------+
| BFGSminblockSiz | OT_INTEGER      | 300             | Block size      |
| e               |                 |                 | parameter used  |
|                 |                 |                 | by certain BFGS |
|                 |                 |                 | methods         |
+-----------------+-----------------+-----------------+-----------------+
| BFGSrestart     | OT_INTEGER      | 50              | Restart BFGS    |
|                 |                 |                 | update after    |
|                 |                 |                 | this many       |
|                 |                 |                 | iterations      |
+-----------------+-----------------+-----------------+-----------------+
| BettsFactor     | OT_REAL         | 2.100           | Update factor   |
|                 |                 |                 | for Betts'      |
|                 |                 |                 | Hessian         |
|                 |                 |                 | regularisation  |
+-----------------+-----------------+-----------------+-----------------+
| BettsPoint      | OT_REAL         | 1               | Smallest        |
|                 |                 |                 | eigenvalue of   |
|                 |                 |                 | the regularised |
|                 |                 |                 | Hessian         |
+-----------------+-----------------+-----------------+-----------------+
| BoundTolFac     | OT_REAL         | 1000            | Factor in       |
|                 |                 |                 | determining     |
|                 |                 |                 | active          |
|                 |                 |                 | constraints by  |
|                 |                 |                 | KKT             |
+-----------------+-----------------+-----------------+-----------------+
| CheckFJ         | OT_REAL         | 1.000e+12       | Upper bound     |
|                 |                 |                 | used by Fritz-  |
|                 |                 |                 | John heuristic  |
+-----------------+-----------------+-----------------+-----------------+
| CheckStructureD | OT_BOOLEAN      | True            | Enable          |
| F               |                 |                 | structural      |
|                 |                 |                 | checking of DF  |
+-----------------+-----------------+-----------------+-----------------+
| CheckStructureD | OT_BOOLEAN      | True            | Enable          |
| G               |                 |                 | structural      |
|                 |                 |                 | checking of DG  |
+-----------------+-----------------+-----------------+-----------------+
| CheckStructureH | OT_BOOLEAN      | True            | Enable          |
| M               |                 |                 | structural      |
|                 |                 |                 | checking of HM  |
+-----------------+-----------------+-----------------+-----------------+
| CorStepBettsSum | OT_REAL         | 0.500           | (experimental)  |
+-----------------+-----------------+-----------------+-----------------+
| CorStepConStop  | OT_REAL         | 0.000           | (experimental)  |
+-----------------+-----------------+-----------------+-----------------+
| CorStepConvio   | OT_REAL         | 1               | (experimental)  |
+-----------------+-----------------+-----------------+-----------------+
| CorStepMaxIter  | OT_INTEGER      | 50              | (experimental)  |
+-----------------+-----------------+-----------------+-----------------+
| CorStepMethod   | OT_INTEGER      | 0               | (experimental)  |
+-----------------+-----------------+-----------------+-----------------+
| CorStepMode     | OT_INTEGER      | 1               | (experimental)  |
+-----------------+-----------------+-----------------+-----------------+
| CorStepPFactor  | OT_REAL         | 1               | (experimental)  |
+-----------------+-----------------+-----------------+-----------------+
| CorStepPMax     | OT_REAL         | 1000000         | (experimental)  |
+-----------------+-----------------+-----------------+-----------------+
| CorStepRecovery | OT_BOOLEAN      | False           | (experimental)  |
| DX              |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| CurvBCond       | OT_REAL         | 0.020           | Block BFGS      |
|                 |                 |                 | curvature       |
|                 |                 |                 | condition bound |
+-----------------+-----------------+-----------------+-----------------+
| CurvBFac        | OT_REAL         | 0.300           | Block BFGS      |
|                 |                 |                 | curvature       |
|                 |                 |                 | condition       |
|                 |                 |                 | regularisation  |
|                 |                 |                 | factor          |
+-----------------+-----------------+-----------------+-----------------+
| CurvCond        | OT_REAL         | 0.020           | BFGS Curvature  |
|                 |                 |                 | condition bound |
+-----------------+-----------------+-----------------+-----------------+
| CurvFac         | OT_REAL         | 0.300           | BFGS curvature  |
|                 |                 |                 | condition       |
|                 |                 |                 | regularisation  |
|                 |                 |                 | factor          |
+-----------------+-----------------+-----------------+-----------------+
| DebugMarker05   | OT_INTEGER      | 42              | Debug marker.   |
|                 |                 |                 | Used to find    |
|                 |                 |                 | memory alignmen |
|                 |                 |                 | t/padding       |
|                 |                 |                 | issues          |
+-----------------+-----------------+-----------------+-----------------+
| DebugMarker06   | OT_INTEGER      | 42              | Debug marker.   |
|                 |                 |                 | Used to find    |
|                 |                 |                 | memory alignmen |
|                 |                 |                 | t/padding       |
|                 |                 |                 | issues          |
+-----------------+-----------------+-----------------+-----------------+
| FGtogether      | OT_BOOLEAN      | False           | F and G cannot  |
|                 |                 |                 | be evaluated    |
|                 |                 |                 | separately      |
+-----------------+-----------------+-----------------+-----------------+
| FJandND         | OT_BOOLEAN      | False           | Enable Fritz-   |
|                 |                 |                 | John and non-   |
|                 |                 |                 | differentiable  |
|                 |                 |                 | check           |
|                 |                 |                 | heuristics      |
+-----------------+-----------------+-----------------+-----------------+
| FeasibleDual    | OT_BOOLEAN      | False           | Activate dual   |
|                 |                 |                 | feasibility     |
|                 |                 |                 | mode            |
+-----------------+-----------------+-----------------+-----------------+
| FeasibleInit    | OT_BOOLEAN      | False           | Activate        |
|                 |                 |                 | initial         |
|                 |                 |                 | feasibility     |
|                 |                 |                 | mode            |
+-----------------+-----------------+-----------------+-----------------+
| FeasibleInitTol | OT_REAL         | 0.001           | Feasibility     |
|                 |                 |                 | tolerance for   |
|                 |                 |                 | no-objective    |
|                 |                 |                 | feasible mode   |
+-----------------+-----------------+-----------------+-----------------+
| FeasibleOnly    | OT_BOOLEAN      | False           | Activate        |
|                 |                 |                 | feasible-only   |
|                 |                 |                 | mode            |
+-----------------+-----------------+-----------------+-----------------+
| FidifEps        | OT_REAL         | 0.000           | Finite          |
|                 |                 |                 | difference      |
|                 |                 |                 | perturbation    |
+-----------------+-----------------+-----------------+-----------------+
| FidifHM         | OT_BOOLEAN      | False           | Approximate     |
|                 |                 |                 | Hessian by      |
|                 |                 |                 | finite          |
|                 |                 |                 | differences     |
|                 |                 |                 | (otherwise      |
|                 |                 |                 | BFGS)           |
+-----------------+-----------------+-----------------+-----------------+
| FilterBisecAlph | OT_BOOLEAN      | True            | Filter          |
| a               |                 |                 | heuristic to    |
|                 |                 |                 | save Armijo     |
|                 |                 |                 | iterations      |
+-----------------+-----------------+-----------------+-----------------+
| FilterGammaCV   | OT_REAL         | 0.000           | Constraint      |
|                 |                 |                 | violation       |
|                 |                 |                 | decrease factor |
|                 |                 |                 | in Filter       |
|                 |                 |                 | acceptance      |
|                 |                 |                 | check           |
+-----------------+-----------------+-----------------+-----------------+
| FilterGammaF    | OT_REAL         | 0.000           | Objective       |
|                 |                 |                 | decrease factor |
|                 |                 |                 | in Filter       |
|                 |                 |                 | acceptance      |
|                 |                 |                 | check           |
+-----------------+-----------------+-----------------+-----------------+
| FilterIntersecA | OT_BOOLEAN      | True            | Filter          |
| lpha            |                 |                 | heuristic to    |
|                 |                 |                 | save Armijo     |
|                 |                 |                 | iterations      |
+-----------------+-----------------+-----------------+-----------------+
| FirstDifCentral | OT_BOOLEAN      | True            | Use central     |
|                 |                 |                 | finite          |
|                 |                 |                 | difference      |
|                 |                 |                 | quotient for    |
|                 |                 |                 | first           |
|                 |                 |                 | derivatives     |
+-----------------+-----------------+-----------------+-----------------+
| FocusOnFeas     | OT_BOOLEAN      | True            | Enable Focus-   |
|                 |                 |                 | on-Feasibility  |
|                 |                 |                 | mode            |
+-----------------+-----------------+-----------------+-----------------+
| FocusOnFeasFact | OT_REAL         | 1.360           | Factor in       |
| or              |                 |                 | Focus-on-       |
|                 |                 |                 | Feasibility     |
|                 |                 |                 | mode            |
+-----------------+-----------------+-----------------+-----------------+
| GammaAlpha      | OT_REAL         | 0.050           | Safety factor   |
|                 |                 |                 | for alphamin    |
|                 |                 |                 | calculation by  |
|                 |                 |                 | Filter          |
+-----------------+-----------------+-----------------+-----------------+
| GroupMethod     | OT_INTEGER      | 1               | Select method   |
|                 |                 |                 | to determine    |
|                 |                 |                 | graph colouring |
|                 |                 |                 | groups          |
+-----------------+-----------------+-----------------+-----------------+
| IgnoreFilterCri | OT_BOOLEAN      | False           | Activate        |
| t               |                 |                 | accelerating    |
|                 |                 |                 | heuristics for  |
|                 |                 |                 | Filter          |
+-----------------+-----------------+-----------------+-----------------+
| IncBettsTau     | OT_REAL         | 2               | Increase factor |
|                 |                 |                 | for Betts'      |
|                 |                 |                 | update          |
|                 |                 |                 | dampening term  |
+-----------------+-----------------+-----------------+-----------------+
| IncBettsTauMore | OT_REAL         | 100             | Larger increase |
|                 |                 |                 | factor for      |
|                 |                 |                 | Betts' update   |
|                 |                 |                 | dampening term  |
+-----------------+-----------------+-----------------+-----------------+
| IncreaseIWS     | OT_REAL         | 1               | Increase factor |
|                 |                 |                 | for estimated   |
|                 |                 |                 | integer         |
|                 |                 |                 | workspace       |
|                 |                 |                 | requirement     |
+-----------------+-----------------+-----------------+-----------------+
| IncreaseRWS     | OT_REAL         | 1               | Increase factor |
|                 |                 |                 | for estimated   |
|                 |                 |                 | real workspace  |
|                 |                 |                 | requirement     |
+-----------------+-----------------+-----------------+-----------------+
| Infty           | OT_REAL         | 1.000e+20       | Upper bound for |
|                 |                 |                 | numbers to be   |
|                 |                 |                 | regarded as     |
|                 |                 |                 | finite          |
+-----------------+-----------------+-----------------+-----------------+
| InftyUnbounded  | OT_REAL         | 1.000e+20       | Tolerance for   |
|                 |                 |                 | unboundedness   |
|                 |                 |                 | detection       |
|                 |                 |                 | heuristic       |
+-----------------+-----------------+-----------------+-----------------+
| InitialLMest    | OT_BOOLEAN      | True            | Enable initial  |
|                 |                 |                 | Lagrange        |
|                 |                 |                 | multiplier      |
|                 |                 |                 | estimate        |
+-----------------+-----------------+-----------------+-----------------+
| KeepAcceptableS | OT_BOOLEAN      | True            | Save acceptable |
| ol              |                 |                 | solutions as    |
|                 |                 |                 | fallback        |
+-----------------+-----------------+-----------------+-----------------+
| LMestQPipComTol | OT_REAL         | 0.003           | IP              |
|                 |                 |                 | complementarity |
|                 |                 |                 | tolerance in    |
|                 |                 |                 | initial         |
|                 |                 |                 | multiplier      |
|                 |                 |                 | estimate        |
+-----------------+-----------------+-----------------+-----------------+
| LMestQPipResTol | OT_REAL         | 1               | IP residual     |
|                 |                 |                 | tolerance in    |
|                 |                 |                 | initial         |
|                 |                 |                 | multiplier      |
|                 |                 |                 | estimate        |
+-----------------+-----------------+-----------------+-----------------+
| LinMult         | OT_BOOLEAN      | False           | Control         |
|                 |                 |                 | Lagrange        |
|                 |                 |                 | multiplier      |
|                 |                 |                 | update          |
+-----------------+-----------------+-----------------+-----------------+
| LogLevel        | OT_INTEGER      | 0               | Enable XML      |
|                 |                 |                 | logfiles and    |
|                 |                 |                 | writing         |
|                 |                 |                 | interval        |
+-----------------+-----------------+-----------------+-----------------+
| LogResult       | OT_INTEGER      | 0               | Enable XML      |
|                 |                 |                 | result logging  |
|                 |                 |                 | and detail      |
|                 |                 |                 | level           |
+-----------------+-----------------+-----------------+-----------------+
| LowPassAlphaF   | OT_REAL         | 0.950           | Lowpass-filter  |
|                 |                 |                 | update factor   |
|                 |                 |                 | for objective   |
|                 |                 |                 | values          |
+-----------------+-----------------+-----------------+-----------------+
| LowPassAlphaG   | OT_REAL         | 0.950           | Lowpass-filter  |
|                 |                 |                 | update factor   |
|                 |                 |                 | for constraint  |
|                 |                 |                 | values          |
+-----------------+-----------------+-----------------+-----------------+
| LowPassAlphaMer | OT_REAL         | 0.100           | Lowpass-filter  |
| it              |                 |                 | update factor   |
|                 |                 |                 | for merit       |
|                 |                 |                 | function values |
+-----------------+-----------------+-----------------+-----------------+
| LowPassFilter   | OT_BOOLEAN      | True            | Enable lowpass- |
|                 |                 |                 | filter          |
|                 |                 |                 | termination     |
|                 |                 |                 | criterion       |
+-----------------+-----------------+-----------------+-----------------+
| MA97blas3       | OT_BOOLEAN      | False           | Use BLAS level  |
|                 |                 |                 | 3 (dgemm) in    |
|                 |                 |                 | MA97            |
+-----------------+-----------------+-----------------+-----------------+
| MA97mf          | OT_BOOLEAN      | False           | Use             |
|                 |                 |                 | multifrontal-   |
|                 |                 |                 | style forward   |
|                 |                 |                 | solve of MA97   |
+-----------------+-----------------+-----------------+-----------------+
| MA97nemin       | OT_INTEGER      | 8               | Node            |
|                 |                 |                 | amalgation,     |
|                 |                 |                 | controls        |
|                 |                 |                 | merging in      |
|                 |                 |                 | elimination     |
|                 |                 |                 | tree by MA97    |
+-----------------+-----------------+-----------------+-----------------+
| MA97ordering    | OT_INTEGER      | 5               | Ordering used   |
|                 |                 |                 | by MA97         |
+-----------------+-----------------+-----------------+-----------------+
| MA97print       | OT_INTEGER      | -1              | Print level     |
|                 |                 |                 | used by MA97    |
+-----------------+-----------------+-----------------+-----------------+
| MA97scaling     | OT_INTEGER      | 0               | Scaling used by |
|                 |                 |                 | MA97            |
+-----------------+-----------------+-----------------+-----------------+
| MA97small       | OT_REAL         | 0.000           | Any pivot whose |
|                 |                 |                 | modulus is less |
|                 |                 |                 | than this is    |
|                 |                 |                 | treated as zero |
|                 |                 |                 | by MA97         |
+-----------------+-----------------+-----------------+-----------------+
| MA97u           | OT_REAL         | 0.010           | Relative pivot  |
|                 |                 |                 | tolerance of    |
|                 |                 |                 | MA97            |
+-----------------+-----------------+-----------------+-----------------+
| MatrixCC        | OT_BOOLEAN      | False           | Not to be       |
|                 |                 |                 | included into a |
|                 |                 |                 | parameter file! |
+-----------------+-----------------+-----------------+-----------------+
| MaxCalls        | OT_INTEGER      | 2.147e+09       | Upper bound to  |
|                 |                 |                 | Reverse         |
|                 |                 |                 | Communication   |
|                 |                 |                 | calls           |
+-----------------+-----------------+-----------------+-----------------+
| MaxForce        | OT_INTEGER      | 1000            | Maximum number  |
|                 |                 |                 | of Force        |
|                 |                 |                 | recovery        |
|                 |                 |                 | strategy steps  |
+-----------------+-----------------+-----------------+-----------------+
| MaxGPart        | OT_INTEGER      | 1               | (experimental)  |
+-----------------+-----------------+-----------------+-----------------+
| MaxIter         | OT_INTEGER      | 500             | Upper bound on  |
|                 |                 |                 | major           |
|                 |                 |                 | iterations      |
+-----------------+-----------------+-----------------+-----------------+
| MaxLScounter    | OT_INTEGER      | 3               | Control         |
|                 |                 |                 | activation of   |
|                 |                 |                 | Filter          |
|                 |                 |                 | acceleration    |
|                 |                 |                 | heuristics      |
+-----------------+-----------------+-----------------+-----------------+
| MaxNorm         | OT_BOOLEAN      | True            | Select max-norm |
|                 |                 |                 | instead of      |
|                 |                 |                 | 1-norm in       |
|                 |                 |                 | Filter          |
+-----------------+-----------------+-----------------+-----------------+
| MeritFunction   | OT_INTEGER      | 4               | Select merit    |
|                 |                 |                 | function and    |
|                 |                 |                 | penalty update  |
|                 |                 |                 | [0, 3..5]       |
+-----------------+-----------------+-----------------+-----------------+
| MeritGradTol    | OT_REAL         | 0.000           | Threshold of    |
|                 |                 |                 | meritfunction   |
|                 |                 |                 | gradient for    |
|                 |                 |                 | increasing      |
|                 |                 |                 | Hessian         |
|                 |                 |                 | regularisation  |
+-----------------+-----------------+-----------------+-----------------+
| MinBettsTau     | OT_REAL         | 0.000           | Lower bound for |
|                 |                 |                 | Betts' update   |
|                 |                 |                 | dampening term  |
+-----------------+-----------------+-----------------+-----------------+
| MoreRelax       | OT_BOOLEAN      | False           | Introduce one   |
|                 |                 |                 | relaxation      |
|                 |                 |                 | variable for    |
|                 |                 |                 | every           |
|                 |                 |                 | constraint      |
+-----------------+-----------------+-----------------+-----------------+
| NLPmethod       | OT_INTEGER      | 1               | Select (1)      |
|                 |                 |                 | Meritfunction   |
|                 |                 |                 | or (3) Filter   |
|                 |                 |                 | globalisation   |
+-----------------+-----------------+-----------------+-----------------+
| NLPprint        | OT_INTEGER      | 2               | NLP print level |
|                 |                 |                 | [-1..4]         |
+-----------------+-----------------+-----------------+-----------------+
| PairMethod      | OT_INTEGER      | 1               | Select method   |
|                 |                 |                 | to determine    |
|                 |                 |                 | graph colouring |
|                 |                 |                 | pairgroups      |
+-----------------+-----------------+-----------------+-----------------+
| PenUpdEpsBar    | OT_REAL         | 0.900           | Penalty update  |
|                 |                 |                 | parameter       |
|                 |                 |                 | factor for      |
|                 |                 |                 | MeritFunction = |
|                 |                 |                 | 3               |
+-----------------+-----------------+-----------------+-----------------+
| PenUpdEpsKFac   | OT_REAL         | 2               | Penalty update  |
|                 |                 |                 | parameter       |
|                 |                 |                 | factor for      |
|                 |                 |                 | MeritFunction = |
|                 |                 |                 | 4               |
+-----------------+-----------------+-----------------+-----------------+
| PenUpdEpsKSeque | OT_INTEGER      | 2               | Penalty update  |
| nce             |                 |                 | parameter       |
+-----------------+-----------------+-----------------+-----------------+
| PenUpdMaxDeltaK | OT_REAL         | 11              | Max penalty for |
|                 |                 |                 | MeritFunction = |
|                 |                 |                 | 4               |
+-----------------+-----------------+-----------------+-----------------+
| PenUpdMaxFac    | OT_REAL         | 100000000       | Max factor for  |
|                 |                 |                 | increasing      |
|                 |                 |                 | penalty for     |
|                 |                 |                 | MeritFunction = |
|                 |                 |                 | 4               |
+-----------------+-----------------+-----------------+-----------------+
| PenUpdRBar      | OT_REAL         | 2               | Penalty update  |
|                 |                 |                 | parameter for   |
|                 |                 |                 | MeritFunction = |
|                 |                 |                 | 3               |
+-----------------+-----------------+-----------------+-----------------+
| PrecisionF      | OT_REAL         | 0.000           | (currently      |
|                 |                 |                 | unused)         |
|                 |                 |                 | Relative        |
|                 |                 |                 | precision of    |
|                 |                 |                 | objective       |
+-----------------+-----------------+-----------------+-----------------+
| PrecisionG      | OT_REAL         | 0.000           | (currently      |
|                 |                 |                 | unused)         |
|                 |                 |                 | Relative        |
|                 |                 |                 | precision of    |
|                 |                 |                 | constraints     |
+-----------------+-----------------+-----------------+-----------------+
| QPscaleParam    | OT_REAL         | 0               | (currently      |
|                 |                 |                 | unused) Scaling |
|                 |                 |                 | factor for QP   |
+-----------------+-----------------+-----------------+-----------------+
| QuadraticProble | OT_BOOLEAN      | False           | Not to be       |
| m               |                 |                 | included into a |
|                 |                 |                 | parameter file! |
+-----------------+-----------------+-----------------+-----------------+
| ReduceBettsTau  | OT_REAL         | 0.300           | Decrease factor |
|                 |                 |                 | for Betts'      |
|                 |                 |                 | update          |
|                 |                 |                 | dampening term  |
+-----------------+-----------------+-----------------+-----------------+
| RefineFeasibili | OT_INTEGER      | 0               | 0 -             |
| ty              |                 |                 | Deactivated, 1  |
|                 |                 |                 | - After first   |
|                 |                 |                 | feasible        |
|                 |                 |                 | iterate, 2 -    |
|                 |                 |                 | Always on,      |
|                 |                 |                 | Activates       |
|                 |                 |                 | iterative       |
|                 |                 |                 | refinement due  |
|                 |                 |                 | to perturbation |
|                 |                 |                 | in constraints  |
|                 |                 |                 | using           |
|                 |                 |                 | parametric      |
|                 |                 |                 | sensitivities   |
+-----------------+-----------------+-----------------+-----------------+
| RefineMaxHMReg  | OT_REAL         | 1000            | Maximum allowed |
|                 |                 |                 | regularisation  |
|                 |                 |                 | of the hessian  |
|                 |                 |                 | CAUTION         |
|                 |                 |                 | absolute value  |
+-----------------+-----------------+-----------------+-----------------+
| RefineMaxRelax  | OT_REAL         | 0.750           | Maximum allowed |
|                 |                 |                 | relaxation to   |
|                 |                 |                 | apply           |
|                 |                 |                 | feasibility     |
|                 |                 |                 | refinement      |
+-----------------+-----------------+-----------------+-----------------+
| RefineOnlyOnAlp | OT_BOOLEAN      | True            | Activates new   |
| ha              |                 |                 | iterative       |
|                 |                 |                 | refinement of   |
|                 |                 |                 | constraints     |
|                 |                 |                 | only when       |
|                 |                 |                 | Armijo alpha    |
|                 |                 |                 | equals one      |
+-----------------+-----------------+-----------------+-----------------+
| RefineStartTol  | OT_REAL         | 0.000           | Start tolerance |
|                 |                 |                 | for successful  |
|                 |                 |                 | termination of  |
|                 |                 |                 | iterative       |
|                 |                 |                 | refinement due  |
|                 |                 |                 | to perturbation |
|                 |                 |                 | in constraints  |
+-----------------+-----------------+-----------------+-----------------+
| RegStrategy     | OT_INTEGER      | 1               | Select Hessian  |
|                 |                 |                 | regularisation  |
|                 |                 |                 | strategy in     |
|                 |                 |                 | Filter          |
+-----------------+-----------------+-----------------+-----------------+
| ReinitFilter    | OT_BOOLEAN      | False           | Enables Filter- |
|                 |                 |                 | reinitialisatio |
|                 |                 |                 | n accelerating  |
|                 |                 |                 | heuristic       |
+-----------------+-----------------+-----------------+-----------------+
| RelaxMaxDelta   | OT_REAL         | 0.920           | Upper bound for |
|                 |                 |                 | accepting the   |
|                 |                 |                 | constraint      |
|                 |                 |                 | relaxation      |
|                 |                 |                 | variable        |
+-----------------+-----------------+-----------------+-----------------+
| RelaxMaxPen     | OT_REAL         | 50000000        | Upper bound on  |
|                 |                 |                 | the constraint  |
|                 |                 |                 | relaxation      |
|                 |                 |                 | penalty         |
+-----------------+-----------------+-----------------+-----------------+
| RelaxRho        | OT_REAL         | 6               | Update factor   |
|                 |                 |                 | for the         |
|                 |                 |                 | constraint      |
|                 |                 |                 | relaxation      |
|                 |                 |                 | penalty         |
+-----------------+-----------------+-----------------+-----------------+
| RelaxStart      | OT_REAL         | 1               | Initial value   |
|                 |                 |                 | of the          |
|                 |                 |                 | constraint      |
|                 |                 |                 | relaxation      |
|                 |                 |                 | penalty         |
+-----------------+-----------------+-----------------+-----------------+
| RestUntilFeas   | OT_BOOLEAN      | False           | Do restoration  |
|                 |                 |                 | until a         |
|                 |                 |                 | feasible        |
|                 |                 |                 | solution is     |
|                 |                 |                 | found           |
+-----------------+-----------------+-----------------+-----------------+
| ScaleConIter    | OT_BOOLEAN      | False           | Scale           |
|                 |                 |                 | constraints in  |
|                 |                 |                 | every iteration |
+-----------------+-----------------+-----------------+-----------------+
| ScaleFacObj     | OT_REAL         | 10              | Value to scale  |
|                 |                 |                 | large objective |
|                 |                 |                 | functions to    |
+-----------------+-----------------+-----------------+-----------------+
| ScaleFacQP      | OT_REAL         | 10              | Upper bound on  |
|                 |                 |                 | resulting       |
|                 |                 |                 | matrix norm for |
|                 |                 |                 | QP scaling      |
+-----------------+-----------------+-----------------+-----------------+
| ScaledFD        | OT_BOOLEAN      | True            | Use a scaled    |
|                 |                 |                 | perturbation    |
|                 |                 |                 | for finite      |
|                 |                 |                 | differences     |
+-----------------+-----------------+-----------------+-----------------+
| ScaledKKT       | OT_BOOLEAN      | True            | Scale KKT       |
|                 |                 |                 | conditions      |
+-----------------+-----------------+-----------------+-----------------+
| ScaledObj       | OT_BOOLEAN      | True            | Scale the       |
|                 |                 |                 | objective       |
|                 |                 |                 | function        |
+-----------------+-----------------+-----------------+-----------------+
| ScaledQP        | OT_BOOLEAN      | True            | Scale some      |
|                 |                 |                 | matrices handed |
|                 |                 |                 | to the QP       |
+-----------------+-----------------+-----------------+-----------------+
| StartBettsTau   | OT_REAL         | 0.100           | Initial value   |
|                 |                 |                 | for Betts'      |
|                 |                 |                 | update          |
|                 |                 |                 | dampening term  |
+-----------------+-----------------+-----------------+-----------------+
| SteffensenOnRef | OT_BOOLEAN      | False           | Use Steffensen  |
| ine             |                 |                 | Extrapolation   |
|                 |                 |                 | during          |
|                 |                 |                 | Feasibility     |
|                 |                 |                 | Refinement      |
+-----------------+-----------------+-----------------+-----------------+
| SwitchingDelta  | OT_REAL         | 0.010           | Filter          |
|                 |                 |                 | switching       |
|                 |                 |                 | condition       |
|                 |                 |                 | parameter       |
+-----------------+-----------------+-----------------+-----------------+
| SwitchingSCV    | OT_REAL         | 1.100           | Filter          |
|                 |                 |                 | switching       |
|                 |                 |                 | condition       |
|                 |                 |                 | parameter       |
+-----------------+-----------------+-----------------+-----------------+
| SwitchingSF     | OT_REAL         | 2.300           | Filter          |
|                 |                 |                 | switching       |
|                 |                 |                 | condition       |
|                 |                 |                 | parameter       |
+-----------------+-----------------+-----------------+-----------------+
| TakeQPSol       | OT_BOOLEAN      | False           | Evaluate QP     |
|                 |                 |                 | search          |
|                 |                 |                 | direction       |
|                 |                 |                 | regardless of   |
|                 |                 |                 | convergence     |
+-----------------+-----------------+-----------------+-----------------+
| Timeout         | OT_REAL         | 300             | Timeout in      |
|                 |                 |                 | seconds         |
+-----------------+-----------------+-----------------+-----------------+
| TolComp         | OT_REAL         | 0.001           | Complementarity |
|                 |                 |                 | tolerance       |
+-----------------+-----------------+-----------------+-----------------+
| TolFeas         | OT_REAL         | 0.000           | Feasibility     |
|                 |                 |                 | tolerance       |
+-----------------+-----------------+-----------------+-----------------+
| TolOpti         | OT_REAL         | 0.000           | Optimality      |
|                 |                 |                 | tolerance       |
+-----------------+-----------------+-----------------+-----------------+
| TolWeakActive   | OT_REAL         | 1               | (experimental)  |
+-----------------+-----------------+-----------------+-----------------+
| TooBig          | OT_BOOLEAN      | True            | Enable too-big  |
|                 |                 |                 | termination     |
|                 |                 |                 | heuristics      |
+-----------------+-----------------+-----------------+-----------------+
| TooBigCV        | OT_REAL         | 1.000e+25       | Upper bound on  |
|                 |                 |                 | constraint      |
|                 |                 |                 | violation for   |
|                 |                 |                 | too-big         |
|                 |                 |                 | heuristic       |
+-----------------+-----------------+-----------------+-----------------+
| TooBigKKT       | OT_REAL         | 1.000e+30       | Upper bound on  |
|                 |                 |                 | KKT values for  |
|                 |                 |                 | too-big         |
|                 |                 |                 | heuristic       |
+-----------------+-----------------+-----------------+-----------------+
| UpdateMu        | OT_BOOLEAN      | True            | Activates       |
|                 |                 |                 | update of       |
|                 |                 |                 | lagrange        |
|                 |                 |                 | multipliers     |
|                 |                 |                 | during          |
|                 |                 |                 | correction step |
+-----------------+-----------------+-----------------+-----------------+
| UseZen          | OT_BOOLEAN      | False           | Run Zen module  |
|                 |                 |                 | after           |
|                 |                 |                 | successful      |
|                 |                 |                 | termination     |
+-----------------+-----------------+-----------------+-----------------+
| UserDF          | OT_BOOLEAN      | True            | Objective       |
|                 |                 |                 | gradient values |
|                 |                 |                 | supplied by     |
|                 |                 |                 | caller          |
+-----------------+-----------------+-----------------+-----------------+
| UserDG          | OT_BOOLEAN      | True            | Jacobian values |
|                 |                 |                 | supplied by     |
|                 |                 |                 | caller          |
+-----------------+-----------------+-----------------+-----------------+
| UserHM          | OT_BOOLEAN      | True            | Hessian values  |
|                 |                 |                 | supplied by     |
|                 |                 |                 | caller          |
+-----------------+-----------------+-----------------+-----------------+
| UserHMstructure | OT_INTEGER      | 2               | Enable          |
|                 |                 |                 | automatic       |
|                 |                 |                 | Hessian         |
|                 |                 |                 | structure       |
|                 |                 |                 | generation or   |
|                 |                 |                 | checking        |
+-----------------+-----------------+-----------------+-----------------+
| UserZenDGp      | OT_BOOLEAN      | False           | Hessian values  |
|                 |                 |                 | supplied by     |
|                 |                 |                 | caller          |
+-----------------+-----------------+-----------------+-----------------+
| UserZenDLp      | OT_BOOLEAN      | False           | Gradient values |
|                 |                 |                 | supplied by     |
|                 |                 |                 | caller          |
+-----------------+-----------------+-----------------+-----------------+
| UserZenDLpp     | OT_BOOLEAN      | False           | Hessian values  |
|                 |                 |                 | supplied by     |
|                 |                 |                 | caller          |
+-----------------+-----------------+-----------------+-----------------+
| UserZenDLxp     | OT_BOOLEAN      | False           | Hessian values  |
|                 |                 |                 | supplied by     |
|                 |                 |                 | caller          |
+-----------------+-----------------+-----------------+-----------------+
| WeakActiveSet   | OT_BOOLEAN      | False           | (experimental)  |
+-----------------+-----------------+-----------------+-----------------+
| ZenCheckMaxPert | OT_BOOLEAN      | False           | Check maximum   |
|                 |                 |                 | of secure       |
|                 |                 |                 | perturbation    |
|                 |                 |                 | when updating   |
|                 |                 |                 | solution        |
+-----------------+-----------------+-----------------+-----------------+
| ZenFDnewMethod  | OT_BOOLEAN      | True            |                 |
+-----------------+-----------------+-----------------+-----------------+
| ZenRenewLU      | OT_BOOLEAN      | False           | false: use LU   |
|                 |                 |                 | from last QP    |
|                 |                 |                 | step; true:     |
|                 |                 |                 | renew LU        |
|                 |                 |                 | decomposition.  |
+-----------------+-----------------+-----------------+-----------------+
| eps             | OT_REAL         | 0.000           | Machine epsilon |
+-----------------+-----------------+-----------------+-----------------+
| internalParChan | OT_INTEGER      | 0               | Counter for     |
| ged             |                 |                 | changed         |
|                 |                 |                 | parameters.     |
|                 |                 |                 | Internal use    |
|                 |                 |                 | only.           |
+-----------------+-----------------+-----------------+-----------------+
| print_time      | OT_BOOLEAN      | True            | Print           |
|                 |                 |                 | information     |
|                 |                 |                 | about execution |
|                 |                 |                 | time            |
+-----------------+-----------------+-----------------+-----------------+
| qp_ipBarrier    | OT_REAL         | 7.800           | IP barrier      |
|                 |                 |                 | parameter.      |
+-----------------+-----------------+-----------------+-----------------+
| qp_ipComTol     | OT_REAL         | 0.000           | IP              |
|                 |                 |                 | complementarity |
|                 |                 |                 | tolerance.      |
+-----------------+-----------------+-----------------+-----------------+
| qp_ipFracBound  | OT_REAL         | 0.880           | IP fraction-to- |
|                 |                 |                 | the-boundary    |
|                 |                 |                 | parameter.      |
+-----------------+-----------------+-----------------+-----------------+
| qp_ipLsMethod   | OT_STRING       | None            | Select the      |
|                 |                 |                 | direct linear   |
|                 |                 |                 | solver used by  |
|                 |                 |                 | the IP method.  |
+-----------------+-----------------+-----------------+-----------------+
| qp_ipMinAlpha   | OT_REAL         | 0.000           | IP line search  |
|                 |                 |                 | minimum step    |
|                 |                 |                 | size.           |
+-----------------+-----------------+-----------------+-----------------+
| qp_ipRelaxDiv   | OT_REAL         | 2               | The relaxation  |
|                 |                 |                 | term is divided |
|                 |                 |                 | by this value   |
|                 |                 |                 | if successful.  |
+-----------------+-----------------+-----------------+-----------------+
| qp_ipRelaxMax   | OT_REAL         | 0.000           | Maximum         |
|                 |                 |                 | relaxation      |
|                 |                 |                 | value.          |
+-----------------+-----------------+-----------------+-----------------+
| qp_ipRelaxMin   | OT_REAL         | 0.000           | Mimimum         |
|                 |                 |                 | relaxation      |
|                 |                 |                 | value.          |
+-----------------+-----------------+-----------------+-----------------+
| qp_ipRelaxMult  | OT_REAL         | 10              | The relaxation  |
|                 |                 |                 | term is         |
|                 |                 |                 | multiplied by   |
|                 |                 |                 | this value if   |
|                 |                 |                 | unsuccessful.   |
+-----------------+-----------------+-----------------+-----------------+
| qp_ipResTol     | OT_REAL         | 0.000           | IP residuals    |
|                 |                 |                 | tolerance.      |
+-----------------+-----------------+-----------------+-----------------+
| qp_ipTryRelax   | OT_BOOLEAN      | True            | Enable          |
|                 |                 |                 | relaxation      |
|                 |                 |                 | strategy when   |
|                 |                 |                 | encountering an |
|                 |                 |                 | error.          |
+-----------------+-----------------+-----------------+-----------------+
| qp_lsItMaxIter  | OT_INTEGER      | 1000            | Maximum number  |
|                 |                 |                 | of iterations   |
|                 |                 |                 | of the          |
|                 |                 |                 | iterative       |
|                 |                 |                 | linear solvers. |
+-----------------+-----------------+-----------------+-----------------+
| qp_lsItMethod   | OT_STRING       | None            | Select the      |
|                 |                 |                 | iterative       |
|                 |                 |                 | linear solver.  |
+-----------------+-----------------+-----------------+-----------------+
| qp_lsItPrecondM | OT_STRING       | None            | Select          |
| ethod           |                 |                 | preconditioner  |
|                 |                 |                 | for the         |
|                 |                 |                 | iterative       |
|                 |                 |                 | linear solver.  |
+-----------------+-----------------+-----------------+-----------------+
| qp_lsRefineMaxI | OT_INTEGER      | 10              | Maximum number  |
| ter             |                 |                 | of iterative    |
|                 |                 |                 | refinement      |
|                 |                 |                 | steps of the    |
|                 |                 |                 | direct linear   |
|                 |                 |                 | solvers.        |
+-----------------+-----------------+-----------------+-----------------+
| qp_lsScale      | OT_BOOLEAN      | True            | Enables scaling |
|                 |                 |                 | on linear       |
|                 |                 |                 | solver level.   |
+-----------------+-----------------+-----------------+-----------------+
| qp_lsTol        | OT_REAL         | 0.000           | Tolerance for   |
|                 |                 |                 | the linear      |
|                 |                 |                 | solver.         |
+-----------------+-----------------+-----------------+-----------------+
| qp_lsTrySimple  | OT_BOOLEAN      | False           | Some matrices   |
|                 |                 |                 | can be solved   |
|                 |                 |                 | without calling |
|                 |                 |                 | a linear        |
|                 |                 |                 | equation solver |
|                 |                 |                 | .Currently only |
|                 |                 |                 | diagonal        |
|                 |                 |                 | matrices are    |
|                 |                 |                 | supported.Non-  |
|                 |                 |                 | diagonal        |
|                 |                 |                 | matrices will   |
|                 |                 |                 | besolved with   |
|                 |                 |                 | the chosen      |
|                 |                 |                 | linear equation |
|                 |                 |                 | solver.         |
+-----------------+-----------------+-----------------+-----------------+
| qp_maxIter      | OT_INTEGER      | 80              | Imposes an      |
|                 |                 |                 | upper limit on  |
|                 |                 |                 | the number of   |
|                 |                 |                 | minor solver    |
|                 |                 |                 | iterations,     |
|                 |                 |                 | i.e. for the    |
|                 |                 |                 | quadratic       |
|                 |                 |                 | subproblem      |
|                 |                 |                 | solver.If the   |
|                 |                 |                 | limit is        |
|                 |                 |                 | reached before  |
|                 |                 |                 | convergence,    |
|                 |                 |                 | WORHP will      |
|                 |                 |                 | activate QP     |
|                 |                 |                 | recovery        |
|                 |                 |                 | strategies to   |
|                 |                 |                 | prevent a       |
|                 |                 |                 | solver          |
|                 |                 |                 | breakdown.      |
+-----------------+-----------------+-----------------+-----------------+
| qp_method       | OT_STRING       | None            | Select the      |
|                 |                 |                 | solution method |
|                 |                 |                 | used by the QP  |
|                 |                 |                 | solver.         |
+-----------------+-----------------+-----------------+-----------------+
| qp_nsnBeta      | OT_REAL         | 0.900           | NSN stepsize    |
|                 |                 |                 | decrease        |
|                 |                 |                 | factor.         |
+-----------------+-----------------+-----------------+-----------------+
| qp_nsnGradStep  | OT_BOOLEAN      | True            | Enable gradient |
|                 |                 |                 | steps in the    |
|                 |                 |                 | NSN method.     |
+-----------------+-----------------+-----------------+-----------------+
| qp_nsnKKT       | OT_REAL         | 0.000           | NSN KKT         |
|                 |                 |                 | tolerance.      |
+-----------------+-----------------+-----------------+-----------------+
| qp_nsnLsMethod  | OT_STRING       | None            | Select the      |
|                 |                 |                 | direct linear   |
|                 |                 |                 | solver used by  |
|                 |                 |                 | the NSN method. |
+-----------------+-----------------+-----------------+-----------------+
| qp_nsnMinAlpha  | OT_REAL         | 0.000           | NSN line search |
|                 |                 |                 | minimum step    |
|                 |                 |                 | size.           |
+-----------------+-----------------+-----------------+-----------------+
| qp_nsnSigma     | OT_REAL         | 0.010           | NSN line search |
|                 |                 |                 | slope           |
|                 |                 |                 | parameter.      |
+-----------------+-----------------+-----------------+-----------------+
| qp_printLevel   | OT_STRING       | None            | Controls the    |
|                 |                 |                 | amount of QP    |
|                 |                 |                 | solver output.  |
+-----------------+-----------------+-----------------+-----------------+
| qp_scaleIntern  | OT_BOOLEAN      | False           | Enable scaling  |
|                 |                 |                 | on QP level.    |
+-----------------+-----------------+-----------------+-----------------+
| qp_strict       | OT_BOOLEAN      | True            | Use strict      |
|                 |                 |                 | termination     |
|                 |                 |                 | criteria in IP  |
|                 |                 |                 | method.         |
+-----------------+-----------------+-----------------+-----------------+

>List of available monitors

+-------------+
|     Id      |
+=============+
| eval_f      |
+-------------+
| eval_g      |
+-------------+
| eval_grad_f |
+-------------+
| eval_h      |
+-------------+
| eval_jac_g  |
+-------------+

>List of available stats

+--------------------+
|         Id         |
+====================+
| iter_count         |
+--------------------+
| iteration          |
+--------------------+
| iterations         |
+--------------------+
| n_eval_f           |
+--------------------+
| n_eval_g           |
+--------------------+
| n_eval_grad_f      |
+--------------------+
| n_eval_h           |
+--------------------+
| n_eval_jac_g       |
+--------------------+
| return_code        |
+--------------------+
| return_status      |
+--------------------+
| t_callback_fun     |
+--------------------+
| t_callback_prepare |
+--------------------+
| t_eval_f           |
+--------------------+
| t_eval_g           |
+--------------------+
| t_eval_grad_f      |
+--------------------+
| t_eval_h           |
+--------------------+
| t_eval_jac_g       |
+--------------------+
| t_mainloop         |
+--------------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

scpgen
------



A structure-exploiting sequential quadratic programming (to be come
sequential convex programming) method for nonlinear programming.

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| beta            | OT_REAL         | 0.800           | Line-search     |
|                 |                 |                 | parameter,      |
|                 |                 |                 | restoration     |
|                 |                 |                 | factor of       |
|                 |                 |                 | stepsize        |
+-----------------+-----------------+-----------------+-----------------+
| c1              | OT_REAL         | 0.000           | Armijo          |
|                 |                 |                 | condition,      |
|                 |                 |                 | coefficient of  |
|                 |                 |                 | decrease in     |
|                 |                 |                 | merit           |
+-----------------+-----------------+-----------------+-----------------+
| codegen         | OT_BOOLEAN      | false           | C-code          |
|                 |                 |                 | generation      |
+-----------------+-----------------+-----------------+-----------------+
| hessian_approxi | OT_STRING       | \"exact\"         | gauss-          |
| mation          |                 |                 | newton|exact    |
+-----------------+-----------------+-----------------+-----------------+
| lbfgs_memory    | OT_INTEGER      | 10              | Size of L-BFGS  |
|                 |                 |                 | memory.         |
+-----------------+-----------------+-----------------+-----------------+
| max_iter        | OT_INTEGER      | 50              | Maximum number  |
|                 |                 |                 | of SQP          |
|                 |                 |                 | iterations      |
+-----------------+-----------------+-----------------+-----------------+
| max_iter_ls     | OT_INTEGER      | 1               | Maximum number  |
|                 |                 |                 | of linesearch   |
|                 |                 |                 | iterations      |
+-----------------+-----------------+-----------------+-----------------+
| merit_memsize   | OT_INTEGER      | 4               | Size of memory  |
|                 |                 |                 | to store        |
|                 |                 |                 | history of      |
|                 |                 |                 | merit function  |
|                 |                 |                 | values          |
+-----------------+-----------------+-----------------+-----------------+
| merit_start     | OT_REAL         | 0.000           | Lower bound for |
|                 |                 |                 | the merit       |
|                 |                 |                 | function        |
|                 |                 |                 | parameter       |
+-----------------+-----------------+-----------------+-----------------+
| name_x          | OT_STRINGVECTOR | GenericType()   | Names of the    |
|                 |                 |                 | variables.      |
+-----------------+-----------------+-----------------+-----------------+
| print_header    | OT_BOOLEAN      | true            | Print the       |
|                 |                 |                 | header with     |
|                 |                 |                 | problem         |
|                 |                 |                 | statistics      |
+-----------------+-----------------+-----------------+-----------------+
| print_time      | OT_BOOLEAN      | true            | Print           |
|                 |                 |                 | information     |
|                 |                 |                 | about execution |
|                 |                 |                 | time            |
+-----------------+-----------------+-----------------+-----------------+
| print_x         | OT_INTEGERVECTO | GenericType()   | Which variables |
|                 | R               |                 | to print.       |
+-----------------+-----------------+-----------------+-----------------+
| qp_solver       | OT_STRING       | GenericType()   | The QP solver   |
|                 |                 |                 | to be used by   |
|                 |                 |                 | the SQP method  |
+-----------------+-----------------+-----------------+-----------------+
| qp_solver_optio | OT_DICT         | GenericType()   | Options to be   |
| ns              |                 |                 | passed to the   |
|                 |                 |                 | QP solver       |
+-----------------+-----------------+-----------------+-----------------+
| reg_threshold   | OT_REAL         | 0.000           | Threshold for   |
|                 |                 |                 | the             |
|                 |                 |                 | regularization. |
+-----------------+-----------------+-----------------+-----------------+
| regularize      | OT_BOOLEAN      | false           | Automatic       |
|                 |                 |                 | regularization  |
|                 |                 |                 | of Lagrange     |
|                 |                 |                 | Hessian.        |
+-----------------+-----------------+-----------------+-----------------+
| tol_du          | OT_REAL         | 0.000           | Stopping        |
|                 |                 |                 | criterion for   |
|                 |                 |                 | dual            |
|                 |                 |                 | infeasability   |
+-----------------+-----------------+-----------------+-----------------+
| tol_pr          | OT_REAL         | 0.000           | Stopping        |
|                 |                 |                 | criterion for   |
|                 |                 |                 | primal          |
|                 |                 |                 | infeasibility   |
+-----------------+-----------------+-----------------+-----------------+
| tol_pr_step     | OT_REAL         | 0.000           | Stopping        |
|                 |                 |                 | criterion for   |
|                 |                 |                 | the step size   |
+-----------------+-----------------+-----------------+-----------------+
| tol_reg         | OT_REAL         | 0.000           | Stopping        |
|                 |                 |                 | criterion for   |
|                 |                 |                 | regularization  |
+-----------------+-----------------+-----------------+-----------------+

>List of available monitors

+-------------+
|     Id      |
+=============+
| dx          |
+-------------+
| eval_f      |
+-------------+
| eval_g      |
+-------------+
| eval_grad_f |
+-------------+
| eval_h      |
+-------------+
| eval_jac_g  |
+-------------+
| qp          |
+-------------+

>List of available stats

+------------+
|     Id     |
+============+
| iter_count |
+------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

sqpmethod
---------



A textbook SQPMethod

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| beta            | OT_REAL         | 0.800           | Line-search     |
|                 |                 |                 | parameter,      |
|                 |                 |                 | restoration     |
|                 |                 |                 | factor of       |
|                 |                 |                 | stepsize        |
+-----------------+-----------------+-----------------+-----------------+
| c1              | OT_REAL         | 0.000           | Armijo          |
|                 |                 |                 | condition,      |
|                 |                 |                 | coefficient of  |
|                 |                 |                 | decrease in     |
|                 |                 |                 | merit           |
+-----------------+-----------------+-----------------+-----------------+
| hessian_approxi | OT_STRING       | \"exact\"         | limited-        |
| mation          |                 |                 | memory|exact    |
+-----------------+-----------------+-----------------+-----------------+
| lbfgs_memory    | OT_INTEGER      | 10              | Size of L-BFGS  |
|                 |                 |                 | memory.         |
+-----------------+-----------------+-----------------+-----------------+
| max_iter        | OT_INTEGER      | 50              | Maximum number  |
|                 |                 |                 | of SQP          |
|                 |                 |                 | iterations      |
+-----------------+-----------------+-----------------+-----------------+
| max_iter_ls     | OT_INTEGER      | 3               | Maximum number  |
|                 |                 |                 | of linesearch   |
|                 |                 |                 | iterations      |
+-----------------+-----------------+-----------------+-----------------+
| merit_memory    | OT_INTEGER      | 4               | Size of memory  |
|                 |                 |                 | to store        |
|                 |                 |                 | history of      |
|                 |                 |                 | merit function  |
|                 |                 |                 | values          |
+-----------------+-----------------+-----------------+-----------------+
| min_step_size   | OT_REAL         | 0.000           | The size (inf-  |
|                 |                 |                 | norm) of the    |
|                 |                 |                 | step size       |
|                 |                 |                 | should not      |
|                 |                 |                 | become smaller  |
|                 |                 |                 | than this.      |
+-----------------+-----------------+-----------------+-----------------+
| print_header    | OT_BOOLEAN      | true            | Print the       |
|                 |                 |                 | header with     |
|                 |                 |                 | problem         |
|                 |                 |                 | statistics      |
+-----------------+-----------------+-----------------+-----------------+
| print_time      | OT_BOOLEAN      | true            | Print           |
|                 |                 |                 | information     |
|                 |                 |                 | about execution |
|                 |                 |                 | time            |
+-----------------+-----------------+-----------------+-----------------+
| qp_solver       | OT_STRING       | GenericType()   | The QP solver   |
|                 |                 |                 | to be used by   |
|                 |                 |                 | the SQP method  |
+-----------------+-----------------+-----------------+-----------------+
| qp_solver_optio | OT_DICT         | GenericType()   | Options to be   |
| ns              |                 |                 | passed to the   |
|                 |                 |                 | QP solver       |
+-----------------+-----------------+-----------------+-----------------+
| regularize      | OT_BOOLEAN      | false           | Automatic       |
|                 |                 |                 | regularization  |
|                 |                 |                 | of Lagrange     |
|                 |                 |                 | Hessian.        |
+-----------------+-----------------+-----------------+-----------------+
| tol_du          | OT_REAL         | 0.000           | Stopping        |
|                 |                 |                 | criterion for   |
|                 |                 |                 | dual            |
|                 |                 |                 | infeasability   |
+-----------------+-----------------+-----------------+-----------------+
| tol_pr          | OT_REAL         | 0.000           | Stopping        |
|                 |                 |                 | criterion for   |
|                 |                 |                 | primal          |
|                 |                 |                 | infeasibility   |
+-----------------+-----------------+-----------------+-----------------+

>List of available monitors

+-------------+
|     Id      |
+=============+
| bfgs        |
+-------------+
| dx          |
+-------------+
| eval_f      |
+-------------+
| eval_g      |
+-------------+
| eval_grad_f |
+-------------+
| eval_h      |
+-------------+
| eval_jac_g  |
+-------------+
| qp          |
+-------------+

>List of available stats

+--------------------+
|         Id         |
+====================+
| iter_count         |
+--------------------+
| iteration          |
+--------------------+
| iterations         |
+--------------------+
| n_eval_f           |
+--------------------+
| n_eval_g           |
+--------------------+
| n_eval_grad_f      |
+--------------------+
| n_eval_h           |
+--------------------+
| n_eval_jac_g       |
+--------------------+
| return_status      |
+--------------------+
| t_callback_fun     |
+--------------------+
| t_callback_prepare |
+--------------------+
| t_eval_f           |
+--------------------+
| t_eval_g           |
+--------------------+
| t_eval_grad_f      |
+--------------------+
| t_eval_h           |
+--------------------+
| t_eval_jac_g       |
+--------------------+
| t_mainloop         |
+--------------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

stabilizedsqp
-------------



Stabilized Sequential Quadratic Programming method.

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| TReta1          | OT_REAL         | 0.800           | Required        |
|                 |                 |                 | predicted /     |
|                 |                 |                 | actual decrease |
|                 |                 |                 | for TR increase |
+-----------------+-----------------+-----------------+-----------------+
| TReta2          | OT_REAL         | 0.200           | Required        |
|                 |                 |                 | predicted /     |
|                 |                 |                 | actual decrease |
|                 |                 |                 | for TR decrease |
+-----------------+-----------------+-----------------+-----------------+
| alphaMin        | OT_REAL         | 0.001           | Used to check   |
|                 |                 |                 | whether to      |
|                 |                 |                 | increase rho.   |
+-----------------+-----------------+-----------------+-----------------+
| beta            | OT_REAL         | 0.500           | Line-search     |
|                 |                 |                 | parameter,      |
|                 |                 |                 | restoration     |
|                 |                 |                 | factor of       |
|                 |                 |                 | stepsize        |
+-----------------+-----------------+-----------------+-----------------+
| c1              | OT_REAL         | 0.001           | Armijo          |
|                 |                 |                 | condition,      |
|                 |                 |                 | coefficient of  |
|                 |                 |                 | decrease in     |
|                 |                 |                 | merit           |
+-----------------+-----------------+-----------------+-----------------+
| dvMax0          | OT_REAL         | 100             | Parameter used  |
|                 |                 |                 | to defined the  |
|                 |                 |                 | max step        |
|                 |                 |                 | length.         |
+-----------------+-----------------+-----------------+-----------------+
| eps_active      | OT_REAL         | 0.000           | Threshold for   |
|                 |                 |                 | the epsilon-    |
|                 |                 |                 | active set.     |
+-----------------+-----------------+-----------------+-----------------+
| gamma1          | OT_REAL         | 2               | Trust region    |
|                 |                 |                 | increase        |
|                 |                 |                 | parameter       |
+-----------------+-----------------+-----------------+-----------------+
| gamma2          | OT_REAL         | 1               | Trust region    |
|                 |                 |                 | update          |
|                 |                 |                 | parameter       |
+-----------------+-----------------+-----------------+-----------------+
| gamma3          | OT_REAL         | 1               | Trust region    |
|                 |                 |                 | decrease        |
|                 |                 |                 | parameter       |
+-----------------+-----------------+-----------------+-----------------+
| hessian_approxi | OT_STRING       | \"exact\"         | limited-        |
| mation          |                 |                 | memory|exact    |
+-----------------+-----------------+-----------------+-----------------+
| lbfgs_memory    | OT_INTEGER      | 10              | Size of L-BFGS  |
|                 |                 |                 | memory.         |
+-----------------+-----------------+-----------------+-----------------+
| max_iter        | OT_INTEGER      | 100             | Maximum number  |
|                 |                 |                 | of SQP          |
|                 |                 |                 | iterations      |
+-----------------+-----------------+-----------------+-----------------+
| max_iter_ls     | OT_INTEGER      | 20              | Maximum number  |
|                 |                 |                 | of linesearch   |
|                 |                 |                 | iterations      |
+-----------------+-----------------+-----------------+-----------------+
| max_time        | OT_REAL         | 1.000e+12       | Timeout         |
+-----------------+-----------------+-----------------+-----------------+
| merit_memory    | OT_INTEGER      | 4               | Size of memory  |
|                 |                 |                 | to store        |
|                 |                 |                 | history of      |
|                 |                 |                 | merit function  |
|                 |                 |                 | values          |
+-----------------+-----------------+-----------------+-----------------+
| min_step_size   | OT_REAL         | 0.000           | The size (inf-  |
|                 |                 |                 | norm) of the    |
|                 |                 |                 | step size       |
|                 |                 |                 | should not      |
|                 |                 |                 | become smaller  |
|                 |                 |                 | than this.      |
+-----------------+-----------------+-----------------+-----------------+
| muR0            | OT_REAL         | 0.000           | Initial choice  |
|                 |                 |                 | of              |
|                 |                 |                 | regularization  |
|                 |                 |                 | parameter       |
+-----------------+-----------------+-----------------+-----------------+
| nu              | OT_REAL         | 1               | Parameter for   |
|                 |                 |                 | primal-dual     |
|                 |                 |                 | augmented       |
|                 |                 |                 | Lagrangian.     |
+-----------------+-----------------+-----------------+-----------------+
| phiWeight       | OT_REAL         | 0.000           | Weight used in  |
|                 |                 |                 | pseudo-filter.  |
+-----------------+-----------------+-----------------+-----------------+
| print_header    | OT_BOOLEAN      | true            | Print the       |
|                 |                 |                 | header with     |
|                 |                 |                 | problem         |
|                 |                 |                 | statistics      |
+-----------------+-----------------+-----------------+-----------------+
| regularize      | OT_BOOLEAN      | false           | Automatic       |
|                 |                 |                 | regularization  |
|                 |                 |                 | of Lagrange     |
|                 |                 |                 | Hessian.        |
+-----------------+-----------------+-----------------+-----------------+
| stabilized_qp_s | OT_STRING       | GenericType()   | The Stabilized  |
| olver           |                 |                 | QP solver to be |
|                 |                 |                 | used by the SQP |
|                 |                 |                 | method          |
+-----------------+-----------------+-----------------+-----------------+
| stabilized_qp_s | OT_DICT         | GenericType()   | Options to be   |
| olver_options   |                 |                 | passed to the   |
|                 |                 |                 | Stabilized QP   |
|                 |                 |                 | solver          |
+-----------------+-----------------+-----------------+-----------------+
| tau0            | OT_REAL         | 0.010           | Initial         |
|                 |                 |                 | parameter for   |
|                 |                 |                 | the merit       |
|                 |                 |                 | function        |
|                 |                 |                 | optimality      |
|                 |                 |                 | threshold.      |
+-----------------+-----------------+-----------------+-----------------+
| tol_du          | OT_REAL         | 0.000           | Stopping        |
|                 |                 |                 | criterion for   |
|                 |                 |                 | dual            |
|                 |                 |                 | infeasability   |
+-----------------+-----------------+-----------------+-----------------+
| tol_pr          | OT_REAL         | 0.000           | Stopping        |
|                 |                 |                 | criterion for   |
|                 |                 |                 | primal          |
|                 |                 |                 | infeasibility   |
+-----------------+-----------------+-----------------+-----------------+
| yEinitial       | OT_STRING       | \"simple\"        | Initial         |
|                 |                 |                 | multiplier.     |
|                 |                 |                 | Simple (all     |
|                 |                 |                 | zero) or least  |
|                 |                 |                 | (LSQ).          |
+-----------------+-----------------+-----------------+-----------------+

>List of available monitors

+-------------+
|     Id      |
+=============+
| dx          |
+-------------+
| eval_f      |
+-------------+
| eval_g      |
+-------------+
| eval_grad_f |
+-------------+
| eval_h      |
+-------------+
| eval_jac_g  |
+-------------+
| qp          |
+-------------+

>List of available stats

+---------------+
|      Id       |
+===============+
| iter_count    |
+---------------+
| return_status |
+---------------+

--------------------------------------------------------------------------------



Joel Andersson
Diagrams
--------



C++ includes: nlp_solver.hpp ";

%feature("docstring") casadi::NlpSolver::repr "

Print a representation of the object.

";

%feature("docstring") casadi::NlpSolver::spEvaluate "[INTERNAL]  Propagate
the sparsity pattern through a set of directional.

derivatives forward or backward (for usage, see the example
propagating_sparsity.cpp)

";

%feature("docstring") casadi::NlpSolver::setOptionByEnumValue "[INTERNAL]
Set a certain option by giving an enum value.

";

%feature("docstring") casadi::NlpSolver::isNull "

Is a null pointer?

";

%feature("docstring") casadi::NlpSolver::name "

Name of the function.

";

%feature("docstring") casadi::NlpSolver::setOptionsFromFile "

Read options from parameter xml.

";

%feature("docstring") casadi::NlpSolver::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::NlpSolver::checkInputs "[INTERNAL]  Check if
the numerical values of the supplied bounds make sense.

";

%feature("docstring") casadi::NlpSolver::nnz_out "

Get of number of output nonzeros For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::NlpSolver::gradF "

Access the objective gradient function>Input scheme: casadi::GradFInput
(GRADF_NUM_IN = 2) [gradFIn]

+-----------+-------+---------------------+
| Full name | Short |     Description     |
+===========+=======+=====================+
| GRADF_X   | x     | Decision variable . |
+-----------+-------+---------------------+
| GRADF_P   | p     | Fixed parameter .   |
+-----------+-------+---------------------+

";

%feature("docstring") casadi::NlpSolver::sparsity_in "

Get sparsity of a given input.

";

%feature("docstring") casadi::NlpSolver::dictionary "

Get the dictionary.

";

%feature("docstring") casadi::NlpSolver::free_sx "

Get all the free variables of the function.

";

%feature("docstring") casadi::NlpSolver::call "

Evaluate the function symbolically or numerically.

";

%feature("docstring") casadi::NlpSolver::sx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::NlpSolver::index_out "

Find the index for a string describing a particular entry of an output
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::NlpSolver::n_out "

Get the number of function outputs.

";

%feature("docstring") casadi::NlpSolver::n_in "

Get the number of function inputs.

";

%feature("docstring") casadi::NlpSolver::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::NlpSolver::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::NlpSolver::evaluate "

Evaluate.

";

%feature("docstring") casadi::NlpSolver::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::NlpSolver::size1_out "

Get output dimension.

";

%feature("docstring") casadi::NlpSolver::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::NlpSolver::name_in "

>  [str] Function.name_in() const 
------------------------------------------------------------------------

Get input scheme.

>  str Function.name_in(int ind) const 
------------------------------------------------------------------------

Get input scheme name by index.

";

%feature("docstring") casadi::NlpSolver::size2_in "

Get input dimension.

";

%feature("docstring") casadi::NlpSolver::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::NlpSolver::size_in "

Get input dimension.

";

%feature("docstring") casadi::NlpSolver::getOptionEnumValue "[INTERNAL]
Get the enum value corresponding to th certain option.

";

%feature("docstring") casadi::NlpSolver::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::NlpSolver::sz_arg "[INTERNAL]  Get required
length of arg field.

";

%feature("docstring") casadi::NlpSolver::callDerivative "[INTERNAL]
Evaluate the function symbolically or numerically with directional
derivatives The first two arguments are the nondifferentiated inputs and
results of the evaluation, the next two arguments are a set of forward
directional seeds and the resulting forward directional derivatives, the
length of the vector being the number of forward directions. The next two
arguments are a set of adjoint directional seeds and the resulting adjoint
directional derivatives, the length of the vector being the number of
adjoint directions.

";

%feature("docstring") casadi::NlpSolver::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::NlpSolver::getReportConstraints "";

%feature("docstring") casadi::NlpSolver::is_a "

Check if the function is of a particular type Optionally check if name
matches one of the base classes (default true)

";

%feature("docstring") casadi::NlpSolver::nnz_in "

Get of number of input nonzeros For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::NlpSolver::spInit "[INTERNAL]  Reset the
sparsity propagation.

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::NlpSolver::getStat "

Get a single statistic obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::NlpSolver::mapaccum "

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
  



";

%feature("docstring") casadi::NlpSolver::hessLag "

Access the Jacobian of the constraint function.

>Input scheme: casadi::HessLagInput (HESSLAG_NUM_IN = 4) [hessLagIn]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| HESSLAG_X              | x                      | Decision variable .    |
+------------------------+------------------------+------------------------+
| HESSLAG_P              | p                      | Fixed parameter .      |
+------------------------+------------------------+------------------------+
| HESSLAG_LAM_F          | lam_f                  | Multiplier for f. Just |
|                        |                        | a scalar factor for    |
|                        |                        | the objective that the |
|                        |                        | NLP solver might use   |
|                        |                        | to scale the           |
|                        |                        | objective.             |
+------------------------+------------------------+------------------------+
| HESSLAG_LAM_G          | lam_g                  | Multiplier for g .     |
+------------------------+------------------------+------------------------+

>Output scheme: casadi::HessLagOutput (HESSLAG_NUM_OUT = 5) [hessLagOut]

+----------------+--------+------------------------------------------------+
|   Full name    | Short  |                  Description                   |
+================+========+================================================+
| HESSLAG_HESS   | hess   | Hessian of the Lagrangian .                    |
+----------------+--------+------------------------------------------------+
| HESSLAG_F      | f      | Objective function .                           |
+----------------+--------+------------------------------------------------+
| HESSLAG_G      | g      | Constraint function .                          |
+----------------+--------+------------------------------------------------+
| HESSLAG_GRAD_X | grad_x | Gradient of the Lagrangian with respect to x . |
+----------------+--------+------------------------------------------------+
| HESSLAG_GRAD_P | grad_p | Gradient of the Lagrangian with respect to p . |
+----------------+--------+------------------------------------------------+

";

%feature("docstring") casadi::NlpSolver::index_in "

Find the index for a string describing a particular entry of an input
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::NlpSolver::setJacSparsity "

Generate the sparsity of a Jacobian block

";


// File: classcasadi_1_1NonZeros.xml
%feature("docstring") casadi::NonZeros::NonZeros "

Constructor.

";

%feature("docstring") casadi::NonZeros "

Access to a set of nonzeros.

NonZeros class for Matrix NonZeros is the return type for operator[] of the
Matrix class, it allows access to the value as well as changing the parent
object Joel Andersson

C++ includes: nonzeros.hpp ";


// File: classcasadi_1_1Norm.xml


// File: classcasadi_1_1Norm1.xml


// File: classcasadi_1_1Norm2.xml


// File: classcasadi_1_1NormF.xml


// File: classcasadi_1_1NormInf.xml


// File: classcasadi_1_1Nullspace.xml


/*  Simple Getters & Setters  */

/*  Advanced Getters  */

/*  Option Functionality  */ %feature("docstring")
casadi::Nullspace::size1_in "

Get input dimension.

";

%feature("docstring") casadi::Nullspace::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::Nullspace::getAtomicOutput "

Get the (integer) output argument of an atomic operation.

";

%feature("docstring") casadi::Nullspace::setInput "

>  void IOInterface< Function  >.setInput(T val, int iind=0)
------------------------------------------------------------------------

Set an input by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.setInput(T val, str iname)
------------------------------------------------------------------------

Set an input by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::Nullspace::derForward "

Get a function that calculates nfwd forward derivatives.

Returns a function with n_in + n_out +nfwd*n_in inputs and nfwd*n_out
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nfwd*n_in inputs correspond to forward seeds, one direction at a time The
nfwd*n_out outputs correspond to forward sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::Nullspace::spCanEvaluate "[INTERNAL]  Is the
class able to propagate seeds through the algorithm?

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::Nullspace::name_out "

>  [str] Function.name_out() const 
------------------------------------------------------------------------

Get output scheme.

>  str Function.name_out(int ind) const 
------------------------------------------------------------------------

Get output scheme name by index.

";

%feature("docstring") casadi::Nullspace::jacobian "

Generate a Jacobian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. If
compact is set to true, only the nonzeros of the input and output
expressions are considered. If symmetric is set to true, the Jacobian being
calculated is known to be symmetric (usually a Hessian), which can be
exploited by the algorithm.

The generated Jacobian has one more output than the calling function
corresponding to the Jacobian and the same number of inputs.

";

%feature("docstring") casadi::Nullspace::size2_out "

Get output dimension.

";

%feature("docstring") casadi::Nullspace::sz_res "[INTERNAL]  Get required
length of res field.

";

%feature("docstring") casadi::Nullspace::getWorkSize "

Get the length of the work vector.

";

%feature("docstring") casadi::Nullspace::call "

Evaluate the function symbolically or numerically.

";

%feature("docstring") casadi::Nullspace::gradient "

Generate a gradient function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the output must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::Nullspace::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

";

%feature("docstring") casadi::Nullspace::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::Nullspace::sx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::Nullspace::getOutput "

>  DMatrix  IOInterface< Function  >.getOutput(int oind=0) const
------------------------------------------------------------------------

Get an output by index.

Parameters:
-----------

oind:  index within the range [0..n_out()-1]

>  DMatrix  IOInterface< Function  >.getOutput(str oname) const
------------------------------------------------------------------------

Get an output by name.

Parameters:
-----------

oname:  output name. Only allowed when an output scheme is set.

>  void IOInterface< Function  >.getOutput(T val, int oind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.getOutput(T val, str oname)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::Nullspace::jacSparsity "

Get, if necessary generate, the sparsity of a Jacobian block

";

%feature("docstring") casadi::Nullspace::addMonitor "

Add modules to be monitored.

";

%feature("docstring") casadi::Nullspace::dictionary "

Get the dictionary.

";

%feature("docstring") casadi::Nullspace::hessian "

Generate a Hessian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The generated Hessian has two more outputs than the calling function
corresponding to the Hessian and the gradients.

";

%feature("docstring") casadi::Nullspace::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::Nullspace::index_in "

Find the index for a string describing a particular entry of an input
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::Nullspace::setOptionByAllowedIndex "[INTERNAL]  Set a certain option by giving its index into the allowed
values.

";

%feature("docstring") casadi::Nullspace::getAtomicInputReal "

Get the floating point output argument of an atomic operation.

";

%feature("docstring") casadi::Nullspace::evaluate "

Evaluate.

";

%feature("docstring") casadi::Nullspace::getAlgorithmSize "

Get the number of atomic operations.

";

%feature("docstring") casadi::Nullspace::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::Nullspace::getSanitizedName "

get function name with all non alphanumeric characters converted to '_'

";

%feature("docstring") casadi::Nullspace::nnz_out "

Get of number of output nonzeros For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::Nullspace::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::Nullspace::nnz_in "

Get of number of input nonzeros For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::Nullspace::numel_out "

Get of number of output elements For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::Nullspace::callReverse "

Create call to (cached) derivative function, reverse mode.

";

%feature("docstring") casadi::Nullspace::size_in "

Get input dimension.

";

%feature("docstring") casadi::Nullspace::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::Nullspace::description_in "

Get input scheme description by index.

";

%feature("docstring") casadi::Nullspace::callForward "

Create call to (cached) derivative function, forward mode.

";

%feature("docstring") casadi::Nullspace::callDerivative "[INTERNAL]
Evaluate the function symbolically or numerically with directional
derivatives The first two arguments are the nondifferentiated inputs and
results of the evaluation, the next two arguments are a set of forward
directional seeds and the resulting forward directional derivatives, the
length of the vector being the number of forward directions. The next two
arguments are a set of adjoint directional seeds and the resulting adjoint
directional derivatives, the length of the vector being the number of
adjoint directions.

";

%feature("docstring") casadi::Nullspace::description_out "

Get output scheme description by index.

";

%feature("docstring") casadi::Nullspace::sz_iw "[INTERNAL]  Get required
length of iw field.

";

%feature("docstring") casadi::Nullspace::size2_in "

Get input dimension.

";

%feature("docstring") casadi::Nullspace::sparsity_in "

Get sparsity of a given input.

";

%feature("docstring") casadi::Nullspace::setJacobian "

Set the Jacobian function of output oind with respect to input iind NOTE:
Does not take ownership, only weak references to the Jacobians are kept
internally

";

%feature("docstring") casadi::Nullspace::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::Nullspace::tangent "

Generate a tangent function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the input must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::Nullspace::getOptionAllowedIndex "[INTERNAL]
Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::Nullspace::generate "

Export / Generate C code for the function.

";

%feature("docstring") casadi::Nullspace::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::Nullspace::free_sx "

Get all the free variables of the function.

";

%feature("docstring") casadi::Nullspace::setOutput "

>  void IOInterface< Function  >.setOutput(T val, int oind=0)
------------------------------------------------------------------------

Set an output by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.setOutput(T val, str oname)
------------------------------------------------------------------------

Set an output by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::Nullspace::getInput "

>  DMatrix  IOInterface< Function  >.getInput(int iind=0) const
------------------------------------------------------------------------

Get an input by index.

Parameters:
-----------

iind:  index within the range [0..n_in()-1]

>  DMatrix  IOInterface< Function  >.getInput(str iname) const
------------------------------------------------------------------------

Get an input by name.

Parameters:
-----------

iname:  input name. Only allowed when an input scheme is set.

>  void IOInterface< Function  >.getInput(T val, int iind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.getInput(T val, str iname)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::Nullspace::sz_arg "[INTERNAL]  Get required
length of arg field.

";

%feature("docstring") casadi::Nullspace::n_out "

Get the number of function outputs.

";

%feature("docstring") casadi::Nullspace::countNodes "

Number of nodes in the algorithm.

";

%feature("docstring") casadi::Nullspace::mx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::Nullspace::setOptionByEnumValue "[INTERNAL]
Set a certain option by giving an enum value.

";

%feature("docstring") casadi::Nullspace::setJacSparsity "

Generate the sparsity of a Jacobian block

";

%feature("docstring") casadi::Nullspace::n_in "

Get the number of function inputs.

";

%feature("docstring") casadi::Nullspace::setFullJacobian "

Set the Jacobian of all the input nonzeros with respect to all output
nonzeros NOTE: Does not take ownership, only weak references to the Jacobian
are kept internally

";

%feature("docstring") casadi::Nullspace::getAtomicInput "

Get the (integer) input arguments of an atomic operation.

";

%feature("docstring") casadi::Nullspace::checkInputs "[INTERNAL]  Check if
the numerical values of the supplied bounds make sense.

";

%feature("docstring") casadi::Nullspace::spInit "[INTERNAL]  Reset the
sparsity propagation.

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::Nullspace::name "

Name of the function.

";

%feature("docstring") casadi::Nullspace::index_out "

Find the index for a string describing a particular entry of an output
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::Nullspace::default_in "

Get default input value (NOTE: constant reference)

";

%feature("docstring") casadi::Nullspace::getStat "

Get a single statistic obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::Nullspace::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::Nullspace::derivative "

Get a function that calculates nfwd forward derivatives and nadj adjoint
derivatives Legacy function: Use derForward and derReverse instead.

Returns a function with (1+nfwd)*n_in+nadj*n_out inputs and (1+nfwd)*n_out +
nadj*n_in outputs. The first n_in inputs correspond to nondifferentiated
inputs. The next nfwd*n_in inputs correspond to forward seeds, one direction
at a time and the last nadj*n_out inputs correspond to adjoint seeds, one
direction at a time. The first n_out outputs correspond to nondifferentiated
outputs. The next nfwd*n_out outputs correspond to forward sensitivities,
one direction at a time and the last nadj*n_in outputs corresponds to
adjoint sensitivities, one direction at a time.

(n_in = n_in(), n_out = n_out())

";

%feature("docstring") casadi::Nullspace::print "

Print a description of the object.

";

%feature("docstring") casadi::Nullspace::removeMonitor "

Remove modules to be monitored.

";

%feature("docstring") casadi::Nullspace::mx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::Nullspace::printDimensions "

Print dimensions of inputs and outputs.

";

%feature("docstring") casadi::Nullspace "

Base class for nullspace construction.

Constructs a basis for the null-space of a fat matrix A. i.e. finds Z such
that AZ = 0 holds.

The nullspace is also known as the orthogonal complement of the rowspace of
a matrix.

It is assumed that the matrix A is of full rank.

Implementations are not required to construct an orthogonal or orthonormal
basis Joris Gillis

>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| ad_weight    | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | derivative c |              |
|              |              |              | alculation.W |              |
|              |              |              | hen there is |              |
|              |              |              | an option of |              |
|              |              |              | either using |              |
|              |              |              | forward or   |              |
|              |              |              | reverse mode |              |
|              |              |              | directional  |              |
|              |              |              | derivatives, |              |
|              |              |              | the          |              |
|              |              |              | condition ad |              |
|              |              |              | _weight*nf<= |              |
|              |              |              | (1-ad_weight |              |
|              |              |              | )*na is used |              |
|              |              |              | where nf and |              |
|              |              |              | na are       |              |
|              |              |              | estimates of |              |
|              |              |              | the number   |              |
|              |              |              | of forward/r |              |
|              |              |              | everse mode  |              |
|              |              |              | directional  |              |
|              |              |              | derivatives  |              |
|              |              |              | needed. By   |              |
|              |              |              | default,     |              |
|              |              |              | ad_weight is |              |
|              |              |              | calculated a |              |
|              |              |              | utomatically |              |
|              |              |              | , but this   |              |
|              |              |              | can be       |              |
|              |              |              | overridden   |              |
|              |              |              | by setting   |              |
|              |              |              | this option. |              |
|              |              |              | In           |              |
|              |              |              | particular,  |              |
|              |              |              | 0 means      |              |
|              |              |              | forcing      |              |
|              |              |              | forward mode |              |
|              |              |              | and 1        |              |
|              |              |              | forcing      |              |
|              |              |              | reverse      |              |
|              |              |              | mode. Leave  |              |
|              |              |              | unset for    |              |
|              |              |              | (class       |              |
|              |              |              | specific)    |              |
|              |              |              | heuristics.  |              |
+--------------+--------------+--------------+--------------+--------------+
| ad_weight_sp | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | sparsity     |              |
|              |              |              | pattern      |              |
|              |              |              | calculation  |              |
|              |              |              | calculation. |              |
|              |              |              | Overrides    |              |
|              |              |              | default      |              |
|              |              |              | behavior.    |              |
|              |              |              | Set to 0 and |              |
|              |              |              | 1 to force   |              |
|              |              |              | forward and  |              |
|              |              |              | reverse mode |              |
|              |              |              | respectively |              |
|              |              |              | . Cf. option |              |
|              |              |              | \"ad_weight\". |              |
+--------------+--------------+--------------+--------------+--------------+
| compiler     | OT_STRING    | \"clang\"      | Just-in-time | casadi::Func |
|              |              |              | compiler     | tionInternal |
|              |              |              | plugin to be |              |
|              |              |              | used.        |              |
+--------------+--------------+--------------+--------------+--------------+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode     |
|              |              |              | according to |              |
|              |              |              | a given      |              |
|              |              |              | recipe (low- |              |
|              |              |              | level)       |              |
+--------------+--------------+--------------+--------------+--------------+
| dense        | OT_BOOLEAN   | true         | Indicates    | casadi::Null |
|              |              |              | that dense   | spaceInterna |
|              |              |              | matrices can | l            |
|              |              |              | be assumed   |              |
+--------------+--------------+--------------+--------------+--------------+
| gather_stats | OT_BOOLEAN   | false        | Flag to      | casadi::Func |
|              |              |              | indicate     | tionInternal |
|              |              |              | whether      |              |
|              |              |              | statistics   |              |
|              |              |              | must be      |              |
|              |              |              | gathered     |              |
+--------------+--------------+--------------+--------------+--------------+
| input_scheme | OT_STRINGVEC | GenericType( | Custom input | casadi::Func |
|              | TOR          | )            | scheme       | tionInternal |
+--------------+--------------+--------------+--------------+--------------+
| inputs_check | OT_BOOLEAN   | true         | Throw        | casadi::Func |
|              |              |              | exceptions   | tionInternal |
|              |              |              | when the     |              |
|              |              |              | numerical    |              |
|              |              |              | values of    |              |
|              |              |              | the inputs   |              |
|              |              |              | don't make   |              |
|              |              |              | sense        |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_penalty  | OT_REAL      | 2            | When         | casadi::Func |
|              |              |              | requested    | tionInternal |
|              |              |              | for a number |              |
|              |              |              | of forward/r |              |
|              |              |              | everse       |              |
|              |              |              | directions,  |              |
|              |              |              | it may be    |              |
|              |              |              | cheaper to   |              |
|              |              |              | compute      |              |
|              |              |              | first the    |              |
|              |              |              | full         |              |
|              |              |              | jacobian and |              |
|              |              |              | then         |              |
|              |              |              | multiply     |              |
|              |              |              | with seeds,  |              |
|              |              |              | rather than  |              |
|              |              |              | obtain the   |              |
|              |              |              | requested    |              |
|              |              |              | directions   |              |
|              |              |              | in a straigh |              |
|              |              |              | tforward     |              |
|              |              |              | manner.      |              |
|              |              |              | Casadi uses  |              |
|              |              |              | a heuristic  |              |
|              |              |              | to decide    |              |
|              |              |              | which is     |              |
|              |              |              | cheaper. A   |              |
|              |              |              | high value   |              |
|              |              |              | of 'jac_pena |              |
|              |              |              | lty' makes   |              |
|              |              |              | it less      |              |
|              |              |              | likely for   |              |
|              |              |              | the heurstic |              |
|              |              |              | to chose the |              |
|              |              |              | full         |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy.    |              |
|              |              |              | The special  |              |
|              |              |              | value -1     |              |
|              |              |              | indicates    |              |
|              |              |              | never to use |              |
|              |              |              | the full     |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy     |              |
+--------------+--------------+--------------+--------------+--------------+
| jit          | OT_BOOLEAN   | false        | Use just-in- | casadi::Func |
|              |              |              | time         | tionInternal |
|              |              |              | compiler to  |              |
|              |              |              | speed up the |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| jit_options  | OT_DICT      | GenericType( | Options to   | casadi::Func |
|              |              | )            | be passed to | tionInternal |
|              |              |              | the jit      |              |
|              |              |              | compiler.    |              |
+--------------+--------------+--------------+--------------+--------------+
| monitor      | OT_STRINGVEC | GenericType( | Monitors to  | casadi::Func |
|              | TOR          | )            | be activated | tionInternal |
|              |              |              | (inputs|outp |              |
|              |              |              | uts)         |              |
+--------------+--------------+--------------+--------------+--------------+
| output_schem | OT_STRINGVEC | GenericType( | Custom       | casadi::Func |
| e            | TOR          | )            | output       | tionInternal |
|              |              |              | scheme       |              |
+--------------+--------------+--------------+--------------+--------------+
| regularity_c | OT_BOOLEAN   | true         | Throw        | casadi::Func |
| heck         |              |              | exceptions   | tionInternal |
|              |              |              | when NaN or  |              |
|              |              |              | Inf appears  |              |
|              |              |              | during       |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| user_data    | OT_VOIDPTR   | GenericType( | A user-      | casadi::Func |
|              |              | )            | defined      | tionInternal |
|              |              |              | field that   |              |
|              |              |              | can be used  |              |
|              |              |              | to identify  |              |
|              |              |              | the function |              |
|              |              |              | or pass      |              |
|              |              |              | additional   |              |
|              |              |              | information  |              |
+--------------+--------------+--------------+--------------+--------------+
| verbose      | OT_BOOLEAN   | false        | Verbose      | casadi::Func |
|              |              |              | evaluation   | tionInternal |
|              |              |              | for          |              |
|              |              |              | debugging    |              |
+--------------+--------------+--------------+--------------+--------------+

Diagrams
--------



C++ includes: nullspace.hpp ";

%feature("docstring") casadi::Nullspace::type_name "

Get type name.

";

%feature("docstring") casadi::Nullspace::sz_w "[INTERNAL]  Get required
length of w field.

";

%feature("docstring") casadi::Nullspace::name_in "

>  [str] Function.name_in() const 
------------------------------------------------------------------------

Get input scheme.

>  str Function.name_in(int ind) const 
------------------------------------------------------------------------

Get input scheme name by index.

";

%feature("docstring") casadi::Nullspace::getAtomicOperation "

Get an atomic operation operator index.

";

%feature("docstring") casadi::Nullspace::size_out "

Get output dimension.

";

%feature("docstring") casadi::Nullspace::getStats "

Get all statistics obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::Nullspace::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::Nullspace::isNull "

Is a null pointer?

";

%feature("docstring") casadi::Nullspace::sparsity_out "

Get sparsity of a given output.

";

%feature("docstring") casadi::Nullspace::fullJacobian "

Generate a Jacobian function of all the inputs elements with respect to all
the output elements).

";

%feature("docstring") casadi::Nullspace::getOptionEnumValue "[INTERNAL]
Get the enum value corresponding to th certain option.

";

%feature("docstring") casadi::Nullspace::numel_in "

Get of number of input elements For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::Nullspace::setDerReverse "

Set a function that calculates nadj adjoint derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::Nullspace::free_mx "

Get all the free variables of the function.

";

%feature("docstring") casadi::Nullspace::spEvaluate "[INTERNAL]  Propagate
the sparsity pattern through a set of directional.

derivatives forward or backward (for usage, see the example
propagating_sparsity.cpp)

";

%feature("docstring") casadi::Nullspace::repr "

Print a representation of the object.

";

%feature("docstring") casadi::Nullspace::Nullspace "

Default constructor.

";

%feature("docstring") casadi::Nullspace::mapaccum "

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
  



";

%feature("docstring") casadi::Nullspace::derReverse "

Get a function that calculates nadj adjoint derivatives.

Returns a function with n_in + n_out +nadj*n_out inputs and nadj*n_in
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nadj*n_out inputs correspond to adjoint seeds, one direction at a time The
nadj*n_in outputs correspond to adjoint sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

(n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::Nullspace::generateLiftingFunctions "

Extract the functions needed for the Lifted Newton method.

";

%feature("docstring") casadi::Nullspace::sx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::Nullspace::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::Nullspace::is_a "

Check if the function is of a particular type Optionally check if name
matches one of the base classes (default true)

";

%feature("docstring") casadi::Nullspace::map "

>  [[MX] ] Function.map([[MX ] ] arg, str parallelization=\"serial\")

>  [MX] Function.map([MX ] arg, str parallelization=\"serial\")
------------------------------------------------------------------------

Evaluate symbolically in parallel (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

>  Function Function.map(str name, int N, Dict options=Dict()) const 
------------------------------------------------------------------------

Create a mapped version of this function.

Suppose the function has a signature of:

::

     f: (a, p) -> ( s )
  



The the mapaccumulated version has the signature:

::

     F: (A, P) -> (S )
  
      with
          a: horzcat([a0, a1, ..., a_(N-1)])
          p: horzcat([p0, p1, ..., p_(N-1)])
          s: horzcat([s0, s1, ..., s_(N-1)])
      and
          s0 <- f(a0, p0)
          s1 <- f(a1, p1)
          ...
          s_(N-1) <- f(a_(N-1), p_(N-1))
  



>  Function Function.map(str name, int n, [bool ] repeat_in, [bool ] repeat_out, Dict opts=Dict()) const 
------------------------------------------------------------------------

Generic map.

";

%feature("docstring") casadi::Nullspace::getOptionTypeName "

Get the type name of a certain option.

";

%feature("docstring") casadi::Nullspace::setDerForward "

Set a function that calculates nfwd forward derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::Nullspace::size1_out "

Get output dimension.

";


// File: classcasadi_1_1OldCollocationIntegrator.xml


// File: classcasadi_1_1OneSX.xml


// File: classcasadi_1_1OptionsFunctionality.xml


/*  Option Functionality  */ %feature("docstring")
casadi::OptionsFunctionality::setOptionByAllowedIndex " [INTERNAL]  Set a
certain option by giving its index into the allowed values.

";

%feature("docstring") casadi::OptionsFunctionality::print "

Print a description of the object.

";

%feature("docstring") casadi::OptionsFunctionality::getOptionAllowedIndex "[INTERNAL]  Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::OptionsFunctionality::dictionary "

Get the dictionary.

";

%feature("docstring") casadi::OptionsFunctionality::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::OptionsFunctionality::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::OptionsFunctionality::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::OptionsFunctionality::printPtr "[INTERNAL]
Print the pointer to the internal class

";

%feature("docstring") casadi::OptionsFunctionality::OptionsFunctionality "

Default constructor.

";

%feature("docstring") casadi::OptionsFunctionality::isNull "

Is a null pointer?

";

%feature("docstring") casadi::OptionsFunctionality::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::OptionsFunctionality "

Provides options setting/getting functionality.

Gives a derived class the ability to set and retrieve options in a
convenient way. It also contains error checking, making sure that the option
exists and that the value type is correct.

A derived class should add option names, types and default values to the
corresponding vectors.

Joel Andersson

>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode     |
|              |              |              | according to |              |
|              |              |              | a given      |              |
|              |              |              | recipe (low- |              |
|              |              |              | level)       |              |
+--------------+--------------+--------------+--------------+--------------+

Diagrams
--------



C++ includes: options_functionality.hpp ";

%feature("docstring") casadi::OptionsFunctionality::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::OptionsFunctionality::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::OptionsFunctionality::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::OptionsFunctionality::getOptionEnumValue "[INTERNAL]  Get the enum value corresponding to th certain option.

";

%feature("docstring") casadi::OptionsFunctionality::setOptionByEnumValue "[INTERNAL]  Set a certain option by giving an enum value.

";

%feature("docstring") casadi::OptionsFunctionality::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::OptionsFunctionality::repr "

Print a representation of the object.

";

%feature("docstring") casadi::OptionsFunctionality::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::OptionsFunctionality::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::OptionsFunctionality::getOptionTypeName "

Get the type name of a certain option.

";

%feature("docstring") casadi::OptionsFunctionality::~OptionsFunctionality "

Destructor.

";


// File: classcasadi_1_1Polynomial.xml
%feature("docstring") casadi::Polynomial "

Helper class for differentiating and integrating polynomials.

Joel Andersson

C++ includes: polynomial.hpp ";

%feature("docstring") casadi::Polynomial::derivative "

Create a new polynomial for the derivative.

";

%feature("docstring") casadi::Polynomial::Polynomial "

>  Polynomial(real_t scalar=1)
------------------------------------------------------------------------

Construct a constant polynomial.

>  Polynomial(real_t p0, real_t p1)
------------------------------------------------------------------------

Construct a linear polynomial.

>  Polynomial(real_t p0, real_t p1, real_t p2)
------------------------------------------------------------------------

Construct a quadratic polynomial.

>  Polynomial(real_t p0, real_t p1, real_t p2, real_t p3)
------------------------------------------------------------------------

Construct a cubic polynomial.

>  Polynomial([T ] coeff)
------------------------------------------------------------------------

Construct from a vector of polynomial coefficients.

";

%feature("docstring") casadi::Polynomial::print "

Print a description of the object.

";

%feature("docstring") casadi::Polynomial::anti_derivative "

Create a new polynomial for the anti-derivative (primitive function)

";

%feature("docstring") casadi::Polynomial::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::Polynomial::degree "

Degree of the polynomial.

";

%feature("docstring") casadi::Polynomial::toScalar "

Get scalar value (error if degree()!=0)

";

%feature("docstring") casadi::Polynomial::trim "

Remove excess zeros.

";

%feature("docstring") casadi::Polynomial::repr "

Print a representation of the object.

";

%feature("docstring") casadi::Polynomial::getRepresentation "

Return a string with a representation (for SWIG)

";


// File: classcasadi_1_1PrintableObject.xml
%feature("docstring") casadi::PrintableObject::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") friendwrap_repr "

Return a string with a representation of the object, cf. repr(Object) in
Python.

";

%feature("docstring") friendwrap_str "

Return a string with a description of the object, cf. str(Object) in Python.

";

%feature("docstring") casadi::PrintableObject "

Base class for objects that have a natural string representation.

Joel Andersson

C++ includes: printable_object.hpp ";

%feature("docstring") casadi::PrintableObject::getDescription "

Return a string with a description (for SWIG)

";


// File: classcasadi_1_1Project.xml


// File: classcasadi_1_1QpSolver.xml


/*  Simple Getters & Setters  */

/*  Advanced Getters  */

/*  Option Functionality  */ %feature("docstring")
casadi::QpSolver::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::QpSolver::printDimensions "

Print dimensions of inputs and outputs.

";

%feature("docstring") casadi::QpSolver::generateLiftingFunctions "

Extract the functions needed for the Lifted Newton method.

";

%feature("docstring") casadi::QpSolver::setOutput "

>  void IOInterface< Function  >.setOutput(T val, int oind=0)
------------------------------------------------------------------------

Set an output by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.setOutput(T val, str oname)
------------------------------------------------------------------------

Set an output by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::QpSolver::jacSparsity "

Get, if necessary generate, the sparsity of a Jacobian block

";

%feature("docstring") casadi::QpSolver::derReverse "

Get a function that calculates nadj adjoint derivatives.

Returns a function with n_in + n_out +nadj*n_out inputs and nadj*n_in
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nadj*n_out inputs correspond to adjoint seeds, one direction at a time The
nadj*n_in outputs correspond to adjoint sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

(n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::QpSolver::mx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::QpSolver::generateNativeCode "

Generate native code in the interfaced language for debugging

";

%feature("docstring") casadi::QpSolver::spEvaluate "[INTERNAL]  Propagate
the sparsity pattern through a set of directional.

derivatives forward or backward (for usage, see the example
propagating_sparsity.cpp)

";

%feature("docstring") casadi::QpSolver::size1_out "

Get output dimension.

";

%feature("docstring") casadi::QpSolver::numel_out "

Get of number of output elements For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::QpSolver::spCanEvaluate "[INTERNAL]  Is the
class able to propagate seeds through the algorithm?

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::QpSolver::jacobian "

Generate a Jacobian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. If
compact is set to true, only the nonzeros of the input and output
expressions are considered. If symmetric is set to true, the Jacobian being
calculated is known to be symmetric (usually a Hessian), which can be
exploited by the algorithm.

The generated Jacobian has one more output than the calling function
corresponding to the Jacobian and the same number of inputs.

";

%feature("docstring") casadi::QpSolver::mx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::QpSolver::size1_in "

Get input dimension.

";

%feature("docstring") casadi::QpSolver::sz_iw "[INTERNAL]  Get required
length of iw field.

";

%feature("docstring") casadi::QpSolver::generate "

Export / Generate C code for the function.

";

%feature("docstring") casadi::QpSolver::getStat "

Get a single statistic obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::QpSolver::setInput "

>  void IOInterface< Function  >.setInput(T val, int iind=0)
------------------------------------------------------------------------

Set an input by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.setInput(T val, str iname)
------------------------------------------------------------------------

Set an input by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::QpSolver::getAtomicInput "

Get the (integer) input arguments of an atomic operation.

";

%feature("docstring") casadi::QpSolver::isNull "

Is a null pointer?

";

%feature("docstring") casadi::QpSolver::callForward "

Create call to (cached) derivative function, forward mode.

";

%feature("docstring") casadi::QpSolver::name "

Name of the function.

";

%feature("docstring") casadi::QpSolver::size2_in "

Get input dimension.

";

%feature("docstring") casadi::QpSolver::sparsity_out "

Get sparsity of a given output.

";

%feature("docstring") casadi::QpSolver::getWorkSize "

Get the length of the work vector.

";

%feature("docstring") casadi::QpSolver::dictionary "

Get the dictionary.

";

%feature("docstring") casadi::QpSolver::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::QpSolver::setOptionByAllowedIndex "[INTERNAL]
Set a certain option by giving its index into the allowed values.

";

%feature("docstring") casadi::QpSolver::getOptionEnumValue "[INTERNAL]  Get
the enum value corresponding to th certain option.

";

%feature("docstring") casadi::QpSolver::QpSolver "

>  QpSolver()
------------------------------------------------------------------------

Default constructor.

>  QpSolver(str name, str solver, const std.map< str, Sparsity > &st, Dict opts=Dict())
------------------------------------------------------------------------

Constructor (new syntax, includes initialization)

Parameters:
-----------

name:

Name of a solver. It might be one of:

- cplex

- ooqp

- qpoases

- sqic

- nlp

Note: some of the plugins in this list might not be available on your
system. Also, there might be extra plugins available to you that are not
listed here. You can obtain their documentation with
QpSolver.doc(\"myextraplugin\")

st:

Problem structure.>Struct scheme: casadi::QPStruct ( = 2) []

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| QP_STRUCT_H            |                        | The square matrix H:   |
|                        |                        | sparse, (n x n). Only  |
|                        |                        | the lower triangular   |
|                        |                        | part is actually used. |
|                        |                        | The matrix is assumed  |
|                        |                        | to be symmetrical.     |
+------------------------+------------------------+------------------------+
| QP_STRUCT_A            |                        | The matrix A: sparse,  |
|                        |                        | (nc x n) - product     |
|                        |                        | with x must be dense.  |
+------------------------+------------------------+------------------------+

";

%feature("docstring") casadi::QpSolver::sz_res "[INTERNAL]  Get required
length of res field.

";

%feature("docstring") casadi::QpSolver::name_out "

>  [str] Function.name_out() const 
------------------------------------------------------------------------

Get output scheme.

>  str Function.name_out(int ind) const 
------------------------------------------------------------------------

Get output scheme name by index.

";

%feature("docstring") casadi::QpSolver::sparsity_in "

Get sparsity of a given input.

";

%feature("docstring") casadi::QpSolver::gradient "

Generate a gradient function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the output must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::QpSolver::is_a "

Check if the function is of a particular type Optionally check if name
matches one of the base classes (default true)

";

%feature("docstring") casadi::QpSolver::setFullJacobian "

Set the Jacobian of all the input nonzeros with respect to all output
nonzeros NOTE: Does not take ownership, only weak references to the Jacobian
are kept internally

";

%feature("docstring") casadi::QpSolver::index_out "

Find the index for a string describing a particular entry of an output
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::QpSolver::name_in "

>  [str] Function.name_in() const 
------------------------------------------------------------------------

Get input scheme.

>  str Function.name_in(int ind) const 
------------------------------------------------------------------------

Get input scheme name by index.

";

%feature("docstring") casadi::QpSolver::getOptionTypeName "

Get the type name of a certain option.

";

%feature("docstring") casadi::QpSolver::getSanitizedName "

get function name with all non alphanumeric characters converted to '_'

";

%feature("docstring") casadi::QpSolver::description_out "

Get output scheme description by index.

";

%feature("docstring") casadi::QpSolver::type_name "

Get type name.

";

%feature("docstring") casadi::QpSolver "

QpSolver.

Solves the following strictly convex problem:



::

  min          1/2 x' H x + g' x
   x
  
  subject to
              LBA <= A x <= UBA
              LBX <= x   <= UBX
  
      with :
        H sparse (n x n) positive definite
        g dense  (n x 1)
  
      n: number of decision variables (x)
      nc: number of constraints (A)



If H is not positive-definite, the solver should throw an error.

General information
===================



>Input scheme: casadi::QpSolverInput (QP_SOLVER_NUM_IN = 9) [qpIn]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| QP_SOLVER_H            | h                      | The square matrix H:   |
|                        |                        | sparse, (n x n). Only  |
|                        |                        | the lower triangular   |
|                        |                        | part is actually used. |
|                        |                        | The matrix is assumed  |
|                        |                        | to be symmetrical.     |
+------------------------+------------------------+------------------------+
| QP_SOLVER_G            | g                      | The vector g: dense,   |
|                        |                        | (n x 1) .              |
+------------------------+------------------------+------------------------+
| QP_SOLVER_A            | a                      | The matrix A: sparse,  |
|                        |                        | (nc x n) - product     |
|                        |                        | with x must be dense.  |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+
| QP_SOLVER_LBA          | lba                    | dense, (nc x 1)        |
+------------------------+------------------------+------------------------+
| QP_SOLVER_UBA          | uba                    | dense, (nc x 1)        |
+------------------------+------------------------+------------------------+
| QP_SOLVER_LBX          | lbx                    | dense, (n x 1)         |
+------------------------+------------------------+------------------------+
| QP_SOLVER_UBX          | ubx                    | dense, (n x 1)         |
+------------------------+------------------------+------------------------+
| QP_SOLVER_X0           | x0                     | dense, (n x 1)         |
+------------------------+------------------------+------------------------+
| QP_SOLVER_LAM_X0       | lam_x0                 | dense                  |
+------------------------+------------------------+------------------------+

>Output scheme: casadi::QpSolverOutput (QP_SOLVER_NUM_OUT = 4) [qpOut]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| QP_SOLVER_X            | x                      | The primal solution .  |
+------------------------+------------------------+------------------------+
| QP_SOLVER_COST         | cost                   | The optimal cost .     |
+------------------------+------------------------+------------------------+
| QP_SOLVER_LAM_A        | lam_a                  | The dual solution      |
|                        |                        | corresponding to       |
|                        |                        | linear bounds .        |
+------------------------+------------------------+------------------------+
| QP_SOLVER_LAM_X        | lam_x                  | The dual solution      |
|                        |                        | corresponding to       |
|                        |                        | simple bounds .        |
+------------------------+------------------------+------------------------+

>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| ad_weight    | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | derivative c |              |
|              |              |              | alculation.W |              |
|              |              |              | hen there is |              |
|              |              |              | an option of |              |
|              |              |              | either using |              |
|              |              |              | forward or   |              |
|              |              |              | reverse mode |              |
|              |              |              | directional  |              |
|              |              |              | derivatives, |              |
|              |              |              | the          |              |
|              |              |              | condition ad |              |
|              |              |              | _weight*nf<= |              |
|              |              |              | (1-ad_weight |              |
|              |              |              | )*na is used |              |
|              |              |              | where nf and |              |
|              |              |              | na are       |              |
|              |              |              | estimates of |              |
|              |              |              | the number   |              |
|              |              |              | of forward/r |              |
|              |              |              | everse mode  |              |
|              |              |              | directional  |              |
|              |              |              | derivatives  |              |
|              |              |              | needed. By   |              |
|              |              |              | default,     |              |
|              |              |              | ad_weight is |              |
|              |              |              | calculated a |              |
|              |              |              | utomatically |              |
|              |              |              | , but this   |              |
|              |              |              | can be       |              |
|              |              |              | overridden   |              |
|              |              |              | by setting   |              |
|              |              |              | this option. |              |
|              |              |              | In           |              |
|              |              |              | particular,  |              |
|              |              |              | 0 means      |              |
|              |              |              | forcing      |              |
|              |              |              | forward mode |              |
|              |              |              | and 1        |              |
|              |              |              | forcing      |              |
|              |              |              | reverse      |              |
|              |              |              | mode. Leave  |              |
|              |              |              | unset for    |              |
|              |              |              | (class       |              |
|              |              |              | specific)    |              |
|              |              |              | heuristics.  |              |
+--------------+--------------+--------------+--------------+--------------+
| ad_weight_sp | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | sparsity     |              |
|              |              |              | pattern      |              |
|              |              |              | calculation  |              |
|              |              |              | calculation. |              |
|              |              |              | Overrides    |              |
|              |              |              | default      |              |
|              |              |              | behavior.    |              |
|              |              |              | Set to 0 and |              |
|              |              |              | 1 to force   |              |
|              |              |              | forward and  |              |
|              |              |              | reverse mode |              |
|              |              |              | respectively |              |
|              |              |              | . Cf. option |              |
|              |              |              | \"ad_weight\". |              |
+--------------+--------------+--------------+--------------+--------------+
| compiler     | OT_STRING    | \"clang\"      | Just-in-time | casadi::Func |
|              |              |              | compiler     | tionInternal |
|              |              |              | plugin to be |              |
|              |              |              | used.        |              |
+--------------+--------------+--------------+--------------+--------------+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode   c |
|              |              |              | according to | asadi::QpSol |
|              |              |              | a given      | verInternal  |
|              |              |              | recipe (low- |              |
|              |              |              | level)  (lp) |              |
+--------------+--------------+--------------+--------------+--------------+
| gather_stats | OT_BOOLEAN   | false        | Flag to      | casadi::Func |
|              |              |              | indicate     | tionInternal |
|              |              |              | whether      |              |
|              |              |              | statistics   |              |
|              |              |              | must be      |              |
|              |              |              | gathered     |              |
+--------------+--------------+--------------+--------------+--------------+
| input_scheme | OT_STRINGVEC | GenericType( | Custom input | casadi::Func |
|              | TOR          | )            | scheme       | tionInternal |
+--------------+--------------+--------------+--------------+--------------+
| inputs_check | OT_BOOLEAN   | true         | Throw        | casadi::Func |
|              |              |              | exceptions   | tionInternal |
|              |              |              | when the     |              |
|              |              |              | numerical    |              |
|              |              |              | values of    |              |
|              |              |              | the inputs   |              |
|              |              |              | don't make   |              |
|              |              |              | sense        |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_penalty  | OT_REAL      | 2            | When         | casadi::Func |
|              |              |              | requested    | tionInternal |
|              |              |              | for a number |              |
|              |              |              | of forward/r |              |
|              |              |              | everse       |              |
|              |              |              | directions,  |              |
|              |              |              | it may be    |              |
|              |              |              | cheaper to   |              |
|              |              |              | compute      |              |
|              |              |              | first the    |              |
|              |              |              | full         |              |
|              |              |              | jacobian and |              |
|              |              |              | then         |              |
|              |              |              | multiply     |              |
|              |              |              | with seeds,  |              |
|              |              |              | rather than  |              |
|              |              |              | obtain the   |              |
|              |              |              | requested    |              |
|              |              |              | directions   |              |
|              |              |              | in a straigh |              |
|              |              |              | tforward     |              |
|              |              |              | manner.      |              |
|              |              |              | Casadi uses  |              |
|              |              |              | a heuristic  |              |
|              |              |              | to decide    |              |
|              |              |              | which is     |              |
|              |              |              | cheaper. A   |              |
|              |              |              | high value   |              |
|              |              |              | of 'jac_pena |              |
|              |              |              | lty' makes   |              |
|              |              |              | it less      |              |
|              |              |              | likely for   |              |
|              |              |              | the heurstic |              |
|              |              |              | to chose the |              |
|              |              |              | full         |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy.    |              |
|              |              |              | The special  |              |
|              |              |              | value -1     |              |
|              |              |              | indicates    |              |
|              |              |              | never to use |              |
|              |              |              | the full     |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy     |              |
+--------------+--------------+--------------+--------------+--------------+
| jit          | OT_BOOLEAN   | false        | Use just-in- | casadi::Func |
|              |              |              | time         | tionInternal |
|              |              |              | compiler to  |              |
|              |              |              | speed up the |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| jit_options  | OT_DICT      | GenericType( | Options to   | casadi::Func |
|              |              | )            | be passed to | tionInternal |
|              |              |              | the jit      |              |
|              |              |              | compiler.    |              |
+--------------+--------------+--------------+--------------+--------------+
| monitor      | OT_STRINGVEC | GenericType( | Monitors to  | casadi::Func |
|              | TOR          | )            | be activated | tionInternal |
|              |              |              | (inputs|outp |              |
|              |              |              | uts)         |              |
+--------------+--------------+--------------+--------------+--------------+
| output_schem | OT_STRINGVEC | GenericType( | Custom       | casadi::Func |
| e            | TOR          | )            | output       | tionInternal |
|              |              |              | scheme       |              |
+--------------+--------------+--------------+--------------+--------------+
| regularity_c | OT_BOOLEAN   | true         | Throw        | casadi::Func |
| heck         |              |              | exceptions   | tionInternal |
|              |              |              | when NaN or  |              |
|              |              |              | Inf appears  |              |
|              |              |              | during       |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| user_data    | OT_VOIDPTR   | GenericType( | A user-      | casadi::Func |
|              |              | )            | defined      | tionInternal |
|              |              |              | field that   |              |
|              |              |              | can be used  |              |
|              |              |              | to identify  |              |
|              |              |              | the function |              |
|              |              |              | or pass      |              |
|              |              |              | additional   |              |
|              |              |              | information  |              |
+--------------+--------------+--------------+--------------+--------------+
| verbose      | OT_BOOLEAN   | false        | Verbose      | casadi::Func |
|              |              |              | evaluation   | tionInternal |
|              |              |              | for          |              |
|              |              |              | debugging    |              |
+--------------+--------------+--------------+--------------+--------------+

List of plugins
===============



- cplex

- ooqp

- qpoases

- sqic

- nlp

Note: some of the plugins in this list might not be available on your
system. Also, there might be extra plugins available to you that are not
listed here. You can obtain their documentation with
QpSolver.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

cplex
-----



Interface to Cplex solver for sparse Quadratic Programs

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| barrier_maxiter | OT_INTEGER      | 2.100e+09       | Maximum number  |
|                 |                 |                 | of barrier      |
|                 |                 |                 | iterations.     |
+-----------------+-----------------+-----------------+-----------------+
| convex          | OT_BOOLEAN      | true            | Indicates if    |
|                 |                 |                 | the QP is       |
|                 |                 |                 | convex or not   |
|                 |                 |                 | (affects only   |
|                 |                 |                 | the barrier     |
|                 |                 |                 | method).        |
+-----------------+-----------------+-----------------+-----------------+
| dep_check       | OT_STRING       | \"off\"           | Detect          |
|                 |                 |                 | redundant       |
|                 |                 |                 | constraints. (a |
|                 |                 |                 | utomatic:-1|off |
|                 |                 |                 | :0|begin:1|end: |
|                 |                 |                 | 2|both:3)       |
+-----------------+-----------------+-----------------+-----------------+
| dump_filename   | OT_STRING       | \"qp.dat\"        | The filename to |
|                 |                 |                 | dump to.        |
+-----------------+-----------------+-----------------+-----------------+
| dump_to_file    | OT_BOOLEAN      | false           | Dumps QP to     |
|                 |                 |                 | file in CPLEX   |
|                 |                 |                 | format.         |
+-----------------+-----------------+-----------------+-----------------+
| qp_method       | OT_STRING       | \"automatic\"     | Determines      |
|                 |                 |                 | which CPLEX     |
|                 |                 |                 | algorithm to    |
|                 |                 |                 | use. (automatic |
|                 |                 |                 | |primal_simplex |
|                 |                 |                 | |dual_simplex|n |
|                 |                 |                 | etwork|barrier| |
|                 |                 |                 | sifting|concurr |
|                 |                 |                 | ent|crossover)  |
+-----------------+-----------------+-----------------+-----------------+
| simplex_maxiter | OT_INTEGER      | 2.100e+09       | Maximum number  |
|                 |                 |                 | of simplex      |
|                 |                 |                 | iterations.     |
+-----------------+-----------------+-----------------+-----------------+
| tol             | OT_REAL         | 0.000           | Tolerance of    |
|                 |                 |                 | solver          |
+-----------------+-----------------+-----------------+-----------------+
| warm_start      | OT_BOOLEAN      | false           | Use warm start  |
|                 |                 |                 | with simplex    |
|                 |                 |                 | methods         |
|                 |                 |                 | (affects only   |
|                 |                 |                 | the simplex     |
|                 |                 |                 | methods).       |
+-----------------+-----------------+-----------------+-----------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

ooqp
----



Interface to the OOQP Solver for quadratic programming The current
implementation assumes that OOQP is configured with the MA27 sparse linear
solver.

NOTE: when doing multiple calls to evaluate(), check if you need to
reInit();

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| artol           | OT_REAL         | 0.000           | tolerance as    |
|                 |                 |                 | provided with   |
|                 |                 |                 | setArTol to     |
|                 |                 |                 | OOQP            |
+-----------------+-----------------+-----------------+-----------------+
| mutol           | OT_REAL         | 0.000           | tolerance as    |
|                 |                 |                 | provided with   |
|                 |                 |                 | setMuTol to     |
|                 |                 |                 | OOQP            |
+-----------------+-----------------+-----------------+-----------------+
| print_level     | OT_INTEGER      | 0               | Print level.    |
|                 |                 |                 | OOQP listens to |
|                 |                 |                 | print_level 0,  |
|                 |                 |                 | 10 and 100      |
+-----------------+-----------------+-----------------+-----------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

qpoases
-------



Interface to QPOases Solver for quadratic programming

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| CPUtime         | OT_REAL         | None            | The maximum     |
|                 |                 |                 | allowed CPU     |
|                 |                 |                 | time in seconds |
|                 |                 |                 | for the whole   |
|                 |                 |                 | initialisation  |
|                 |                 |                 | (and the        |
|                 |                 |                 | actually        |
|                 |                 |                 | required one on |
|                 |                 |                 | output).        |
|                 |                 |                 | Disabled if     |
|                 |                 |                 | unset.          |
+-----------------+-----------------+-----------------+-----------------+
| boundRelaxation | OT_REAL         | 10000           | Initial         |
|                 |                 |                 | relaxation of   |
|                 |                 |                 | bounds to start |
|                 |                 |                 | homotopy and    |
|                 |                 |                 | initial value   |
|                 |                 |                 | for far bounds. |
+-----------------+-----------------+-----------------+-----------------+
| boundTolerance  | OT_REAL         | 0.000           | If upper and    |
|                 |                 |                 | lower bounds    |
|                 |                 |                 | differ less     |
|                 |                 |                 | than this       |
|                 |                 |                 | tolerance, they |
|                 |                 |                 | are regarded    |
|                 |                 |                 | equal, i.e. as  |
|                 |                 |                 | equality        |
|                 |                 |                 | constraint.     |
+-----------------+-----------------+-----------------+-----------------+
| enableCholeskyR | OT_INTEGER      | 0               | Specifies the   |
| efactorisation  |                 |                 | frequency of a  |
|                 |                 |                 | full re-        |
|                 |                 |                 | factorisation   |
|                 |                 |                 | of projected    |
|                 |                 |                 | Hessian matrix: |
|                 |                 |                 | 0: turns them   |
|                 |                 |                 | off, 1: uses    |
|                 |                 |                 | them at each    |
|                 |                 |                 | iteration etc.  |
+-----------------+-----------------+-----------------+-----------------+
| enableDriftCorr | OT_INTEGER      | 1               | Specifies the   |
| ection          |                 |                 | frequency of    |
|                 |                 |                 | drift           |
|                 |                 |                 | corrections: 0: |
|                 |                 |                 | turns them off. |
+-----------------+-----------------+-----------------+-----------------+
| enableEqualitie | OT_BOOLEAN      | False           | Specifies       |
| s               |                 |                 | whether         |
|                 |                 |                 | equalities      |
|                 |                 |                 | should be       |
|                 |                 |                 | treated as      |
|                 |                 |                 | always active   |
|                 |                 |                 | (True) or not   |
|                 |                 |                 | (False)         |
+-----------------+-----------------+-----------------+-----------------+
| enableFarBounds | OT_BOOLEAN      | True            | Enables the use |
|                 |                 |                 | of far bounds.  |
+-----------------+-----------------+-----------------+-----------------+
| enableFlippingB | OT_BOOLEAN      | True            | Enables the use |
| ounds           |                 |                 | of flipping     |
|                 |                 |                 | bounds.         |
+-----------------+-----------------+-----------------+-----------------+
| enableFullLITes | OT_BOOLEAN      | False           | Enables         |
| ts              |                 |                 | condition-      |
|                 |                 |                 | hardened (but   |
|                 |                 |                 | more expensive) |
|                 |                 |                 | LI test.        |
+-----------------+-----------------+-----------------+-----------------+
| enableNZCTests  | OT_BOOLEAN      | True            | Enables nonzero |
|                 |                 |                 | curvature       |
|                 |                 |                 | tests.          |
+-----------------+-----------------+-----------------+-----------------+
| enableRamping   | OT_BOOLEAN      | True            | Enables         |
|                 |                 |                 | ramping.        |
+-----------------+-----------------+-----------------+-----------------+
| enableRegularis | OT_BOOLEAN      | False           | Enables         |
| ation           |                 |                 | automatic       |
|                 |                 |                 | Hessian         |
|                 |                 |                 | regularisation. |
+-----------------+-----------------+-----------------+-----------------+
| epsDen          | OT_REAL         | 0.000           | Denominator     |
|                 |                 |                 | tolerance for   |
|                 |                 |                 | ratio tests.    |
+-----------------+-----------------+-----------------+-----------------+
| epsFlipping     | OT_REAL         | 0.000           | Tolerance of    |
|                 |                 |                 | squared         |
|                 |                 |                 | Cholesky        |
|                 |                 |                 | diagonal factor |
|                 |                 |                 | which triggers  |
|                 |                 |                 | flipping bound. |
+-----------------+-----------------+-----------------+-----------------+
| epsIterRef      | OT_REAL         | 0.000           | Early           |
|                 |                 |                 | termination     |
|                 |                 |                 | tolerance for   |
|                 |                 |                 | iterative       |
|                 |                 |                 | refinement.     |
+-----------------+-----------------+-----------------+-----------------+
| epsLITests      | OT_REAL         | 0.000           | Tolerance for   |
|                 |                 |                 | linear          |
|                 |                 |                 | independence    |
|                 |                 |                 | tests.          |
+-----------------+-----------------+-----------------+-----------------+
| epsNZCTests     | OT_REAL         | 0.000           | Tolerance for   |
|                 |                 |                 | nonzero         |
|                 |                 |                 | curvature       |
|                 |                 |                 | tests.          |
+-----------------+-----------------+-----------------+-----------------+
| epsNum          | OT_REAL         | -0.000          | Numerator       |
|                 |                 |                 | tolerance for   |
|                 |                 |                 | ratio tests.    |
+-----------------+-----------------+-----------------+-----------------+
| epsRegularisati | OT_REAL         | 0.000           | Scaling factor  |
| on              |                 |                 | of identity     |
|                 |                 |                 | matrix used for |
|                 |                 |                 | Hessian         |
|                 |                 |                 | regularisation. |
+-----------------+-----------------+-----------------+-----------------+
| finalRamping    | OT_REAL         | 1               | Final value for |
|                 |                 |                 | ramping         |
|                 |                 |                 | strategy.       |
+-----------------+-----------------+-----------------+-----------------+
| growFarBounds   | OT_REAL         | 1000            | Factor to grow  |
|                 |                 |                 | far bounds.     |
+-----------------+-----------------+-----------------+-----------------+
| initialFarBound | OT_REAL         | 1000000         | Initial size    |
| s               |                 |                 | for far bounds. |
+-----------------+-----------------+-----------------+-----------------+
| initialRamping  | OT_REAL         | 0.500           | Start value for |
|                 |                 |                 | ramping         |
|                 |                 |                 | strategy.       |
+-----------------+-----------------+-----------------+-----------------+
| initialStatusBo | OT_STRING       | lower           | Initial status  |
| unds            |                 |                 | of bounds at    |
|                 |                 |                 | first           |
|                 |                 |                 | iteration.      |
+-----------------+-----------------+-----------------+-----------------+
| maxDualJump     | OT_REAL         | 100000000       | Maximum allowed |
|                 |                 |                 | jump in dual    |
|                 |                 |                 | variables in    |
|                 |                 |                 | linear          |
|                 |                 |                 | independence    |
|                 |                 |                 | tests.          |
+-----------------+-----------------+-----------------+-----------------+
| maxPrimalJump   | OT_REAL         | 100000000       | Maximum allowed |
|                 |                 |                 | jump in primal  |
|                 |                 |                 | variables in    |
|                 |                 |                 | nonzero         |
|                 |                 |                 | curvature       |
|                 |                 |                 | tests.          |
+-----------------+-----------------+-----------------+-----------------+
| nWSR            | OT_INTEGER      | None            | The maximum     |
|                 |                 |                 | number of       |
|                 |                 |                 | working set     |
|                 |                 |                 | recalculations  |
|                 |                 |                 | to be performed |
|                 |                 |                 | during the      |
|                 |                 |                 | initial         |
|                 |                 |                 | homotopy.       |
|                 |                 |                 | Default is 5(nx |
|                 |                 |                 | + nc)           |
+-----------------+-----------------+-----------------+-----------------+
| numRefinementSt | OT_INTEGER      | 1               | Maximum number  |
| eps             |                 |                 | of iterative    |
|                 |                 |                 | refinement      |
|                 |                 |                 | steps.          |
+-----------------+-----------------+-----------------+-----------------+
| numRegularisati | OT_INTEGER      | 0               | Maximum number  |
| onSteps         |                 |                 | of successive   |
|                 |                 |                 | regularisation  |
|                 |                 |                 | steps.          |
+-----------------+-----------------+-----------------+-----------------+
| printLevel      | OT_STRING       | medium          | Defines the     |
|                 |                 |                 | amount of text  |
|                 |                 |                 | output during   |
|                 |                 |                 | QP solution,    |
|                 |                 |                 | see Section 5.7 |
+-----------------+-----------------+-----------------+-----------------+
| terminationTole | OT_REAL         | 0.000           | Relative        |
| rance           |                 |                 | termination     |
|                 |                 |                 | tolerance to    |
|                 |                 |                 | stop homotopy.  |
+-----------------+-----------------+-----------------+-----------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

sqic
----



Interface to the SQIC solver for quadratic programming

>List of available options

+----+------+---------+-------------+
| Id | Type | Default | Description |
+====+======+=========+=============+
+----+------+---------+-------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

nlp
---



Solve QPs using an NlpSolver

>List of available options

+----+------+---------+-------------+
| Id | Type | Default | Description |
+====+======+=========+=============+
+----+------+---------+-------------+

>List of available stats

+------------------+
|        Id        |
+==================+
| nlp_solver_stats |
+------------------+

--------------------------------------------------------------------------------



Joel Andersson
Diagrams
--------



C++ includes: qp_solver.hpp ";

%feature("docstring") casadi::QpSolver::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::QpSolver::size2_out "

Get output dimension.

";

%feature("docstring") casadi::QpSolver::n_in "

Get the number of function inputs.

";

%feature("docstring") casadi::QpSolver::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

";

%feature("docstring") casadi::QpSolver::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::QpSolver::getInput "

>  DMatrix  IOInterface< Function  >.getInput(int iind=0) const
------------------------------------------------------------------------

Get an input by index.

Parameters:
-----------

iind:  index within the range [0..n_in()-1]

>  DMatrix  IOInterface< Function  >.getInput(str iname) const
------------------------------------------------------------------------

Get an input by name.

Parameters:
-----------

iname:  input name. Only allowed when an input scheme is set.

>  void IOInterface< Function  >.getInput(T val, int iind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.getInput(T val, str iname)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::QpSolver::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::QpSolver::removeMonitor "

Remove modules to be monitored.

";

%feature("docstring") casadi::QpSolver::fullJacobian "

Generate a Jacobian function of all the inputs elements with respect to all
the output elements).

";

%feature("docstring") casadi::QpSolver::sz_arg "[INTERNAL]  Get required
length of arg field.

";

%feature("docstring") casadi::QpSolver::hessian "

Generate a Hessian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The generated Hessian has two more outputs than the calling function
corresponding to the Hessian and the gradients.

";

%feature("docstring") casadi::QpSolver::getAtomicInputReal "

Get the floating point output argument of an atomic operation.

";

%feature("docstring") casadi::QpSolver::addMonitor "

Add modules to be monitored.

";

%feature("docstring") casadi::QpSolver::free_sx "

Get all the free variables of the function.

";

%feature("docstring") casadi::QpSolver::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::QpSolver::size_out "

Get output dimension.

";

%feature("docstring") casadi::QpSolver::default_in "

Get default input value (NOTE: constant reference)

";

%feature("docstring") casadi::QpSolver::size_in "

Get input dimension.

";

%feature("docstring") casadi::QpSolver::call "

Evaluate the function symbolically or numerically.

";

%feature("docstring") casadi::QpSolver::nnz_in "

Get of number of input nonzeros For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::QpSolver::callReverse "

Create call to (cached) derivative function, reverse mode.

";

%feature("docstring") casadi::QpSolver::sx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::QpSolver::getAtomicOutput "

Get the (integer) output argument of an atomic operation.

";

%feature("docstring") casadi::QpSolver::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::QpSolver::repr "

Print a representation of the object.

";

%feature("docstring") casadi::QpSolver::map "

>  [[MX] ] Function.map([[MX ] ] arg, str parallelization=\"serial\")

>  [MX] Function.map([MX ] arg, str parallelization=\"serial\")
------------------------------------------------------------------------

Evaluate symbolically in parallel (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

>  Function Function.map(str name, int N, Dict options=Dict()) const 
------------------------------------------------------------------------

Create a mapped version of this function.

Suppose the function has a signature of:

::

     f: (a, p) -> ( s )
  



The the mapaccumulated version has the signature:

::

     F: (A, P) -> (S )
  
      with
          a: horzcat([a0, a1, ..., a_(N-1)])
          p: horzcat([p0, p1, ..., p_(N-1)])
          s: horzcat([s0, s1, ..., s_(N-1)])
      and
          s0 <- f(a0, p0)
          s1 <- f(a1, p1)
          ...
          s_(N-1) <- f(a_(N-1), p_(N-1))
  



>  Function Function.map(str name, int n, [bool ] repeat_in, [bool ] repeat_out, Dict opts=Dict()) const 
------------------------------------------------------------------------

Generic map.

";

%feature("docstring") casadi::QpSolver::sx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::QpSolver::sz_w "[INTERNAL]  Get required
length of w field.

";

%feature("docstring") casadi::QpSolver::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::QpSolver::getOptionAllowedIndex "[INTERNAL]
Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::QpSolver::description_in "

Get input scheme description by index.

";

%feature("docstring") casadi::QpSolver::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::QpSolver::callDerivative "[INTERNAL]
Evaluate the function symbolically or numerically with directional
derivatives The first two arguments are the nondifferentiated inputs and
results of the evaluation, the next two arguments are a set of forward
directional seeds and the resulting forward directional derivatives, the
length of the vector being the number of forward directions. The next two
arguments are a set of adjoint directional seeds and the resulting adjoint
directional derivatives, the length of the vector being the number of
adjoint directions.

";

%feature("docstring") casadi::QpSolver::getOutput "

>  DMatrix  IOInterface< Function  >.getOutput(int oind=0) const
------------------------------------------------------------------------

Get an output by index.

Parameters:
-----------

oind:  index within the range [0..n_out()-1]

>  DMatrix  IOInterface< Function  >.getOutput(str oname) const
------------------------------------------------------------------------

Get an output by name.

Parameters:
-----------

oname:  output name. Only allowed when an output scheme is set.

>  void IOInterface< Function  >.getOutput(T val, int oind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.getOutput(T val, str oname)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::QpSolver::getAtomicOperation "

Get an atomic operation operator index.

";

%feature("docstring") casadi::QpSolver::spInit "[INTERNAL]  Reset the
sparsity propagation.

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::QpSolver::setOptionByEnumValue "[INTERNAL]
Set a certain option by giving an enum value.

";

%feature("docstring") casadi::QpSolver::setJacSparsity "

Generate the sparsity of a Jacobian block

";

%feature("docstring") casadi::QpSolver::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::QpSolver::tangent "

Generate a tangent function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the input must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::QpSolver::print "

Print a description of the object.

";

%feature("docstring") casadi::QpSolver::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::QpSolver::derForward "

Get a function that calculates nfwd forward derivatives.

Returns a function with n_in + n_out +nfwd*n_in inputs and nfwd*n_out
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nfwd*n_in inputs correspond to forward seeds, one direction at a time The
nfwd*n_out outputs correspond to forward sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::QpSolver::setJacobian "

Set the Jacobian function of output oind with respect to input iind NOTE:
Does not take ownership, only weak references to the Jacobians are kept
internally

";

%feature("docstring") casadi::QpSolver::getAlgorithmSize "

Get the number of atomic operations.

";

%feature("docstring") casadi::QpSolver::nnz_out "

Get of number of output nonzeros For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::QpSolver::evaluate "

Evaluate.

";

%feature("docstring") casadi::QpSolver::setDerForward "

Set a function that calculates nfwd forward derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::QpSolver::checkInputs "[INTERNAL]  Check if
the numerical values of the supplied bounds make sense.

";

%feature("docstring") casadi::QpSolver::getStats "

Get all statistics obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::QpSolver::mapaccum "

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
  



";

%feature("docstring") casadi::QpSolver::index_in "

Find the index for a string describing a particular entry of an input
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::QpSolver::derivative "

Get a function that calculates nfwd forward derivatives and nadj adjoint
derivatives Legacy function: Use derForward and derReverse instead.

Returns a function with (1+nfwd)*n_in+nadj*n_out inputs and (1+nfwd)*n_out +
nadj*n_in outputs. The first n_in inputs correspond to nondifferentiated
inputs. The next nfwd*n_in inputs correspond to forward seeds, one direction
at a time and the last nadj*n_out inputs correspond to adjoint seeds, one
direction at a time. The first n_out outputs correspond to nondifferentiated
outputs. The next nfwd*n_out outputs correspond to forward sensitivities,
one direction at a time and the last nadj*n_in outputs corresponds to
adjoint sensitivities, one direction at a time.

(n_in = n_in(), n_out = n_out())

";

%feature("docstring") casadi::QpSolver::n_out "

Get the number of function outputs.

";

%feature("docstring") casadi::QpSolver::countNodes "

Number of nodes in the algorithm.

";

%feature("docstring") casadi::QpSolver::numel_in "

Get of number of input elements For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::QpSolver::free_mx "

Get all the free variables of the function.

";

%feature("docstring") casadi::QpSolver::setDerReverse "

Set a function that calculates nadj adjoint derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";


// File: classcasadi_1_1QpToImplicit.xml


// File: classcasadi_1_1QpToNlp.xml


// File: classcasadi_1_1RealtypeSX.xml


// File: classcasadi_1_1Reshape.xml


// File: classcasadi_1_1RkIntegrator.xml


// File: classcasadi_1_1Scpgen.xml


// File: classcasadi_1_1SetNonzeros.xml


// File: classcasadi_1_1SetNonzerosSlice.xml


// File: classcasadi_1_1SetNonzerosSlice2.xml


// File: classcasadi_1_1SetNonzerosVector.xml


// File: classcasadi_1_1SharedObject.xml
%feature("docstring") casadi::SharedObject::isNull "

Is a null pointer?

";

%feature("docstring") casadi::SharedObject::print "

Print a description of the object.

";

%feature("docstring") casadi::SharedObject::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::SharedObject::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::SharedObject "

SharedObject implements a reference counting framework similar for efficient
and easily-maintained memory management.

To use the class, both the SharedObject class (the public class), and the
SharedObjectNode class (the internal class) must be inherited from. It can
be done in two different files and together with memory management, this
approach provides a clear distinction of which methods of the class are to
be considered \"public\", i.e. methods for public use that can be considered
to remain over time with small changes, and the internal memory.

When interfacing a software, which typically includes including some header
file, this is best done only in the file where the internal class is
defined, to avoid polluting the global namespace and other side effects.

The default constructor always means creating a null pointer to an internal
class only. To allocate an internal class (this works only when the internal
class isn't abstract), use the constructor with arguments.

The copy constructor and the assignment operator perform shallow copies
only, to make a deep copy you must use the clone method explicitly. This
will give a shared pointer instance.

In an inheritance hierarchy, you can cast down automatically, e.g.
(SXFunction is a child class of Function): SXFunction derived(...); Function
base = derived;

To cast up, use the shared_cast template function, which works analogously
to dynamic_cast, static_cast, const_cast etc, e.g.: SXFunction derived(...);
Function base = derived; SXFunction derived_from_base =
shared_cast<SXFunction>(base);

A failed shared_cast will result in a null pointer (cf. dynamic_cast)

Joel Andersson

C++ includes: shared_object.hpp ";

%feature("docstring") casadi::SharedObject::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::SharedObject::repr "

Print a representation of the object.

";

%feature("docstring") casadi::SharedObject::getDescription "

Return a string with a description (for SWIG)

";


// File: classcasadi_1_1ShellCompiler.xml


// File: classcasadi_1_1SimplifiedExternal.xml


// File: classcasadi_1_1Simulator.xml


/*  Simple Getters & Setters  */

/*  Advanced Getters  */

/*  Option Functionality  */ %feature("docstring")
casadi::Simulator::getOptionEnumValue " [INTERNAL]  Get the enum value
corresponding to th certain option.

";

%feature("docstring") casadi::Simulator::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::Simulator::fullJacobian "

Generate a Jacobian function of all the inputs elements with respect to all
the output elements).

";

%feature("docstring") casadi::Simulator::index_in "

Find the index for a string describing a particular entry of an input
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::Simulator::getAtomicInput "

Get the (integer) input arguments of an atomic operation.

";

%feature("docstring") casadi::Simulator::map "

>  [[MX] ] Function.map([[MX ] ] arg, str parallelization=\"serial\")

>  [MX] Function.map([MX ] arg, str parallelization=\"serial\")
------------------------------------------------------------------------

Evaluate symbolically in parallel (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

>  Function Function.map(str name, int N, Dict options=Dict()) const 
------------------------------------------------------------------------

Create a mapped version of this function.

Suppose the function has a signature of:

::

     f: (a, p) -> ( s )
  



The the mapaccumulated version has the signature:

::

     F: (A, P) -> (S )
  
      with
          a: horzcat([a0, a1, ..., a_(N-1)])
          p: horzcat([p0, p1, ..., p_(N-1)])
          s: horzcat([s0, s1, ..., s_(N-1)])
      and
          s0 <- f(a0, p0)
          s1 <- f(a1, p1)
          ...
          s_(N-1) <- f(a_(N-1), p_(N-1))
  



>  Function Function.map(str name, int n, [bool ] repeat_in, [bool ] repeat_out, Dict opts=Dict()) const 
------------------------------------------------------------------------

Generic map.

";

%feature("docstring") casadi::Simulator::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::Simulator::evaluate "

Evaluate.

";

%feature("docstring") casadi::Simulator::sx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::Simulator::sz_iw "[INTERNAL]  Get required
length of iw field.

";

%feature("docstring") casadi::Simulator::jacSparsity "

Get, if necessary generate, the sparsity of a Jacobian block

";

%feature("docstring") casadi::Simulator::getStat "

Get a single statistic obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::Simulator::getOptionAllowedIndex "[INTERNAL]
Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::Simulator::description_in "

Get input scheme description by index.

";

%feature("docstring") casadi::Simulator::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::Simulator::getInput "

>  DMatrix  IOInterface< Function  >.getInput(int iind=0) const
------------------------------------------------------------------------

Get an input by index.

Parameters:
-----------

iind:  index within the range [0..n_in()-1]

>  DMatrix  IOInterface< Function  >.getInput(str iname) const
------------------------------------------------------------------------

Get an input by name.

Parameters:
-----------

iname:  input name. Only allowed when an input scheme is set.

>  void IOInterface< Function  >.getInput(T val, int iind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.getInput(T val, str iname)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::Simulator::is_a "

Check if the function is of a particular type Optionally check if name
matches one of the base classes (default true)

";

%feature("docstring") casadi::Simulator::name_in "

>  [str] Function.name_in() const 
------------------------------------------------------------------------

Get input scheme.

>  str Function.name_in(int ind) const 
------------------------------------------------------------------------

Get input scheme name by index.

";

%feature("docstring") casadi::Simulator::size1_out "

Get output dimension.

";

%feature("docstring") casadi::Simulator "

Integrator class.

An \"simulator\" integrates an IVP, stopping at a (fixed) number of grid
points and evaluates a set of output functions at these points. The internal
stepsizes of the integrator need not coincide with the gridpoints.

Simulator is an casadi::Function mapping from casadi::IntegratorInput to n.
\\\\

The output function needs to be a mapping from casadi::DAEInput to n. The
default output has n=1 and the output is the (vectorized) differential state
for each time step.

Joel Andersson

>Input scheme: casadi::IntegratorInput (INTEGRATOR_NUM_IN = 6) [integratorIn]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| INTEGRATOR_X0          | x0                     | Differential state at  |
|                        |                        | the initial time .     |
+------------------------+------------------------+------------------------+
| INTEGRATOR_P           | p                      | Parameters .           |
+------------------------+------------------------+------------------------+
| INTEGRATOR_Z0          | z0                     | Initial guess for the  |
|                        |                        | algebraic variable .   |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RX0         | rx0                    | Backward differential  |
|                        |                        | state at the final     |
|                        |                        | time .                 |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RP          | rp                     | Backward parameter     |
|                        |                        | vector .               |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RZ0         | rz0                    | Initial guess for the  |
|                        |                        | backwards algebraic    |
|                        |                        | variable .             |
+------------------------+------------------------+------------------------+

>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| ad_weight    | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | derivative c |              |
|              |              |              | alculation.W |              |
|              |              |              | hen there is |              |
|              |              |              | an option of |              |
|              |              |              | either using |              |
|              |              |              | forward or   |              |
|              |              |              | reverse mode |              |
|              |              |              | directional  |              |
|              |              |              | derivatives, |              |
|              |              |              | the          |              |
|              |              |              | condition ad |              |
|              |              |              | _weight*nf<= |              |
|              |              |              | (1-ad_weight |              |
|              |              |              | )*na is used |              |
|              |              |              | where nf and |              |
|              |              |              | na are       |              |
|              |              |              | estimates of |              |
|              |              |              | the number   |              |
|              |              |              | of forward/r |              |
|              |              |              | everse mode  |              |
|              |              |              | directional  |              |
|              |              |              | derivatives  |              |
|              |              |              | needed. By   |              |
|              |              |              | default,     |              |
|              |              |              | ad_weight is |              |
|              |              |              | calculated a |              |
|              |              |              | utomatically |              |
|              |              |              | , but this   |              |
|              |              |              | can be       |              |
|              |              |              | overridden   |              |
|              |              |              | by setting   |              |
|              |              |              | this option. |              |
|              |              |              | In           |              |
|              |              |              | particular,  |              |
|              |              |              | 0 means      |              |
|              |              |              | forcing      |              |
|              |              |              | forward mode |              |
|              |              |              | and 1        |              |
|              |              |              | forcing      |              |
|              |              |              | reverse      |              |
|              |              |              | mode. Leave  |              |
|              |              |              | unset for    |              |
|              |              |              | (class       |              |
|              |              |              | specific)    |              |
|              |              |              | heuristics.  |              |
+--------------+--------------+--------------+--------------+--------------+
| ad_weight_sp | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | sparsity     |              |
|              |              |              | pattern      |              |
|              |              |              | calculation  |              |
|              |              |              | calculation. |              |
|              |              |              | Overrides    |              |
|              |              |              | default      |              |
|              |              |              | behavior.    |              |
|              |              |              | Set to 0 and |              |
|              |              |              | 1 to force   |              |
|              |              |              | forward and  |              |
|              |              |              | reverse mode |              |
|              |              |              | respectively |              |
|              |              |              | . Cf. option |              |
|              |              |              | \"ad_weight\". |              |
+--------------+--------------+--------------+--------------+--------------+
| compiler     | OT_STRING    | \"clang\"      | Just-in-time | casadi::Func |
|              |              |              | compiler     | tionInternal |
|              |              |              | plugin to be |              |
|              |              |              | used.        |              |
+--------------+--------------+--------------+--------------+--------------+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode     |
|              |              |              | according to |              |
|              |              |              | a given      |              |
|              |              |              | recipe (low- |              |
|              |              |              | level)       |              |
+--------------+--------------+--------------+--------------+--------------+
| gather_stats | OT_BOOLEAN   | false        | Flag to      | casadi::Func |
|              |              |              | indicate     | tionInternal |
|              |              |              | whether      |              |
|              |              |              | statistics   |              |
|              |              |              | must be      |              |
|              |              |              | gathered     |              |
+--------------+--------------+--------------+--------------+--------------+
| input_scheme | OT_STRINGVEC | GenericType( | Custom input | casadi::Func |
|              | TOR          | )            | scheme       | tionInternal |
+--------------+--------------+--------------+--------------+--------------+
| inputs_check | OT_BOOLEAN   | true         | Throw        | casadi::Func |
|              |              |              | exceptions   | tionInternal |
|              |              |              | when the     |              |
|              |              |              | numerical    |              |
|              |              |              | values of    |              |
|              |              |              | the inputs   |              |
|              |              |              | don't make   |              |
|              |              |              | sense        |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_penalty  | OT_REAL      | 2            | When         | casadi::Func |
|              |              |              | requested    | tionInternal |
|              |              |              | for a number |              |
|              |              |              | of forward/r |              |
|              |              |              | everse       |              |
|              |              |              | directions,  |              |
|              |              |              | it may be    |              |
|              |              |              | cheaper to   |              |
|              |              |              | compute      |              |
|              |              |              | first the    |              |
|              |              |              | full         |              |
|              |              |              | jacobian and |              |
|              |              |              | then         |              |
|              |              |              | multiply     |              |
|              |              |              | with seeds,  |              |
|              |              |              | rather than  |              |
|              |              |              | obtain the   |              |
|              |              |              | requested    |              |
|              |              |              | directions   |              |
|              |              |              | in a straigh |              |
|              |              |              | tforward     |              |
|              |              |              | manner.      |              |
|              |              |              | Casadi uses  |              |
|              |              |              | a heuristic  |              |
|              |              |              | to decide    |              |
|              |              |              | which is     |              |
|              |              |              | cheaper. A   |              |
|              |              |              | high value   |              |
|              |              |              | of 'jac_pena |              |
|              |              |              | lty' makes   |              |
|              |              |              | it less      |              |
|              |              |              | likely for   |              |
|              |              |              | the heurstic |              |
|              |              |              | to chose the |              |
|              |              |              | full         |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy.    |              |
|              |              |              | The special  |              |
|              |              |              | value -1     |              |
|              |              |              | indicates    |              |
|              |              |              | never to use |              |
|              |              |              | the full     |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy     |              |
+--------------+--------------+--------------+--------------+--------------+
| jit          | OT_BOOLEAN   | false        | Use just-in- | casadi::Func |
|              |              |              | time         | tionInternal |
|              |              |              | compiler to  |              |
|              |              |              | speed up the |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| jit_options  | OT_DICT      | GenericType( | Options to   | casadi::Func |
|              |              | )            | be passed to | tionInternal |
|              |              |              | the jit      |              |
|              |              |              | compiler.    |              |
+--------------+--------------+--------------+--------------+--------------+
| monitor      | OT_STRINGVEC | GenericType( | Monitors to  | casadi::Func |
|              | TOR          | )            | be activated | tionInternal |
|              |              |              | (inputs|outp | casadi::Simu |
|              |              |              | uts)  (initi | latorInterna |
|              |              |              | al|step)     | l            |
+--------------+--------------+--------------+--------------+--------------+
| output_schem | OT_STRINGVEC | GenericType( | Custom       | casadi::Func |
| e            | TOR          | )            | output       | tionInternal |
|              |              |              | scheme       |              |
+--------------+--------------+--------------+--------------+--------------+
| regularity_c | OT_BOOLEAN   | true         | Throw        | casadi::Func |
| heck         |              |              | exceptions   | tionInternal |
|              |              |              | when NaN or  |              |
|              |              |              | Inf appears  |              |
|              |              |              | during       |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| user_data    | OT_VOIDPTR   | GenericType( | A user-      | casadi::Func |
|              |              | )            | defined      | tionInternal |
|              |              |              | field that   |              |
|              |              |              | can be used  |              |
|              |              |              | to identify  |              |
|              |              |              | the function |              |
|              |              |              | or pass      |              |
|              |              |              | additional   |              |
|              |              |              | information  |              |
+--------------+--------------+--------------+--------------+--------------+
| verbose      | OT_BOOLEAN   | false        | Verbose      | casadi::Func |
|              |              |              | evaluation   | tionInternal |
|              |              |              | for          |              |
|              |              |              | debugging    |              |
+--------------+--------------+--------------+--------------+--------------+

>List of available monitors

+---------+---------------------------+
|   Id    |          Used in          |
+=========+===========================+
| initial | casadi::SimulatorInternal |
+---------+---------------------------+
| inputs  | casadi::FunctionInternal  |
+---------+---------------------------+
| outputs | casadi::FunctionInternal  |
+---------+---------------------------+
| step    | casadi::SimulatorInternal |
+---------+---------------------------+

Diagrams
--------



C++ includes: simulator.hpp ";

%feature("docstring") casadi::Simulator::getSanitizedName "

get function name with all non alphanumeric characters converted to '_'

";

%feature("docstring") casadi::Simulator::tangent "

Generate a tangent function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the input must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::Simulator::getAtomicInputReal "

Get the floating point output argument of an atomic operation.

";

%feature("docstring") casadi::Simulator::setJacobian "

Set the Jacobian function of output oind with respect to input iind NOTE:
Does not take ownership, only weak references to the Jacobians are kept
internally

";

%feature("docstring") casadi::Simulator::generate "

Export / Generate C code for the function.

";

%feature("docstring") casadi::Simulator::getAtomicOutput "

Get the (integer) output argument of an atomic operation.

";

%feature("docstring") casadi::Simulator::dictionary "

Get the dictionary.

";

%feature("docstring") casadi::Simulator::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::Simulator::size_in "

Get input dimension.

";

%feature("docstring") casadi::Simulator::callDerivative "[INTERNAL]
Evaluate the function symbolically or numerically with directional
derivatives The first two arguments are the nondifferentiated inputs and
results of the evaluation, the next two arguments are a set of forward
directional seeds and the resulting forward directional derivatives, the
length of the vector being the number of forward directions. The next two
arguments are a set of adjoint directional seeds and the resulting adjoint
directional derivatives, the length of the vector being the number of
adjoint directions.

";

%feature("docstring") casadi::Simulator::size2_out "

Get output dimension.

";

%feature("docstring") casadi::Simulator::setInput "

>  void IOInterface< Function  >.setInput(T val, int iind=0)
------------------------------------------------------------------------

Set an input by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.setInput(T val, str iname)
------------------------------------------------------------------------

Set an input by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::Simulator::callForward "

Create call to (cached) derivative function, forward mode.

";

%feature("docstring") casadi::Simulator::sparsity_in "

Get sparsity of a given input.

";

%feature("docstring") casadi::Simulator::generateLiftingFunctions "

Extract the functions needed for the Lifted Newton method.

";

%feature("docstring") casadi::Simulator::spCanEvaluate "[INTERNAL]  Is the
class able to propagate seeds through the algorithm?

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::Simulator::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::Simulator::removeMonitor "

Remove modules to be monitored.

";

%feature("docstring") casadi::Simulator::getStats "

Get all statistics obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::Simulator::sz_arg "[INTERNAL]  Get required
length of arg field.

";

%feature("docstring") casadi::Simulator::size2_in "

Get input dimension.

";

%feature("docstring") casadi::Simulator::getWorkSize "

Get the length of the work vector.

";

%feature("docstring") casadi::Simulator::sz_w "[INTERNAL]  Get required
length of w field.

";

%feature("docstring") casadi::Simulator::sx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::Simulator::default_in "

Get default input value (NOTE: constant reference)

";

%feature("docstring") casadi::Simulator::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::Simulator::isNull "

Is a null pointer?

";

%feature("docstring") casadi::Simulator::repr "

Print a representation of the object.

";

%feature("docstring") casadi::Simulator::gradient "

Generate a gradient function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the output must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::Simulator::nnz_in "

Get of number of input nonzeros For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::Simulator::getOptionTypeName "

Get the type name of a certain option.

";

%feature("docstring") casadi::Simulator::sz_res "[INTERNAL]  Get required
length of res field.

";

%feature("docstring") casadi::Simulator::name "

Name of the function.

";

%feature("docstring") casadi::Simulator::checkInputs "[INTERNAL]  Check if
the numerical values of the supplied bounds make sense.

";

%feature("docstring") casadi::Simulator::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::Simulator::jacobian "

Generate a Jacobian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. If
compact is set to true, only the nonzeros of the input and output
expressions are considered. If symmetric is set to true, the Jacobian being
calculated is known to be symmetric (usually a Hessian), which can be
exploited by the algorithm.

The generated Jacobian has one more output than the calling function
corresponding to the Jacobian and the same number of inputs.

";

%feature("docstring") casadi::Simulator::size1_in "

Get input dimension.

";

%feature("docstring") casadi::Simulator::setFullJacobian "

Set the Jacobian of all the input nonzeros with respect to all output
nonzeros NOTE: Does not take ownership, only weak references to the Jacobian
are kept internally

";

%feature("docstring") casadi::Simulator::sparsity_out "

Get sparsity of a given output.

";

%feature("docstring") casadi::Simulator::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::Simulator::setOptionByAllowedIndex "[INTERNAL]  Set a certain option by giving its index into the allowed
values.

";

%feature("docstring") casadi::Simulator::spEvaluate "[INTERNAL]  Propagate
the sparsity pattern through a set of directional.

derivatives forward or backward (for usage, see the example
propagating_sparsity.cpp)

";

%feature("docstring") casadi::Simulator::derForward "

Get a function that calculates nfwd forward derivatives.

Returns a function with n_in + n_out +nfwd*n_in inputs and nfwd*n_out
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nfwd*n_in inputs correspond to forward seeds, one direction at a time The
nfwd*n_out outputs correspond to forward sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::Simulator::setDerReverse "

Set a function that calculates nadj adjoint derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::Simulator::getAtomicOperation "

Get an atomic operation operator index.

";

%feature("docstring") casadi::Simulator::getAlgorithmSize "

Get the number of atomic operations.

";

%feature("docstring") casadi::Simulator::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::Simulator::getOutput "

>  DMatrix  IOInterface< Function  >.getOutput(int oind=0) const
------------------------------------------------------------------------

Get an output by index.

Parameters:
-----------

oind:  index within the range [0..n_out()-1]

>  DMatrix  IOInterface< Function  >.getOutput(str oname) const
------------------------------------------------------------------------

Get an output by name.

Parameters:
-----------

oname:  output name. Only allowed when an output scheme is set.

>  void IOInterface< Function  >.getOutput(T val, int oind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.getOutput(T val, str oname)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::Simulator::setOptionByEnumValue "[INTERNAL]
Set a certain option by giving an enum value.

";

%feature("docstring") casadi::Simulator::callReverse "

Create call to (cached) derivative function, reverse mode.

";

%feature("docstring") casadi::Simulator::numel_out "

Get of number of output elements For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::Simulator::derReverse "

Get a function that calculates nadj adjoint derivatives.

Returns a function with n_in + n_out +nadj*n_out inputs and nadj*n_in
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nadj*n_out inputs correspond to adjoint seeds, one direction at a time The
nadj*n_in outputs correspond to adjoint sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

(n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::Simulator::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::Simulator::printDimensions "

Print dimensions of inputs and outputs.

";

%feature("docstring") casadi::Simulator::Simulator "

>  Simulator()
------------------------------------------------------------------------

Default constructor.

>  Simulator(str name, Integrator integrator, Function output_fcn, DMatrix grid, Dict opts=Dict())
------------------------------------------------------------------------

Constructor.

Parameters:
-----------

output_fcn:  output function which maps to n outputs. (new syntax, includes
initialization)

>Input scheme: casadi::DAEInput (DAE_NUM_IN = 4) [daeIn]

+-----------+-------+----------------------------+
| Full name | Short |        Description         |
+===========+=======+============================+
| DAE_X     | x     | Differential state .       |
+-----------+-------+----------------------------+
| DAE_Z     | z     | Algebraic state .          |
+-----------+-------+----------------------------+
| DAE_P     | p     | Parameter .                |
+-----------+-------+----------------------------+
| DAE_T     | t     | Explicit time dependence . |
+-----------+-------+----------------------------+

>  Simulator(str name, Integrator integrator, DMatrix grid, Dict opts=Dict())
------------------------------------------------------------------------

Output function equal to the state (new syntax, includes initialization)

";

%feature("docstring") casadi::Simulator::name_out "

>  [str] Function.name_out() const 
------------------------------------------------------------------------

Get output scheme.

>  str Function.name_out(int ind) const 
------------------------------------------------------------------------

Get output scheme name by index.

";

%feature("docstring") casadi::Simulator::n_in "

Get the number of function inputs.

";

%feature("docstring") casadi::Simulator::mx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::Simulator::countNodes "

Number of nodes in the algorithm.

";

%feature("docstring") casadi::Simulator::free_sx "

Get all the free variables of the function.

";

%feature("docstring") casadi::Simulator::mapaccum "

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
  



";

%feature("docstring") casadi::Simulator::derivative "

Get a function that calculates nfwd forward derivatives and nadj adjoint
derivatives Legacy function: Use derForward and derReverse instead.

Returns a function with (1+nfwd)*n_in+nadj*n_out inputs and (1+nfwd)*n_out +
nadj*n_in outputs. The first n_in inputs correspond to nondifferentiated
inputs. The next nfwd*n_in inputs correspond to forward seeds, one direction
at a time and the last nadj*n_out inputs correspond to adjoint seeds, one
direction at a time. The first n_out outputs correspond to nondifferentiated
outputs. The next nfwd*n_out outputs correspond to forward sensitivities,
one direction at a time and the last nadj*n_in outputs corresponds to
adjoint sensitivities, one direction at a time.

(n_in = n_in(), n_out = n_out())

";

%feature("docstring") casadi::Simulator::numel_in "

Get of number of input elements For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::Simulator::n_out "

Get the number of function outputs.

";

%feature("docstring") casadi::Simulator::setDerForward "

Set a function that calculates nfwd forward derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::Simulator::spInit "[INTERNAL]  Reset the
sparsity propagation.

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::Simulator::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::Simulator::index_out "

Find the index for a string describing a particular entry of an output
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::Simulator::call "

Evaluate the function symbolically or numerically.

";

%feature("docstring") casadi::Simulator::free_mx "

Get all the free variables of the function.

";

%feature("docstring") casadi::Simulator::size_out "

Get output dimension.

";

%feature("docstring") casadi::Simulator::addMonitor "

Add modules to be monitored.

";

%feature("docstring") casadi::Simulator::description_out "

Get output scheme description by index.

";

%feature("docstring") casadi::Simulator::type_name "

Get type name.

";

%feature("docstring") casadi::Simulator::setOutput "

>  void IOInterface< Function  >.setOutput(T val, int oind=0)
------------------------------------------------------------------------

Set an output by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.setOutput(T val, str oname)
------------------------------------------------------------------------

Set an output by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::Simulator::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

";

%feature("docstring") casadi::Simulator::nnz_out "

Get of number of output nonzeros For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::Simulator::setJacSparsity "

Generate the sparsity of a Jacobian block

";

%feature("docstring") casadi::Simulator::hessian "

Generate a Hessian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The generated Hessian has two more outputs than the calling function
corresponding to the Hessian and the gradients.

";

%feature("docstring") casadi::Simulator::mx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::Simulator::print "

Print a description of the object.

";


// File: classcasadi_1_1Slice.xml
%feature("docstring") casadi::Slice::isscalar "

Is the slice a scalar.

";

%feature("docstring") casadi::Slice::getAll "

>  [int] Slice.getAll(int len, bool ind1=false) const 
------------------------------------------------------------------------

Get a vector of indices.

>  [int] Slice.getAll(Slice outer, int len) const 
------------------------------------------------------------------------

Get a vector of indices (nested slice)

";

%feature("docstring") casadi::Slice::print "

Print a description of the object.

";

%feature("docstring") casadi::Slice "

Class representing a Slice.

Note that Python or Octave do not need to use this class. They can just use
slicing utility from the host language ( M[0:6] in Python, M(1:7) )

C++ includes: slice.hpp ";

%feature("docstring") casadi::Slice::toScalar "

Get scalar (if isscalar)

";

%feature("docstring") casadi::Slice::Slice "

>  Slice()
------------------------------------------------------------------------

Default constructor - all elements.

>  Slice(int i, bool ind1=false)
------------------------------------------------------------------------

A single element (explicit to avoid ambiguity with IMatrix overload.

>  Slice(int start, int stop, int step=1)
------------------------------------------------------------------------

A slice.

>  Slice([int ] v, bool ind1=false)
------------------------------------------------------------------------

Construct from an index vector (requires isSlice(v) to be true)

>  Slice([int ] v, Slice &outer)
------------------------------------------------------------------------

Construct nested slices from an index vector (requires isSlice2(v) to be
true)

";

%feature("docstring") casadi::Slice::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::Slice::repr "

Print a representation of the object.

";

%feature("docstring") casadi::Slice::getRepresentation "

Return a string with a representation (for SWIG)

";


// File: classsnoptProblem.xml
%feature("docstring") snoptProblem::setParameter "[INTERNAL] ";

%feature("docstring") snoptProblem::setLog "[INTERNAL] ";

%feature("docstring") snoptProblem::getIntParameter "[INTERNAL] ";

%feature("docstring") snoptProblem::setSpecsFile "[INTERNAL] ";

%feature("docstring") snoptProblem::setIntParameter "[INTERNAL] ";

%feature("docstring") snoptProblem::getRealParameter "[INTERNAL] ";

%feature("docstring") snoptProblem::setPrintFile "[INTERNAL] ";

%feature("docstring") snoptProblem::getParameter "[INTERNAL] ";

%feature("docstring") snoptProblem::setUserR "[INTERNAL] ";

%feature("docstring") snoptProblem::setSTOP "[INTERNAL] ";

%feature("docstring") snoptProblem::setRealParameter "[INTERNAL] ";

%feature("docstring") snoptProblem::setUserI "[INTERNAL] ";

%feature("docstring") snoptProblem::setProbName "[INTERNAL] ";

%feature("docstring") snoptProblem::solve "[INTERNAL] ";

%feature("docstring") snoptProblem::setUserspace "[INTERNAL] ";

%feature("docstring") snoptProblem "[INTERNAL] C++ includes:
snoptProblem.hpp ";


// File: classsnoptProblemA.xml
%feature("docstring") snoptProblemA::getIntParameter "[INTERNAL] ";

%feature("docstring") snoptProblemA::computeJac "[INTERNAL] ";

%feature("docstring") snoptProblemA::snoptProblemA "[INTERNAL] ";

%feature("docstring") snoptProblemA::getRealParameter "[INTERNAL] ";

%feature("docstring") snoptProblemA::setSpecsFile "[INTERNAL] ";

%feature("docstring") snoptProblemA::setRealParameter "[INTERNAL] ";

%feature("docstring") snoptProblemA::setWorkspace "[INTERNAL] ";

%feature("docstring") snoptProblemA::setNeG "[INTERNAL] ";

%feature("docstring") snoptProblemA::setNeA "[INTERNAL] ";

%feature("docstring") snoptProblemA::setSTOP "[INTERNAL] ";

%feature("docstring") snoptProblemA::setUserspace "[INTERNAL] ";

%feature("docstring") snoptProblemA::setIntParameter "[INTERNAL] ";

%feature("docstring") snoptProblemA::setG "[INTERNAL] ";

%feature("docstring") snoptProblemA::setF "[INTERNAL] ";

%feature("docstring") snoptProblemA::~snoptProblemA "[INTERNAL] ";

%feature("docstring") snoptProblemA::setPrintFile "[INTERNAL] ";

%feature("docstring") snoptProblemA::setA "[INTERNAL] ";

%feature("docstring") snoptProblemA::solve "[INTERNAL] ";

%feature("docstring") snoptProblemA "[INTERNAL] C++ includes:
snoptProblem.hpp ";

%feature("docstring") snoptProblemA::setX "[INTERNAL] ";

%feature("docstring") snoptProblemA::setUserFun "[INTERNAL] ";

%feature("docstring") snoptProblemA::getParameter "[INTERNAL] ";

%feature("docstring") snoptProblemA::setProblemSize "[INTERNAL] ";

%feature("docstring") snoptProblemA::setParameter "[INTERNAL] ";

%feature("docstring") snoptProblemA::setLog "[INTERNAL] ";

%feature("docstring") snoptProblemA::setUserI "[INTERNAL] ";

%feature("docstring") snoptProblemA::setObjective "[INTERNAL] ";

%feature("docstring") snoptProblemA::setUserR "[INTERNAL] ";

%feature("docstring") snoptProblemA::setProbName "[INTERNAL] ";


// File: classsnoptProblemB.xml
%feature("docstring") snoptProblemB::setPrintFile "[INTERNAL] ";

%feature("docstring") snoptProblemB::getRealParameter "[INTERNAL] ";

%feature("docstring") snoptProblemB::setFuncon "[INTERNAL] ";

%feature("docstring") snoptProblemB::setParameter "[INTERNAL] ";

%feature("docstring") snoptProblemB::setObjective "[INTERNAL] ";

%feature("docstring") snoptProblemB::setWorkspace "[INTERNAL] ";

%feature("docstring") snoptProblemB::setSTOP "[INTERNAL] ";

%feature("docstring") snoptProblemB::getIntParameter "[INTERNAL] ";

%feature("docstring") snoptProblemB::setFunobj "[INTERNAL] ";

%feature("docstring") snoptProblemB::setUserspace "[INTERNAL] ";

%feature("docstring") snoptProblemB::snoptProblemB "[INTERNAL] ";

%feature("docstring") snoptProblemB::getParameter "[INTERNAL] ";

%feature("docstring") snoptProblemB::setUserFun "[INTERNAL] ";

%feature("docstring") snoptProblemB::setIntParameter "[INTERNAL] ";

%feature("docstring") snoptProblemB::setLog "[INTERNAL] ";

%feature("docstring") snoptProblemB::setUserI "[INTERNAL] ";

%feature("docstring") snoptProblemB::setJ "[INTERNAL] ";

%feature("docstring") snoptProblemB::setProbName "[INTERNAL] ";

%feature("docstring") snoptProblemB::setUserR "[INTERNAL] ";

%feature("docstring") snoptProblemB "[INTERNAL] C++ includes:
snoptProblem.hpp ";

%feature("docstring") snoptProblemB::setProblemSize "[INTERNAL] ";

%feature("docstring") snoptProblemB::setX "[INTERNAL] ";

%feature("docstring") snoptProblemB::setSpecsFile "[INTERNAL] ";

%feature("docstring") snoptProblemB::setRealParameter "[INTERNAL] ";

%feature("docstring") snoptProblemB::solve "[INTERNAL] ";

%feature("docstring") snoptProblemB::~snoptProblemB "[INTERNAL] ";


// File: classsnoptProblemC.xml
%feature("docstring") snoptProblemC::setPrintFile "[INTERNAL] ";

%feature("docstring") snoptProblemC "[INTERNAL] C++ includes:
snoptProblem.hpp ";

%feature("docstring") snoptProblemC::getRealParameter "[INTERNAL] ";

%feature("docstring") snoptProblemC::setRealParameter "[INTERNAL] ";

%feature("docstring") snoptProblemC::solve "[INTERNAL] ";

%feature("docstring") snoptProblemC::setParameter "[INTERNAL] ";

%feature("docstring") snoptProblemC::setUserR "[INTERNAL] ";

%feature("docstring") snoptProblemC::getIntParameter "[INTERNAL] ";

%feature("docstring") snoptProblemC::getParameter "[INTERNAL] ";

%feature("docstring") snoptProblemC::snoptProblemC "[INTERNAL] ";

%feature("docstring") snoptProblemC::setUserspace "[INTERNAL] ";

%feature("docstring") snoptProblemC::setUserI "[INTERNAL] ";

%feature("docstring") snoptProblemC::setWorkspace "[INTERNAL] ";

%feature("docstring") snoptProblemC::setJ "[INTERNAL] ";

%feature("docstring") snoptProblemC::setProbName "[INTERNAL] ";

%feature("docstring") snoptProblemC::setSpecsFile "[INTERNAL] ";

%feature("docstring") snoptProblemC::setSTOP "[INTERNAL] ";

%feature("docstring") snoptProblemC::setX "[INTERNAL] ";

%feature("docstring") snoptProblemC::setUserFun "[INTERNAL] ";

%feature("docstring") snoptProblemC::~snoptProblemC "[INTERNAL] ";

%feature("docstring") snoptProblemC::setObjective "[INTERNAL] ";

%feature("docstring") snoptProblemC::setLog "[INTERNAL] ";

%feature("docstring") snoptProblemC::setProblemSize "[INTERNAL] ";

%feature("docstring") snoptProblemC::setIntParameter "[INTERNAL] ";


// File: classcasadi_1_1Solve.xml


// File: classcasadi_1_1SparseStorage.xml
%feature("docstring") casadi::SparseStorage::hasNZ "[INTERNAL]  Returns
true if the matrix has a non-zero at location rr, cc.

";

%feature("docstring") casadi::SparseStorage::elem "[INTERNAL]  get a
reference to an element

";

%feature("docstring") casadi::SparseStorage::clear "[INTERNAL] ";

%feature("docstring") casadi::SparseStorage::reserve "[INTERNAL] ";

%feature("docstring") casadi::SparseStorage::toScalar "[INTERNAL]  Convert
to scalar type.

";

%feature("docstring") casadi::SparseStorage::resize "[INTERNAL] ";

%feature("docstring") casadi::SparseStorage::sparsityRef "[INTERNAL]
Access the sparsity, make a copy if there are multiple references to it.

";

%feature("docstring") casadi::SparseStorage "[INTERNAL] C++ includes:
sparse_storage.hpp ";

%feature("docstring") casadi::SparseStorage::data "

>  [DataType ]  SparseStorage< DataType >.data()
------------------------------------------------------------------------
[INTERNAL] 
Access the non-zero elements.

>  [DataType ]  SparseStorage< DataType >.data() const 
------------------------------------------------------------------------
[INTERNAL] 
Const access the non-zero elements.

";

%feature("docstring") casadi::SparseStorage::sparsity "[INTERNAL]  Const
access the sparsity - reference to data member.

";

%feature("docstring") casadi::SparseStorage::SparseStorage "

>  SparseStorage(Sparsity sparsity, DataType val=DataType(0))
------------------------------------------------------------------------
[INTERNAL] 
Sparse matrix with a given sparsity

>  SparseStorage()
------------------------------------------------------------------------
[INTERNAL] 
constructors

empty 0-by-0 matrix constructor

>  SparseStorage(const SparseStorage< DataType > &m)
------------------------------------------------------------------------
[INTERNAL] 
Copy constructor.

";


// File: classcasadi_1_1Sparsity.xml


/*  Check if two sparsity patterns are identical  */

/*  Size and element counting  */ %feature("docstring")
casadi::Sparsity::enlargeRows "

Enlarge the matrix along the first dimension (i.e. insert rows)

";

%feature("docstring") casadi::Sparsity "

General sparsity class.

The storage format is a compressed column storage (CCS) format.  In this
format, the structural non-zero elements are stored in column-major order,
starting from the upper left corner of the matrix and ending in the lower
right corner.

In addition to the dimension ( size1(), size2()), (i.e. the number of rows
and the number of columns respectively), there are also two vectors of
integers:

\"colind\" [length size2()+1], which contains the index to the first non-
zero element on or after the corresponding column. All the non-zero elements
of a particular i are thus the elements with index el that fulfills:
colind[i] <= el < colind[i+1].

\"row\" [same length as the number of non-zero elements, nnz()] The rows for
each of the structural non-zeros.

Note that with this format, it is cheap to loop over all the non-zero
elements of a particular column, at constant time per element, but expensive
to jump to access a location (i, j).

If the matrix is dense, i.e. length(row) == size1()*size2(), the format
reduces to standard dense column major format, which allows access to an
arbitrary element in constant time.

Since the object is reference counted (it inherits from SharedObject),
several matrices are allowed to share the same sparsity pattern.

The implementations of some methods of this class has been taken from the
CSparse package and modified to use C++ standard library and CasADi data
structures.

See:   Matrix

Joel Andersson

C++ includes: sparsity.hpp ";

%feature("docstring") casadi::Sparsity::largest_first "

Order the columns by decreasing degree.

";

%feature("docstring") casadi::Sparsity::dim "

Get the dimension as a string.

";

%feature("docstring") casadi::Sparsity::issymmetric "

Is symmetric?

";

%feature("docstring") casadi::Sparsity::addNZ "

Get the index of a non-zero element Add the element if it does not exist and
copy object if it's not unique.

";

%feature("docstring") casadi::Sparsity::rowsSequential "

Do the rows appear sequentially on each column.

Parameters:
-----------

strictly:  if true, then do not allow multiple entries

";

%feature("docstring") casadi::Sparsity::get_diag "

Get the diagonal of the matrix/create a diagonal matrix (mapping will
contain the nonzero mapping) When the input is square, the diagonal elements
are returned. If the input is vector-like, a diagonal matrix is constructed
with it.

";

%feature("docstring") casadi::Sparsity::isdiag "

Is diagonal?

";

%feature("docstring") casadi::Sparsity::enlargeColumns "

Enlarge the matrix along the second dimension (i.e. insert columns)

";

%feature("docstring") casadi::Sparsity::hash "";

%feature("docstring") casadi::Sparsity::nnz_upper "

Number of non-zeros in the upper triangular half, i.e. the number of
elements (i, j) with j>=i.

";

%feature("docstring") casadi::Sparsity::getLowerNZ "

Get nonzeros in lower triangular part.

";

%feature("docstring") casadi::Sparsity::find "

Get the location of all non-zero elements as they would appear in a Dense
matrix A : DenseMatrix 4 x 3 B : SparseMatrix 4 x 3 , 5 structural non-
zeros.

k = A.find() A[k] will contain the elements of A that are non-zero in B

";

%feature("docstring") casadi::Sparsity::getUpperNZ "

Get nonzeros in upper triangular part.

";

%feature("docstring") casadi::Sparsity::spy_matlab "

Generate a script for Matlab or Octave which visualizes the sparsity using
the spy command.

";

%feature("docstring") casadi::Sparsity::isdense "

Is dense?

";

%feature("docstring") casadi::Sparsity::isscalar "

Is scalar?

";

%feature("docstring") casadi::Sparsity::repr "

Print a representation of the object.

";

%feature("docstring") casadi::Sparsity::bw_lower "

Lower half-bandwidth.

";

%feature("docstring") casadi::Sparsity::T "

Transpose the matrix.

";

%feature("docstring") casadi::Sparsity::nnz_diag "

Number of non-zeros on the diagonal, i.e. the number of elements (i, j) with
j==i.

";

%feature("docstring") casadi::Sparsity::Sparsity "

>  Sparsity(int dummy=0)
------------------------------------------------------------------------

Default constructor.

>  Sparsity(int nrow, int ncol)
------------------------------------------------------------------------

Pattern with all structural zeros.

>  Sparsity(int nrow, int ncol, [int ] colind, [int ] row)
------------------------------------------------------------------------

Construct from sparsity pattern vectors given in compressed column storage
format.

";

%feature("docstring") casadi::Sparsity::isTranspose "

Check if the sparsity is the transpose of another.

";

%feature("docstring") casadi::Sparsity::colind "

Get a reference to the colindex of column cc (see class description)

";

%feature("docstring") casadi::Sparsity::elimination_tree "

Calculate the elimination tree See Direct Methods for Sparse Linear Systems
by Davis (2006). If the parameter ata is false, the algorithm is equivalent
to Matlab's etree(A), except that the indices are zero- based. If ata is
true, the algorithm is equivalent to Matlab's etree(A, 'row').

";

%feature("docstring") casadi::Sparsity::getCCS "

Get the sparsity in compressed column storage (CCS) format.

";

%feature("docstring") casadi::Sparsity::numel "

The total number of elements, including structural zeros, i.e.
size2()*size1()

See:   nnz()

";

%feature("docstring") casadi::Sparsity::getCRS "

Get the sparsity in compressed row storage (CRS) format.

";

%feature("docstring") casadi::Sparsity::isempty "

Check if the sparsity is empty.

A sparsity is considered empty if one of the dimensions is zero (or
optionally both dimensions)

";

%feature("docstring") casadi::Sparsity::unite "

Union of two sparsity patterns.

";

%feature("docstring") casadi::Sparsity::nnz_lower "

Number of non-zeros in the lower triangular half, i.e. the number of
elements (i, j) with j<=i.

";

%feature("docstring") casadi::Sparsity::getTriplet "

Get the sparsity in sparse triplet format.

";

%feature("docstring") casadi::Sparsity::appendColumns "

Append another sparsity patten horizontally.

";

%feature("docstring") casadi::Sparsity::issingular "

Check whether the sparsity-pattern indicates structural singularity.

";

%feature("docstring") casadi::Sparsity::isReshape "

Check if the sparsity is a reshape of another.

";

%feature("docstring") casadi::Sparsity::removeDuplicates "

Remove duplicate entries.

The same indices will be removed from the mapping vector, which must have
the same length as the number of nonzeros

";

%feature("docstring") casadi::Sparsity::makeDense "

Make a patten dense.

";

%feature("docstring") casadi::Sparsity::pattern_inverse "

Take the inverse of a sparsity pattern; flip zeros and non-zeros.

";

%feature("docstring") casadi::Sparsity::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::Sparsity::iscolumn "

Check if the pattern is a column vector (i.e. size2()==1)

";

%feature("docstring") casadi::Sparsity::dense "

Create a dense rectangular sparsity pattern.

";

%feature("docstring") casadi::Sparsity::diag "

Create diagonal sparsity pattern.

";

%feature("docstring") casadi::Sparsity::sub "

>  Sparsity Sparsity.sub([int ] rr, [int ] cc,[int ] output_mapping, bool ind1=false) const 
------------------------------------------------------------------------

Get a submatrix.

Returns the sparsity of the submatrix, with a mapping such that submatrix[k]
= originalmatrix[mapping[k]]

>  Sparsity Sparsity.sub([int ] rr, Sparsity sp,[int ] output_mapping, bool ind1=false) const 
------------------------------------------------------------------------

Get a set of elements.

Returns the sparsity of the corresponding elements, with a mapping such that
submatrix[k] = originalmatrix[mapping[k]]

";

%feature("docstring") casadi::Sparsity::unit "

Create the sparsity pattern for a unit vector of length n and a nonzero on
position el.

";

%feature("docstring") casadi::Sparsity::append "

Append another sparsity patten vertically (NOTE: only efficient if vector)

";

%feature("docstring") casadi::Sparsity::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::Sparsity::depth_first_search "

Depth-first search on the adjacency graph of the sparsity See Direct Methods
for Sparse Linear Systems by Davis (2006).

";

%feature("docstring") casadi::Sparsity::combine "

Combine two sparsity patterns Returns the new sparsity pattern as well as a
mapping with the same length as the number of non-zero elements The mapping
matrix contains the arguments for each nonzero, the first bit indicates if
the first argument is nonzero, the second bit indicates if the second
argument is nonzero (note that none of, one of or both of the arguments can
be nonzero)

";

%feature("docstring") casadi::Sparsity::nnz "

Get the number of (structural) non-zeros.

See:   numel()

";

%feature("docstring") casadi::Sparsity::intersect "

Intersection of two sparsity patterns Returns the new sparsity pattern as
well as a mapping with the same length as the number of non-zero elements
The value is 1 if the non-zero comes from the first (i.e. this) object, 2 if
it is from the second and 3 (i.e. 1 | 2) if from both.

";

%feature("docstring") casadi::Sparsity::row "

Get the row of a non-zero element.

";

%feature("docstring") casadi::Sparsity::print "

Print a description of the object.

";

%feature("docstring") casadi::Sparsity::get_colind "

Get the column index for each column Together with the row-vector, one
obtains the sparsity pattern in the column compressed format.

";

%feature("docstring") casadi::Sparsity::get_row "

Get the row for each non-zero entry Together with the column-vector, this
vector gives the sparsity of the matrix in sparse triplet format, and
together with the colind vector, one obtains the sparsity in column
compressed format.

";

%feature("docstring") casadi::Sparsity::issquare "

Is square?

";

%feature("docstring") casadi::Sparsity::transpose "

Transpose the matrix and get the reordering of the non-zero entries.

Parameters:
-----------

mapping:  the non-zeros of the original matrix for each non-zero of the new
matrix

";

%feature("docstring") casadi::Sparsity::size "

Get the shape.

";

%feature("docstring") casadi::Sparsity::star_coloring "

Perform a star coloring of a symmetric matrix: A greedy distance-2 coloring
algorithm Algorithm 4.1 in What Color Is Your Jacobian? Graph Coloring for
Computing Derivatives A. H. GEBREMEDHIN, F. MANNE, A. POTHEN SIAM Rev.,
47(4), 629705 (2006)

Ordering options: None (0), largest first (1)

";

%feature("docstring") casadi::Sparsity::hasNZ "

Returns true if the pattern has a non-zero at location rr, cc.

";

%feature("docstring") casadi::Sparsity::print_compact "

Print a compact description of the sparsity pattern.

";

%feature("docstring") casadi::Sparsity::sanity_check "

Check if the dimensions and colind, row vectors are compatible.

Parameters:
-----------

complete:  set to true to also check elementwise throws an error as possible
result

";

%feature("docstring") casadi::Sparsity::compressed "

Create from a single vector containing the pattern in compressed column
storage format: The format: The first two entries are the number of rows
(nrow) and columns (ncol) The next ncol+1 entries are the column offsets
(colind). Note that the last element, colind[ncol], gives the number of
nonzeros The last colind[ncol] entries are the row indices

";

%feature("docstring") casadi::Sparsity::uni_coloring "

Perform a unidirectional coloring: A greedy distance-2 coloring algorithm
(Algorithm 3.1 in A. H. GEBREMEDHIN, F. MANNE, A. POTHEN)

";

%feature("docstring") casadi::Sparsity::isNull "

Is a null pointer?

";

%feature("docstring") casadi::Sparsity::enlarge "

Enlarge matrix Make the matrix larger by inserting empty rows and columns,
keeping the existing non-zeros.

For the matrices A to B A(m, n) length(jj)=m , length(ii)=n B(nrow, ncol)

A=enlarge(m, n, ii, jj) makes sure that

B[jj, ii] == A

";

%feature("docstring") casadi::Sparsity::compress "

Compress a sparsity pattern.

";

%feature("docstring") casadi::Sparsity::istriu "

Is upper triangular?

";

%feature("docstring") casadi::Sparsity::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::Sparsity::isvector "

Check if the pattern is a row or column vector.

";

%feature("docstring") casadi::Sparsity::getNZ "

>  int Sparsity.getNZ(int rr, int cc) const 
------------------------------------------------------------------------

Get the index of an existing non-zero element return -1 if the element does
not exist.

>  [int] Sparsity.getNZ([int ] rr, [int ] cc) const 
------------------------------------------------------------------------

Get a set of non-zero element return -1 if the element does not exist.

>  void Sparsity.getNZ([int ] INOUT) const 
------------------------------------------------------------------------

Get the nonzero index for a set of elements The index vector is used both
for input and outputs and must be sorted by increasing nonzero index, i.e.
column-wise. Elements not found in the sparsity pattern are set to -1.

";

%feature("docstring") casadi::Sparsity::size2 "

Get the number of columns.

";

%feature("docstring") casadi::Sparsity::get_col "

Get the column for each non-zero entry Together with the row-vector, this
vector gives the sparsity of the matrix in sparse triplet format, i.e. the
column and row for each non-zero elements.

";

%feature("docstring") casadi::Sparsity::size1 "

Get the number of rows.

";

%feature("docstring") casadi::Sparsity::star_coloring2 "

Perform a star coloring of a symmetric matrix: A new greedy distance-2
coloring algorithm Algorithm 4.1 in NEW ACYCLIC AND STAR COLORING ALGORITHMS
WITH APPLICATION TO COMPUTING HESSIANS A. H. GEBREMEDHIN, A. TARAFDAR, F.
MANNE, A. POTHEN SIAM J. SCI. COMPUT. Vol. 29, No. 3, pp. 10421072 (2007)

Ordering options: None (0), largest first (1)

";

%feature("docstring") casadi::Sparsity::dulmage_mendelsohn "

Compute the Dulmage-Mendelsohn decomposition See Direct Methods for Sparse
Linear Systems by Davis (2006).

Dulmage-Mendelsohn will try to bring your matrix into lower block-
triangular (LBT) form. It will not care about the distance of off- diagonal
elements to the diagonal: there is no guarantee you will get a block-
diagonal matrix if you supply a randomly permuted block- diagonal matrix.

If your matrix is symmetrical, this method is of limited use; permutation
can make it non-symmetric.

See:   strongly_connected_components

";

%feature("docstring") casadi::Sparsity::istril "

Is lower triangular?

";

%feature("docstring") casadi::Sparsity::spy "

Print a textual representation of sparsity.

";

%feature("docstring") casadi::Sparsity::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::Sparsity::strongly_connected_components "

Find the strongly connected components of the bigraph defined by the
sparsity pattern of a square matrix.

See Direct Methods for Sparse Linear Systems by Davis (2006). Returns:
Number of components

Offset for each components (length: 1 + number of components)

Indices for each components, component i has indices index[offset[i]], ...,
index[offset[i+1]]

In the case that the matrix is symmetric, the result has a particular
interpretation: Given a symmetric matrix A and n =
A.strongly_connected_components(p, r)

=> A[p, p] will appear block-diagonal with n blocks and with the indices of
the block boundaries to be found in r.

";

%feature("docstring") casadi::Sparsity::pmult "

Permute rows and/or columns Multiply the sparsity with a permutation matrix
from the left and/or from the right P * A * trans(P), A * trans(P) or A *
trans(P) with P defined by an index vector containing the row for each col.
As an alternative, P can be transposed (inverted).

";

%feature("docstring") casadi::Sparsity::isrow "

Check if the pattern is a row vector (i.e. size1()==1)

";

%feature("docstring") casadi::Sparsity::is_equal "";

%feature("docstring") casadi::Sparsity::erase "

>  [int] Sparsity.erase([int ] rr, [int ] cc, bool ind1=false)
------------------------------------------------------------------------

Erase rows and/or columns of a matrix.

>  [int] Sparsity.erase([int ] rr, bool ind1=false)
------------------------------------------------------------------------

Erase elements of a matrix.

";

%feature("docstring") casadi::Sparsity::bw_upper "

Upper half-bandwidth.

";

%feature("docstring") casadi::Sparsity::scalar "

Create a scalar sparsity pattern.

";

%feature("docstring") casadi::Sparsity::resize "

Resize.

";


// File: classcasadi_1_1SparsityInterface.xml
%feature("docstring") friendwrap_diagsplit "

>  [MatType ] diagsplit(MatType x, [int ] output_offset1, [int ] output_offset2)
------------------------------------------------------------------------

split diagonally, retaining square matrices

Parameters:
-----------

output_offset1:  List of all start locations (row) for each group the last
matrix will run to the end.

output_offset2:  List of all start locations (row) for each group the last
matrix will run to the end.

diagcat(diagsplit(x, ...)) = x

>  [MatType ] diagsplit(MatType x, [int ] output_offset)
------------------------------------------------------------------------

split diagonally, retaining square matrices

Parameters:
-----------

output_offset:  List of all start locations for each group the last matrix
will run to the end.

diagcat(diagsplit(x, ...)) = x

>  [MatType ] diagsplit(MatType x, int incr=1)
------------------------------------------------------------------------

split diagonally, retaining groups of square matrices

Parameters:
-----------

incr:  Size of each matrix

diagsplit(diagsplit(x, ...)) = x

>  [MatType ] diagsplit(MatType x, int incr1, int incr2)
------------------------------------------------------------------------

split diagonally, retaining fixed-sized matrices

Parameters:
-----------

incr1:  Row dimension of each matrix

incr2:  Column dimension of each matrix

diagsplit(diagsplit(x, ...)) = x

";

%feature("docstring") friendwrap_triu "

Get the upper triangular part of a matrix.

";

%feature("docstring") friendwrap_mac "

Multiply-accumulate operation Matrix product of two matrices (X and Y),
adding the result to a third matrix Z. The result has the same sparsity
pattern as C meaning that other entries of (X*Y) are ignored. The operation
is equivalent to: Z+mul(X,Y).project(Z.sparsity()).

";

%feature("docstring") friendwrap_transpose "

Transpose.

";

%feature("docstring") friendwrap_tril "

Get the lower triangular part of a matrix.

";

%feature("docstring") friendwrap_offset "

Helper function, get offsets corresponding to a vector of matrices.

";

%feature("docstring") friendwrap_vec "

make a vector Reshapes/vectorizes the matrix such that the shape becomes
(expr.numel(), 1). Columns are stacked on top of each other. Same as
reshape(expr, expr.numel(), 1)

a c b d  turns into

a b c d

";

%feature("docstring") friendwrap_horzcat "

>  MatType horzcat([MatType ] v)
------------------------------------------------------------------------

Concatenate a list of matrices horizontally Alternative terminology:
horizontal stack, hstack, horizontal append, [a b].

horzcat(horzsplit(x, ...)) = x

>  MatType horzcat(MatType x, MatType y)
------------------------------------------------------------------------

Concatenate horizontally, two matrices.

>  MatType horzcat(MatType x, MatType y, MatType z)
------------------------------------------------------------------------

Concatenate horizontally, three matrices.

>  MatType horzcat(MatType x, MatType y, MatType z, MatType w)
------------------------------------------------------------------------

Concatenate horizontally, four matrices.

";

%feature("docstring") friendwrap_vecNZ "

Returns a flattened version of the matrix, preserving only nonzeros.

";

%feature("docstring") casadi::SparsityInterface "

Sparsity interface class.

This is a common base class for GenericMatrix (i.e. MX and Matrix<>) and
Sparsity, introducing a uniform syntax and implementing common functionality
using the curiously recurring template pattern (CRTP) idiom. Joel Andersson

C++ includes: sparsity_interface.hpp ";

%feature("docstring") friendwrap_horzsplit "

>  [MatType ] horzsplit(MatType v, [int ] offset)
------------------------------------------------------------------------

split horizontally, retaining groups of columns

Parameters:
-----------

offset:  List of all start columns for each group the last column group will
run to the end.

horzcat(horzsplit(x, ...)) = x

>  [MatType ] horzsplit(MatType v, int incr=1)
------------------------------------------------------------------------

split horizontally, retaining fixed-sized groups of columns

Parameters:
-----------

incr:  Size of each group of columns

horzcat(horzsplit(x, ...)) = x

";

%feature("docstring") friendwrap_veccat "

concatenate vertically while vectorizing all arguments with vec

";

%feature("docstring") friendwrap_blocksplit "

>  [[MatType ] ] blocksplit(MatType x, [int ] vert_offset, [int ] horz_offset)
------------------------------------------------------------------------

chop up into blocks

Parameters:
-----------

vert_offset:  Defines the boundaries of the block rows

horz_offset:  Defines the boundaries of the block columns

blockcat(blocksplit(x,..., ...)) = x

>  [[MatType ] ] blocksplit(MatType x, int vert_incr=1, int horz_incr=1)
------------------------------------------------------------------------

chop up into blocks

Parameters:
-----------

vert_incr:  Defines the increment for block boundaries in row dimension

horz_incr:  Defines the increment for block boundaries in column dimension

blockcat(blocksplit(x,..., ...)) = x

";

%feature("docstring") friendwrap_repmat "

Repeat matrix A n times vertically and m times horizontally.

";

%feature("docstring") friendwrap_vertcat "

>  MatType vertcat([MatType ] v)
------------------------------------------------------------------------

Concatenate a list of matrices vertically Alternative terminology: vertical
stack, vstack, vertical append, [a;b].

vertcat(vertsplit(x, ...)) = x

>  MatType vertcat(MatType x, MatType y)
------------------------------------------------------------------------

Concatenate vertically, two matrices.

>  MatType vertcat(MatType x, MatType y, MatType z)
------------------------------------------------------------------------

Concatenate vertically, three matrices.

>  MatType vertcat(MatType x, MatType y, MatType z, MatType w)
------------------------------------------------------------------------

Concatenate vertically, four matrices.

";

%feature("docstring") friendwrap_sprank "

Obtain the structural rank of a sparsity-pattern.

";

%feature("docstring") friendwrap_kron "

Kronecker tensor product.

Creates a block matrix in which each element (i, j) is a_ij*b

";

%feature("docstring") friendwrap_reshape "

>  MatType reshape(MatType a, int nrow, int ncol)
------------------------------------------------------------------------

Returns a reshaped version of the matrix.

>  MatType reshape(MatType a,(int,int) rc)
------------------------------------------------------------------------

Returns a reshaped version of the matrix, dimensions as a vector.

>  MatType reshape(MatType a, Sparsity sp)
------------------------------------------------------------------------

Reshape the matrix.

";

%feature("docstring") friendwrap_norm_0_mul "

0-norm (nonzero count) of a Matrix-matrix product

";

%feature("docstring") friendwrap_diagcat "

>  MatType diagcat([MatType ] A)
------------------------------------------------------------------------

Construct a matrix with given block on the diagonal.

>  MatType diagcat(MatType x, MatType y)
------------------------------------------------------------------------

Concatenate along diagonal, two matrices.

>  MatType diagcat(MatType x, MatType y, MatType z)
------------------------------------------------------------------------

Concatenate along diagonal, three matrices.

>  MatType diagcat(MatType x, MatType y, MatType z, MatType w)
------------------------------------------------------------------------

Concatenate along diagonal, four matrices.

";

%feature("docstring") friendwrap_vertsplit "

>  [MatType ] vertsplit(MatType v, [int ] offset)
------------------------------------------------------------------------

split vertically, retaining groups of rows

*

Parameters:
-----------

output_offset:  List of all start rows for each group the last row group
will run to the end.

vertcat(vertsplit(x, ...)) = x

>  [MatType ] vertsplit(MatType v, int incr=1)
------------------------------------------------------------------------

split vertically, retaining fixed-sized groups of rows

Parameters:
-----------

incr:  Size of each group of rows

vertcat(vertsplit(x, ...)) = x



::

  >>> print vertsplit(SX.sym(\"a\",4))
  [SX(a_0), SX(a_1), SX(a_2), SX(a_3)]
  





::

  >>> print vertsplit(SX.sym(\"a\",4),2)
  [SX([a_0, a_1]), SX([a_2, a_3])]
  



If the number of rows is not a multiple of incr, the last entry returned
will have a size smaller than incr.



::

  >>> print vertsplit(DMatrix([0,1,2,3,4]),2)
  [DMatrix([0, 1]), DMatrix([2, 3]), DMatrix(4)]
  



";

%feature("docstring") friendwrap_mul "

>  MatType mul(MatType X, MatType Y)
------------------------------------------------------------------------

Matrix product of two matrices.

>  MatType mul([MatType ] args)
------------------------------------------------------------------------

Matrix product of n matrices.

";

%feature("docstring") friendwrap_blockcat "

>  MatType blockcat([[MatType ] ] v)
------------------------------------------------------------------------

Construct a matrix from a list of list of blocks.

>  MatType blockcat(MatType A, MatType B, MatType C, MatType D)
------------------------------------------------------------------------

Construct a matrix from 4 blocks.

";


// File: classcasadi_1_1Split.xml


// File: classcasadi_1_1Sqpmethod.xml


// File: classcasadi_1_1StabilizedQpSolver.xml


/*  Simple Getters & Setters  */

/*  Advanced Getters  */

/*  Option Functionality  */ %feature("docstring")
casadi::StabilizedQpSolver::n_in "

Get the number of function inputs.

";

%feature("docstring") casadi::StabilizedQpSolver::mapsum "

Evaluate symbolically in parallel and sum (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

";

%feature("docstring") casadi::StabilizedQpSolver::derForward "

Get a function that calculates nfwd forward derivatives.

Returns a function with n_in + n_out +nfwd*n_in inputs and nfwd*n_out
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nfwd*n_in inputs correspond to forward seeds, one direction at a time The
nfwd*n_out outputs correspond to forward sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::StabilizedQpSolver::printDimensions "

Print dimensions of inputs and outputs.

";

%feature("docstring") casadi::StabilizedQpSolver::getOptionDescription "

Get the description of a certain option.

";

%feature("docstring") casadi::StabilizedQpSolver::getAtomicInputReal "

Get the floating point output argument of an atomic operation.

";

%feature("docstring") casadi::StabilizedQpSolver::nnz_out "

Get of number of output nonzeros For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::StabilizedQpSolver::jacobian "

Generate a Jacobian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. If
compact is set to true, only the nonzeros of the input and output
expressions are considered. If symmetric is set to true, the Jacobian being
calculated is known to be symmetric (usually a Hessian), which can be
exploited by the algorithm.

The generated Jacobian has one more output than the calling function
corresponding to the Jacobian and the same number of inputs.

";

%feature("docstring") casadi::StabilizedQpSolver::index_out "

Find the index for a string describing a particular entry of an output
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::StabilizedQpSolver::setFullJacobian "

Set the Jacobian of all the input nonzeros with respect to all output
nonzeros NOTE: Does not take ownership, only weak references to the Jacobian
are kept internally

";

%feature("docstring") casadi::StabilizedQpSolver::getInput "

>  DMatrix  IOInterface< Function  >.getInput(int iind=0) const
------------------------------------------------------------------------

Get an input by index.

Parameters:
-----------

iind:  index within the range [0..n_in()-1]

>  DMatrix  IOInterface< Function  >.getInput(str iname) const
------------------------------------------------------------------------

Get an input by name.

Parameters:
-----------

iname:  input name. Only allowed when an input scheme is set.

>  void IOInterface< Function  >.getInput(T val, int iind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.getInput(T val, str iname)
------------------------------------------------------------------------
[INTERNAL] 
Get an input by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::StabilizedQpSolver::setOptionByAllowedIndex "[INTERNAL]  Set a certain option by giving its index into the allowed
values.

";

%feature("docstring") casadi::StabilizedQpSolver::gradient "

Generate a gradient function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the output must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::StabilizedQpSolver::derReverse "

Get a function that calculates nadj adjoint derivatives.

Returns a function with n_in + n_out +nadj*n_out inputs and nadj*n_in
outputs. The first n_in inputs correspond to nondifferentiated inputs. The
next n_out inputs correspond to nondifferentiated outputs. and the last
nadj*n_out inputs correspond to adjoint seeds, one direction at a time The
nadj*n_in outputs correspond to adjoint sensitivities, one direction at a
time. * (n_in = n_in(), n_out = n_out())

(n_in = n_in(), n_out = n_out())

The functions returned are cached, meaning that if called multiple timed
with the same value, then multiple references to the same function will be
returned.

";

%feature("docstring") casadi::StabilizedQpSolver::mx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::StabilizedQpSolver::size1_out "

Get output dimension.

";

%feature("docstring") casadi::StabilizedQpSolver::setDerReverse "

Set a function that calculates nadj adjoint derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::StabilizedQpSolver::evaluate "

Evaluate.

";

%feature("docstring") casadi::StabilizedQpSolver::getOptionType "

Get the type of a certain option.

";

%feature("docstring") casadi::StabilizedQpSolver::jacSparsity "

Get, if necessary generate, the sparsity of a Jacobian block

";

%feature("docstring") casadi::StabilizedQpSolver::setJacobian "

Set the Jacobian function of output oind with respect to input iind NOTE:
Does not take ownership, only weak references to the Jacobians are kept
internally

";

%feature("docstring") casadi::StabilizedQpSolver::generateLiftingFunctions "

Extract the functions needed for the Lifted Newton method.

";

%feature("docstring") casadi::StabilizedQpSolver::setOutput "

>  void IOInterface< Function  >.setOutput(T val, int oind=0)
------------------------------------------------------------------------

Set an output by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.setOutput(T val, str oname)
------------------------------------------------------------------------

Set an output by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::StabilizedQpSolver::removeMonitor "

Remove modules to be monitored.

";

%feature("docstring") casadi::StabilizedQpSolver::getOptionAllowed "

Get the allowed values of a certain option.

";

%feature("docstring") casadi::StabilizedQpSolver::getStat "

Get a single statistic obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::StabilizedQpSolver::setOptionByEnumValue "[INTERNAL]  Set a certain option by giving an enum value.

";

%feature("docstring") casadi::StabilizedQpSolver::call "

Evaluate the function symbolically or numerically.

";

%feature("docstring") casadi::StabilizedQpSolver::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::StabilizedQpSolver::nnz_in "

Get of number of input nonzeros For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::StabilizedQpSolver::getOptionDefault "

Get the default of a certain option.

";

%feature("docstring") casadi::StabilizedQpSolver::countNodes "

Number of nodes in the algorithm.

";

%feature("docstring") casadi::StabilizedQpSolver::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::StabilizedQpSolver::fullJacobian "

Generate a Jacobian function of all the inputs elements with respect to all
the output elements).

";

%feature("docstring") casadi::StabilizedQpSolver::printPtr "[INTERNAL]
Print the pointer to the internal class

";

%feature("docstring") casadi::StabilizedQpSolver "

StabilizedQpSolver.

Solves the following strictly convex problem:



::

  min          1/2 x' H x + g' x
  x
  
  subject to
  LBA <= A x <= UBA
  LBX <= x   <= UBX
  
  with :
  H sparse (n x n) positive definite
  g dense  (n x 1)
  
  n: number of decision variables (x)
  nc: number of constraints (A)



If H is not positive-definite, the solver should throw an error.

General information
===================



>Input scheme: casadi::StabilizedQpSolverInput (STABILIZED_QP_SOLVER_NUM_IN = 12) [stabilizedQpIn]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| STABILIZED_QP_SOLVER_H | h                      | The square matrix H:   |
|                        |                        | sparse, (n x n). Only  |
|                        |                        | the lower triangular   |
|                        |                        | part is actually used. |
|                        |                        | The matrix is assumed  |
|                        |                        | to be symmetrical.     |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_G | g                      | The vector g: dense,   |
|                        |                        | (n x 1) .              |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_A | a                      | The matrix A: sparse,  |
|                        |                        | (nc x n) - product     |
|                        |                        | with x must be dense.  |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_L | lba                    | dense, (nc x 1)        |
| BA                     |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_U | uba                    | dense, (nc x 1)        |
| BA                     |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_L | lbx                    | dense, (n x 1)         |
| BX                     |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_U | ubx                    | dense, (n x 1)         |
| BX                     |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_X | x0                     | dense, (n x 1)         |
| 0                      |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_L | lam_x0                 | dense                  |
| AM_X0                  |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_M | muR                    | dense (1 x 1)          |
| UR                     |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_M | muE                    | dense (nc x 1)         |
| UE                     |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_M | mu                     | dense (nc x 1)         |
| U                      |                        |                        |
+------------------------+------------------------+------------------------+

>Output scheme: casadi::QpSolverOutput (QP_SOLVER_NUM_OUT = 4) [qpOut]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| QP_SOLVER_X            | x                      | The primal solution .  |
+------------------------+------------------------+------------------------+
| QP_SOLVER_COST         | cost                   | The optimal cost .     |
+------------------------+------------------------+------------------------+
| QP_SOLVER_LAM_A        | lam_a                  | The dual solution      |
|                        |                        | corresponding to       |
|                        |                        | linear bounds .        |
+------------------------+------------------------+------------------------+
| QP_SOLVER_LAM_X        | lam_x                  | The dual solution      |
|                        |                        | corresponding to       |
|                        |                        | simple bounds .        |
+------------------------+------------------------+------------------------+

>List of available options

+--------------+--------------+--------------+--------------+--------------+
|      Id      |     Type     |   Default    | Description  |   Used in    |
+==============+==============+==============+==============+==============+
| ad_weight    | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | derivative c |              |
|              |              |              | alculation.W |              |
|              |              |              | hen there is |              |
|              |              |              | an option of |              |
|              |              |              | either using |              |
|              |              |              | forward or   |              |
|              |              |              | reverse mode |              |
|              |              |              | directional  |              |
|              |              |              | derivatives, |              |
|              |              |              | the          |              |
|              |              |              | condition ad |              |
|              |              |              | _weight*nf<= |              |
|              |              |              | (1-ad_weight |              |
|              |              |              | )*na is used |              |
|              |              |              | where nf and |              |
|              |              |              | na are       |              |
|              |              |              | estimates of |              |
|              |              |              | the number   |              |
|              |              |              | of forward/r |              |
|              |              |              | everse mode  |              |
|              |              |              | directional  |              |
|              |              |              | derivatives  |              |
|              |              |              | needed. By   |              |
|              |              |              | default,     |              |
|              |              |              | ad_weight is |              |
|              |              |              | calculated a |              |
|              |              |              | utomatically |              |
|              |              |              | , but this   |              |
|              |              |              | can be       |              |
|              |              |              | overridden   |              |
|              |              |              | by setting   |              |
|              |              |              | this option. |              |
|              |              |              | In           |              |
|              |              |              | particular,  |              |
|              |              |              | 0 means      |              |
|              |              |              | forcing      |              |
|              |              |              | forward mode |              |
|              |              |              | and 1        |              |
|              |              |              | forcing      |              |
|              |              |              | reverse      |              |
|              |              |              | mode. Leave  |              |
|              |              |              | unset for    |              |
|              |              |              | (class       |              |
|              |              |              | specific)    |              |
|              |              |              | heuristics.  |              |
+--------------+--------------+--------------+--------------+--------------+
| ad_weight_sp | OT_REAL      | GenericType( | Weighting    | casadi::Func |
|              |              | )            | factor for   | tionInternal |
|              |              |              | sparsity     |              |
|              |              |              | pattern      |              |
|              |              |              | calculation  |              |
|              |              |              | calculation. |              |
|              |              |              | Overrides    |              |
|              |              |              | default      |              |
|              |              |              | behavior.    |              |
|              |              |              | Set to 0 and |              |
|              |              |              | 1 to force   |              |
|              |              |              | forward and  |              |
|              |              |              | reverse mode |              |
|              |              |              | respectively |              |
|              |              |              | . Cf. option |              |
|              |              |              | \"ad_weight\". |              |
+--------------+--------------+--------------+--------------+--------------+
| compiler     | OT_STRING    | \"clang\"      | Just-in-time | casadi::Func |
|              |              |              | compiler     | tionInternal |
|              |              |              | plugin to be |              |
|              |              |              | used.        |              |
+--------------+--------------+--------------+--------------+--------------+
| defaults_rec | OT_STRINGVEC | GenericType( | Changes      | casadi::Opti |
| ipes         | TOR          | )            | default      | onsFunctiona |
|              |              |              | options      | lityNode   c |
|              |              |              | according to | asadi::Stabi |
|              |              |              | a given      | lizedQpSolve |
|              |              |              | recipe (low- | rInternal    |
|              |              |              | level)  (qp) |              |
+--------------+--------------+--------------+--------------+--------------+
| gather_stats | OT_BOOLEAN   | false        | Flag to      | casadi::Func |
|              |              |              | indicate     | tionInternal |
|              |              |              | whether      |              |
|              |              |              | statistics   |              |
|              |              |              | must be      |              |
|              |              |              | gathered     |              |
+--------------+--------------+--------------+--------------+--------------+
| input_scheme | OT_STRINGVEC | GenericType( | Custom input | casadi::Func |
|              | TOR          | )            | scheme       | tionInternal |
+--------------+--------------+--------------+--------------+--------------+
| inputs_check | OT_BOOLEAN   | true         | Throw        | casadi::Func |
|              |              |              | exceptions   | tionInternal |
|              |              |              | when the     |              |
|              |              |              | numerical    |              |
|              |              |              | values of    |              |
|              |              |              | the inputs   |              |
|              |              |              | don't make   |              |
|              |              |              | sense        |              |
+--------------+--------------+--------------+--------------+--------------+
| jac_penalty  | OT_REAL      | 2            | When         | casadi::Func |
|              |              |              | requested    | tionInternal |
|              |              |              | for a number |              |
|              |              |              | of forward/r |              |
|              |              |              | everse       |              |
|              |              |              | directions,  |              |
|              |              |              | it may be    |              |
|              |              |              | cheaper to   |              |
|              |              |              | compute      |              |
|              |              |              | first the    |              |
|              |              |              | full         |              |
|              |              |              | jacobian and |              |
|              |              |              | then         |              |
|              |              |              | multiply     |              |
|              |              |              | with seeds,  |              |
|              |              |              | rather than  |              |
|              |              |              | obtain the   |              |
|              |              |              | requested    |              |
|              |              |              | directions   |              |
|              |              |              | in a straigh |              |
|              |              |              | tforward     |              |
|              |              |              | manner.      |              |
|              |              |              | Casadi uses  |              |
|              |              |              | a heuristic  |              |
|              |              |              | to decide    |              |
|              |              |              | which is     |              |
|              |              |              | cheaper. A   |              |
|              |              |              | high value   |              |
|              |              |              | of 'jac_pena |              |
|              |              |              | lty' makes   |              |
|              |              |              | it less      |              |
|              |              |              | likely for   |              |
|              |              |              | the heurstic |              |
|              |              |              | to chose the |              |
|              |              |              | full         |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy.    |              |
|              |              |              | The special  |              |
|              |              |              | value -1     |              |
|              |              |              | indicates    |              |
|              |              |              | never to use |              |
|              |              |              | the full     |              |
|              |              |              | Jacobian     |              |
|              |              |              | strategy     |              |
+--------------+--------------+--------------+--------------+--------------+
| jit          | OT_BOOLEAN   | false        | Use just-in- | casadi::Func |
|              |              |              | time         | tionInternal |
|              |              |              | compiler to  |              |
|              |              |              | speed up the |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| jit_options  | OT_DICT      | GenericType( | Options to   | casadi::Func |
|              |              | )            | be passed to | tionInternal |
|              |              |              | the jit      |              |
|              |              |              | compiler.    |              |
+--------------+--------------+--------------+--------------+--------------+
| monitor      | OT_STRINGVEC | GenericType( | Monitors to  | casadi::Func |
|              | TOR          | )            | be activated | tionInternal |
|              |              |              | (inputs|outp |              |
|              |              |              | uts)         |              |
+--------------+--------------+--------------+--------------+--------------+
| output_schem | OT_STRINGVEC | GenericType( | Custom       | casadi::Func |
| e            | TOR          | )            | output       | tionInternal |
|              |              |              | scheme       |              |
+--------------+--------------+--------------+--------------+--------------+
| regularity_c | OT_BOOLEAN   | true         | Throw        | casadi::Func |
| heck         |              |              | exceptions   | tionInternal |
|              |              |              | when NaN or  |              |
|              |              |              | Inf appears  |              |
|              |              |              | during       |              |
|              |              |              | evaluation   |              |
+--------------+--------------+--------------+--------------+--------------+
| user_data    | OT_VOIDPTR   | GenericType( | A user-      | casadi::Func |
|              |              | )            | defined      | tionInternal |
|              |              |              | field that   |              |
|              |              |              | can be used  |              |
|              |              |              | to identify  |              |
|              |              |              | the function |              |
|              |              |              | or pass      |              |
|              |              |              | additional   |              |
|              |              |              | information  |              |
+--------------+--------------+--------------+--------------+--------------+
| verbose      | OT_BOOLEAN   | false        | Verbose      | casadi::Func |
|              |              |              | evaluation   | tionInternal |
|              |              |              | for          |              |
|              |              |              | debugging    |              |
+--------------+--------------+--------------+--------------+--------------+

List of plugins
===============



- sqic

- qp

Note: some of the plugins in this list might not be available on your
system. Also, there might be extra plugins available to you that are not
listed here. You can obtain their documentation with
StabilizedQpSolver.doc(\"myextraplugin\")



--------------------------------------------------------------------------------

sqic
----



Interface to SQIC

>List of available options

+----+------+---------+-------------+
| Id | Type | Default | Description |
+====+======+=========+=============+
+----+------+---------+-------------+

--------------------------------------------------------------------------------



--------------------------------------------------------------------------------

qp --



Solved a stabilized QP using a standard QP solver

>List of available options

+-----------------+-----------------+-----------------+-----------------+
|       Id        |      Type       |     Default     |   Description   |
+=================+=================+=================+=================+
| qp_solver       | OT_STRING       | GenericType()   | The QP solver   |
|                 |                 |                 | used to solve   |
|                 |                 |                 | the stabilized  |
|                 |                 |                 | QPs.            |
+-----------------+-----------------+-----------------+-----------------+
| qp_solver_optio | OT_DICT         | GenericType()   | Options to be   |
| ns              |                 |                 | passed to the   |
|                 |                 |                 | QP solver       |
|                 |                 |                 | instance        |
+-----------------+-----------------+-----------------+-----------------+

>List of available stats

+-----------------+
|       Id        |
+=================+
| qp_solver_stats |
+-----------------+

--------------------------------------------------------------------------------



Joel Andersson
Diagrams
--------



C++ includes: stabilized_qp_solver.hpp ";

%feature("docstring") casadi::StabilizedQpSolver::getAtomicOperation "

Get an atomic operation operator index.

";

%feature("docstring") casadi::StabilizedQpSolver::name "

Name of the function.

";

%feature("docstring") casadi::StabilizedQpSolver::addMonitor "

Add modules to be monitored.

";

%feature("docstring") casadi::StabilizedQpSolver::sz_res "[INTERNAL]  Get
required length of res field.

";

%feature("docstring") casadi::StabilizedQpSolver::copyOptions "

Copy all options from another object.

";

%feature("docstring") casadi::StabilizedQpSolver::default_in "

Get default input value (NOTE: constant reference)

";

%feature("docstring") casadi::StabilizedQpSolver::size1_in "

Get input dimension.

";

%feature("docstring") casadi::StabilizedQpSolver::getSanitizedName "

get function name with all non alphanumeric characters converted to '_'

";

%feature("docstring") casadi::StabilizedQpSolver::generate "

Export / Generate C code for the function.

";

%feature("docstring") casadi::StabilizedQpSolver::callForward "

Create call to (cached) derivative function, forward mode.

";

%feature("docstring") casadi::StabilizedQpSolver::spInit "[INTERNAL]  Reset
the sparsity propagation.

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::StabilizedQpSolver::dictionary "

Get the dictionary.

";

%feature("docstring") casadi::StabilizedQpSolver::size2_out "

Get output dimension.

";

%feature("docstring") casadi::StabilizedQpSolver::sparsity_in "

Get sparsity of a given input.

";

%feature("docstring") casadi::StabilizedQpSolver::setJacSparsity "

Generate the sparsity of a Jacobian block

";

%feature("docstring") casadi::StabilizedQpSolver::n_out "

Get the number of function outputs.

";

%feature("docstring") casadi::StabilizedQpSolver::sz_arg "[INTERNAL]  Get
required length of arg field.

";

%feature("docstring") casadi::StabilizedQpSolver::sz_w "[INTERNAL]  Get
required length of w field.

";

%feature("docstring") casadi::StabilizedQpSolver::isNull "

Is a null pointer?

";

%feature("docstring") casadi::StabilizedQpSolver::getOutput "

>  DMatrix  IOInterface< Function  >.getOutput(int oind=0) const
------------------------------------------------------------------------

Get an output by index.

Parameters:
-----------

oind:  index within the range [0..n_out()-1]

>  DMatrix  IOInterface< Function  >.getOutput(str oname) const
------------------------------------------------------------------------

Get an output by name.

Parameters:
-----------

oname:  output name. Only allowed when an output scheme is set.

>  void IOInterface< Function  >.getOutput(T val, int oind=0)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by index.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oind:  index within the range [0..n_out()-1]

>  void IOInterface< Function  >.getOutput(T val, str oname)
------------------------------------------------------------------------
[INTERNAL] 
Get an output by name.

Parameters:
-----------

val:  can be double&, std::vector<double>&, Matrix<double>&, double *

oname:  output name. Only allowed when an output scheme is set.

";

%feature("docstring") casadi::StabilizedQpSolver::free_sx "

Get all the free variables of the function.

";

%feature("docstring") casadi::StabilizedQpSolver::derivative "

Get a function that calculates nfwd forward derivatives and nadj adjoint
derivatives Legacy function: Use derForward and derReverse instead.

Returns a function with (1+nfwd)*n_in+nadj*n_out inputs and (1+nfwd)*n_out +
nadj*n_in outputs. The first n_in inputs correspond to nondifferentiated
inputs. The next nfwd*n_in inputs correspond to forward seeds, one direction
at a time and the last nadj*n_out inputs correspond to adjoint seeds, one
direction at a time. The first n_out outputs correspond to nondifferentiated
outputs. The next nfwd*n_out outputs correspond to forward sensitivities,
one direction at a time and the last nadj*n_in outputs corresponds to
adjoint sensitivities, one direction at a time.

(n_in = n_in(), n_out = n_out())

";

%feature("docstring") casadi::StabilizedQpSolver::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::StabilizedQpSolver::getAtomicInput "

Get the (integer) input arguments of an atomic operation.

";

%feature("docstring") casadi::StabilizedQpSolver::mx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::StabilizedQpSolver::tangent "

Generate a tangent function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The default behavior of this class is defined by the derived class. Note
that the input must be scalar. In other cases, use the Jacobian instead.

";

%feature("docstring") casadi::StabilizedQpSolver::hessian "

Generate a Hessian function of output oind with respect to input iind.

Parameters:
-----------

iind:  The index of the input

oind:  The index of the output

The generated Hessian has two more outputs than the calling function
corresponding to the Hessian and the gradients.

";

%feature("docstring") casadi::StabilizedQpSolver::name_out "

>  [str] Function.name_out() const 
------------------------------------------------------------------------

Get output scheme.

>  str Function.name_out(int ind) const 
------------------------------------------------------------------------

Get output scheme name by index.

";

%feature("docstring") casadi::StabilizedQpSolver::numel_in "

Get of number of input elements For a particular input or for all for all of
the inputs.

";

%feature("docstring") casadi::StabilizedQpSolver::getOptionEnumValue "[INTERNAL]  Get the enum value corresponding to th certain option.

";

%feature("docstring") casadi::StabilizedQpSolver::callDerivative "[INTERNAL]  Evaluate the function symbolically or numerically with
directional derivatives The first two arguments are the nondifferentiated
inputs and results of the evaluation, the next two arguments are a set of
forward directional seeds and the resulting forward directional derivatives,
the length of the vector being the number of forward directions. The next
two arguments are a set of adjoint directional seeds and the resulting
adjoint directional derivatives, the length of the vector being the number
of adjoint directions.

";

%feature("docstring") casadi::StabilizedQpSolver::sparsity_out "

Get sparsity of a given output.

";

%feature("docstring") casadi::StabilizedQpSolver::StabilizedQpSolver "

>  StabilizedQpSolver()
------------------------------------------------------------------------

Default constructor.

>  StabilizedQpSolver(str name, str solver, const std.map< str, Sparsity > &st, Dict opts=Dict())
------------------------------------------------------------------------

Constructor (new syntax, includes initialization)

Parameters:
-----------

solver:

Name of a solver. It might be one of:

- sqic

- qp

Note: some of the plugins in this list might not be available on your
system. Also, there might be extra plugins available to you that are not
listed here. You can obtain their documentation with
StabilizedQpSolver.doc(\"myextraplugin\")

st:  Problem structure

>Struct scheme: casadi::QPStruct ( = 2) []

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| QP_STRUCT_H            |                        | The square matrix H:   |
|                        |                        | sparse, (n x n). Only  |
|                        |                        | the lower triangular   |
|                        |                        | part is actually used. |
|                        |                        | The matrix is assumed  |
|                        |                        | to be symmetrical.     |
+------------------------+------------------------+------------------------+
| QP_STRUCT_A            |                        | The matrix A: sparse,  |
|                        |                        | (nc x n) - product     |
|                        |                        | with x must be dense.  |
+------------------------+------------------------+------------------------+

";

%feature("docstring") casadi::StabilizedQpSolver::numel_out "

Get of number of output elements For a particular output or for all for all
of the outputs.

";

%feature("docstring") casadi::StabilizedQpSolver::printOptions "

Print options to a stream.

";

%feature("docstring") casadi::StabilizedQpSolver::getAtomicOutput "

Get the (integer) output argument of an atomic operation.

";

%feature("docstring") casadi::StabilizedQpSolver::setInput "

>  void IOInterface< Function  >.setInput(T val, int iind=0)
------------------------------------------------------------------------

Set an input by index.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iind:  index within the range [0..n_in()-1]

>  void IOInterface< Function  >.setInput(T val, str iname)
------------------------------------------------------------------------

Set an input by name.

Parameters:
-----------

val:  can be double, const std::vector<double>&, const Matrix<double>&,
double *

iname:  input name. Only allowed when an input scheme is set.

";

%feature("docstring") casadi::StabilizedQpSolver::size_out "

Get output dimension.

";

%feature("docstring") casadi::StabilizedQpSolver::description_out "

Get output scheme description by index.

";

%feature("docstring") casadi::StabilizedQpSolver::free_mx "

Get all the free variables of the function.

";

%feature("docstring") casadi::StabilizedQpSolver::sx_in "

Get symbolic primitives equivalent to the input expressions.

";

%feature("docstring") casadi::StabilizedQpSolver::sz_iw "[INTERNAL]  Get
required length of iw field.

";

%feature("docstring") casadi::StabilizedQpSolver::is_a "

Check if the function is of a particular type Optionally check if name
matches one of the base classes (default true)

";

%feature("docstring") casadi::StabilizedQpSolver::map "

>  [[MX] ] Function.map([[MX ] ] arg, str parallelization=\"serial\")

>  [MX] Function.map([MX ] arg, str parallelization=\"serial\")
------------------------------------------------------------------------

Evaluate symbolically in parallel (matrix graph)

Parameters:
-----------

parallelization:  Type of parallelization used: expand|serial|openmp

>  Function Function.map(str name, int N, Dict options=Dict()) const 
------------------------------------------------------------------------

Create a mapped version of this function.

Suppose the function has a signature of:

::

     f: (a, p) -> ( s )
  



The the mapaccumulated version has the signature:

::

     F: (A, P) -> (S )
  
      with
          a: horzcat([a0, a1, ..., a_(N-1)])
          p: horzcat([p0, p1, ..., p_(N-1)])
          s: horzcat([s0, s1, ..., s_(N-1)])
      and
          s0 <- f(a0, p0)
          s1 <- f(a1, p1)
          ...
          s_(N-1) <- f(a_(N-1), p_(N-1))
  



>  Function Function.map(str name, int n, [bool ] repeat_in, [bool ] repeat_out, Dict opts=Dict()) const 
------------------------------------------------------------------------

Generic map.

";

%feature("docstring") casadi::StabilizedQpSolver::description_in "

Get input scheme description by index.

";

%feature("docstring") casadi::StabilizedQpSolver::repr "

Print a representation of the object.

";

%feature("docstring") casadi::StabilizedQpSolver::size2_in "

Get input dimension.

";

%feature("docstring") casadi::StabilizedQpSolver::index_in "

Find the index for a string describing a particular entry of an input
scheme.

example: schemeEntry(\"x_opt\") -> returns NLP_SOLVER_X if FunctionInternal
adheres to SCHEME_NLPINput

";

%feature("docstring") casadi::StabilizedQpSolver::size_in "

Get input dimension.

";

%feature("docstring") casadi::StabilizedQpSolver::sx_out "

Get symbolic primitives equivalent to the output expressions.

";

%feature("docstring") casadi::StabilizedQpSolver::type_name "

Get type name.

";

%feature("docstring") casadi::StabilizedQpSolver::callReverse "

Create call to (cached) derivative function, reverse mode.

";

%feature("docstring") casadi::StabilizedQpSolver::getWorkSize "

Get the length of the work vector.

";

%feature("docstring") casadi::StabilizedQpSolver::checkInputs "[INTERNAL]
Check if the numerical values of the supplied bounds make sense.

";

%feature("docstring") casadi::StabilizedQpSolver::getStats "

Get all statistics obtained at the end of the last evaluate call.

";

%feature("docstring") casadi::StabilizedQpSolver::spEvaluate "[INTERNAL]
Propagate the sparsity pattern through a set of directional.

derivatives forward or backward (for usage, see the example
propagating_sparsity.cpp)

";

%feature("docstring") casadi::StabilizedQpSolver::getAlgorithmSize "

Get the number of atomic operations.

";

%feature("docstring") casadi::StabilizedQpSolver::print "

Print a description of the object.

";

%feature("docstring") casadi::StabilizedQpSolver::spCanEvaluate "[INTERNAL]
Is the class able to propagate seeds through the algorithm?

(for usage, see the example propagating_sparsity.cpp)

";

%feature("docstring") casadi::StabilizedQpSolver::mapaccum "

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
  



";

%feature("docstring") casadi::StabilizedQpSolver::setDerForward "

Set a function that calculates nfwd forward derivatives NOTE: Does not take
ownership, only weak references to the derivatives are kept internally.

";

%feature("docstring") casadi::StabilizedQpSolver::generateNativeCode "

Generate native code in the interfaced language for debugging

";

%feature("docstring") casadi::StabilizedQpSolver::getOptionAllowedIndex "[INTERNAL]  Get the index into allowed options of a certain option.

";

%feature("docstring") casadi::StabilizedQpSolver::name_in "

>  [str] Function.name_in() const 
------------------------------------------------------------------------

Get input scheme.

>  str Function.name_in(int ind) const 
------------------------------------------------------------------------

Get input scheme name by index.

";

%feature("docstring") casadi::StabilizedQpSolver::getOptionNames "

Get a list of all option names.

";

%feature("docstring") casadi::StabilizedQpSolver::getOptionTypeName "

Get the type name of a certain option.

";


// File: classcasadi_1_1StabilizedQpToQp.xml


// File: classcasadi_1_1StabilizedSqp.xml


// File: classcasadi_1_1Logger_1_1Stream.xml
%feature("docstring") casadi::Logger::Stream "C++ includes:
casadi_logger.hpp ";

%feature("docstring") casadi::Logger::Stream::Stream "";


// File: classcasadi_1_1Logger_1_1Streambuf.xml
%feature("docstring") casadi::Logger::Streambuf "C++ includes:
casadi_logger.hpp ";

%feature("docstring") casadi::Logger::Streambuf::Streambuf "";


// File: classcasadi_1_1SubAssign.xml


// File: classcasadi_1_1SubIndex.xml
%feature("docstring") casadi::SubIndex "

SubIndex class for Matrix Same as the above class but for single argument
return for operator() Joel Andersson

C++ includes: submatrix.hpp ";

%feature("docstring") casadi::SubIndex::SubIndex "

Constructor.

";


// File: classcasadi_1_1SubMatrix.xml
%feature("docstring") casadi::SubMatrix "

SubMatrix class for Matrix SubMatrix is the return type for operator() of
the Matrix class, it allows access to the value as well as changing the
parent object Joel Andersson

C++ includes: submatrix.hpp ";

%feature("docstring") casadi::SubMatrix::SubMatrix "

Constructor.

";


// File: classcasadi_1_1SubRef.xml


// File: classcasadi_1_1SymbolicMX.xml


// File: classcasadi_1_1SymbolicQr.xml


// File: classcasadi_1_1SymbolicSX.xml


// File: classcasadi_1_1Transpose.xml


// File: classcasadi_1_1UnaryMX.xml


// File: classcasadi_1_1UnarySX.xml


// File: classcasadi_1_1Vertcat.xml


// File: classcasadi_1_1Vertsplit.xml


// File: classcasadi_1_1WeakRef.xml
%feature("docstring") casadi::WeakRef::isNull "[INTERNAL]  Is a null
pointer?

";

%feature("docstring") casadi::WeakRef "[INTERNAL]  Weak reference type A
weak reference to a SharedObject.

Joel Andersson

C++ includes: weak_ref.hpp ";

%feature("docstring") casadi::WeakRef::shared "[INTERNAL]  Get a shared
(owning) reference.

";

%feature("docstring") casadi::WeakRef::__hash__ "[INTERNAL]  Returns a
number that is unique for a given Node. If the Object does not point to any
node, \"0\" is returned.

";

%feature("docstring") casadi::WeakRef::print "[INTERNAL]  Print a
description of the object.

";

%feature("docstring") casadi::WeakRef::getDescription "[INTERNAL]  Return a
string with a description (for SWIG)

";

%feature("docstring") casadi::WeakRef::WeakRef "

>  WeakRef(int dummy=0)
------------------------------------------------------------------------
[INTERNAL] 
Default constructor.

>  WeakRef(SharedObject shared)
------------------------------------------------------------------------
[INTERNAL] 
Construct from a shared object (also implicit type conversion)

";

%feature("docstring") casadi::WeakRef::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::WeakRef::getRepresentation "[INTERNAL]
Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::WeakRef::repr "[INTERNAL]  Print a
representation of the object.

";

%feature("docstring") casadi::WeakRef::alive "[INTERNAL]  Check if alive.

";


// File: classcasadi_1_1Wrapper.xml
%feature("docstring") casadi::Wrapper::checkDimensions "[INTERNAL]  Check
the dimensions of the internal function after initialization.

";

%feature("docstring") casadi::Wrapper::evaluate "[INTERNAL]  Evaluate the
internal function and make it external.

";

%feature("docstring") casadi::Wrapper "[INTERNAL]  A helper class for a
Function that wrap another Function.

Joris Gillis

C++ includes: wrapper.hpp ";


// File: classcasadi_1_1XmlFile.xml
%feature("docstring") casadi::XmlFile "

XML parser Can be used for parsing XML files into CasADi data structures.

Joel Andersson

C++ includes: xml_file.hpp ";

%feature("docstring") casadi::XmlFile::getDescription "

Return a string with a description (for SWIG)

";

%feature("docstring") casadi::XmlFile::print "

Print a description of the object.

";

%feature("docstring") casadi::XmlFile::printPtr "[INTERNAL]  Print the
pointer to the internal class

";

%feature("docstring") casadi::XmlFile::getRepresentation "

Return a string with a representation (for SWIG)

";

%feature("docstring") casadi::XmlFile::__hash__ "

Returns a number that is unique for a given Node. If the Object does not
point to any node, \"0\" is returned.

";

%feature("docstring") casadi::XmlFile::XmlFile "";

%feature("docstring") casadi::XmlFile::~XmlFile "";

%feature("docstring") casadi::XmlFile::isNull "

Is a null pointer?

";

%feature("docstring") casadi::XmlFile::repr "

Print a representation of the object.

";


// File: classcasadi_1_1ZeroByZero.xml


// File: classcasadi_1_1ZeroSX.xml


// File: namespacecasadi.xml
%feature("docstring") casadi::matrixName< double > "
Get typename.

";

%feature("docstring") casadi::complement "

Returns the list of all i in [0, size[ not found in supplied list.

The supplied vector may contain duplicates and may be non-monotonous The
supplied vector will be checked for bounds The result vector is guaranteed
to be monotonously increasing

";

%feature("docstring") casadi::inBounds "

>  bool inBounds([T ] v, int upper)
------------------------------------------------------------------------

Check if for each element of v holds: v_i < upper.

>  bool inBounds([T ] v, int lower, int upper)
------------------------------------------------------------------------

Check if for each element of v holds: lower <= v_i < upper.

";

%feature("docstring") casadi::getSchemeSize "";

%feature("docstring") casadi::swapIndices "

swap inner and outer indices of list of lists



::

  * [[apple0,apple1,...],[pear0,pear1,...]] ->
  *   [[apple0,pear0],[apple1,pear1],...]
  * 



";

%feature("docstring") casadi::isNon_increasing "

Check if the vector is non-increasing.

";

%feature("docstring") casadi::dlaqge_ "[INTERNAL]  Equilibrate the system.

";

%feature("docstring") casadi::iszero "[INTERNAL]  Check if entry is zero
(false negative allowed)

";

%feature("docstring") casadi::profileWriteEntry "[INTERNAL] ";

%feature("docstring") casadi::qpOut "

Output arguments of an QP Solver

>Output scheme: casadi::QpSolverOutput (QP_SOLVER_NUM_OUT = 4) [qpOut]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| QP_SOLVER_X            | x                      | The primal solution .  |
+------------------------+------------------------+------------------------+
| QP_SOLVER_COST         | cost                   | The optimal cost .     |
+------------------------+------------------------+------------------------+
| QP_SOLVER_LAM_A        | lam_a                  | The dual solution      |
|                        |                        | corresponding to       |
|                        |                        | linear bounds .        |
+------------------------+------------------------+------------------------+
| QP_SOLVER_LAM_X        | lam_x                  | The dual solution      |
|                        |                        | corresponding to       |
|                        |                        | simple bounds .        |
+------------------------+------------------------+------------------------+

";

%feature("docstring") casadi::ProfilingType< ProfilingData_SOURCE > "
[INTERNAL] ";

%feature("docstring") casadi::daeIn "

Input arguments of an ODE/DAE function

>Input scheme: casadi::DAEInput (DAE_NUM_IN = 4) [daeIn]

+-----------+-------+----------------------------+
| Full name | Short |        Description         |
+===========+=======+============================+
| DAE_X     | x     | Differential state .       |
+-----------+-------+----------------------------+
| DAE_Z     | z     | Algebraic state .          |
+-----------+-------+----------------------------+
| DAE_P     | p     | Parameter .                |
+-----------+-------+----------------------------+
| DAE_T     | t     | Explicit time dependence . |
+-----------+-------+----------------------------+

";

%feature("docstring") casadi::integratorIn "

Input arguments of an integrator

>Input scheme: casadi::IntegratorInput (INTEGRATOR_NUM_IN = 6) [integratorIn]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| INTEGRATOR_X0          | x0                     | Differential state at  |
|                        |                        | the initial time .     |
+------------------------+------------------------+------------------------+
| INTEGRATOR_P           | p                      | Parameters .           |
+------------------------+------------------------+------------------------+
| INTEGRATOR_Z0          | z0                     | Initial guess for the  |
|                        |                        | algebraic variable .   |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RX0         | rx0                    | Backward differential  |
|                        |                        | state at the final     |
|                        |                        | time .                 |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RP          | rp                     | Backward parameter     |
|                        |                        | vector .               |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RZ0         | rz0                    | Initial guess for the  |
|                        |                        | backwards algebraic    |
|                        |                        | variable .             |
+------------------------+------------------------+------------------------+

";

%feature("docstring") casadi::profileWriteSourceLine "[INTERNAL] ";

%feature("docstring") casadi::isDecreasing "

Check if the vector is strictly decreasing.

";

%feature("docstring") casadi::casadi_to_sparse "[INTERNAL]  Convert dense
to sparse.

";

%feature("docstring") casadi::profileWriteName "[INTERNAL] ";

%feature("docstring") casadi::ProfilingType< ProfilingData_EXIT > "
[INTERNAL] ";

%feature("docstring") casadi::casadi_fill_n "[INTERNAL]  FILL: x <- alpha.

";

%feature("docstring") casadi::casadi_scal "[INTERNAL]  SCAL: x <- alpha*x.

";

%feature("docstring") casadi::getSchemeEntryDoc "";

%feature("docstring") casadi::ProfilingType< ProfilingData_ENTRY > "
[INTERNAL] ";

%feature("docstring") casadi::timerPlusEq "[INTERNAL] ";

%feature("docstring") casadi::hash_combine "

>  void hash_combine(std.size_t &seed, T v)

>  void hash_combine(std.size_t &seed, [int ] v)
------------------------------------------------------------------------
[INTERNAL] 
Generate a hash value incrementally (function taken from boost)

>  void hash_combine(std.size_t &seed, const int *v, int sz)
------------------------------------------------------------------------
[INTERNAL] 
Generate a hash value incrementally, array.

";

%feature("docstring") casadi::casadi_swap "[INTERNAL]  SWAP: x <-> y.

";

%feature("docstring") casadi::getSchemeName "";

%feature("docstring") casadi::profileWriteTime "[INTERNAL] ";

%feature("docstring") casadi::hasNegative "

Check if the vector has negative entries.

";

%feature("docstring") casadi::rdaeIn "

Input arguments of an ODE/DAE backward integration function

>Input scheme: casadi::RDAEInput (RDAE_NUM_IN = 7) [rdaeIn]

+-----------+-------+-------------------------------+
| Full name | Short |          Description          |
+===========+=======+===============================+
| RDAE_RX   | rx    | Backward differential state . |
+-----------+-------+-------------------------------+
| RDAE_RZ   | rz    | Backward algebraic state .    |
+-----------+-------+-------------------------------+
| RDAE_RP   | rp    | Backward parameter vector .   |
+-----------+-------+-------------------------------+
| RDAE_X    | x     | Forward differential state .  |
+-----------+-------+-------------------------------+
| RDAE_Z    | z     | Forward algebraic state .     |
+-----------+-------+-------------------------------+
| RDAE_P    | p     | Parameter vector .            |
+-----------+-------+-------------------------------+
| RDAE_T    | t     | Explicit time dependence .    |
+-----------+-------+-------------------------------+

";

%feature("docstring") casadi::check_exposed "[INTERNAL] ";

%feature("docstring") casadi::getTimerTime "[INTERNAL] ";

%feature("docstring") casadi::read_matlab "

>  void read_matlab(std.istream &stream,[T ] v)
------------------------------------------------------------------------

Read vector, matlab style.

>  void read_matlab(std.ifstream &file,[[T ] ] v)
------------------------------------------------------------------------

Read matrix, matlab style.

";

%feature("docstring") casadi::write_matlab "

>  void write_matlab(std.ostream &stream, [T ] v)
------------------------------------------------------------------------

Print vector, matlab style.

>  void write_matlab(std.ostream &stream, [[T ] ] v)
------------------------------------------------------------------------

Print matrix, matlab style.

";

%feature("docstring") casadi::hash_sparsity "

>  std.size_t hash_sparsity(int nrow, int ncol, [int ] colind, [int ] row)
------------------------------------------------------------------------
[INTERNAL] 
Hash a sparsity pattern.

>  std.size_t hash_sparsity(int nrow, int ncol, const int *colind, const int *row)
------------------------------------------------------------------------
[INTERNAL] 
";

%feature("docstring") casadi::isStrictlyMonotone "

Check if the vector is strictly monotone.

";

%feature("docstring") casadi::casadi_to_dense_tr "[INTERNAL]  Convert
sparse to transposed dense.

";

%feature("docstring") casadi::IOScheme "";

%feature("docstring") casadi::isRegular "

Checks if vector does not contain NaN or Inf.

";

%feature("docstring") casadi::jacGIn "

Input arguments of an NLP Jacobian function

>Input scheme: casadi::JacGInput (JACG_NUM_IN = 2) [jacGIn]

+-----------+-------+---------------------+
| Full name | Short |     Description     |
+===========+=======+=====================+
| JACG_X    | x     | Decision variable . |
+-----------+-------+---------------------+
| JACG_P    | p     | Fixed parameter .   |
+-----------+-------+---------------------+

";

%feature("docstring") casadi::rdaeOut "

Output arguments of an ODE/DAE backward integration function

>Output scheme: casadi::RDAEOutput (RDAE_NUM_OUT = 3) [rdaeOut]

+-----------+-------+-------------------------------------------+
| Full name | Short |                Description                |
+===========+=======+===========================================+
| RDAE_ODE  | ode   | Right hand side of ODE. .                 |
+-----------+-------+-------------------------------------------+
| RDAE_ALG  | alg   | Right hand side of algebraic equations. . |
+-----------+-------+-------------------------------------------+
| RDAE_QUAD | quad  | Right hand side of quadratures. .         |
+-----------+-------+-------------------------------------------+

";

%feature("docstring") casadi::matrixName< int > "

Get typename.

";

%feature("docstring") casadi::dgeequ_ "[INTERNAL]  Calculate col and row
scaling.

";

%feature("docstring") casadi::lookupvector "

Returns a vector for quickly looking up entries of supplied list.

lookupvector[i]!=-1 <=> v contains i v[lookupvector[i]] == i <=> v contains
i

Duplicates are treated by looking up last occurrence

";

%feature("docstring") casadi::ProfilingType< ProfilingData_IO > " [INTERNAL]
";

%feature("docstring") casadi::dormqr_ "[INTERNAL]  Multiply right hand side
with Q-transpose (lapack)

";

%feature("docstring") casadi::getRealTime "[INTERNAL]  Returns the real
time, in seconds, or -1.0 if an error occurred.

Time is measured since an arbitrary and OS-dependent start time. The
returned real time is only useful for computing an elapsed time between two
calls to this function.

";

%feature("docstring") casadi::integratorOut "

Output arguments of an integrator

>Output scheme: casadi::IntegratorOutput (INTEGRATOR_NUM_OUT = 6) [integratorOut]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| INTEGRATOR_XF          | xf                     | Differential state at  |
|                        |                        | the final time .       |
+------------------------+------------------------+------------------------+
| INTEGRATOR_QF          | qf                     | Quadrature state at    |
|                        |                        | the final time .       |
+------------------------+------------------------+------------------------+
| INTEGRATOR_ZF          | zf                     | Algebraic variable at  |
|                        |                        | the final time .       |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RXF         | rxf                    | Backward differential  |
|                        |                        | state at the initial   |
|                        |                        | time .                 |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RQF         | rqf                    | Backward quadrature    |
|                        |                        | state at the initial   |
|                        |                        | time .                 |
+------------------------+------------------------+------------------------+
| INTEGRATOR_RZF         | rzf                    | Backward algebraic     |
|                        |                        | variable at the        |
|                        |                        | initial time .         |
+------------------------+------------------------+------------------------+

";

%feature("docstring") casadi::ProfilingType< ProfilingData_NAME > "
[INTERNAL] ";

%feature("docstring") casadi::operation_checker "[INTERNAL] ";

%feature("docstring") casadi::controldaeIn "

Input arguments of an ODE/DAE function

";

%feature("docstring") casadi::hessLagOut "

Output arguments of an NLP Hessian function

>Output scheme: casadi::HessLagOutput (HESSLAG_NUM_OUT = 5) [hessLagOut]

+----------------+--------+------------------------------------------------+
|   Full name    | Short  |                  Description                   |
+================+========+================================================+
| HESSLAG_HESS   | hess   | Hessian of the Lagrangian .                    |
+----------------+--------+------------------------------------------------+
| HESSLAG_F      | f      | Objective function .                           |
+----------------+--------+------------------------------------------------+
| HESSLAG_G      | g      | Constraint function .                          |
+----------------+--------+------------------------------------------------+
| HESSLAG_GRAD_X | grad_x | Gradient of the Lagrangian with respect to x . |
+----------------+--------+------------------------------------------------+
| HESSLAG_GRAD_P | grad_p | Gradient of the Lagrangian with respect to p . |
+----------------+--------+------------------------------------------------+

";

%feature("docstring") casadi::casadi_mm_sparse_t "[INTERNAL]  Sparse
matrix-matrix multiplication, first factor transposed: z <- z + trans(x)*y.

";

%feature("docstring") casadi::isNonDecreasing "

Check if the vector is non-decreasing.

";

%feature("docstring") casadi::profileWriteBare "[INTERNAL] ";

%feature("docstring") casadi::casadi_mv_t "[INTERNAL]  Sparse matrix-vector
multiplication, first factor transposed: z <- z + trans(x)*y.

";

%feature("docstring") casadi::nlpIn "

Input arguments of an NLP function

>Input scheme: casadi::NLPInput (NL_NUM_IN = 2) [nlpIn]

+-----------+-------+---------------------+
| Full name | Short |     Description     |
+===========+=======+=====================+
| NL_X      | x     | Decision variable . |
+-----------+-------+---------------------+
| NL_P      | p     | Fixed parameter .   |
+-----------+-------+---------------------+

";

%feature("docstring") casadi::isIncreasing "

Check if the vector is strictly increasing.

";

%feature("docstring") casadi::gradFIn "

Input arguments of an NLP objective gradient function

>Input scheme: casadi::GradFInput (GRADF_NUM_IN = 2) [gradFIn]

+-----------+-------+---------------------+
| Full name | Short |     Description     |
+===========+=======+=====================+
| GRADF_X   | x     | Decision variable . |
+-----------+-------+---------------------+
| GRADF_P   | p     | Fixed parameter .   |
+-----------+-------+---------------------+

";

%feature("docstring") casadi::simpleIRK "

Construct an implicit Runge-Kutta integrator using a collocation scheme The
constructed function has three inputs, corresponding to initial state (x0),
parameter (p) and integration time (h) and one output, corresponding to
final state (xf).

Parameters:
-----------

f:  ODE function with two inputs (x and p) and one output (xdot)

N:  Number of integrator steps

order:  Order of interpolating polynomials

scheme:  Collocation scheme, as excepted by collocationPoints function.

";

%feature("docstring") casadi::ptrVec "[INTERNAL]  Convenience function,
convert vectors to vectors of pointers.

";

%feature("docstring") casadi::casadi_axpy "[INTERNAL]  AXPY: y <- a*x + y.

";

%feature("docstring") casadi::isMonotone "

Check if the vector is monotone.

";

%feature("docstring") casadi::jacGOut "

Output arguments of an NLP Jacobian function

>Output scheme: casadi::JacGOutput (JACG_NUM_OUT = 3) [jacGOut]

+-----------+-------+-------------------------------+
| Full name | Short |          Description          |
+===========+=======+===============================+
| JACG_JAC  | jac   | Jacobian of the constraints . |
+-----------+-------+-------------------------------+
| JACG_F    | f     | Objective function .          |
+-----------+-------+-------------------------------+
| JACG_G    | g     | Constraint function .         |
+-----------+-------+-------------------------------+

";

%feature("docstring") casadi::casadi_iamax "[INTERNAL]  IAMAX: index
corresponding to the entry with the largest absolute value.

";

%feature("docstring") casadi::getSchemeEntryNames "";

%feature("docstring") casadi::hessLagIn "

Input arguments of an NLP Hessian function

>Input scheme: casadi::HessLagInput (HESSLAG_NUM_IN = 4) [hessLagIn]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| HESSLAG_X              | x                      | Decision variable .    |
+------------------------+------------------------+------------------------+
| HESSLAG_P              | p                      | Fixed parameter .      |
+------------------------+------------------------+------------------------+
| HESSLAG_LAM_F          | lam_f                  | Multiplier for f. Just |
|                        |                        | a scalar factor for    |
|                        |                        | the objective that the |
|                        |                        | NLP solver might use   |
|                        |                        | to scale the           |
|                        |                        | objective.             |
+------------------------+------------------------+------------------------+
| HESSLAG_LAM_G          | lam_g                  | Multiplier for g .     |
+------------------------+------------------------+------------------------+

";

%feature("docstring") casadi::daeOut "

Output arguments of an DAE function

>Output scheme: casadi::DAEOutput (DAE_NUM_OUT = 3) [daeOut]

+-----------+-------+--------------------------------------------+
| Full name | Short |                Description                 |
+===========+=======+============================================+
| DAE_ODE   | ode   | Right hand side of the implicit ODE .      |
+-----------+-------+--------------------------------------------+
| DAE_ALG   | alg   | Right hand side of algebraic equations .   |
+-----------+-------+--------------------------------------------+
| DAE_QUAD  | quad  | Right hand side of quadratures equations . |
+-----------+-------+--------------------------------------------+

";

%feature("docstring") casadi::getSchemeEntryEnum "";

%feature("docstring") casadi::collocationInterpolators "

Obtain collocation interpolating matrices.

Parameters:
-----------

tau_root:  location of collocation points, as obtained from
collocationPoints

C:  interpolating coefficients to obtain derivatives Length: order+1, order
+ 1



::

dX/dt @collPoint(j) ~ Sum_i C[j][i]*X@collPoint(i)



Parameters:
-----------

D:  interpolating coefficients to obtain end state Length: order+1

";

%feature("docstring") casadi::matrixName "

Get typename.

";

%feature("docstring") casadi::casadi_to_sparse_tr "[INTERNAL]  Convert
transposed dense to sparse.

";

%feature("docstring") casadi::casadi_nrm2 "[INTERNAL]  NRM2: ||x||_2 ->
return.

";

%feature("docstring") casadi::nlpSolverIn "

Input arguments of an NLP Solver

>Input scheme: casadi::NlpSolverInput (NLP_SOLVER_NUM_IN = 8) [nlpSolverIn]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| NLP_SOLVER_X0          | x0                     | Decision variables,    |
|                        |                        | initial guess (nx x 1) |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_P           | p                      | Value of fixed         |
|                        |                        | parameters (np x 1) .  |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LBX         | lbx                    | Decision variables     |
|                        |                        | lower bound (nx x 1),  |
|                        |                        | default -inf .         |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_UBX         | ubx                    | Decision variables     |
|                        |                        | upper bound (nx x 1),  |
|                        |                        | default +inf .         |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LBG         | lbg                    | Constraints lower      |
|                        |                        | bound (ng x 1),        |
|                        |                        | default -inf .         |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_UBG         | ubg                    | Constraints upper      |
|                        |                        | bound (ng x 1),        |
|                        |                        | default +inf .         |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_X0      | lam_x0                 | Lagrange multipliers   |
|                        |                        | for bounds on X,       |
|                        |                        | initial guess (nx x 1) |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_G0      | lam_g0                 | Lagrange multipliers   |
|                        |                        | for bounds on G,       |
|                        |                        | initial guess (ng x 1) |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+

";

%feature("docstring") casadi::hash_value "[INTERNAL]  Hash value of an
integer.

";

%feature("docstring") casadi::casadi_copy_n "[INTERNAL]  COPY: y <-x.

";

%feature("docstring") casadi::collocationPoints "

Obtain collocation points of specific order and scheme.

Parameters:
-----------

scheme:  'radau' or 'legendre'

";

%feature("docstring") casadi::getSchemeEntryEnumName "";

%feature("docstring") casadi::simpleIntegrator "

Simplified wrapper for the Integrator class Constructs an integrator using
the same syntax as simpleRK and simpleIRK. The constructed function has
three inputs, corresponding to initial state (x0), parameter (p) and
integration time (h) and one output, corresponding to final state (xf).

Parameters:
-----------

f:  ODE function with two inputs (x and p) and one output (xdot)

N:  Number of integrator steps

order:  Order of interpolating polynomials

scheme:  Collocation scheme, as excepted by collocationPoints function.

";

%feature("docstring") casadi::casadi_to_dense "[INTERNAL]  Convert sparse
to dense.

";

%feature("docstring") casadi::controlsimulatorIn "

Input arguments of a control simulator

";

%feature("docstring") casadi::dgetrf_ "[INTERNAL]  LU-Factorize dense
matrix (lapack)

";

%feature("docstring") casadi::casadi_project "[INTERNAL]  Sparse copy: y <-
x, w work vector (length >= number of rows)

";

%feature("docstring") casadi::profileWriteExit "[INTERNAL] ";

%feature("docstring") casadi::casadi_inner_prod "[INTERNAL]  Inner product.

";

%feature("docstring") casadi::diffTimers "[INTERNAL] ";

%feature("docstring") casadi::nlpSolverOut "

Output arguments of an NLP Solver

>Output scheme: casadi::NlpSolverOutput (NLP_SOLVER_NUM_OUT = 6) [nlpSolverOut]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| NLP_SOLVER_X           | x                      | Decision variables at  |
|                        |                        | the optimal solution   |
|                        |                        | (nx x 1) .             |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_F           | f                      | Cost function value at |
|                        |                        | the optimal solution   |
|                        |                        | (1 x 1) .              |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_G           | g                      | Constraints function   |
|                        |                        | at the optimal         |
|                        |                        | solution (ng x 1) .    |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_X       | lam_x                  | Lagrange multipliers   |
|                        |                        | for bounds on X at the |
|                        |                        | solution (nx x 1) .    |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_G       | lam_g                  | Lagrange multipliers   |
|                        |                        | for bounds on G at the |
|                        |                        | solution (ng x 1) .    |
+------------------------+------------------------+------------------------+
| NLP_SOLVER_LAM_P       | lam_p                  | Lagrange multipliers   |
|                        |                        | for bounds on P at the |
|                        |                        | solution (np x 1) .    |
+------------------------+------------------------+------------------------+

";

%feature("docstring") casadi::profileWrite "[INTERNAL] ";

%feature("docstring") casadi::dtrsm_ "[INTERNAL]   Solve upper triangular
system (lapack)

";

%feature("docstring") casadi::casadi_quad_form "[INTERNAL]  Calculates
inner_prod(x, mul(A, x))

Calculates inner_prod(x, mul(A, x)) without memory allocation.

";

%feature("docstring") casadi::linsolIn "

Input arguments of a linear solver

>Input scheme: casadi::LinsolInput (LINSOL_NUM_IN = 2) [linsolIn]

+-----------+-------+------------------------------------------------+
| Full name | Short |                  Description                   |
+===========+=======+================================================+
| LINSOL_A  | A     | The square matrix A: sparse, (n x n). .        |
+-----------+-------+------------------------------------------------+
| LINSOL_B  | B     | The right-hand-side matrix b: dense, (n x m) . |
+-----------+-------+------------------------------------------------+

";

%feature("docstring") casadi::casadi_mm_sparse "[INTERNAL]  Sparse matrix-
matrix multiplication: z <- z + x*y.

";

%feature("docstring") casadi::linsolOut "

Output arguments of a linear solver

>Output scheme: casadi::LinsolOutput (LINSOL_NUM_OUT = 1) [linsolOut]

+-----------+-------+----------------------------------------------+
| Full name | Short |                 Description                  |
+===========+=======+==============================================+
| LINSOL_X  | X     | Solution to the linear system of equations . |
+-----------+-------+----------------------------------------------+

";

%feature("docstring") casadi::diffToDict "[INTERNAL] ";

%feature("docstring") casadi::casadi_mv "[INTERNAL]  Sparse matrix-vector
multiplication: z <- z + x*y.

";

%feature("docstring") casadi::dgetrs_ "[INTERNAL]   Solve a system of
equation using an LU-factorized matrix (lapack)

";

%feature("docstring") casadi::ProfilingType "[INTERNAL] ";

%feature("docstring") casadi::qpIn "

Input arguments of a QP problem

>Input scheme: casadi::QpSolverInput (QP_SOLVER_NUM_IN = 9) [qpIn]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| QP_SOLVER_H            | h                      | The square matrix H:   |
|                        |                        | sparse, (n x n). Only  |
|                        |                        | the lower triangular   |
|                        |                        | part is actually used. |
|                        |                        | The matrix is assumed  |
|                        |                        | to be symmetrical.     |
+------------------------+------------------------+------------------------+
| QP_SOLVER_G            | g                      | The vector g: dense,   |
|                        |                        | (n x 1) .              |
+------------------------+------------------------+------------------------+
| QP_SOLVER_A            | a                      | The matrix A: sparse,  |
|                        |                        | (nc x n) - product     |
|                        |                        | with x must be dense.  |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+
| QP_SOLVER_LBA          | lba                    | dense, (nc x 1)        |
+------------------------+------------------------+------------------------+
| QP_SOLVER_UBA          | uba                    | dense, (nc x 1)        |
+------------------------+------------------------+------------------------+
| QP_SOLVER_LBX          | lbx                    | dense, (n x 1)         |
+------------------------+------------------------+------------------------+
| QP_SOLVER_UBX          | ubx                    | dense, (n x 1)         |
+------------------------+------------------------+------------------------+
| QP_SOLVER_X0           | x0                     | dense, (n x 1)         |
+------------------------+------------------------+------------------------+
| QP_SOLVER_LAM_X0       | lam_x0                 | dense                  |
+------------------------+------------------------+------------------------+

";

%feature("docstring") casadi::ProfilingType< ProfilingData_TIMELINE > "
[INTERNAL] ";

%feature("docstring") casadi::casadi_trans "[INTERNAL]  TRANS: y <-
trans(x)

";

%feature("docstring") casadi::stabilizedQpIn "

Input arguments of a QP problem

>Input scheme: casadi::StabilizedQpSolverInput (STABILIZED_QP_SOLVER_NUM_IN = 12) [stabilizedQpIn]

+------------------------+------------------------+------------------------+
|       Full name        |         Short          |      Description       |
+========================+========================+========================+
| STABILIZED_QP_SOLVER_H | h                      | The square matrix H:   |
|                        |                        | sparse, (n x n). Only  |
|                        |                        | the lower triangular   |
|                        |                        | part is actually used. |
|                        |                        | The matrix is assumed  |
|                        |                        | to be symmetrical.     |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_G | g                      | The vector g: dense,   |
|                        |                        | (n x 1) .              |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_A | a                      | The matrix A: sparse,  |
|                        |                        | (nc x n) - product     |
|                        |                        | with x must be dense.  |
|                        |                        | .                      |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_L | lba                    | dense, (nc x 1)        |
| BA                     |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_U | uba                    | dense, (nc x 1)        |
| BA                     |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_L | lbx                    | dense, (n x 1)         |
| BX                     |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_U | ubx                    | dense, (n x 1)         |
| BX                     |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_X | x0                     | dense, (n x 1)         |
| 0                      |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_L | lam_x0                 | dense                  |
| AM_X0                  |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_M | muR                    | dense (1 x 1)          |
| UR                     |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_M | muE                    | dense (nc x 1)         |
| UE                     |                        |                        |
+------------------------+------------------------+------------------------+
| STABILIZED_QP_SOLVER_M | mu                     | dense (nc x 1)         |
| U                      |                        |                        |
+------------------------+------------------------+------------------------+

";

%feature("docstring") casadi::gradFOut "

Output arguments of an NLP objective gradient function

>Output scheme: casadi::GradFOutput (GRADF_NUM_OUT = 3) [gradFOut]

+------------+-------+-------------------------------+
| Full name  | Short |          Description          |
+============+=======+===============================+
| GRADF_GRAD | grad  | Jacobian of the constraints . |
+------------+-------+-------------------------------+
| GRADF_F    | f     | Objective function .          |
+------------+-------+-------------------------------+
| GRADF_G    | g     | Constraint function .         |
+------------+-------+-------------------------------+

";

%feature("docstring") casadi::casadi_asum "[INTERNAL]  ASUM: ||x||_1 ->
return.

";

%feature("docstring") casadi::nlpOut "

Output arguments of an NLP function

>Output scheme: casadi::NLPOutput (NL_NUM_OUT = 2) [nlpOut]

+-----------+-------+-----------------------+
| Full name | Short |      Description      |
+===========+=======+=======================+
| NL_F      | f     | Objective function .  |
+-----------+-------+-----------------------+
| NL_G      | g     | Constraint function . |
+-----------+-------+-----------------------+

";

%feature("docstring") casadi::simpleRK "

Construct an explicit Runge-Kutta integrator The constructed function has
three inputs, corresponding to initial state (x0), parameter (p) and
integration time (h) and one output, corresponding to final state (xf).

Parameters:
-----------

f:  ODE function with two inputs (x and p) and one output (xdot)

N:  Number of integrator steps

order:  Order of interpolating polynomials

";

%feature("docstring") casadi::dgeqrf_ "[INTERNAL]  QR-factorize dense
matrix (lapack)

";

%feature("docstring") casadi::ptrToLong "[INTERNAL] ";

%feature("docstring") casadi::profileWriteSourceLineDep "[INTERNAL] ";

%feature("docstring") casadi::userOut "";

%feature("docstring") casadi::casadi_norm_inf_mul "[INTERNAL]  Inf-norm of
a Matrix-matrix product,*

Parameters:
-----------

dwork:  A real work vector that you must allocate Minimum size: y.size1()

iwork:  A integer work vector that you must allocate Minimum size:
y.size1()+x.size2()+1

";

%feature("docstring") casadi::matrixName< SXElement > " [INTERNAL] ";

%feature("docstring") casadi::getSchemeEntryName "";


// File: namespaceIpopt.xml


// File: namespacestd.xml


// File: chapter1.xml


// File: chapter2.xml


// File: chapter3.xml


// File: chapter4.xml


// File: chapter5.xml


// File: chapter6.xml

