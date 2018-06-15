.. _sec-syntax_differences:

Difference in usage from different languages
============================================


.. include:: defs.rst

General usage
-------------


Example:

.. list-table:: General usage
   :header-rows: 1

   * - Task
     - Python
     - C++
     - MATLAB/Octave
   * - Starting |casadi|
     - .. code-block:: python

            from casadi import *
     - .. code-block:: cpp

            #include <casadi/casadi.hpp>
            using namespace casadi;
     - .. code-block:: octave

            import casadi.*
   * - Printing object
     - .. code-block:: python

            print(A)
     - .. code-block:: cpp

            std::cout << A;
     - .. code-block:: octave

            disp(A)

   * - Printing with type information
     - ``A <ENTER>`` (interactive), ``print(repr(A))``
     - .. code-block:: cpp

            std::cout << repr(A);
     - ``A <ENTER>`` (interactive), ``disp(repr(A))``

   * - Get (extended) representation, ``more=false`` by default
     - ``A.str(more)``
     - ``str(A, more)``
     - ``str(A, more)``

   * - Calling a class function
     - .. code-block:: python

            SX.zeros(3,4)
     - .. code-block:: cpp

            SX::zeros(3,4)
     - .. code-block:: octave

            SX.zeros(3,4)

   * - Creating a dictionary (e.g. for options) 
     - .. code-block:: python

            d = dict(opt1=opt1)
            #or
            d = {'opt1':opt1}
            #or
            d = {}; a['opt1'] = opt1
     - .. code-block:: cpp

            Dict a;
            a["opt1"] = opt1;
     - .. code-block:: octave

            a = struct;
            a.opt1 = opt1;
            
            #or
            a = struct('opt1',opt1);

   * - Creating a symbol
     - .. code-block:: python

            MX.sym("x",2,2)
     - .. code-block:: cpp

            MX::sym("x",2,2)
     - .. code-block:: octave

            MX.sym('x',2,2)

   * - Creating a function
     - .. code-block:: python

            Function("f",[x,y],[x+y])
     - .. code-block:: cpp

            Function("f",{x,y},{x+y})
     - .. code-block:: octave

            Function('f',{x,y},{x+y})


   * - Calling a function (form 1)
     - .. code-block:: python

            z=f(x,y)
     - .. code-block:: cpp

            z = f({x,y})
     - .. code-block:: octave

            z=f(x,y)

   * - Calling a function (form 2)
     - .. code-block:: python

            res = f(x=x,y=y)
            res["z"]
     - .. code-block:: cpp

            auto res =\
            f({{"x",x},{"y",y}});
            res["z"]
     - .. code-block:: octave

            res=f('x',x,'y',y)
            res.z

   * - Create an NLP solver
     - .. code-block:: python

            nlp = {"x":x,"f":f}
            nlpsol("S","ipopt",nlp)
     - .. code-block:: cpp

            MXDict nlp =\
            {{"x",x},{"f",f}};
            nlpsol("S","ipopt",nlp);
     - .. code-block:: octave

            nlp=struct('x',x,'f',f);
            nlpsol('S','ipopt',nlp);
 


List of operations
------------------

The following is a list of the most important operations. Operations that differ between the different
languages are marked with a star (*). This list is neither complete, nor does it show all the variants of
each operation. Further information is available in the API documentation.



.. list-table:: General usage
   :header-rows: 1

   * - Task
     - Python
     - C++
     - MATLAB/Octave
   * - Addition, subtraction
     - .. code-block:: python

            x+y, x-y, -x
     - .. code-block:: cpp

            x+y, x-y, -x
     - .. code-block:: octave

            x+y, x-y, -x
   * - Elementwise multiplication, division
     - .. code-block:: python

            x*y, x/y
     - .. code-block:: cpp

            x*y, x/y
     - .. code-block:: octave

            x.*y, x./y
   * -  Matrix multiplication
     - .. code-block:: python3

            x @ y
            # or before Py3.4:
            mtimes(x,y)
     - .. code-block:: cpp

            mtimes(x,y)
     - .. code-block:: octave

            x*y
   * - Linear system solve
     - .. code-block:: python3

            solve(A,b)
     - .. code-block:: cpp

            solve(A,b)
     - .. code-block:: octave

            A\b
   * - Natural exponential function and logarithm
     - .. code-block:: python

            exp(x), log(x)
     - .. code-block:: cpp

            exp(x), log(x)
     - .. code-block:: octave

            exp(x), log(x)
   * - Exponentiation
     - .. code-block:: python

            x**y
     - .. code-block:: cpp

            pow(x,y)
     - .. code-block:: octave

            x^y, x.^y

   * - Square root
     - .. code-block:: python

            sqrt(x)
     - .. code-block:: cpp

            sqrt(x)
     - .. code-block:: octave

            sqrt(x)
   * -  Trigonometric functions
     - .. code-block:: python

            sin(x), cos(x), tan(x)
     - .. code-block:: cpp

            sin(x), cos(x), tan(x)
     - .. code-block:: octave

            sin(x), cos(x), tan(x)
   * - Inverse trigonometric
     - .. code-block:: python

            asin(x), acos(x)
     - .. code-block:: cpp

            asin(x), acos(x)
     - .. code-block:: octave

            asin(x), acos(x)
   * - Two argument arctangent
     - .. code-block:: python

            atan2(x, y)
     - .. code-block:: cpp

            atan2(x, y)
     - .. code-block:: octave

            atan2(x, y)
   * - Hyperbolic functions
     - .. code-block:: python

            sinh(x), cosh(x), tanh(x)
     - .. code-block:: cpp

            sinh(x), cosh(x), tanh(x)
     - .. code-block:: octave

            sinh(x), cosh(x), tanh(x)
   * - Inverse hyperbolic
     - .. code-block:: python

            asinh(x), acosh(x)
     - .. code-block:: cpp

            asinh(x), acosh(x)
     - .. code-block:: octave

            asinh(x), acosh(x)

   * - Inequalities
     - .. code-block:: python

            a<b, a<=b, a>b, a>=b
     - .. code-block:: cpp

            a<b, a<=b, a>b, a>=b
     - .. code-block:: octave

            a<b, a<=b, a>b, a>=b
   * - (*Not*) equal to
     - .. code-block:: python

            a==b, a!=b
     - .. code-block:: cpp

            a==b, a!=b
     - .. code-block:: octave

            a==b, a~=b
   * - Logical and
     - .. code-block:: python

            logic_and(a, b)
     - .. code-block:: cpp

            a && b
     - .. code-block:: octave

            a & b
   * - Logical or
     - .. code-block:: python

            logic_or(a, b)
     - .. code-block:: cpp

            a || b
     - .. code-block:: octave

            a | b
   * - Logical not
     - .. code-block:: python

            logic_not(a)
     - .. code-block:: cpp

            !a
     - .. code-block:: octave

            ~a
   * - Round to integer
     - .. code-block:: python

            floor(x), ceil(x)
     - .. code-block:: cpp

            floor(x), ceil(x)
     - .. code-block:: octave

            floor(x), ceil(x)
   * - Modulus after division
     - .. code-block:: python

            fmod(x, y)
     - .. code-block:: cpp

            fmod(x, y)
     - .. code-block:: octave

            mod(x, y)
   * - Modulus after division
     - .. code-block:: python

            fabs(x)
     - .. code-block:: cpp

            fabs(x)
     - .. code-block:: octave

            abs(x)
   * - Norm
     - .. code-block:: python

            norm_1(x)
     - .. code-block:: cpp

            norm_1(x)
     - .. code-block:: octave

            norm(x,1)
   * - Sign function
     - .. code-block:: python

            sign(x)
     - .. code-block:: cpp

            sign(x)
     - .. code-block:: octave

            sign(x)
   * - (Inverse) error function
     - .. code-block:: python

            erf(x), erfinv(x)
     - .. code-block:: cpp

            erf(x), erfinv(x)
     - .. code-block:: octave

            erf(x), erfinv(x)
   * - Elementwise min and max
     - .. code-block:: python

            fmin(x, y), fmax(x, y)
     - .. code-block:: cpp

            fmin(x, y), fmax(x, y)
     - .. code-block:: octave

            min(x, y), max(x, y)
   * - Global min and max
     - .. code-block:: python

            mmin(x), mmax(x)
     - .. code-block:: cpp

            mmin(x), mmax(x)
     - .. code-block:: octave

            mmin(x), mmax(x)
   * - Index of first nonzero
     - .. code-block:: python

            find(x)
     - .. code-block:: cpp

            find(x)
     - .. code-block:: octave

            find(x)
   * -  If-then-else
     - .. code-block:: python

            if_else(c, x, y)
     - .. code-block:: cpp

            if_else(c, x, y)
     - .. code-block:: octave

            if_else(c, x, y)
   * - Transpose
     - .. code-block:: python

            A.T
     - .. code-block:: cpp

            A.T()
     - .. code-block:: octave

            A',A.'
   * - Inner product
     - .. code-block:: python

            dot(x, y)
     - .. code-block:: cpp

            dot(x, y)
     - .. code-block:: octave

            dot(x, y)
   * - Horizontal/vertical concatenation
     - .. code-block:: python

            horzcat(x, y)
            vertcat(x, y)
            hcat([x, y])
            vcat([x, y])
     - .. code-block:: cpp

            horzcat(x, y)
            vertcat(x, y)
            MX::horzcat({x,y})
            MX::vertcat({x,y})
     - .. code-block:: octave

            [x,y]
            [x; y]
            c = {x,y};
            [c{:}]
            vertcat(c{:})
   * - Horizontal/vertical split (inverse of concatenation)
     - .. code-block:: python

            vertsplit(x)
            horzsplit(x)
     - .. code-block:: cpp

            vertsplit(x)
            horzsplit(x)
     - .. code-block:: octave

            vertsplit(x)
            horzsplit(x)
   * - Element access
     - .. code-block:: python

            # 0-based
            A[i,j]
            A[i]
     - .. code-block:: cpp

            // 0-based
            A(i,j)
            A(i)
     - .. code-block:: octave

            %  1-based
            A(i,j)
            A(i)

   * - Element assignment
     - .. code-block:: python

            # 0-based
            A[i,j] = b
            A[i] = b
     - .. code-block:: cpp

            // 0-based
            A(i,j) = b
            A(i) = b
     - .. code-block:: octave

            %  1-based
            A(i,j) = b
            A(i) = b
   * - Nonzero access
     - .. code-block:: python

            # 0-based
            A.nz[k]
     - .. code-block:: cpp

            // 0-based
            A.nz(k)
     - .. code-block:: octave

            (currently unsupported)
   * - Nonzero assignment
     - .. code-block:: python

            # 0-based
            A.nz[k] = b
     - .. code-block:: cpp

            // 0-based
            A.nz(k) = b
     - .. code-block:: octave

            (currently unsupported)
   * - Project to a different sparsity
     - .. code-block:: python

            project(x, s)
     - .. code-block:: cpp

            project(x, s)
     - .. code-block:: octave

            project(x, s)

