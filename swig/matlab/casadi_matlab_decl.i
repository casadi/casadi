%typemap(argout, noblock=1) std::ostream &stream {
  *$1 << "\n" << std::flush;
}

// Abstract mixins avoid the Octave 10+ SwigRef diamond; before their %include
%feature("abstract") casadi::PrintableCommon;
%feature("abstract") casadi::GenericExpressionCommon;
%feature("abstract") casadi::GenericMatrixCommon;
%feature("abstract") casadi::SparsityInterfaceCommon;
%feature("abstract") casadi::MatrixCommon;
// Applied to the template, not to GenSharedObject
%feature("abstract") casadi::GenericShared<casadi::SharedObject, casadi::SharedObjectInternal>;
%feature("abstract") casadi::SharedObject;
%feature("abstract") casadi::IndexAbstraction;

%typemap(in, doc="index", pystub_in="int", noblock=1) casadi_index {
  if (!swig::assign($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input $argnum to type ' index '.");
  if ($1==0) SWIG_exception_fail(SWIG_TypeError,"Index starts at 1, got index '0'.");
  if ($1>=1) $1--;
}

%rename(uminus) operator-;
// Pattern intersection is elementwise; * and mtimes give the pattern of the matrix product
%rename(times) casadi::Sparsity::operator*;
// inv is the matrix inverse in MATLAB; elementwise inverse is 1./x
%ignore casadi::MX::inv() const;
// Would hide the einstein method of GenericMatrixCommon for MX arguments
%ignore casadi::MX::einstein;
%feature("varargin","1") casadi_vertcat;
%feature("varargin","1") casadi_horzcat;
%feature("varargin","1") casadi_diagcat;
%feature("varargin","1") casadi_veccat;
%feature("optionalunpack","1") size;

// Raise an error if "this" not correct
%typemap(check, noblock=1) SWIGTYPE *self %{
if (!$1) {
  SWIG_Error(SWIG_RuntimeError, "Invalid 'self' object");
  SWIG_fail;
 }
%}

// Workarounds, pending proper fix
%rename(nonzero) __nonzero__;
%rename(hash) __hash__;

%rename(rem) casadi_mod;

%{
  namespace casadi {
    /// Helper function: Convert ':' to Slice
    inline Slice char2Slice(char ch) {
      casadi_assert_dev(ch==':');
      return Slice();
    }
  } // namespace casadi
%}


namespace casadi{
%extend GenericMatrixCommon {
  %matlabcode %{
    function varargout = spy(self,varargin)
      [varargout{1:nargout}] = spy(sparse(casadi.DM(self.sparsity(),1)),varargin{:});
    end
    function varargout = subsref(self,s)
      if strcmp(s(1).type,'()')
        ret = paren(self, s(1).subs{:});
      elseif strcmp(s(1).type,'{}')
        ret = brace(self, s(1).subs{:});
      else
        [varargout{1:nargout}] = builtin('subsref',self,s);
        return
      end
      % Remaining subscripts act on the result, e.g. x(1,:).dim()
      if numel(s)==1
        varargout{1} = ret;
      else
        [varargout{1:nargout}] = subsref(ret, s(2:end));
      end
    end
    function n = numArgumentsFromSubscript(self,s,context)
      % One output for e.g. disp(x.sparsity); statements keep the default
      if ~exist('OCTAVE_VERSION', 'builtin') && context==matlab.mixin.util.IndexingContext.Statement
        n = builtin('numArgumentsFromSubscript',self,s,context);
      else
        n = 1;
      end
    end
    function self = subsasgn(self,s,v)
      if numel(s)==1 && strcmp(s.type,'()')
        paren_asgn(self, v, s.subs{:});
      elseif numel(s)==1 && strcmp(s.type,'{}')
        brace_asgn(self, v, s.subs{:});
      else
        self = builtin('subsasgn',self,s,v);
      end
    end
    function out = norm(self,varargin)
      narginchk(1,2);
      % 2-norm by default
      if nargin==1
        ind = 2;
      else
        ind = varargin{1};
      end
      % Typecheck
      assert((isnumeric(ind) && isscalar(ind)) || ischar(ind))
      % Pick the right norm
      if isnumeric(ind)
        switch ind
          case 1
            out = norm_1(self);
          case 2
            out = norm_2(self);
          case inf
            out = norm_inf(self);
          otherwise
            error(sprintf('Unknown norm argument: %g', ind))
        end
      else
        switch ind
          case 'fro'
            out = norm_fro(self);
          case 'inf'
            out = norm_inf(self);
          otherwise
            error(sprintf('Unknown norm argument: ''%s''', ind))
        end
      end
    end
    function out = min(varargin)
      narginchk(1,2);
      if nargin==1
        out = mmin(varargin{1});
      else
        out = fmin(varargin{1}, varargin{2});
      end
    end
    function out = max(varargin)
      narginchk(1,2);
      if nargin==1
        out = mmax(varargin{1});
      else
        out = fmax(varargin{1}, varargin{2});
      end
    end
    function b = isrow(self)
      b = is_row(self);
    end
    function b = iscolumn(self)
      b = is_column(self);
    end
    function b = isvector(self)
      b = is_vector(self);
    end
    function b = isscalar(self)
      b = is_scalar(self);
    end
  %}
}
%extend DaeBuilder {
  %matlabcode %{
    function varargout = subsref(self,s)
      if strcmp(s(1).type,'()')
        ret = paren(self, s(1).subs{:});
        if numel(s)==1
          varargout{1} = ret;
        else
          [varargout{1:nargout}] = subsref(ret, s(2:end));
        end
      else
        [varargout{1:nargout}] = builtin('subsref',self,s);
      end
    end
    function n = numArgumentsFromSubscript(self,s,context)
      % One output for e.g. fprintf('%g', dae.start_time); statements keep the default
      if ~exist('OCTAVE_VERSION', 'builtin') && context==matlab.mixin.util.IndexingContext.Statement
        n = builtin('numArgumentsFromSubscript',self,s,context);
      else
        n = 1;
      end
    end
  %}
}
%extend Function {
  %matlabcode %{
    function varargout = subsref(self,s)
      if numel(s)==1 && strcmp(s.type,'()')
        [varargout{1:nargout}]= paren(self, s.subs{:});
      else
        [varargout{1:nargout}] = builtin('subsref',self,s);
      end
   end
   function varargout = paren(self, varargin)
      if nargin==1 || (nargin>=2 && ischar(varargin{1}))
        % Named inputs: return struct
        assert(nargout<2, 'Syntax error');
        assert(mod(nargin,2)==1, 'Syntax error');
        arg = struct;
        for i=1:2:nargin-1
          assert(ischar(varargin{i}), 'Syntax error');
          arg.(varargin{i}) = varargin{i+1};
        end
        res = self.call(arg);
        varargout{1} = res;
      else
        % Ordered inputs: return variable number of outputs
        res = self.call(varargin);
        assert(nargout<=numel(res), 'Too many outputs');
        for i=1:max(min(1,numel(res)),nargout)
          varargout{i} = res{i};
        end
      end
    end
  %}
 }
}

  %rename(internal_variable) casadi::Opti::variable;
  %rename(internal_parameter) casadi::Opti::parameter;
  %rename(internal_subject_to) casadi::Opti::subject_to;
  %extend casadi::Opti {
    %matlabcode %{
      function out = variable(self, varargin)
        st = dbstack('-completenames',1);
        if length(st)>0
          meta = struct('stacktrace', {num2cell(st)});
        else
          meta = struct;
        end
        out = self.internal_variable(varargin{:});
        self.update_user_dict(out, meta);
      end
      function out = parameter(self, varargin)
        st = dbstack('-completenames',1);
        if length(st)>0
          meta = struct('stacktrace', {num2cell(st)});
        else
          meta = struct;
        end
        out = self.internal_parameter(varargin{:});
        self.update_user_dict(out, meta);
      end
      function [] = subject_to(self, varargin)
        if length(varargin)==0
          self.internal_subject_to();
          return;
        end
        st = dbstack('-completenames',1);
        if length(st)>0
          meta = struct('stacktrace', {num2cell(st)});
        else
          meta = struct;
        end
        self.internal_subject_to(varargin{:});
        self.update_user_dict(varargin{1}, meta);
      end
    %}
  }
// solver() without arguments is the deprecated MATLAB name of casadi_solver()
%rename(casadi_solver) casadi::OptiAdvanced::casadi_solver;
%extend casadi::OptiAdvanced {
  %matlabcode %{
    function varargout = solver(self, varargin)
      if isempty(varargin)
        warning('casadi:deprecated', ...
          'OptiAdvanced.solver() is deprecated, use casadi_solver() instead.');
        [varargout{1:nargout}] = self.casadi_solver();
      else
        [varargout{1:nargout}] = solver@casadi.Opti(self, varargin{:});
      end
    end
  %}
}

%rename("internal_pop_type") casadi::DeserializerBase::pop_type;
%rename("%(regex:/(SERIALIZED_\w+)/internal_\\1/)s", regextarget=1, fullname=1) "casadi::SerializerBase::SERIALIZED_\w+";
