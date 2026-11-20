namespace casadi{


%extend Matrix<double> {

  %matlabcode %{
     function s = saveobj(obj)
        try
            s.serialization = obj.serialize();
        catch exception
            warning(['Serializing of CasADi DM failed:' getReport(exception) ]);
            s = struct;
        end
     end
  %}
  %matlabcode_static %{
     function obj = loadobj(s)
        try
          if isstruct(s)
             obj = casadi.DM.deserialize(s.serialization);
          else
             obj = s;
          end
        catch exception
            warning(['Serializing of CasADi DM failed:' getReport(exception) ]);
            s = struct;
        end
     end
  %}
}

%extend Sparsity {
  %matlabcode %{
     function s = saveobj(obj)
        try
            s.serialization = obj.serialize();
        catch exception
            warning(['Serializing of CasADi Sparsity failed:' getReport(exception) ]);
            s = struct;
        end
     end
  %}
  %matlabcode_static %{
     function obj = loadobj(s)
        try
          if isstruct(s)
             obj = casadi.Sparsity.deserialize(s.serialization);
          else
             obj = s;
          end
        catch exception
            warning(['Serializing of CasADi Sparsity failed:' getReport(exception) ]);
            s = struct;
        end
     end
  %}

  // Transpose using the A' syntax in addition to A.'
  Sparsity ctranspose() const { return $self->T();}
}


%extend Function {

  %matlabcode %{
     function s = saveobj(obj)
        try
            s.serialization = obj.serialize();
        catch exception
            warning(['Serializing of CasADi Function failed:' getReport(exception) ]);
            s = struct;
        end
     end
  %}
  %matlabcode_static %{
     function obj = loadobj(s)
        try
          if isstruct(s)
             obj = casadi.Function.deserialize(s.serialization);
          else
             obj = s;
          end
        catch exception
            warning(['Serializing of CasADi Function failed:' getReport(exception) ]);
            s = struct;
        end
     end
  %}

}

} // namespace casadi

%extend casadi::DeserializerBase {
  %matlabcode %{
    function out = unpack(self)
      type = casadi.SerializerBase.type_to_string(self.internal_pop_type);
      out = self.(['blind_unpack_' type]);
    end
  %}
}

%extend casadi::Opti {
  %matlabcode %{
    function [] = callback(self, varargin)
      casadi.OptiCallbackHelper.callback_setup(self, varargin{:})
    end
  %}
}

%{
#ifdef HAVE_OCTAVE
  // Mandatory as of Octave 10
  extern "C" const int __octave_mex_soversion__ = 1;
#endif
%}

namespace casadi{
%extend PrintableCommon {
  %matlabcode %{
    function s = repr(self)
      s = [self.type_name() '(' self.str() ')'];
    end
  %}
}
} // namespace casadi
