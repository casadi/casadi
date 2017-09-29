classdef OptiCallbackHelper < casadi.OptiCallback

  properties
    callback
  end  
        
  methods
    function self = OptiCallbackHelper(callback)
      self@casadi.OptiCallback();
      self.callback = callback;
    end
    function [] = call(self, i)
      self.callback(i);
    end    
  end
  
  methods(Static)
    function [] = callback_setup(self, varargin)
      if length(varargin)==1
        fh = varargin{1};
      else
        self.callback_class();
        drawnow
        pause(0.0001)
        return
      end
      persistent cb;
      if isempty(cb)
        cb = {};
      end
      cb = {cb{:} casadi.OptiCallbackHelper(fh)};
      self.callback_class(cb{end});
    end
  end
end
