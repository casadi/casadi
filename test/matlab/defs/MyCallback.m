classdef MyCallback < casadi.Callback
  methods
    function self = MyCallback(name)
      self@casadi.Callback();
      construct(self, name);
    end
    function [res] = eval(self, arg)
      res{1} = arg{1}^2;
    end
  end
end

