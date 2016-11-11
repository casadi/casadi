function ret = returntypes(type,f)
  % Wraps a function, applying 'full' or 'sparse' on each output
  %
  % Accepted values for type:
  %   'full'        -- Apply full
  %   'sparse'      -- Apply sparse
  %   'full|sparse' -- Apply full when dense, sparse otherwise
  % Type may be a cell array, indicating a type for each individual output
  ret = @(varargin) returntypes_call(type,f,varargin{:});
end
