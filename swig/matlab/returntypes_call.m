function varargout = returntypes_call(type,f, varargin)
  % Call the supplied function, applying 'full' or 'sparse' on each output
  %
  % Accepted values for type:
  %   'full'        -- Apply full
  %   'sparse'      -- Apply sparse
  %   'full|sparse' -- Apply full when dense, sparse otherwise
  % Type may be a cell array, indicating a type for each individual output
  [varargout{1:nargout}] = f(varargin{:});
  if ~iscell(type)
    typecell = cell(1,numel(varargout));
    for i=1:numel(varargout)
      typecell{i} = type;
    end
    type = typecell;
  end
  for i=1:numel(varargout)
    t = type{i};
    if strcmp(t,'full')
      varargout{i} = full(varargout{i});
    elseif strcmp(t,'sparse')
      varargout{i} = sparse(varargout{i});
    elseif strcmp(t,'full|sparse')
      if numel(varargout{i})==nnz(varargout{i})
        varargout{i} = full(varargout{i});
      else
        varargout{i} = sparse(varargout{i});
      end
    else
      error('Choose full, sparse or full|sparse');
    end
  end
end
