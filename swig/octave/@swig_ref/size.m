function s=size(x,varargin)
  s = [ x.size1() x.size2() ];
  if nargin>1
    s=s(varargin{1});
  end
end
