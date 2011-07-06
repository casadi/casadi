casadi

a = zeros(2,3)
b = DMatrix(2,3,1)
c = symbolic("x",2,3)
d = MX("x",2,3)


nums = {a,b};
syms = {c,d};

for i=1:2
  for j=1:2
    n = nums{i}
    s = syms{j}
    [n s]
    [s n]
    [n;s]
    [s;n]
    [s;s]
  end
end

