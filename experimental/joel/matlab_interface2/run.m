system('gcc -fPIC -shared -Os base.cpp -o base.dylib')
system('gcc -fPIC -shared -Os derived.cpp -o derived.dylib')
mex base_wrap.cpp base.dylib
mex derived_wrap.cpp derived.dylib base.dylib

b = derived_wrap()
disp b
base_wrap(b)


