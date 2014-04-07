system('gcc -fPIC -shared -Os base.cpp -o base.dylib')
system('gcc -fPIC -shared -Os derived.cpp -o derived.dylib')
mex 'base_wrap.cpp'
mex 'derived_wrap.cpp'

