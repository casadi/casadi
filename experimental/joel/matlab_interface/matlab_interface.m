root = '/home/janderss/dev/casadi/trunk/experimental/joel/matlab_interface/';
libname = 'libcasadi_matlab'
shrlib = [root 'build/' libname '.so'];
hfile = [root 'matlab_interface.h'];

if(libisloaded(libname))
    unloadlibrary(libname);
end
[notfound,warnings] =  loadlibrary(shrlib,hfile,'alias',libname);

calllib(libname,'test',4)

a = 5;
calllib(libname,'test_mex',2)
calllib(libname,'test_mex',ones(2,3))
calllib(libname,'test_mex',2.3)
calllib(libname,'test_mex',[1,2,4])
libfunctions(libname)
