root = '/home/janderss/dev/casadi/trunk/experimental/joel/matlab_interface/';
libname = 'libcasadi_matlab'
shrlib = [root 'build/' libname '.so'];
hfile = [root 'matlab_interface.h'];

if(libisloaded(libname))
    unloadlibrary(libname);
end
[notfound,warnings] =  loadlibrary(shrlib,hfile,'alias',libname);

%  calllib(libname,'test',4)


a = 5;
%  b = calllib(libname,'convert_to_swig',2)
%  c = calllib(libname,'convert_from_swig',b)
%  
%  b = calllib(libname,'convert_to_swig',ones(2,3))
%  c = calllib(libname,'convert_from_swig',b)
%  
%  b = calllib(libname,'convert_to_swig',2.3)
%  c = calllib(libname,'convert_from_swig',b)
%  
%  b = calllib(libname,'convert_to_swig',[1,2,4])
%  c = calllib(libname,'convert_from_swig',b)

%  d = libpointer('int32', 5)
d = libpointer
d2 = swig_ref
b1 = calllib(libname,'convert_to_swig',d)
b2 = calllib(libname,'convert_from_swig',d)
b3 = calllib(libname,'convert_to_swig',d2)




libfunctions(libname)
