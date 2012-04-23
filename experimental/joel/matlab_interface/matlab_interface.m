%
%     This file is part of CasADi.
% 
%     CasADi -- A symbolic framework for dynamic optimization.
%     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
% 
%     CasADi is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 3 of the License, or (at your option) any later version.
% 
%     CasADi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public
%     License along with CasADi; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
% 
% 
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
