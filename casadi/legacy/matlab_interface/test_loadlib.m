shrlib = '/home/janderss/dev/OPTICON-SOFTWARE/build/lib/libcasadi.so';
hfile = '/home/janderss/dev/OPTICON-SOFTWARE/casadi/c_interface/c_interface.h';

if(libisloaded('libcasadi'))
    unloadlibrary('libcasadi');
end
[notfound,warnings] =  loadlibrary(shrlib,hfile,'alias','libcasadi')

% shr_ipopt = '/usr/local/lib/libipopt.so';
% hdr_ipopt  = '/usr/local/include/coin/IpStdCInterface.h';
% if(libisloaded('libipopt'))
%     unloadlibrary('libipopt');
% end
% [notfound,warnings] =  loadlibrary(shr_ipopt,hdr_ipopt,'alias','libipopt')


% shrlib2 = '/home/janderss/dev/OPTICON-SOFTWARE/build/lib/libipopt_interface.so';
% hfile2 = '/home/janderss/dev/OPTICON-SOFTWARE/ipopt_interface/ipopt_solver_c.h';
% if(libisloaded('libcasadi_ipopt'))
%     unloadlibrary('libcasadi_ipopt');
% end
% [notfound,warnings] =  loadlibrary(shrlib2,hfile2,'alias','libcasadi_ipopt')



