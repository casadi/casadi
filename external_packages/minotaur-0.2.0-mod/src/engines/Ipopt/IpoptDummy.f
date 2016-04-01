c
c Ipopt depends on ma57/ma27 or mumps. ma57/ma27 depend on the 
c fortran library that they were compiled with. There is no way Cmake 
c or any other program can find out how Ipopt links to ma57, and how ma57 
c was compiled. To overcome this problem, we add this dummy fortran file, so 
c that CMake links to default fortran compiler. The hope is that the two fortran 
c compilers are the same. 
c See http://www.cmake.org/pipermail/cmake/2011-January/041992.html
c 


      subroutine mntrIpDummy ()
      write(nout,*) 'helloworld'

      return
      end
