#include <iostream>
#include <cstring>
#include <stdio.h>
#include <algorithm>
#include <vector>

#include "wsqic.hpp"

// http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html

int main(int argc, char *argv[])
{
  int n = 15;
  int m = 18;
  std::vector<double> bl(n+m);
  std::vector<double> bu(n+m);

  double inf = 1.0e+20;

  std::fill(bl.begin(),bl.begin()+n,0);

  std::fill(bu.begin(),bu.begin()+n,inf);


  // http://apmonitor.com/wiki/uploads/Apps/hs118.apm

  // LBX
  bl[0] =   8.0;
  bu[0] =  21.0;

  bl[1] =  43.0;
  bu[1] =  57.0;

  bl[2] =   3.0;
  bu[2] =  16.0;

  bu[3] =  90.0;
  bu[4] = 120.0;
  bu[5] =  60.0;
  bu[6] =  90.0;
  bu[7] = 120.0;
  bu[8] =  60.0;
  bu[9] =  90.0;
  bu[10] = 120.0;
  bu[11] =  60.0;
  bu[12] =  90.0;
  bu[13] = 120.0;
  bu[14] =  60.0;

  //
  std::fill(bl.begin()+n,bl.begin()+n+m,0);
  std::fill(bu.begin()+n,bu.begin()+n+m,inf);

  //iObj   = 18  means the linear objective is row 18 in valA(*).
  //The objective row is free.

  int iObj   = 17;
  bl[n+iObj] = -inf;

  // LBG
  bl[n+0]  =  -7.0;
  bu[n+0]  =   6.0;

  bl[n+1]  =  -7.0;
  bu[n+1]  =   6.0;

  bl[n+2]  =  -7.0;
  bu[n+2]  =   6.0;

  bl[n+3]  =  -7.0;
  bu[n+3]  =   6.0;

  bl[n+4]  =  -7.0;
  bu[n+4]  =   7.0;

  bl[n+5]  =  -7.0;
  bu[n+5]  =   7.0;

  bl[n+6]  =  -7.0;
  bu[n+6]  =   7.0;

  bl[n+7]  =  -7.0;
  bu[n+7]  =   7.0;

  bl[n+8]  =  -7.0;
  bu[n+8]  =   6.0;

  bl[n+9] =  -7.0;
  bu[n+9] =   6.0;

  bl[n+10] =  -7.0;
  bu[n+10] =   6.0;

  bl[n+11] =  -7.0;
  bu[n+11] =   6.0;

  bl[n+12] =  60.0;
  bl[n+13] =  50.0;
  bl[n+14] =  70.0;
  bl[n+15] =  85.0;
  bl[n+16] = 100.0;


  std::vector<double> x0(n);

  x0[ 0]  =  20.0;
  x0[ 1]  =  55.0;
  x0[ 2]  =  15.0;
  x0[ 3]  =  20.0;
  x0[ 4]  =  60.0;
  x0[ 5]  =  20.0;
  x0[ 6]  =  20.0;
  x0[ 7]  =  60.0;
  x0[ 8]  =  20.0;
  x0[9]  =  20.0;
  x0[10]  =  60.0;
  x0[11]  =  20.0;
  x0[12]  =  20.0;
  x0[13]  =  60.0;
  x0[14]  =  20.0;

  sqic(&m , &n, &bl[0], &bu[0]);

  return 0;
}
