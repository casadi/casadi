/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <iostream>  
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
#include <vector>  
#include <iostream>  

const int max_num_platforms = 10;
const int max_num_devices = 10;

using namespace std;

int main(int, char**) {  
  cl_uint ret;
  size_t ret_size;

  cl_platform_id platform_id[max_num_platforms];
  cl_uint ret_num_platforms = 0;
  ret = clGetPlatformIDs(max_num_platforms, platform_id, &ret_num_platforms);
  if(ret != CL_SUCCESS) return 1;

  for(int i=0; i<ret_num_platforms && i<max_num_platforms; ++i){
    cout << "Platform " << i << " of " << ret_num_platforms << ":" << endl;
    cl_device_id device_id[max_num_devices];
    cl_uint ret_num_devices = 0;
    ret = clGetDeviceIDs(platform_id[i], CL_DEVICE_TYPE_DEFAULT, max_num_devices, device_id, &ret_num_devices);
    if(ret != CL_SUCCESS) return 1;
    
    for(int j=0; j<ret_num_devices && j<max_num_devices ; ++j){
      cout << "Device " << j << " of " << ret_num_devices << ":" << endl;
      
      // Separator
      cout << ">>>>>>>>>>>" << endl;
    
      // Device name
      char name[256];
      ret = clGetDeviceInfo(device_id[j],CL_DEVICE_NAME,sizeof(name),&name,&ret_size);
      if(ret != CL_SUCCESS) return 1;
      cout << "Name " << name << endl; 

      // Print type
      cl_device_type type;
      ret = clGetDeviceInfo(device_id[j],CL_DEVICE_TYPE,sizeof(type),&type,&ret_size);
      if(ret != CL_SUCCESS) return 1;
      cout << "Device type: ";
      switch(type){
        case CL_DEVICE_TYPE_CPU: cout << "CPU"; break;
        case CL_DEVICE_TYPE_GPU: cout << "GPU"; break;
        case CL_DEVICE_TYPE_ACCELERATOR: cout << "ACCELERATOR"; break;
        case CL_DEVICE_TYPE_DEFAULT: cout << "DEFAULT"; break;
        default: return 1;
      }
      cout << endl;
      
      // Separator
      cout << "<<<<<<<<<<<" << endl;
    }
  }
  return 0;
}
