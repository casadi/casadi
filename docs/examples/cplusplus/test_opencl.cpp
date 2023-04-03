/*
 *    MIT No Attribution
 *
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
 *
 *    Permission is hereby granted, free of charge, to any person obtaining a copy of this
 *    software and associated documentation files (the "Software"), to deal in the Software
 *    without restriction, including without limitation the rights to use, copy, modify,
 *    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 *    permit persons to whom the Software is furnished to do so.
 *
 *    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 *    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 *    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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

int main(int, char**) {  
  cl_uint ret;
  size_t ret_size;

  cl_platform_id platform_id[max_num_platforms];
  cl_uint ret_num_platforms = 0;
  ret = clGetPlatformIDs(max_num_platforms, platform_id, &ret_num_platforms);
  if(ret != CL_SUCCESS) return 1;

  for(int i=0; i<ret_num_platforms && i<max_num_platforms; ++i){
    std::cout << "Platform " << i << " of " << ret_num_platforms << ":" << std::endl;
    cl_device_id device_id[max_num_devices];
    cl_uint ret_num_devices = 0;
    ret = clGetDeviceIDs(platform_id[i], CL_DEVICE_TYPE_DEFAULT, max_num_devices, device_id, &ret_num_devices);
    if(ret != CL_SUCCESS) return 1;
    
    for(int j=0; j<ret_num_devices && j<max_num_devices ; ++j){
      std::cout << "Device " << j << " of " << ret_num_devices << ":" << std::endl;
      
      // Separator
      std::cout << ">>>>>>>>>>>" << std::endl;
    
      // Device name
      char name[256];
      ret = clGetDeviceInfo(device_id[j],CL_DEVICE_NAME,sizeof(name),&name,&ret_size);
      if(ret != CL_SUCCESS) return 1;
      std::cout << "Name " << name << std::endl; 

      // Print type
      cl_device_type type;
      ret = clGetDeviceInfo(device_id[j],CL_DEVICE_TYPE,sizeof(type),&type,&ret_size);
      if(ret != CL_SUCCESS) return 1;
      std::cout << "Device type: ";
      switch(type){
        case CL_DEVICE_TYPE_CPU: cout << "CPU"; break;
        case CL_DEVICE_TYPE_GPU: cout << "GPU"; break;
        case CL_DEVICE_TYPE_ACCELERATOR: cout << "ACCELERATOR"; break;
        case CL_DEVICE_TYPE_DEFAULT: cout << "DEFAULT"; break;
        default: return 1;
      }
      std::cout << std::endl;
      
      // Separator
      std::cout << "<<<<<<<<<<<" << std::endl;
    }
  }
  return 0;
}
