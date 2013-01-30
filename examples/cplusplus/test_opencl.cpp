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
  ret = clGetPlatformIDs(1, platform_id, &ret_num_platforms);
  if(ret != CL_SUCCESS) return 1;
  cout << ret_num_platforms << " platforms" << endl; 

  for(int i=0; i<ret_num_platforms; ++i){
    cout << "Platform " << i << ":" << endl;
    cl_device_id device_id[max_num_devices];
    cl_uint ret_num_devices = 0;
    ret = clGetDeviceIDs(platform_id[i], CL_DEVICE_TYPE_DEFAULT, 1, device_id, &ret_num_devices);
    if(ret != CL_SUCCESS) return 1;
    cout << ret_num_devices << " devices" << endl; 
    
    for(int j=0; j<ret_num_devices; ++j){
      cout << "Device " << j << ":" << endl;
      
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
    }
  }
  return 0;
}
