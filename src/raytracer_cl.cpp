
#include <OpenCL/cl.h>
#include <sys/stat.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "camera.hpp"
#include "image.hpp"
#include "material/dielectric.hpp"
#include "material/lambertian.hpp"
#include "material/material.hpp"
#include "material/metal.hpp"
#include "ray.hpp"
#include "shape/hitable.hpp"
#include "shape/shape_list.hpp"
#include "shape/sphere.hpp"
#include "stb/stb_image.h"
#include "stb/stb_image_write.h"
#include "utils.hpp"
#include "vector.hpp"

using raytracer::Vec3;
using raytracer::Ray;
using raytracer::Image;
using raytracer::Sphere;
using raytracer::Hitable;
using raytracer::HitRecord;
using raytracer::ShapeList;
using raytracer::Camera;
using raytracer::Lambertian;
using raytracer::Metal;
using raytracer::Material;
using raytracer::Dielectric;

namespace {
constexpr char kFilename[] = "cl/raytracer.cl";
constexpr int kRow = 100;
constexpr int kCol = 200;
constexpr int kLocalSize = 100;
constexpr int kColorSize = kRow * kCol * kLocalSize * 3;
const char kOutputDir[] = "output";
}  // namespace

void Write(std::string filename, const raytracer::Image &image) {
  mkdir(kOutputDir, 0775);
  stbi_write_bmp((std::string(kOutputDir) + "/" + filename).c_str(), image.GetWidth(), image.GetHeight(), sizeof(raytracer::Color3D), image.GetPixels());
}

void Wait(cl_command_queue cmdQueue) {
  cl_event wait;
  cl_int status;
  
  status = clEnqueueMarker(cmdQueue, &wait);
  if (status != CL_SUCCESS) {
    std::cout << "clEnqueueMarker failed " << status << std::endl;
    exit(-1);
  }
  status = clWaitForEvents(1, &wait);
  if (status != CL_SUCCESS) {
    std::cout << "clWaitForEvents failed " << status << std::endl;
    exit(-1);
  }
}

int main(int argc, char **args) {
#ifdef _OPENMP
  double start = omp_get_wtime();
#else
  clock_t start = clock();
#endif
  if (argc < 3) {
    std::cout << "Too few arguments" << std::endl;
  }

  std::ifstream cl_file(kFilename);
  
  if (cl_file.fail()) {
    std::cout << "file cannot be opened." << kFilename << std::endl;
    exit(-1);
  }
  cl_int status;

  cl_platform_id platform;
  status = clGetPlatformIDs(1, &platform, NULL); 

  if (status != CL_SUCCESS) {
    std::cout << "clGetPlatformIDs failed " << status << std::endl;
    exit(-1);
  }

  cl_device_id device;
  status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
  if (status != CL_SUCCESS) {
    std::cout << "clGetDeviceIDs failed " << status << std::endl;
    exit(-1);
  }

  // allocate host memory
  float *h_results = new float[kColorSize];
  size_t data_size = kColorSize * sizeof(float);

  cl_int n_sphere = 488;
  cl_float4 *h_spheres = new cl_float4[n_sphere];
  // h_spheres[0] = {{0.0, 0.0, -1.0, 0.5}};
  // h_spheres[1] = {{0.0, -100.5, -1.0, 100.0}};
  // h_spheres[2] = {{1.0, 0.0, -1.0, 0.5}};
  // h_spheres[3] = {{-1.0, 0.0, -1.0, 0.5}};
  // h_spheres[4] = {{-1.0, 0.0, -1.0, -0.45}};
  cl_int *h_types = new cl_int[n_sphere];
  // h_types[0] = 0;
  // h_types[1] = 0;
  // h_types[2] = 1;
  // h_types[3] = 2;
  // h_types[4] = 2;
  cl_float4 *h_materials = new cl_float4[n_sphere];
  // h_materials[0] = {{0.8, 0.3, 0.3, 0.0}};
  // h_materials[1] = {{0.8, 0.8, 0.0, 0.0}};
  // h_materials[2] = {{0.8, 0.6, 0.2, 0.3}};
  // h_materials[3] = {{1.5, 0.0, 0.0, 0.0}};
  // h_materials[4] = {{1.5, 0.0, 0.0, 0.0}};
  int num = 0;
  h_spheres[num] = {{0.0, -1000.0, 0.0, 1000.0}};
  h_types[num] = 0;
  h_materials[num++] = {{0.5, 0.5, 0.5, 1.0}};
  for (int a = -11; a < 11; ++a) {
    for (int b = -11; b < 11; ++b) {
      float choose_mat = drand48();
      float x = a + 0.9 * drand48();
      float y = 0.2;
      float z = b + 0.9 * drand48();
      float adj[3] = {x - 4.0f, y - 0.2f, z};
      if (sqrt(adj[0] * adj[0] + adj[1] * adj[1] + adj[2] * adj[2]) > 0.9) {
        h_spheres[num] = {{x, y, z, 0.2}};
        if (choose_mat < 0.8) {
          h_types[num] = 0;
          h_materials[num] = {{static_cast<float>(drand48() * drand48()),
                               static_cast<float>(drand48() * drand48()),
                               static_cast<float>(drand48() * drand48()), 1.0}};
        } else if (choose_mat < 0.95) {
          h_types[num] = 1;
          h_materials[num] = {{static_cast<float>(0.5 + (1.0 + drand48())),
                               static_cast<float>(0.5 * (1.0 + drand48())),
                               static_cast<float>(0.5 * (1.0 + drand48())),
                               static_cast<float>(0.5 * drand48())}};
        } else {
          h_types[num] = 2;
          h_materials[num] = {{1.5, 1.0, 1.0, 1.0}};
        }
        ++num;
      }
    }
  }
  h_spheres[num] = {{0.0, 1.0, 0.0, 1.0}};
  h_types[num] = 2;
  h_materials[num++] = {{1.5, 1.0, 1.0, 1.0}};
  h_spheres[num] = {{-4.0, 1.0, 0.0, 1.0}};
  h_types[num] = 0;
  h_materials[num++] = {{0.4, 0.2, 0.1, 1.0}};
  h_spheres[num] = {{4.0, 1.0, 0.0, 1.0}};
  h_types[num] = 1;
  h_materials[num++] = {{0.7, 0.6, 0.5, 0.0}};

  n_sphere = num;

  cl_float3 *h_uvs = new cl_float3[kLocalSize];
  for (int i = 0, end = kLocalSize; i < end; ++i) {
    cl_float p1 = drand48();
    cl_float p2 = drand48();
    cl_float p3 = drand48();
    while (p1 * p1 + p2 * p2 + p3 * p3 >= 1.0) {
      p1 = drand48();
      p2 = drand48();
      p3 = drand48();
    }
    h_uvs[i] = {{p1, p2, p3}};
  }

  size_t sphere_size = n_sphere * sizeof(cl_float4);
  size_t type_size = n_sphere * sizeof(cl_int);
  size_t material_attr_size = n_sphere * sizeof(cl_float4);
  size_t uv_size = sizeof(cl_float3) * kLocalSize;

  // create openCL context
  cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, &status);
  if (status != CL_SUCCESS) {
    std::cout << "clCreateContext failed " << status << std::endl;
    exit(-1);
  }

  // create openCL command queue
  cl_command_queue cmd_queue = clCreateCommandQueue(context, device, 0, &status);
  if (status != CL_SUCCESS) {
    std::cout << "clCreateCommandQueue failed " << status << std::endl;
    exit(-1);
  }

  // allocate device buffer
  cl_mem d_spheres = clCreateBuffer(context, CL_MEM_READ_ONLY, sphere_size, NULL, &status);
  if (status != CL_SUCCESS) {
    std::cout << "clCreateBuffer(0) failed " << status << std::endl;
    exit(-1);
  }
  cl_mem d_types = clCreateBuffer(context, CL_MEM_READ_ONLY, type_size, NULL, &status);
  if (status != CL_SUCCESS) {
    std::cout << "clCreateBuffer(1) failed " << status << std::endl;
    exit(-1);
  }
  cl_mem d_materials = clCreateBuffer(context, CL_MEM_READ_ONLY, material_attr_size, NULL, &status);
  if (status != CL_SUCCESS) {
    std::cout << "clCreateBuffer(2) failed " << status << std::endl;
    exit(-1);
  }
  cl_mem d_results = clCreateBuffer(context, CL_MEM_WRITE_ONLY, data_size, NULL, &status);
  if (status != CL_SUCCESS) {
    std::cout << "clCreateBuffer(3) failed " << status << std::endl;
    exit(-1);
  }
  cl_mem d_uvs = clCreateBuffer(context, CL_MEM_READ_ONLY, uv_size, NULL, &status);
  if (status != CL_SUCCESS) {
    std::cout << "clCreateBuffer(4) failed " << status << std::endl;
    exit(-1);
  }

  // enqueue the command to write the data from the host buffers to the device buffers:
  status = clEnqueueWriteBuffer(cmd_queue, d_spheres, CL_FALSE, 0, sphere_size, h_spheres, 0, NULL, NULL);
  if (status != CL_SUCCESS) {
    std::cout << "clEnqueueWriteBuffer(0) failed " << status << std::endl;
    exit(-1);
  }
  status = clEnqueueWriteBuffer(cmd_queue, d_types, CL_FALSE, 0, type_size, h_types, 0, NULL, NULL);
  if (status != CL_SUCCESS) {
    std::cout << "clEnqueueWriteBuffer(1) failed " << status << std::endl;
    exit(-1);
  }
  status = clEnqueueWriteBuffer(cmd_queue, d_materials, CL_FALSE, 0, material_attr_size, h_materials, 0, NULL, NULL);
  if (status != CL_SUCCESS) {
    std::cout << "clEnqueueWriteBuffer(2) failed " << status << std::endl;
    exit(-1);
  }
  status = clEnqueueWriteBuffer(cmd_queue, d_uvs, CL_FALSE, 0, uv_size, h_uvs, 0, NULL, NULL);
  if (status != CL_SUCCESS) {
    std::cout << "clEnqueueWriteBuffer(3) failed " << status << std::endl;
    exit(-1);
  }

  Wait(cmd_queue);

  // read kernel code from file
  cl_file.seekg (0, cl_file.end);
  int length = cl_file.tellg();
  cl_file.seekg (0, cl_file.beg);

  char *cl_program_text = new char[length + 1];
  cl_file.read(cl_program_text, length);
  cl_program_text[length] = '\0';

  cl_file.close();

  char *strings[1];
  strings[0] = cl_program_text;

  cl_program program = clCreateProgramWithSource(context, 1, (const char**)strings, NULL, &status);
  if (status != CL_SUCCESS) {
    std::cout << "clCreateProgramWithSource failed " << status << std::endl;
    exit(-1);
  }
  delete[] cl_program_text;

  // compile and link the kernel code
  char *options = {""};
  status = clBuildProgram(program, 1, &device, options, NULL, NULL);
  if (status != CL_SUCCESS) {
    size_t size;
    clGetProgramBuildInfo( program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &size );
    cl_char *log = new cl_char[ size ];
    clGetProgramBuildInfo( program, device, CL_PROGRAM_BUILD_LOG, size, log, NULL );
    std::cout << "clBuildProgram failed:" << log << std::endl;
    delete [ ] log;
  }
  // create kernel object
  cl_kernel kernel = clCreateKernel(program, args[1], &status);
  if (status != CL_SUCCESS) {
    std::cout << "clCreateKernel failed " << status << std::endl;
    exit(-1);
  }

  // setup arguments to the kernel object
  status = clSetKernelArg(kernel, 0, sizeof(cl_mem), &d_results);
  if (status != CL_SUCCESS) {
    std::cout << "clSetKernelArg(0) failed " << status << std::endl;
    exit(-1);
  }
  status = clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_spheres);
  if (status != CL_SUCCESS) {
    std::cout << "clSetKernelArg(1) failed " << status << std::endl;
    exit(-1);
  }
  status = clSetKernelArg(kernel, 2, sizeof(cl_int), &n_sphere);
  if (status != CL_SUCCESS) {
    std::cout << "clSetKernelArg(2) failed " << status << std::endl;
    exit(-1);
  }
  status = clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_types);
  if (status != CL_SUCCESS) {
    std::cout << "clSetKernelArg(3) failed " << status << std::endl;
    exit(-1);
  }
  status = clSetKernelArg(kernel, 4, sizeof(cl_mem), &d_materials);
  if (status != CL_SUCCESS) {
    std::cout << "clSetKernelArg(4) failed " << status << std::endl;
    exit(-1);
  }
  status = clSetKernelArg(kernel, 5, sizeof(cl_mem), &d_uvs);

  if (status != CL_SUCCESS) {
    std::cout << "clSetKernelArg(5) failed " << status << std::endl;
    exit(-1);
  }
  // enqueue the kernel object for execution
  size_t global_worker_size[3] = {kCol, kRow, kLocalSize};
  size_t local_worker_size[3] = {kLocalSize, 1, 1};

  Wait(cmd_queue);

  status = clEnqueueNDRangeKernel(cmd_queue, kernel, 3, NULL, global_worker_size, NULL, 0, NULL, NULL);
  if (status != CL_SUCCESS) {
    std::cout << "clEnqueueNDRangeKernel failed " << status << std::endl;
    exit(-1);
  }

  Wait(cmd_queue);

  // read the result
  status = clEnqueueReadBuffer(cmd_queue, d_results, CL_TRUE, 0, data_size, h_results, 0, NULL, NULL);
  if (status != CL_SUCCESS) {
    std::cout << "clEnqueueReadBuffer failed " << status << std::endl;
    exit(-1);
  }

  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(cmd_queue);
  clReleaseMemObject(d_results);
  clReleaseMemObject(d_spheres);

  delete[] h_spheres;
  Image image(kCol, kRow);
  int i = 0;
  for (int y = 0; y < kRow; ++y) {
    for (int x = 0; x < kCol; ++x) {
      float r = 0.0;
      float g = 0.0;
      float b = 0.0;
      for (int s = 0; s < kLocalSize; ++s) {
        r += h_results[i++];
        g += h_results[i++];
        b += h_results[i++];
      }
      r /= static_cast<float>(kLocalSize);
      g /= static_cast<float>(kLocalSize);
      b /= static_cast<float>(kLocalSize);
      image.WritePixel(x, y, sqrt(r), sqrt(g), sqrt(b));
      // image.WritePixel(x, y, r, g, b);
      // std::cout << "(" << x << ", " << y << "): " << r << "," << g << "," << b << " ";
    }
    // std::cout << std::endl;
  }
#ifdef _OPENMP
  double end = omp_get_wtime();
  std::cout << "Done! " << (end - start) * 1000.0 << " msec" << std::endl;
#else
  clock_t end = clock();
  std::cout << "Done: " << (end - start) / float(CLOCKS_PER_SEC) * 1000.0 << " msec" << std::endl;
#endif

  Write(args[2], image);

  return 0;
}