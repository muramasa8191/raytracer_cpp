#include <sys/stat.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <OpenCL/opencl.h>
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
  const char kOutputDir[] = "output";
}
void Write(std::string filename, const raytracer::Image &image) {
  mkdir(kOutputDir, 0775);
  stbi_write_bmp((std::string(kOutputDir) + "/" + filename).c_str(), image.GetWidth(), image.GetHeight(), sizeof(raytracer::Color3D), image.GetPixels());
}

std::string gradation(int nx, int ny, raytracer::Image &image) {
  for (int y = ny-1; y >= 0; --y) {
    for (int x = 0; x < nx; ++x) {
      float r = float(x) / nx;
      float g = float(y) / ny;
      float b = 0.5;
      image.WritePixel(x, y, r, g, b);
    }
  }
  return "gradation.bmp";
}

std::string DraySky(int nx, int ny, Image &image) {
  
  Vec3 low_left(-2.0, -1.0, -1.0);
  Vec3 horizontal(4.0, 0.0, 0.0);
  Vec3 vertical(0.0, 2.0, 0.0);
  Vec3 origin(0.0, 0.0, 0.0);

  for (int y = ny-1; y >= 0; --y) {
    for (int x = 0; x < nx; ++x) {
      float u = float(x) / float(nx);
      float v = float(y) / float(ny);
      Ray ray(origin, low_left + u * horizontal + v * vertical);

      Vec3 unit_dir = raytracer::UnitVector(ray.GetDirection());
      float t = 0.5 * (unit_dir.y() + 1.0);
      Vec3 col = (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0);
      image.WritePixel(x, y, col[0], col[1], col[2]);
    }
  }
  return "sky.bmp";
}

Vec3 Sky(Ray ray) {
  Vec3 unit_dir = raytracer::UnitVector(ray.GetDirection());
  float t = 0.5 * (unit_dir.y() + 1.0);
  return (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0);
}

std::string DrawShpere(int nx, int ny, Image &image) {
  Vec3 low_left(-2.0, -1.0, -1.0);
  Vec3 horizontal(4.0, 0.0, 0.0);
  Vec3 vertical(0.0, 2.0, 0.0);
  Vec3 origin(0.0, 0.0, 0.0);

  Vec3 center(0.0f, 0.0f, -1.0f);
  Sphere sphere(center, 0.5f);
  for (int y = ny-1; y >= 0; --y) {
    for (int x = 0; x < nx; ++x) {
      float u = float(x) / float(nx);
      float v = float(y) / float(ny);
      Ray ray(origin, low_left + u * horizontal + v * vertical);

      Vec3 col;
      float t = sphere.HitF(ray);
      if (t > 0.0) {
        Vec3 N = raytracer::UnitVector(ray.PointAt(t) - Vec3(0.0, 0.0, -1.0));
        col = 0.5 * Vec3(N.x() + 1.0, N.y() + 1.0, N.z() + 1.0);
      } else {
        col = Sky(ray);
      }
      image.WritePixel(x, y, col[0], col[1], col[2]);
    }
  }
  return "sphere.bmp";
}
Vec3 Color(const Ray &ray, const ShapeList &world, int depth) {
  HitRecord rec;
  if (world.Hit(ray, 0.001, MAXFLOAT, rec)) {
    Ray scattered;
    Vec3 attenuation;
    if (depth < 50 && rec.mat->scatter(ray, rec, attenuation, scattered)) {
      return attenuation * Color(scattered, world, depth + 1);
    }
    return Vec3(0.0, 0.0, 0.0);
  }
  return Sky(ray);
}

Vec3 Color3(const Ray &ray, const ShapeList &world) {
  HitRecord rec;
  if (world.Hit(ray, 0.001, MAXFLOAT, rec)) {
    Vec3 target = rec.p + rec.normal + raytracer::RandomInUnitSphere();
    return 0.5 * Color3(Ray(rec.p, target - rec.p), world);
  }
  return Sky(ray);
}

Vec3 Color2(const Ray &ray, const ShapeList &world) {
  HitRecord rec;
  if (world.Hit(ray, 0.0, MAXFLOAT, rec)) {
    return 0.5 * Vec3(rec.normal.x() + 1.0, rec.normal.y() + 1.0, rec.normal.z() + 1.0);
  }
  return Sky(ray);
}

std::string ShapeListTest(int nx, int ny, Image &image) {

  Vec3 low_left(-2.0, -1.0, -1.0);
  Vec3 horizontal(4.0, 0.0, 0.0);
  Vec3 vertical(0.0, 2.0, 0.0);
  Vec3 origin(0.0, 0.0, 0.0);

  std::vector<std::shared_ptr<Hitable>> hitable_list;
  hitable_list.push_back(std::make_shared<Sphere>(Vec3(0.0, 0.0, -1.0), 0.5f));
  hitable_list.push_back(std::make_shared<Sphere>(Vec3(0.0, -100.5, -1.0), 100.0));
  ShapeList world(hitable_list);
  for (int y = ny-1; y >= 0; --y) {
    for (int x = 0; x < nx; ++x) {
      float u = float(x) / float(nx);
      float v = float(y) / float(ny);
      Ray ray(origin, low_left + u * horizontal + v * vertical);

      Vec3 col = Color2(ray, world);
      image.WritePixel(x, y, col[0], col[1], col[2]);
    }
  }
  return "shapeList.bmp";
}

std::string Antialias(int nx, int ny, int ns, Image &image) {
  std::vector<std::shared_ptr<Hitable>> hitable_list;
  hitable_list.push_back(std::make_shared<Sphere>(Vec3(0.0, 0.0, -1.0), 0.5f));
  hitable_list.push_back(std::make_shared<Sphere>(Vec3(0.0, -100.5, -1.0), 100.0));
  ShapeList world(hitable_list);
  Camera cam;
  for (int y = ny-1; y >= 0; --y) {
    for (int x = 0; x < nx; ++x) {
      Vec3 col(0.0, 0.0, 0.0);
      for (int s = 0; s < ns; ++s) {
        float u = float(x + drand48()) / float(nx);
        float v = float(y + drand48()) / float(ny);
        Ray ray = cam.GetRay(u, v);

        col += Color2(ray, world);
      }
      col /= float(ns);
      image.WritePixel(x, y, col[0], col[1], col[2]);
    }
  }
  return "antialias.bmp";
}

std::string Diffuse(int nx, int ny, int ns, Image &image) {
  std::vector<std::shared_ptr<Hitable>> hitable_list;
  hitable_list.push_back(std::make_shared<Sphere>(Vec3(0.0, 0.0, -1.0), 0.5f));
  hitable_list.push_back(std::make_shared<Sphere>(Vec3(0.0, -100.5, -1.0), 100.0));
  ShapeList world(hitable_list);
  Camera cam;
  for (int y = ny-1; y >= 0; --y) {
    for (int x = 0; x < nx; ++x) {
      Vec3 col(0.0, 0.0, 0.0);
      for (int s = 0; s < ns; ++s) {
        float u = float(x + drand48()) / float(nx);
        float v = float(y + drand48()) / float(ny);
        Ray ray = cam.GetRay(u, v);

        col += Color3(ray, world);
      }
      col /= float(ns);
      col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
      image.WritePixel(x, y, col[0], col[1], col[2]);
    }
  }
  return "diffuse.bmp";
}

std::string MaterialTest(int nx, int ny, int ns, Image &image) {
  ShapeList world;
  world.add(std::make_shared<Sphere>(Vec3(0.0, 0.0, -1.0), 0.5f, std::make_shared<Lambertian>(Vec3(0.8, 0.3, 0.3))));
  world.add(std::make_shared<Sphere>(Vec3(0.0, -100.5, -1.0), 100.0, std::make_shared<Lambertian>(Vec3(0.8, 0.8, 0.0))));
  world.add(std::make_shared<Sphere>(Vec3(1.0, 0.0, -1.0), 0.5f, std::make_shared<Metal>(Vec3(0.8, 0.6, 0.2), 1.0)));
  world.add(std::make_shared<Sphere>(Vec3(-1.0, 0.0, -1.0), 0.5f, std::make_shared<Dielectric>(1.5)));
  world.add(std::make_shared<Sphere>(Vec3(-1.0, 0.0, -1.0), -0.45f, std::make_shared<Dielectric>(1.5)));
  Camera cam;
  for (int y = ny-1; y >= 0; --y) {
    for (int x = 0; x < nx; ++x) {
      Vec3 col(0.0, 0.0, 0.0);
      for (int s = 0; s < ns; ++s) {
        float u = float(x + drand48()) / float(nx);
        float v = float(y + drand48()) / float(ny);
        Ray ray = cam.GetRay(u, v);

        col += Color(ray, world, 0);
      }
      col /= float(ns);
      col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
      image.WritePixel(x, y, col[0], col[1], col[2]);
    }
  }
  return "material.bmp";
}

std::string CameraTest(int nx, int ny, int ns, Image &image) {
  ShapeList world;
  world.add(std::make_shared<Sphere>(Vec3(0.0, 0.0, -1.0), 0.5f, std::make_shared<Lambertian>(Vec3(0.8, 0.3, 0.3))));
  world.add(std::make_shared<Sphere>(Vec3(0.0, -100.5, -1.0), 100.0, std::make_shared<Lambertian>(Vec3(0.8, 0.8, 0.0))));
  world.add(std::make_shared<Sphere>(Vec3(1.0, 0.0, -1.0), 0.5f, std::make_shared<Metal>(Vec3(0.8, 0.6, 0.2), 1.0)));
  world.add(std::make_shared<Sphere>(Vec3(-1.0, 0.0, -1.0), 0.5f, std::make_shared<Dielectric>(1.5)));
  world.add(std::make_shared<Sphere>(Vec3(-1.0, 0.0, -1.0), -0.45f, std::make_shared<Dielectric>(1.5)));

  Vec3 lookfrom(3.0, 3.0, 2.0);
  Vec3 lookat(0.0, 0.0, -1.0);
  float dist_to_focus = (lookfrom - lookat).Length();
  float aperture = 2.0;
  Camera cam(lookfrom, lookat, Vec3(0.0, 1.0, 0.0), 20.0, float(nx)/float(ny), aperture, dist_to_focus);
  for (int y = ny-1; y >= 0; --y) {
    for (int x = 0; x < nx; ++x) {
      Vec3 col(0.0, 0.0, 0.0);
      for (int s = 0; s < ns; ++s) {
        float u = float(x + drand48()) / float(nx);
        float v = float(y + drand48()) / float(ny);
        Ray ray = cam.GetRay(u, v);

        col += Color(ray, world, 0);
      }
      col /= float(ns);
      col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
      image.WritePixel(x, y, col[0], col[1], col[2]);
    }
  }
  return "camera.bmp";
}

ShapeList RandomScene() {
  ShapeList world;
  world.add(std::make_shared<Sphere>(Vec3(0.0, -1000.0, 0.0), 1000.0, std::make_shared<Lambertian>(Vec3(0.5, 0.5, 0.5))));
  for (int a = -11; a < 11; ++a) {
    for (int b = -11; b < 11; ++b) {
      float choose_mat = drand48();
      Vec3 center(a + 0.9 * drand48(), 0.2, b + 0.9 * drand48());
      if ((center - Vec3(4.0, 0.2, 0.0)).Length() > 0.9) {
        if (choose_mat < 0.8) {
          world.add(std::make_shared<Sphere>(center, 
                                             0.2,
                                             std::make_shared<Lambertian>(
                                               Vec3(drand48() * drand48(), drand48() * drand48(), drand48() * drand48()))));
         } else if (choose_mat < 0.95) {
           world.add(std::make_shared<Sphere>(center,
                                             0.2,
                                             std::make_shared<Metal>(
                                               Vec3(0.5 * (1.0 + drand48()), 0.5 * (1.0 + drand48()), 0.5 * (1.0 + drand48())), 0.5 * drand48())));
         } else {
           world.add(std::make_shared<Sphere>(center, 0.2, std::make_shared<Dielectric>(1.5)));
         }
        }
      }
    }
  world.add(std::make_shared<Sphere>(Vec3(0.0, 1.0, 0.0), 1.0, std::make_shared<Dielectric>(1.5)));
  world.add(std::make_shared<Sphere>(Vec3(-4.0, 1.0, 0.0), 1.0, std::make_shared<Lambertian>(Vec3(0.4, 0.2, 0.1))));
  world.add(std::make_shared<Sphere>(Vec3(4.0, 1.0, 0.0), 1.0, std::make_shared<Metal>(Vec3(0.7, 0.6, 0.5), 0.0)));

  return world;
}

std::string RandomTest(int nx, int ny, int ns, Image &image) {
  ShapeList world = RandomScene();

  Vec3 lookfrom(11.0, 2.0, 3.0);
  Vec3 lookat(0.0, 0.0, 0.0);
  float dist_to_focus = (lookfrom - lookat).Length();
  float aperture = 0.1;
  Camera cam(lookfrom, lookat, Vec3(0.0, 1.0, 0.0), 20.0, float(nx)/float(ny), aperture, dist_to_focus);
  for (int y = ny-1; y >= 0; --y) {
    for (int x = 0; x < nx; ++x) {
      Vec3 col(0.0, 0.0, 0.0);
      for (int s = 0; s < ns; ++s) {
        float u = float(x + drand48()) / float(nx);
        float v = float(y + drand48()) / float(ny);
        Ray ray = cam.GetRay(u, v);

        col += Color(ray, world, 0);
      }
      col /= float(ns);
      col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
      image.WritePixel(x, y, col[0], col[1], col[2]);
    }
  }
  return "random.bmp";
}

std::string ParallelTest(int nx, int ny, int ns, Image &image) {
#ifndef _OPENMP
    std::cout<< "OpenMP is not available." << std::endl;
    exit(-1);
#endif
  ShapeList world = RandomScene();

#ifdef _OPENMP   
  int numProcessors = omp_get_num_procs( );
  omp_set_num_threads(numProcessors);
  std::cout << "Have " << numProcessors << " processors.";
  fflush(stdout);
#endif
  Vec3 lookfrom(11.0, 2.0, 3.0);
  Vec3 lookat(0.0, 0.0, 0.0);
  float dist_to_focus = (lookfrom - lookat).Length();
  float aperture = 0.1;
  Camera cam(lookfrom, lookat, Vec3(0.0, 1.0, 0.0), 20.0, float(nx)/float(ny), aperture, dist_to_focus);
#ifdef _OPENMP
  #pragma omp parallel for schedule(static, 1)
#endif
  for (int y = ny-1; y >= 0; --y) {
    for (int x = 0; x < nx; ++x) {
      Vec3 col(0.0, 0.0, 0.0);
      for (int s = 0; s < ns; ++s) {
        float u = float(x + drand48()) / float(nx);
        float v = float(y + drand48()) / float(ny);
        Ray ray = cam.GetRay(u, v);

        col += Color(ray, world, 0);
      }
      col /= float(ns);
      col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
      image.WritePixel(x, y, col[0], col[1], col[2]);
    }
    std::cout << "\r" << ny - y << "/" << ny << "                 ";
    fflush(stdout);
  }
  std::cout << std::endl;
  return "parallel.bmp";
}

int main(int argc, char** args) {
#ifdef _OPENMP
  double start = omp_get_wtime();
#else
  clock_t start = clock();
#endif
  int nx = 200;
  int ny = 100;
  int ns = 100;

  std::cout << "Image(" << nx << ", " << ny << ")" << std::endl;
  Image image = Image(nx, ny);

  std::string filename = ParallelTest(nx, ny, ns, image);

#ifdef _OPENMP
  double end = omp_get_wtime();
  std::cout << "Done(openMP): " << (end - start) * 1000.0 << " msec" << std::endl;
#else
  clock_t end = clock();
  std::cout << "Done: " << (end - start) / float(CLOCKS_PER_SEC) * 1000.0 << " msec" << std::endl;
#endif

  Write(filename, image);

  return 0;
}