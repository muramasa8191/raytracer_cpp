#ifndef RAYTRACER_CPP_SPHERE_HPP_
#define RAYTRACER_CPP_SPHERE_HPP_

#include <iostream>
#include <memory>
#include "hitable.hpp"
#include "material/material.hpp"
#include "ray.hpp"
#include "vector.hpp"

namespace raytracer {
class Sphere : public Hitable {
 public:
  Sphere() {}
  virtual ~Sphere() noexcept = default;
  Sphere(Vec3 center, float radius)
   : center_(center)
   , radius_(radius)
   , mat_() { std::cout << "wrong instructor called" << std::endl; }
  Sphere(Vec3 center, float radius, const std::shared_ptr<Material> &mat): center_(center), radius_(radius), mat_(mat) {}

  virtual bool Hit(const Ray &ray, float t_min, float t_max, HitRecord &rec) const;
  float HitF(const Ray &ray);

 private:
  Vec3 center_;
  float radius_;
  std::shared_ptr<Material> mat_;
};
}

#endif  // RAYTRACER_CPP_SPHERE_HPP_