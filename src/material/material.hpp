#ifndef RAYTRACER_CPP_MATERIAL_HPP_
#define RAYTRACER_CPP_MATERIAL_HPP_
#include <memory>

#include "shape/hitable.hpp"
#include "ray.hpp"
#include "vector.hpp"

namespace raytracer {

class Material {
 public:
  virtual bool scatter(const Ray &ray, const HitRecord &rec, Vec3 &attenuation, Ray &scattered) const = 0;
};
}
#endif  // RAYTRACER_CPP_MATERIAL_HPP_