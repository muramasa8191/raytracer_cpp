#ifndef RAYTRACER_CPP_HITABLE_HPP_
#define RAYTRACER_CPP_HITABLE_HPP_

#include <memory>

#include "ray.hpp"
#include "vector.hpp"

namespace raytracer {
class Material;
struct HitRecord {
  float t;
  Vec3 p;
  Vec3 normal;
  std::shared_ptr<Material> mat;
};

class Hitable {
 public:
  virtual bool Hit(const Ray &ray, float t_min, float t_max, HitRecord &rec) const = 0;
};
}  // namespace raytracer
#endif  // RAYTRACER_CPP_HITABLE_HPP_