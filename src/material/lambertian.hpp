#ifndef RAYTRACER_CPP_LAMBERTIAN_HPP_
#define RAYTRACER_CPP_LAMBERTIAN_HPP_

#include <iostream>
#include "material.hpp"
#include "shape/hitable.hpp"
#include "utils.hpp"
#include "vector.hpp"

namespace raytracer {
class Lambertian : public Material {
 public:
  Lambertian(const Vec3 &albedo) : albedo_(albedo) {}
  virtual ~Lambertian() noexcept = default;
  virtual bool scatter(const Ray &ray, const HitRecord &rec, Vec3 &attenuation, Ray &scattered) const override {
    Vec3 target = rec.p + rec.normal + RandomInUnitSphere();
    scattered = Ray(rec.p, target - rec.p);
    attenuation = albedo_;
    return true;
  }
 private:
  Vec3 albedo_;
};
}
#endif  // RAYTRACER_CPP_LAMBERTIAN_HPP_