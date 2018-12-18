#ifndef RAYTRACER_CPP_METAL_HPP_
#define RAYTRACER_CPP_METAL_HPP_

#include <iostream>
#include "material.hpp"
#include "ray.hpp"
#include "utils.hpp"
#include "vector.hpp"

namespace raytracer {
class Metal : public Material {
 public:
  Metal(const Vec3 &albedo, float fuzz) : albedo_(albedo), fuzz_(fuzz) {}
  virtual ~Metal() = default;
  virtual bool scatter(const Ray &ray, const HitRecord &rec, Vec3 &attenuation, Ray &scattered) const override {
    Vec3 reflected = reflect(UnitVector(ray.GetDirection()), rec.normal);
    scattered = Ray(rec.p, reflected + fuzz_ * RandomInUnitSphere());
    attenuation = albedo_;
    return (Dot(scattered.GetDirection(), rec.normal) > 0.0);
  }

 private:
  Vec3 albedo_;
  float fuzz_;
};
}
#endif  // RAYTRACER_CPP_METAL_HPP_