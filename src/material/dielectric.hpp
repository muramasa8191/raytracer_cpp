#ifndef RAYTRACER_CPP_DIELECTRIC_HPP_
#define RAYTRACER_CPP_DIELECTRIC_HPP_

#include <random>

#include "../shape/hitable.hpp"
#include "material.hpp"
#include "ray.hpp"
#include "utils.hpp"
#include "vector.hpp"

namespace raytracer {
namespace {
float schlick(float cosine, float ref_idx) {
  float r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
  r0 = r0 * r0;
  return r0 + (1.0 - r0) * pow(1.0 - cosine, 5);
}
}  // namespace

class Dielectric : public Material {
 public:
  Dielectric(float ri) : ref_idx_(ri) {}
  virtual ~Dielectric() noexcept = default;
  virtual bool scatter(const Ray &ray, const HitRecord &rec, Vec3 &attenuation, Ray &scattered) const override {
    Vec3 outward_normal;
    Vec3 reflected = reflect(ray.GetDirection(), rec.normal);
    float ni_over_nt;
    attenuation = Vec3(1.0, 1.0, 1.0);
    Vec3 refracted;
    float reflect_prob;
    float cosine;
    if (Dot(ray.GetDirection(), rec.normal) > 0.0) {
      outward_normal = -rec.normal;
      ni_over_nt = ref_idx_;
      cosine = ref_idx_ * Dot(ray.GetDirection(), rec.normal) / ray.GetDirection().Length();
    } else {
      outward_normal = rec.normal;
      ni_over_nt = 1.0 / ref_idx_;
      cosine = -Dot(ray.GetDirection(), rec.normal) / ray.GetDirection().Length();
    }
    if (refract(ray.GetDirection(), outward_normal, ni_over_nt, refracted)) {
      reflect_prob = schlick(cosine, ref_idx_);
    } else {
      reflect_prob = 1.0;
    }
    if (drand48() < reflect_prob) {
      scattered = Ray(rec.p, reflected);
    } else {
      scattered = Ray(rec.p, refracted);
    }
    return true;
  }

 private:
  float ref_idx_;
};
}
#endif