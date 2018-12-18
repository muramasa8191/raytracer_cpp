#include "sphere.hpp"
#include "hitable.hpp"

#include "ray.hpp"
#include "vector.hpp"

namespace raytracer {
bool Sphere::Hit(const Ray &ray, float t_min, float t_max, HitRecord &rec) const {
  Vec3 oc = ray.GetOrigin() - center_;

  float a = Dot(ray.GetDirection(), ray.GetDirection());
  float b = Dot(oc, ray.GetDirection());
  float c = Dot(oc, oc) - radius_ * radius_;
  float discriminant = b * b - a * c;

  if (discriminant > 0.0f) {
    float tmp = (-b - sqrt(b * b - a * c)) / a;
    if (t_min < tmp && tmp < t_max) {
      rec.t = tmp;
      rec.p = ray.PointAt(tmp);
      rec.normal = (rec.p - center_) / radius_;
      rec.mat = mat_;
      return true;
    } 
    tmp = (-b + sqrt(b * b - a * c)) / a;
    if (t_min < tmp && tmp < t_max) {
      rec.t = tmp;
      rec.p = ray.PointAt(tmp);
      rec.normal = (rec.p - center_) / radius_;
      rec.mat = mat_;
      return true;
    }
  }
  return false;
}

float Sphere::HitF(const Ray &ray) {
  Vec3 oc = ray.GetOrigin() - center_;

  float a = Dot(ray.GetDirection(), ray.GetDirection());
  float b = 2.0 * Dot(oc, ray.GetDirection());
  float c = Dot(oc, oc) - radius_ * radius_;
  float discriminant = b * b - 4 * a * c;

  if (discriminant < 0.0f) {
    return -1.0;
  } else {
    return (-b - sqrt(discriminant)) / (2.0 * a);
  }
}
}