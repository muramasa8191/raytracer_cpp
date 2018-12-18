#ifndef RAYTRACER_CPP_CAMERA_HPP_
#define RAYTRACER_CPP_CAMERA_HPP_

#include <random>

#include "ray.hpp"
#include "vector.hpp"

namespace raytracer{
namespace {
Vec3 RandomInUnitDisk() {
  Vec3 p;
  do {
    p = 2.0 * Vec3(drand48(), drand48(), 0) - Vec3(1.0, 1.0, 0.0);
  } while (Dot(p, p) >= 1.0);

  return p;
}
}
class Camera {
 public:
  Camera() {
    low_left_ = Vec3(-2.0, -1.0, -1.0);
    horizontal_ = Vec3(4.0, 0.0, 0.0);
    vertical_ = Vec3(0.0, 2.0, 0.0);
    origin_ = Vec3(0.0, 0.0, 0.0);
  }

  Camera(Vec3 lookfrom, Vec3 lookat, Vec3 vup, float vfov, float aspect, float aperture, float focus_dist) {
    lens_radius_ = aperture * 0.5;
    float theta = vfov * M_PI / 180.0;
    float half_height = tan(theta * 0.5);
    float half_width = aspect * half_height;
    origin_ = lookfrom;
    w_ = UnitVector(lookfrom - lookat);
    u_ = UnitVector(Cross(vup, w_));
    v_ = Cross(w_, u_);
    low_left_ = origin_ - half_width * focus_dist * u_ -half_height * focus_dist * v_ - focus_dist * w_;
    horizontal_ = 2.0 * half_width * focus_dist * u_;
    vertical_ = 2.0 * half_height * focus_dist * v_;
  }

  Ray GetRay(float s, float t) {
    Vec3 rd = lens_radius_ * RandomInUnitDisk();
    Vec3 offset = u_ * rd.x() + v_ * rd.y();
    return Ray(origin_ + offset, low_left_ + s * horizontal_ + t * vertical_ - origin_ - offset);
  }

 private:
  Vec3 low_left_;
  Vec3 horizontal_;
  Vec3 vertical_;
  Vec3 origin_;
  Vec3 u_, v_, w_;
  float lens_radius_;
};
}
#endif