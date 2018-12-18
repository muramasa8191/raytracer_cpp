#ifndef RAYTRACER_CPP_RAY_HPP_
#define RAYTRACER_CPP_RAY_HPP_

#include "vector.hpp"

namespace raytracer {
class Ray {
 public:
  Ray() {}
  Ray(const Vec3 &ori, const Vec3 &dir): origin_(ori), direction_(dir) {}
  Vec3 GetOrigin() const { return origin_; }
  Vec3 GetDirection() const { return direction_; }
  Vec3 PointAt(float t) const { return origin_ + t * direction_; }
 private:
  Vec3 origin_;
  Vec3 direction_;
};
}
#endif  // RAYTRACER_CPP_RAY_HPP_