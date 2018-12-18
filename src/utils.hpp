#ifndef RAYTRACER_CPP_UTILS_HPP_
#define RAYTRACER_CPP_UTILS_HPP_

#include <random>
#include "vector.hpp"

namespace raytracer {
Vec3 RandomInUnitSphere() {
  Vec3 p;
  do {
    p = 2.0 * Vec3(drand48(), drand48(), drand48()) - Vec3(1.0, 1.0, 1.0);
  } while (p.SquaredLength() >= 1.0);
  return p;
}

Vec3 reflect(const Vec3 &v, const Vec3 &n) {
  return v - 2 * Dot(v, n) * n;
}

bool refract(const Vec3 &v, const Vec3 &n, float ni_over_nt, Vec3 &refracted) {
  Vec3 uv = UnitVector(v);
  float dt = Dot(uv, n);
  float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
  if (discriminant > 0.0) {
    refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
    return true;
  }
  return false;
}
}
#endif  // RAYTRACER_CPP_UTILS_HPP_