#include "vector.hpp"

namespace raytracer {
inline Vec3 & Vec3::operator+=(const Vec3 &v) {
  val[0] += v[0];
  val[1] += v[1];
  val[2] += v[2];
  
  return *this;
}

inline Vec3 & Vec3::operator-=(const Vec3 &v) {
  val[0] -= v[0];
  val[1] -= v[1];
  val[2] -= v[2];

  return *this;
}

inline Vec3 & Vec3::operator*=(const Vec3 &v) {
  val[0] *= v[0];
  val[1] *= v[1];
  val[2] *= v[2];

  return *this;
}

inline Vec3 & Vec3::operator/=(const Vec3 &v) {
  val[0] /= v[0];
  val[1] /= v[1];
  val[2] /= v[2];

  return *this;
}

inline Vec3 & Vec3::operator*=(const float s) {
  val[0] *= s;
  val[1] *= s;
  val[2] *= s;

  return *this;
}

inline Vec3 & Vec3::operator/=(const float s) {
  val[0] /= s;
  val[1] /= s;
  val[2] /= s;

  return *this;
}

inline void Vec3::MakeUnitVector() {
  float k = 1.0 / Length();
  val[0] *= k;
  val[1] *= k;
  val[2] *= k;
}
}