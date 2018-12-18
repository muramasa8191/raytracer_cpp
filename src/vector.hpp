#ifndef RAYTRACER_CPP_VECTOR_HPP_
#define RAYTRACER_CPP_VECTOR_HPP_

#include "math.h"

namespace raytracer {
class Vec3 {
 public:
  Vec3():val() {}
  Vec3(float x, float y, float z) { val[0] = x, val[1] = y, val[2] = z;}
  inline float x() const { return val[0]; }
  inline float y() const { return val[1]; }
  inline float z() const { return val[2]; }
  inline float r() const { return val[0]; }
  inline float g() const { return val[1]; }
  inline float b() const { return val[2]; }

  inline const Vec3 & operator+() const { return *this; }
  inline Vec3 operator-() const { return Vec3(-val[0], -val[1], -val[2]);}
  inline float operator[] (int i) const { return val[i]; }
  inline float & operator[] (int i) { return val[i]; }

  inline Vec3 & operator+=(const Vec3 &v) {
    val[0] += v[0];
    val[1] += v[1];
    val[2] += v[2];
    
    return *this;
  }

  inline Vec3 & operator-=(const Vec3 &v) {
    val[0] -= v[0];
    val[1] -= v[1];
    val[2] -= v[2];

    return *this;
  }

  inline Vec3 & operator*=(const Vec3 &v) {
    val[0] *= v[0];
    val[1] *= v[1];
    val[2] *= v[2];

    return *this;
  }

  inline Vec3 & operator/=(const Vec3 &v) {
    val[0] /= v[0];
    val[1] /= v[1];
    val[2] /= v[2];

    return *this;
  }

  inline Vec3 & operator*=(const float s) {
    val[0] *= s;
    val[1] *= s;
    val[2] *= s;

    return *this;
  }

  inline Vec3 & operator/=(const float s) {
    val[0] /= s;
    val[1] /= s;
    val[2] /= s;

    return *this;
  }

  inline void MakeUnitVector() {
    float k = 1.0 / Length();
    val[0] *= k;
    val[1] *= k;
    val[2] *= k;
  }

  inline float Length() const { return sqrt(val[0] * val[0] + val[1] * val[1] + val[2] * val[2]); }
  inline float SquaredLength() const { return val[0] * val[0] + val[1] * val[1] + val[2] * val[2];}

 private:
  float val[3];
};

inline Vec3 operator+(const Vec3 &v1, const Vec3 &v2) {
  return Vec3(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
}
inline Vec3 operator-(const Vec3 &v1, const Vec3 &v2) {
  return Vec3(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);
}
inline Vec3 operator*(const Vec3 &v1, const Vec3 &v2) {
  return Vec3(v1[0] * v2[0], v1[1] * v2[1], v1[2] * v2[2]);
}
inline Vec3 operator/(const Vec3 &v1, const Vec3 &v2) {
  return Vec3(v1[0] / v2[0], v1[1] / v2[1], v1[2] / v2[2]);
}
inline Vec3 operator*(float s, const Vec3 &v) {
  return Vec3(s * v[0], s * v[1], s* v[2]);
}
inline Vec3 operator/(const Vec3 &v, float s) {
  return Vec3(v[0] / s, v[1] / s, v[2] / s);
}
inline Vec3 operator*(const Vec3 &v, float s) {
  return Vec3(v[0] * s, v[1] * s, v[2] * s);
}
inline float Dot(const Vec3 &v1, const Vec3 &v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
inline Vec3 Cross(const Vec3 &v1, const Vec3 &v2) {
  return Vec3((v1[1] * v2[2] - v1[2] * v2[1]),
              (-(v1[0] * v2[2] - v1[2] * v2[0])),
              (v1[0] * v2[1] - v1[1] * v2[0]));
}
inline Vec3 UnitVector(Vec3 v) {
  return v / v.Length();
}
}  // namespace raytracer
#endif  // RAYTRACER_CPP_VECTOR_HPP_