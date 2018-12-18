#ifndef RAYTRACER_CPP_SHAPE_LIST_HPP_
#define RAYTRACER_CPP_SHAPE_LIST_HPP_

#include "hitable.hpp"
#include "ray.hpp"

#include <memory>
#include <vector>

namespace raytracer {
class ShapeList : public Hitable {
 public:
  ShapeList() {}
  explicit ShapeList(std::vector<std::shared_ptr<Hitable>> list) : list_(list) {};
  virtual bool Hit(const Ray &ray, float t_min, float t_max, HitRecord &rec) const;
  void add(const std::shared_ptr<Hitable> &hitable_ptr) {list_.push_back(hitable_ptr);}
 private:
  std::vector<std::shared_ptr<Hitable>> list_;
};
}
#endif  // RAYTRACER_CPP_SHAPE_LIST_HPP_