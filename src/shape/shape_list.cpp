#include "shape_list.hpp"

#include <memory>
#include <vector>

#include "hitable.hpp"

namespace raytracer {
bool ShapeList::Hit(const Ray &ray, float t_min, float t_max, HitRecord &rec) const {
  HitRecord tmp_rec;
  bool is_hit = false;
  double closest = t_max;
  for (std::shared_ptr<Hitable> hitable : list_) {
    if (hitable->Hit(ray, t_min, closest, tmp_rec)) {
      is_hit = true;
      closest = tmp_rec.t;
      rec = tmp_rec;
    }
  }
  return is_hit;
}
}