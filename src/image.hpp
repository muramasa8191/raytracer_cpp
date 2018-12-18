#ifndef RAYTRACER_CPP_IMAGE_HPP_
#define RAYTRACER_CPP_IMAGE_HPP_

#include <memory>
#include <string>
#include <vector>

namespace raytracer {
struct Color3D {
  uint8_t r;
  uint8_t g;
  uint8_t b;
};

class Image {
 public:
  Image() {};
  Image(int w, int h): width_(w), height_(h) {
    pixels_.reset(new Color3D[w * h]);
  }
  int GetWidth() const { return width_; }
  int GetHeight() const { return height_; }
  void* GetPixels() const { return pixels_.get(); }
  void WritePixel(int x, int y, float r, float g, float b);

 private:
  int width_;
  int height_;
  std::unique_ptr<Color3D[]> pixels_;
};
}
#endif  // RAYTRACER_CPP_IMAGE_HPP_