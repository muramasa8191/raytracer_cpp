#include "image.hpp"
// #include <iostream>

namespace raytracer {
void Image::WritePixel(int x, int y, float r, float g, float b) {
  int idx = (height_ - y - 1) * width_ + x;
  pixels_[idx].r = static_cast<uint8_t>(r * 255.99);
  pixels_[idx].g = static_cast<uint8_t>(g * 255.99);
  pixels_[idx].b = static_cast<uint8_t>(b * 255.99);
}
}  // namespace raytracer