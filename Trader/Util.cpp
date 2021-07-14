#include "Util.h"

auto CombinationsIndexHelper(std::size_t sample_size, std::size_t length, std::size_t skip) -> std::vector<std::vector<std::size_t>> {
  std::vector<std::vector<std::size_t>> ret;
  if (length == 1) {
    for (std::size_t i = skip; i < sample_size; ++i) {
      std::vector<std::size_t> vec;
      vec.push_back(i);
      ret.push_back(vec);
    }
    return ret;
  }
  for (std::size_t i = skip; i + length <= sample_size; ++i) {
    auto combs = CombinationsIndexHelper(sample_size, length - 1, i + 1);
    for (auto& c : combs) {
      std::vector<std::size_t> vec;
      vec.push_back(i);
      vec.insert(vec.end(), c.begin(), c.end());
      ret.push_back(vec);
    }
  }
  return ret;
}
auto CombinationsIndex(std::size_t sample_size, std::size_t length) -> std::vector<std::vector<std::size_t>> {
  return CombinationsIndexHelper(sample_size, length, 0);
}

// std::array<std::array<double, 4>, 4> OrthoProjection(double left, double right,
//                                                     double bottom, double top,
//                                                     double far, double near) {
//  return {2 / (right - left),
//          0,
//          0,
//          -((right + left) / (right - left)),
//          0,
//          2 / (top - bottom),
//          0,
//          -((top + bottom) / (top - bottom)),
//          0,
//          0,
//          -2 / (far - near),
//          -((far + near) / (far - near)),
//          0,
//          0,
//          0,
//          1};
//}
//
// std::array<std::array<double, 4>, 4> PerspectiveProjection(
//    double left, double right, double bottom, double top, double near,
//    double far) {
//  return {2 * near / (right - left),
//          0,
//          (right + left) / (right - left),
//          0,
//          0,
//          2 * near / (top - bottom),
//          (top + bottom) / (top - bottom),
//          0,
//          0,
//          0,
//          -(far + near) / (far - near),
//          -2 * far * near / (far - near),
//          0,
//          0,
//          -1,
//          0};
//  // return {2 * near / (right - left),
//  //        0,
//  //        0,
//  //        0,
//  //        0,
//  //        2 * near / (top - bottom),
//  //        0,
//  //        0,
//  //        (right + left) / (right - left),
//  //        (top + bottom) / (top - bottom),
//  //        -(far + near) / (far - near),
//  //        -1,
//  //        0,
//  //        0,
//  //        -2 * far * near / (far - near),
//  //        0};
//}
