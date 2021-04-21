#include "Util.h"

//std::array<std::array<double, 4>, 4> OrthoProjection(double left, double right,
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
//std::array<std::array<double, 4>, 4> PerspectiveProjection(
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
