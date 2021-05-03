#pragma once

#include <array>
#include <cmath>
#include <iostream>

#include <boost/math/constants/constants.hpp>

#include "SFINAE.h"

namespace UD {
namespace Tensor {
template <class Type, std::size_t Size>
struct IArray {
  virtual ~IArray() {}
  virtual std::array<Type, Size>& data() = 0;

  auto operator[](const std::size_t& index) { return data()[index]; }
  auto at(const std::size_t& index) { return data().at(index); }

  auto begin() { return data().begin(); }
  auto end() { return data().end(); }
  auto rbegin() { return data().rbegin(); }
  auto rend() { return data().rend(); }
  auto cbegin() { return data().cbegin(); }
  auto cend() { return data().cend(); }
  auto crbegin() { return data().crbegin(); }
  auto crend() { return data().crend(); }

  auto front() { return data().front(); }
  auto back() { return data().back(); }
  auto empty() { return data().empty(); }
  auto size() { return data().size(); }
};

template <class Type, std::size_t Size, std::size_t... Sizes>
struct ITensor : public IArray<ITensor<Type, Sizes...>, Size> {
  using Element = ITensor<Type, Sizes...>;
  virtual ~ITensor() {}
};
template <class Type, std::size_t Size>
struct ITensor<Type, Size> : public IArray<Type, Size> {
  using Element = Type;
  virtual ~ITensor() {}
};
template <class Type, std::size_t Size, std::size_t... Sizes>
struct Tensor : public ITensor<Type, Size, Sizes...> {
  using Element = typename ITensor<Type, Size, Sizes...>::Element;
  std::array<Element, Size> _data;
  virtual ~Tensor() {}
  virtual std::array<Element, Size>& data() override { return _data; }
};

template <class Type, std::size_t N>
using IVector = ITensor<Type, N>;
template <class Type, std::size_t N>
using Vector = Tensor<Type, N>;

template <class Type, std::size_t N, std::size_t M>
using IMatrix = ITensor<Type, N, M>;
template <class Type, std::size_t N, std::size_t M>
using Matrix = Tensor<Type, N, M>;

template <class Type, std::size_t Size, std::size_t... Sizes>
ITensor<Type, Size, Sizes...>& operator+=(
    ITensor<Type, Size, Sizes...>& lhs,
    const ITensor<Type, Size, Sizes...>& rhs) {
  for (std::size_t i = 0; i < Size; ++i) lhs[i] += rhs[i];
  return lhs;
}
template <class Type, std::size_t... Sizes>
ITensor<Type, Sizes...>& operator+(ITensor<Type, Sizes...> lhs,
                                   const ITensor<Type, Sizes...>& rhs) {
  return lhs += rhs;
}

template <class Type, std::size_t Size, std::size_t... Sizes>
ITensor<Type, Size, Sizes...>& operator-=(
    ITensor<Type, Size, Sizes...>& lhs,
    const ITensor<Type, Size, Sizes...>& rhs) {
  for (std::size_t i = 0; i < Size; ++i) lhs[i] -= rhs[i];
  return lhs;
}
template <class Type, std::size_t... Sizes>
ITensor<Type, Sizes...>& operator-(ITensor<Type, Sizes...> lhs,
                                   const ITensor<Type, Sizes...>& rhs) {
  return lhs -= rhs;
}

template <class Type, std::size_t Size, std::size_t... Sizes>
ITensor<Type, Size, Sizes...>& operator*=(ITensor<Type, Size, Sizes...>& lhs,
                                          const Type& rhs) {
  for (std::size_t i = 0; i < Size; ++i) lhs[i] *= rhs;
  return lhs;
}
template <class Type, std::size_t... Sizes>
ITensor<Type, Sizes...> operator*(ITensor<Type, Sizes...> lhs,
                                  const Type& rhs) {
  return lhs *= rhs;
}
template <class Type, std::size_t... Sizes>
ITensor<Type, Sizes...> operator*(const Type& lhs,
                                  ITensor<Type, Sizes...> rhs) {
  return rhs *= lhs;
}

template <class Type, std::size_t Size, std::size_t... Sizes>
ITensor<Type, Size, Sizes...>& operator/=(ITensor<Type, Size, Sizes...>& lhs,
                                          const Type& rhs) {
  for (std::size_t i = 0; i < Size; ++i) lhs[i] /= rhs;
  return lhs;
}
template <class Type, std::size_t Size, std::size_t... Sizes>
ITensor<Type, Size, Sizes...> operator/(ITensor<Type, Size, Sizes...> lhs,
                                        const Type& rhs) {
  return lhs /= rhs;
}
template <class Type, std::size_t... Sizes>
ITensor<Type, Sizes...> operator/(const Type& lhs,
                                  ITensor<Type, Sizes...> rhs) {
  return rhs /= lhs;
}

template <class Type, std::size_t... Sizes>
ITensor<Type, Sizes...> operator-(ITensor<Type, Sizes...> a) {
  return a *= -1;
}
template <class Type, std::size_t Size, std::size_t... Sizes>
std::ostream& operator<<(std::ostream& os,
                         const ITensor<Type, Size, Sizes...>& a) {
  if (Size == 0) return os << "[]";
  os << '[';
  for (std::size_t i = 0; i < Size - 1; ++i) os << a[i] << ',';
  return os << a[Size - 1] << ']';
}

template <class Type, std::size_t... Sizes>
constexpr std::size_t Order(const ITensor<Type, Sizes...>& a) {
  return sizeof...(Sizes);
}
template <class Type, std::size_t S, std::size_t... Ss>
constexpr std::size_t Size(const ITensor<Type, S, Ss...>& a) {
  return S;
}
template <std::size_t Start, std::size_t N, class Type, std::size_t Size,
          typename std::enable_if_t<(((Start + N) <= Size), N > 0)>* = nullptr>
IArray<Type, N> Sub(const IArray<Type, Size>& a) {
  IArray<Type, N> ret{};
  for (std::size_t i = 0; i < N; ++i) ret[i] = a[Start + i];
  return ret;
}
template <std::size_t Start, class Type, std::size_t Size>
IArray<Type, Size - Start> Sub(const IArray<Type, Size>& a) {
  return Sub<Start, Size - Start, Type, Size>(a);
}

// Merge case: lhs convertible to rhs || not rhs convertible to lhs
template <class Type, std::size_t Size1, std::size_t Size2,
          std::size_t... Sizes1, std::size_t... Sizes2,
          typename std::enable_if_t<
              (std::is_convertible_v<ITensor<Type, Sizes1...>,
                                     ITensor<Type, Sizes2...>> ||
               !std::is_convertible_v<ITensor<Type, Sizes2...>,
                                      ITensor<Type, Sizes1...>>)>* = nullptr>
ITensor<Type, Size1 + Size2, Sizes2...> Merge(
    const ITensor<Type, Size1, Sizes1...>& lhs,
    const ITensor<Type, Size2, Sizes2...>& rhs) {
  ITensor<Type, Size1 + Size2, Sizes2...> ret{};
  ret = lhs;
  for (std::size_t i = 0; i < Size2; ++i) ret[Size1 + i] = rhs[i];
  return ret;
}
// Merge case: not lhs convertible to rhs && rhs convertible to lhs
template <class Type, std::size_t Size1, std::size_t Size2,
          std::size_t... Sizes1, std::size_t... Sizes2,
          typename std::enable_if_t<
              (!std::is_convertible_v<ITensor<Type, Sizes1...>,
                                      ITensor<Type, Sizes2...>> &&
               std::is_convertible_v<ITensor<Type, Sizes2...>,
                                     ITensor<Type, Sizes1...>>)>* = nullptr>
ITensor<Type, Size1 + Size2, Sizes1...> Merge(
    const ITensor<Type, Size1, Sizes1...>& lhs,
    const ITensor<Type, Size2, Sizes2...>& rhs) {
  ITensor<Type, Size1 + Size2, Sizes1...> ret{};
  ret = lhs;
  for (std::size_t i = 0; i < Size2; ++i) ret[Size1 + i] = rhs[i];
  return ret;
}
template <class Type, std::size_t Size
          // , typename std::enable_if_t
          // <(
          //     Index>0
          //     && Index < Size-1
          // )>* = nullptr
          >
ITensor<Type, Size - 1> Minor(const ITensor<Type, Size>& a,
                              const std::size_t& index) {
  ITensor<Type, Size - 1> ret{};
  for (std::size_t i = 0; i < index; ++i) ret[i] = a[i];
  for (std::size_t i = index + 1; i < Size; ++i) ret[i - 1] = a[i];
  return ret;
}
template <class... Indices, class Type, std::size_t Size, std::size_t... Sizes,
          typename std::enable_if_t<
              (std::conjunction_v<std::is_same<std::size_t, Indices>...> &&
               sizeof...(Indices) == sizeof...(Sizes))>* = nullptr>
ITensor<Type, Size - 1, (Sizes - 1)...> Minor(
    const ITensor<Type, Size, Sizes...>& a, const std::size_t& index,
    Indices... indices) {
  ITensor<Type, Size - 1, (Sizes - 1)...> ret{};
  for (std::size_t i = 0; i < index; ++i) ret[i] = Minor(a[i], indices...);
  for (std::size_t i = index + 1; i < Size; ++i)
    ret[i - 1] = Minor(a[i], indices...);
  return ret;
}
template <class Type, std::size_t N>
ITensor<Type, N> Transpose(const ITensor<Type, N>& a) {
  return a;
}
template <class Type, std::size_t N, std::size_t M>
ITensor<Type, M, N> Transpose(const ITensor<Type, N, M>& a) {
  ITensor<Type, M, N> ret;
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < M; ++j) {
      ret[j][i] = a[i][j];
    }
  }
  return ret;
}
template <class Type, std::size_t N>
Type Trace(const ITensor<Type, N, N>& matrix) {
  Type ret{0};
  for (std::size_t i = 0; i < N; ++i) {
    ret += matrix[i][i];
  }
  return ret;
}
template <class Type, std::size_t N>
Type RevTrace(const ITensor<Type, N, N>& matrix) {
  Type ret{0};
  for (std::size_t i = 0; i < N; ++i) {
    ret += matrix[i][N - i - 1];
  }
  return ret;
}
template <class Type, std::size_t N, std::size_t M, std::size_t U>
IMatrix<Type, N, U> operator*(const IMatrix<Type, N, M>& lhs,
                              const IMatrix<Type, M, U>& rhs) {
  IMatrix<Type, N, U> ret{};
  IMatrix<Type, U, M> trhs{Transpose(rhs)};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < U; ++j) {
      ret[i][j] = DotProduct(lhs[i], trhs[j]);
    }
  }
  return ret;
}
// Matrix
// Generalized Outer Product
template <class Type, std::size_t N, std::size_t M>
IMatrix<Type, N, M> TensorProduct(const IVector<Type, N>& lhs,
                                  const IVector<Type, M> rhs) {
  return Transpose(lhs) * IMatrix<Type, 1, M>(rhs);
}
// Exterior/Wedge Product - generalized Cross Product
template <class Type, std::size_t N, std::size_t M>
IMatrix<Type, N, M> ExteriorProduct(const IVector<Type, N>& lhs,
                                    const IVector<Type, M> rhs) {
  return TensorProduct(lhs, rhs) - TensorProduct(rhs, lhs);
}
// Inner Product - generalized Dot Product
template <class Type, std::size_t N>
Type InnerProduct(const IVector<Type, N>& lhs, const IVector<Type, N> rhs) {
  return Trace(TensorProduct(lhs, rhs));
}
template <class Type, std::size_t N>
Type Determinant(const IMatrix<Type, N, N>& matrix) {
  Type ret{};
  std::array<IMatrix<Type, N - 1, N - 1>, N> minors{Minors(matrix, 0)};
  for (std::size_t i = 0; i < N; ++i)
    ret += Determinant(minors[i]) * matrix[0][i] * (i % 2 == 1 ? -1 : 1);
  return ret;
}
template <class Type>
Type Determinant(const IMatrix<Type, 2, 2>& matrix) {
  return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
}
template <class Type>
Type Determinant(const IMatrix<Type, 1, 1>& matrix) {
  return matrix[0][0];
}

template <class Type, std::size_t N, std::size_t M>
std::array<std::array<IMatrix<Type, N - 1, M - 1>, M>, N> Minors(
    const IMatrix<Type, N, M>& matrix) {
  std::array<std::array<IMatrix<Type, N - 1, M - 1>, M>, N> ret;
  for (std::size_t i = 0; i < N; ++i) ret[i] = Minors(matrix, i);
  return ret;
}
template <class Type, std::size_t N, std::size_t M>
std::array<IMatrix<Type, N - 1, M - 1>, M> Minors(
    const IMatrix<Type, N, M>& matrix, const std::size_t& index) {
  std::array<IMatrix<Type, N - 1, M - 1>, M> ret;
  for (std::size_t i = 0; i < M; ++i) ret[i] = Minor(matrix, index, i);
  return ret;
}
template <class Type, std::size_t N>
IVector<Type, N> OrthogonalVector(const IMatrix<Type, N - 1, N>& matrix) {
  std::array<IMatrix<Type, N - 1, N - 1>, N> arr{
      Minors(static_cast<IMatrix<Type, N, N>>(matrix), N - 1)};
  IVector<Type, N> ret;
  for (std::size_t i = 0; i < N; ++i) {
    ret[i] = Determinant(arr[i]);
    if (N % 2 == 0) {
      ret[i] *= (i % 2 == 1 ? 1 : -1);
    }
  }
  return ret;
}
// Vector
template <class Type, std::size_t N>
Type DotProduct(const IVector<Type, N>& lhs, const IVector<Type, N>& rhs) {
  Type ret{0};
  for (std::size_t i = 0; i < N; ++i) {
    ret += lhs[i] * rhs[i];
  }
  return ret;
}
template <class Type, std::size_t N>
IMatrix<Type, N, 1> Transpose(const IVector<Type, N>& a) {
  return Transpose(IMatrix<Type, 1, N>(a));
}
template <class Type, std::size_t N>
IVector<Type, N> Normalize(const IVector<Type, N>& a) {
  Type sum{0};
  for (std::size_t i = 0; i < N; ++i) {
    sum += std::pow(a[i], 2);
  }
  return a / std::sqrt(sum);
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
IVector<Type, N> CrossProduct(const IVector<Type, N>& lhs,
                              const IVector<Type, N>& rhs) {
  IMatrix<Type, N, N> m{ExteriorProduct(lhs, rhs)};
  return IVector<Type, N>{m[1][2], m[2][0], m[0][1]};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
IMatrix<Type, N + 1, N + 1> TranslationMatrix(Type Tx, Type Ty, Type Tz) {
  return IMatrix<Type, N + 1, N + 1>{
      {1, 0, 0, Tx}, {0, 1, 0, Ty}, {0, 0, 1, Tz}, {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
IMatrix<Type, N + 1, N + 1> RotationMatrix(Type Rx, Type Ry, Type Rz) {
  return IMatrix<Type, N + 1, N + 1>{{1, 0, 0, 0},
                                     {0, std::cos(Rx), -std::sin(Rx), 0},
                                     {0, std::sin(Rx), std::cos(Rx), 0},
                                     {0, 0, 0, 1}} *
         IMatrix<Type, N + 1, N + 1>{{std::cos(Ry), 0, std::sin(Ry), 0},
                                     {0, 1, 0, 0},
                                     {-std::sin(Ry), 0, std::cos(Ry), 0},
                                     {0, 0, 0, 1}} *
         IMatrix<Type, N + 1, N + 1>{{std::cos(Rz), -std::sin(Rz), 0, 0},
                                     {std::sin(Rz), std::cos(Rz), 0, 0},
                                     {0, 0, 1, 0},
                                     {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
IMatrix<Type, N + 1, N + 1> ScaleMatrix(Type Sx, Type Sy, Type Sz) {
  return IMatrix<Type, N + 1, N + 1>{
      {Sx, 0, 0, 0}, {0, Sy, 0, 0}, {0, 0, Sz, 0}, {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
IMatrix<Type, N + 1, N + 1> ModelMatrix(Type Tx, Type Ty, Type Tz, Type Rx,
                                        Type Ry, Type Rz, Type Sx, Type Sy,
                                        Type Sz) {
  return TranslationMatrix<N>() * RotationMatrix<N>() * ScaleMatrix<N>();
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
IMatrix<Type, N + 1, N + 1> ViewMatrix(IVector<Type, N> position,
                                       IVector<Type, N> target,
                                       IVector<Type, N> up) {
  IVector<Type, 3> view_direction{Normalize(position - target)};
  IVector<Type, 3> view_right{Normalize(CrossProduct(up, view_direction))};
  IVector<Type, 3> view_up{CrossProduct(view_direction, view_right)};

  return IMatrix<Type, N + 1, N + 1>{
             {view_right[0], view_right[1], view_right[2], 0},
             {view_up[0], view_up[1], view_up[2], 0},
             {view_direction[0], view_direction[1], view_direction[2], 0},
             {0, 0, 0, 1}} *
         IMatrix<Type, N + 1, N + 1>{{1, 0, 0, -position[0]},
                                     {0, 1, 0, -position[1]},
                                     {0, 0, 1, -position[2]},
                                     {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
IMatrix<Type, N + 1, N + 1> OrthographicMatrix(Type near_val, Type far_val,
                                               Type aspect, Type fov) {
  Type top = near_val *
             std::tan((boost::math::constants::pi<double>() / 180) * fov / 2);
  Type bottom = -top;
  Type right = top * aspect;
  Type left = -right;

  return IMatrix<Type, N + 1, N + 1>{
      {1 / (right - left), 0, 0, -((right + left) / (right - left))},
      {0, 2 / (top - bottom), 0, -((top + bottom) / (top - bottom))},
      {0, 0, -(2 / (far_val - near_val)),
       -((far_val + near_val) / (far_val - near_val))},
      {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
IMatrix<Type, N + 1, N + 1> PerspectiveMatrix(Type near_val, Type far_val,
                                              Type aspect, Type fov) {
  Type top = near_val *
             std::tan((boost::math::constants::pi<double>() / 180) * fov / 2);
  Type bottom = -top;
  Type right = top * aspect;
  Type left = -right;

  return IMatrix<Type, N + 1, N + 1>{
      {(2 * near_val) / (right - left), 0, (right + left) / (right - left), 0},
      {0, (2 * near_val) / (top - bottom), (top + bottom) / (top - bottom), 0},
      {0, 0, -((far_val + near_val) / (far_val - near_val)),
       -((2 * far_val * near_val) / (far_val - near_val))},
      {0, 0, -1, 0}};
}
template <class Type, std::size_t N>
IMatrix<Type, N, N> IdentityMatrix() {
  IMatrix<Type, N, N> ret{};
  for (std::size_t i = 0; i < N; ++i) {
    ret[i][i] = 1;
  }
  return ret;
}
}  // namespace Tensor

}  // namespace UD