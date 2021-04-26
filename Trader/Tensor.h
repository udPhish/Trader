#pragma once

#include <array>
#include <iostream>
#include <cmath>

#include "SFINAE.h"

namespace UD {
template <class Type, std::size_t... Sizes>
struct ten;
template <class Type, std::size_t Size>
struct ten<Type, Size> {
  using Element = Type;
  std::array<Element, Size> data;
  Element& operator[](const std::size_t& index) { return data[index]; }
  const Element& operator[](const std::size_t& index) const {
    return data[index];
  }

  // Implicit lossless conversions:
  // Convert Elements
  template <std::size_t S, typename std::enable_if_t<(S >= Size)>* = nullptr>
  operator ten<Type, S>() const {
    ten<Type, S> ret{};
    for (std::size_t i = 0; i < Size; ++i) {
      ret[i] = (*this)[i];
    }
    return ret;
  }
  // Strip outer arrays from return
  template <std::size_t S, std::size_t... Ss,
            typename std::enable_if_t<(sizeof...(Ss) > 0)>* =
                nullptr  // && std::is_convertible_v<ten<Type, Size>, ten<Type,
                         // Ss...>>)>* = nullptr
            >
  operator ten<Type, S, Ss...>() const {
    ten<Type, S, Ss...> ret{};
    ret[0] = *this;
    return ret;
  }
  operator std::array<Type, Size>() const { return data; }
};
template <class Type, std::size_t Size, std::size_t... Sizes>
struct ten<Type, Size, Sizes...> {
  using Element = ten<Type, Sizes...>;
  std::array<Element, Size> data;
  Element& operator[](const std::size_t& index) { return data[index]; }
  const Element& operator[](const std::size_t& index) const {
    return data[index];
  }

  // Implicit lossless conversions:
  // Convert Elements
  template <std::size_t S, std::size_t... Ss,
            typename std::enable_if_t<
                (S >= Size && sizeof...(Ss) == sizeof...(Sizes))>* = nullptr>
  operator ten<Type, S, Ss...>() const {
    ten<Type, S, Ss...> ret{};
    for (std::size_t i = 0; i < Size; ++i) {
      ret[i] = (*this)[i];
    }
    return ret;
  }
  // Strip outer arrays from return
  template <
      std::size_t S, std::size_t... Ss,
      typename std::enable_if_t<(sizeof...(Ss) > sizeof...(Sizes))>* = nullptr>
  operator ten<Type, S, Ss...>() const {
    ten<Type, S, Ss...> ret{};
    ret[0] = *this;
    return ret;
  }
  // Strip outer arrays from *this if Size==1
  template <std::size_t S, std::size_t... Ss,
            typename std::enable_if_t<(sizeof...(Ss) < sizeof...(Sizes) &&
                                       Size == 1)>* = nullptr>
  operator ten<Type, S, Ss...>() const {
    return (*this)[0];
  }
  operator std::array<Element, Size>() const { return data; }
};
template <class Type, std::size_t N>
using vec = ten<Type, N>;
template <class Type, std::size_t N, std::size_t M>
using mat = ten<Type, N, M>;

// Tensor
template <class Type, std::size_t Size, std::size_t... Sizes>
ten<Type, Size, Sizes...>& operator+=(ten<Type, Size, Sizes...>& lhs,
                                      const ten<Type, Size, Sizes...>& rhs) {
  for (std::size_t i = 0; i < Size; ++i) lhs[i] += rhs[i];
  return lhs;
}
template <class Type, std::size_t Size, std::size_t... Sizes>
ten<Type, Size, Sizes...>& operator+(ten<Type, Size, Sizes...> lhs,
                                     const ten<Type, Size, Sizes...>& rhs) {
  return lhs += rhs;
}

template <class Type, std::size_t Size, std::size_t... Sizes>
ten<Type, Size, Sizes...>& operator-=(ten<Type, Size, Sizes...>& lhs,
                                      const ten<Type, Size, Sizes...>& rhs) {
  for (std::size_t i = 0; i < Size; ++i) lhs[i] -= rhs[i];
  return lhs;
}
template <class Type, std::size_t Size, std::size_t... Sizes>
ten<Type, Size, Sizes...>& operator-(ten<Type, Size, Sizes...> lhs,
                                     const ten<Type, Size, Sizes...>& rhs) {
  return lhs -= rhs;
}

template <class Type, std::size_t Size, std::size_t... Sizes>
ten<Type, Size, Sizes...>& operator*=(ten<Type, Size, Sizes...>& lhs,
                                      const Type& rhs) {
  for (std::size_t i = 0; i < Size; ++i) lhs[i] *= rhs;
  return lhs;
}
template <class Type, std::size_t Size, std::size_t... Sizes>
ten<Type, Size, Sizes...> operator*(ten<Type, Size, Sizes...> lhs,
                                    const Type& rhs) {
  return lhs *= rhs;
}
template <class Type, std::size_t Size, std::size_t... Sizes>
ten<Type, Size, Sizes...> operator*(const Type& lhs,
                                    ten<Type, Size, Sizes...> rhs) {
  return rhs *= lhs;
}

template <class Type, std::size_t Size, std::size_t... Sizes>
ten<Type, Size, Sizes...>& operator/=(ten<Type, Size, Sizes...>& lhs,
                                      const Type& rhs) {
  for (std::size_t i = 0; i < Size; ++i) lhs[i] /= rhs;
  return lhs;
}
template <class Type, std::size_t Size, std::size_t... Sizes>
ten<Type, Size, Sizes...> operator/(ten<Type, Size, Sizes...> lhs,
                                    const Type& rhs) {
  return lhs /= rhs;
}
template <class Type, std::size_t Size, std::size_t... Sizes>
ten<Type, Size, Sizes...> operator/(const Type& lhs,
                                    ten<Type, Size, Sizes...> rhs) {
  return rhs /= lhs;
}

template <class Type, std::size_t Size, std::size_t... Sizes>
ten<Type, Size, Sizes...>& operator-(const ten<Type, Size, Sizes...>& a) {
  return -1 * a;
}

template <class Type, std::size_t Size, std::size_t... Sizes>
std::ostream& operator<<(std::ostream& os, const ten<Type, Size, Sizes...>& a) {
  os << '[';
  for (std::size_t i = 0; i < Size - 1; ++i) os << a[i] << ',';
  return os << a[Size - 1] << ']';
}

template <class Type, std::size_t Size, std::size_t... Sizes>
constexpr std::size_t Order(const ten<Type, Size, Sizes...>& a) {
  return sizeof...(Sizes) + 1;
}
template <class Type, std::size_t S, std::size_t... Ss>
constexpr std::size_t Size(const ten<Type, S, Ss...>& a) {
  return S;
}

template <std::size_t Start, std::size_t N, class Type, std::size_t Size,
          std::size_t... Sizes,
          typename std::enable_if_t<(((Start + N) <= Size), N > 0)>* = nullptr>
ten<Type, N, Sizes...> Sub(const ten<Type, Size, Sizes...>& a) {
  ten<Type, N, Sizes...> ret{};
  for (std::size_t i = 0; i < N; ++i) ret[i] = a[Start + i];
  return ret;
}
template <std::size_t Start, class Type, std::size_t Size, std::size_t... Sizes>
ten<Type, Size - Start, Sizes...> Sub(const ten<Type, Size, Sizes...>& a) {
  ten<Type, Size - Start, Sizes...> ret{};
  for (std::size_t i = 0; i < Size - Start; ++i) ret[i] = a[Start + i];
  return ret;
}
// Merge case: lhs convertible to rhs || not rhs convertible to lhs
template <
    class Type, std::size_t Size1, std::size_t Size2, std::size_t... Sizes1,
    std::size_t... Sizes2,
    typename std::enable_if_t<
        (std::is_convertible_v<ten<Type, Sizes1...>, ten<Type, Sizes2...>> ||
         !std::is_convertible_v<ten<Type, Sizes2...>, ten<Type, Sizes1...>>)>* =
        nullptr>
ten<Type, Size1 + Size2, Sizes2...> Merge(
    const ten<Type, Size1, Sizes1...>& lhs,
    const ten<Type, Size2, Sizes2...>& rhs) {
  ten<Type, Size1 + Size2, Sizes2...> ret{};
  ret = lhs;
  for (std::size_t i = 0; i < Size2; ++i) ret[Size1 + i] = rhs[i];
  return ret;
}
// Merge case: not lhs convertible to rhs && rhs convertible to lhs
template <
    class Type, std::size_t Size1, std::size_t Size2, std::size_t... Sizes1,
    std::size_t... Sizes2,
    typename std::enable_if_t<
        (!std::is_convertible_v<ten<Type, Sizes1...>, ten<Type, Sizes2...>> &&
         std::is_convertible_v<ten<Type, Sizes2...>, ten<Type, Sizes1...>>)>* =
        nullptr>
ten<Type, Size1 + Size2, Sizes1...> Merge(
    const ten<Type, Size1, Sizes1...>& lhs,
    const ten<Type, Size2, Sizes2...>& rhs) {
  ten<Type, Size1 + Size2, Sizes1...> ret{};
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
ten<Type, Size - 1> Minor(const ten<Type, Size>& a, const std::size_t& index) {
  ten<Type, Size - 1> ret{};
  for (std::size_t i = 0; i < index; ++i) ret[i] = a[i];
  for (std::size_t i = index + 1; i < Size; ++i) ret[i - 1] = a[i];
  return ret;
}

// template<std::size_t Index, class Type, std::size_t Size, typename
// std::enable_if_t<(Index==0)>* = nullptr> ten<Type, Size-1> Minor(const
// ten<Type, Size>& a){return Sub<1>(a);} template<std::size_t Index, class
// Type, std::size_t Size, typename std::enable_if_t<(Index == Size-1)>* =
// nullptr> ten<Type, Size-1> Minor(const ten<Type, Size>& a){return Sub<0,
// Size-1>(a);}

template <class... Indices, class Type, std::size_t Size, std::size_t... Sizes,
          typename std::enable_if_t<
              (std::conjunction_v<std::is_same<std::size_t, Indices>...> &&
               sizeof...(Indices) == sizeof...(Sizes))>* = nullptr>
ten<Type, Size - 1, (Sizes - 1)...> Minor(const ten<Type, Size, Sizes...>& a,
                                          const std::size_t& index,
                                          Indices... indices) {
  ten<Type, Size - 1, (Sizes - 1)...> ret{};
  for (std::size_t i = 0; i < index; ++i) ret[i] = Minor(a[i], indices...);
  for (std::size_t i = index + 1; i < Size; ++i)
    ret[i - 1] = Minor(a[i], indices...);
  return ret;
}

// Matrix
template <class Type, std::size_t N, std::size_t M, std::size_t U>
mat<Type, N, U> operator*(const mat<Type, N, M>& lhs,
                          const mat<Type, M, U>& rhs) {
  mat<Type, N, U> ret{};
  mat<Type, U, M> trhs{Transpose(rhs)};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < U; ++j) {
      ret[i][j] = DotProduct(lhs[i], trhs[j]);
    }
  }
  return ret;
}
template <class Type, std::size_t N, std::size_t M>
mat<Type, M, N> Transpose(const mat<Type, N, M>& a) {
  mat<Type, M, N> ret;
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < M; ++j) {
      ret[j][i] = a[i][j];
    }
  }
  return ret;
}
template <class Type, std::size_t N>
Type Trace(const mat<Type, N, N>& matrix) {
  Type ret{0};
  for (std::size_t i = 0; i < N; ++i) {
    ret += matrix[i][i];
  }
  return ret;
}
template <class Type, std::size_t N>
Type RevTrace(const mat<Type, N, N>& matrix) {
  Type ret{0};
  for (std::size_t i = 0; i < N; ++i) {
    ret += matrix[i][N - i - 1];
  }
  return ret;
}
// Generalized Outer Product
template <class Type, std::size_t N, std::size_t M>
mat<Type, N, M> TensorProduct(const vec<Type, N>& lhs, const vec<Type, M> rhs) {
  return Transpose(lhs) * mat<Type, 1, M>(rhs);
}
// Exterior/Wedge Product - generalized Cross Product
template <class Type, std::size_t N, std::size_t M>
mat<Type, N, M> ExteriorProduct(const vec<Type, N>& lhs,
                                const vec<Type, M> rhs) {
  return TensorProduct(lhs, rhs) - TensorProduct(rhs, lhs);
}
// Inner Product - generalized Dot Product
template <class Type, std::size_t N>
Type InnerProduct(const vec<Type, N>& lhs, const vec<Type, N> rhs) {
  return Trace(TensorProduct(lhs, rhs));
}
template <class Type, std::size_t N>
Type Determinant(const mat<Type, N, N>& matrix) {
  Type ret{};
  std::array<mat<Type, N - 1, N - 1>, N> minors{Minors(matrix, 0)};
  for (std::size_t i = 0; i < N; ++i)
    ret += Determinant(minors[i]) * matrix[0][i] * (i % 2 == 1 ? -1 : 1);
  return ret;
}
template <class Type>
Type Determinant(const mat<Type, 2, 2>& matrix) {
  return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
}
template <class Type>
Type Determinant(const mat<Type, 1, 1>& matrix) {
  return matrix[0][0];
}

template <class Type, std::size_t N, std::size_t M>
std::array<std::array<mat<Type, N - 1, M - 1>, M>, N> Minors(
    const mat<Type, N, M>& matrix) {
  std::array<std::array<mat<Type, N - 1, M - 1>, M>, N> ret;
  for (std::size_t i = 0; i < N; ++i) ret[i] = Minors(matrix, i);
  return ret;
}
template <class Type, std::size_t N, std::size_t M>
std::array<mat<Type, N - 1, M - 1>, M> Minors(const mat<Type, N, M>& matrix,
                                              const std::size_t& index) {
  std::array<mat<Type, N - 1, M - 1>, M> ret;
  for (std::size_t i = 0; i < M; ++i) ret[i] = Minor(matrix, index, i);
  return ret;
}
template <class Type, std::size_t N>
vec<Type, N> OrthogonalVector(const mat<Type, N - 1, N>& matrix) {
  std::array<mat<Type, N - 1, N - 1>, N> arr{
      Minors(static_cast<mat<Type, N, N>>(matrix), N - 1)};
  vec<Type, N> ret;
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
Type DotProduct(const vec<Type, N>& lhs, const vec<Type, N>& rhs) {
  Type ret{0};
  for (std::size_t i = 0; i < N; ++i) {
    ret += lhs[i] * rhs[i];
  }
  return ret;
}
template <class Type, std::size_t N>
mat<Type, N, 1> Transpose(const vec<Type, N>& a) {
  mat<Type, N, 1> ret;
  for (std::size_t i = 0; i < N; ++i) {
    ret[i][0] = a[i];
  }
  return ret;
}
template <class Type, std::size_t N>
vec<Type, N> Normalize(const vec<Type, N>& a) {
  Type sum{0};
  for (std::size_t i = 0; i < N; ++i) {
    sum += std::pow(a[i], 2);
  }
  return a / std::sqrt(sum);
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
vec<Type, N> CrossProduct(const vec<Type, N>& lhs, const vec<Type, N>& rhs) {
  mat<Type, N, N> m{ExteriorProduct(lhs, rhs)};
  return vec<Type, N>{m[1][2], m[2][0], m[0][1]};
}

template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
mat<Type, N + 1, N + 1> TranslationMatrix(Type Tx, Type Ty, Type Tz) {
  return mat<Type, N + 1, N + 1>{
      {1, 0, 0, Tx}, {0, 1, 0, Ty}, {0, 0, 1, Tz}, {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
mat<Type, N + 1, N + 1> RotationMatrix(Type Rx, Type Ry, Type Rz) {
  return mat<Type, N + 1, N + 1>{{1, 0, 0, 0},
                                 {0, std::cos(Rx), -std::sin(Rx), 0},
                                 {0, std::sin(Rx), std::cos(Rx), 0},
                                 {0, 0, 0, 1}} *
         mat<Type, N + 1, N + 1>{{std::cos(Ry), 0, std::sin(Ry), 0},
                                 {0, 1, 0, 0},
                                 {-std::sin(Ry), 0, std::cos(Ry), 0},
                                 {0, 0, 0, 1}} *
         mat<Type, N + 1, N + 1>{{std::cos(Rz), -std::sin(Rz), 0, 0},
                                 {std::sin(Rz), std::cos(Rz), 0, 0},
                                 {0, 0, 1, 0},
                                 {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
mat<Type, N + 1, N + 1> ScaleMatrix(Type Sx, Type Sy, Type Sz) {
  return mat<Type, N + 1, N + 1>{
      {Sx, 0, 0, 0}, {0, Sy, 0, 0}, {0, 0, Sz, 0}, {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
mat<Type, N + 1, N + 1> ModelMatrix(Type Tx, Type Ty, Type Tz, Type Rx, Type Ry,
                                    Type Rz, Type Sx, Type Sy, Type Sz) {
  return TranslationMatrix<N>() * RotationMatrix<N>() * ScaleMatrix<N>();
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
mat<Type, N + 1, N + 1> ViewMatrix(vec<Type, N> position, vec<Type, N> target,
                                   vec<Type, N> up) {
  vec<Type, 3> view_direction{Normalize(position - target)};
  vec<Type, 3> view_right{Normalize(CrossProduct(up, view_direction))};
  vec<Type, 3> view_up{CrossProduct(view_direction, view_right)};

  return mat<Type, N + 1, N + 1>{
             {view_right[0], view_right[1], view_right[2], 0},
             {view_up[0], view_up[1], view_up[2], 0},
             {view_direction[0], view_direction[1], view_direction[2], 0},
             {0, 0, 0, 1}} *
         mat<Type, N + 1, N + 1>{{1, 0, 0, -position[0]},
                                 {0, 1, 0, -position[1]},
                                 {0, 0, 1, -position[2]},
                                 {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
mat<Type, N + 1, N + 1> OrthographicMatrix(Type near_val, Type far_val,
                                           Type aspect, Type fov) {
  Type top = near_val * std::tan((M_PI / 180) * fov / 2);
  Type bottom = -top;
  Type right = top * aspect;
  Type left = -right;

  return mat<Type, N + 1, N + 1>{
      {1 / (right - left), 0, 0, -((right + left) / (right - left))},
      {0, 2 / (top - bottom), 0, -((top + bottom) / (top - bottom))},
      {0, 0, -(2 / (far_val - near_val)),
       -((far_val + near_val) / (far_val - near_val))},
      {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
mat<Type, N + 1, N + 1> PerspectiveMatrix(Type near_val, Type far_val,
                                          Type aspect, Type fov) {
  Type top = near_val * std::tan((M_PI / 180) * fov / 2);
  Type bottom = -top;
  Type right = top * aspect;
  Type left = -right;

  return mat<Type, N + 1, N + 1>{
      {(2 * near_val) / (right - left), 0, (right + left) / (right - left), 0},
      {0, (2 * near_val) / (top - bottom), (top + bottom) / (top - bottom), 0},
      {0, 0, -((far_val + near_val) / (far_val - near_val)),
       -((2 * far_val * near_val) / (far_val - near_val))},
      {0, 0, -1, 0}};
}
template <class Type, std::size_t N>
mat<Type, N, N> IdentityMatrix() {
  mat<Type, N, N> ret{};
  for (std::size_t i = 0; i < N; ++i) {
    ret[i][i] = 1;
  }
  return ret;
}
}  // namespace UD