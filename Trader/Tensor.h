#pragma once

#include <array>
#include <cmath>
#include <iostream>

#include <boost/math/constants/constants.hpp>

#include "SFINAE.h"

namespace UD {
namespace Tensor {
struct Incrementor {
  virtual ~Incrementor() {}
  virtual Incrementor& operator++() = 0;
  virtual Incrementor operator++(int) = 0;
};
template <class _ValueType>
struct Iterator : public Incrementor {
  using ValueType = _ValueType;
  virtual ~Iterator() {}
  virtual ValueType& operator*() = 0;
};
template <class _ContainerType>
struct IteratorWrapper : public Iterator<typename _ContainerType::value_type> {
  using ContainerType = _ContainerType;
  typename ContainerType::iterator _iterator;

  virtual ~IteratorWrapper() {}

  IteratorWrapper& operator++() override {
    ++_iterator;
    return *this;
  }
  IteratorWrapper operator++(int) override {
    _iterator++;
    return *this;
  };
  typename Iterator<typename _ContainerType::value_type>::ValueType& operator*()
      override {
    return *_iterator;
  }
};
template <class _ContainerType>
struct ReverseIteratorWrapper
    : public Iterator<typename _ContainerType::value_type> {
  using ContainerType = _ContainerType;
  typename ContainerType::reverse_iterator _iterator;

  virtual ~ReverseIteratorWrapper() {}

  ReverseIteratorWrapper& operator++() override {
    ++_iterator;
    return *this;
  }
  ReverseIteratorWrapper operator++(int) override {
    _iterator++;
    return *this;
  };
  typename Iterator<typename _ContainerType::value_type>::ValueType& operator*()
      override {
    return *_iterator;
  }
};
template <class _Type>
struct Range {
  using Type = _Type;
  using Iterator = Iterator<Type>;
  virtual ~Range() {}
  virtual std::size_t size() = 0;
  virtual Type& at(const std::size_t& index) = 0;
  Type& operator[](const std::size_t& index) { return this->at(index); }
  virtual Iterator begin() = 0;
  virtual Iterator end() = 0;
  virtual Iterator rbegin() = 0;
  virtual Iterator rend() = 0;
};
template <class _ContainerType>
struct ContainerWrapper : public Range<typename _ContainerType::value_type> {
  using ContainerType = _ContainerType;
  using Type = typename Range<typename ContainerType::value_type>::Type;
  using Iterator = typename Range<typename ContainerType::value_type>::Iterator;
  ContainerType _container;
  ContainerWrapper() {}
  virtual ~ContainerWrapper() {}
  ContainerWrapper(ContainerType container) : _container(container) {}
  operator ContainerType() { return _container; }

  ContainerType& data() { return _container; }

  std::size_t size() { return _container.size(); }
  Type& at(const std::size_t& index) { return _container.at(index); }
  Iterator begin() { return IteratorWrapper{_container.begin()}; }
  Iterator end() { return IteratorWrapper{_container.end()}; }
  Iterator rbegin() { return ReverseIteratorWrapper{_container.rbegin()}; }
  Iterator rend() { return ReverseIteratorWrapper{_container.rend()}; }
};
template <class Type, std::size_t Size>
struct Array : public ContainerWrapper<std::array<Type, Size>> {
  virtual ~Array() {}

  Array() {}
  Array(std::array<Type, Size> array)
      : ContainerWrapper<std::array<Type, Size>>{array} {}
  operator std::array<Type, Size>() {
    return static_cast<ContainerWrapper<std::array<Type, Size>>>(*this);
  }

  auto front() { return this->data().front(); }
  auto back() { return this->data().back(); }
  auto empty() const { return this->data().empty(); }
};

template <class Type, std::size_t Size, std::size_t... Sizes>
struct Tensor : public Array<Tensor<Type, Sizes...>, Size> {
  using Element = Tensor<Type, Sizes...>;
  virtual ~Tensor() {}
};
template <class Type, std::size_t Size>
struct Tensor<Type, Size> : public Array<Type, Size> {
  using Element = Type;
  virtual ~Tensor() {}
};

template <class Type, std::size_t N, std::size_t M>
using Matrix = Tensor<Type, N, M>;
template <class Type, std::size_t N>
using Vector = Tensor<Type, N>;

template <class Type, std::size_t Size, std::size_t... Sizes>
Tensor<Type, Size, Sizes...>& operator+=(
    Tensor<Type, Size, Sizes...>& lhs,
    const Tensor<Type, Size, Sizes...>& rhs) {
  for (std::size_t i = 0; i < Size; ++i) lhs[i] += rhs[i];
  return lhs;
}
template <class Type, std::size_t... Sizes>
Tensor<Type, Sizes...>& operator+(Tensor<Type, Sizes...> lhs,
                                  const Tensor<Type, Sizes...>& rhs) {
  return lhs += rhs;
}

template <class Type, std::size_t Size, std::size_t... Sizes>
Tensor<Type, Size, Sizes...>& operator-=(
    Tensor<Type, Size, Sizes...>& lhs,
    const Tensor<Type, Size, Sizes...>& rhs) {
  for (std::size_t i = 0; i < Size; ++i) lhs[i] -= rhs[i];
  return lhs;
}
template <class Type, std::size_t... Sizes>
Tensor<Type, Sizes...>& operator-(Tensor<Type, Sizes...> lhs,
                                  const Tensor<Type, Sizes...>& rhs) {
  return lhs -= rhs;
}

template <class Type, std::size_t Size, std::size_t... Sizes>
Tensor<Type, Size, Sizes...>& operator*=(Tensor<Type, Size, Sizes...>& lhs,
                                         const Type& rhs) {
  for (std::size_t i = 0; i < Size; ++i) lhs[i] *= rhs;
  return lhs;
}
template <class Type, std::size_t... Sizes>
Tensor<Type, Sizes...> operator*(Tensor<Type, Sizes...> lhs, const Type& rhs) {
  return lhs *= rhs;
}
template <class Type, std::size_t... Sizes>
Tensor<Type, Sizes...> operator*(const Type& lhs, Tensor<Type, Sizes...> rhs) {
  return rhs *= lhs;
}

template <class Type, std::size_t Size, std::size_t... Sizes>
Tensor<Type, Size, Sizes...>& operator/=(Tensor<Type, Size, Sizes...>& lhs,
                                         const Type& rhs) {
  for (std::size_t i = 0; i < Size; ++i) lhs[i] /= rhs;
  return lhs;
}
template <class Type, std::size_t Size, std::size_t... Sizes>
Tensor<Type, Size, Sizes...> operator/(Tensor<Type, Size, Sizes...> lhs,
                                       const Type& rhs) {
  return lhs /= rhs;
}
template <class Type, std::size_t... Sizes>
Tensor<Type, Sizes...> operator/(const Type& lhs, Tensor<Type, Sizes...> rhs) {
  return rhs /= lhs;
}

template <class Type, std::size_t... Sizes>
Tensor<Type, Sizes...> operator-(Tensor<Type, Sizes...> a) {
  return a *= -1;
}
template <class Type, std::size_t Size, std::size_t... Sizes>
std::ostream& operator<<(std::ostream& os,
                         const Tensor<Type, Size, Sizes...>& a) {
  if (Size == 0) return os << "[]";
  os << '[';
  for (std::size_t i = 0; i < Size - 1; ++i) os << a[i] << ',';
  return os << a[Size - 1] << ']';
}

template <class Type, std::size_t... Sizes>
constexpr std::size_t Order(const Tensor<Type, Sizes...>& a) {
  return sizeof...(Sizes);
}
template <class Type, std::size_t S, std::size_t... Ss>
constexpr std::size_t Size(const Tensor<Type, S, Ss...>& a) {
  return S;
}
template <std::size_t Start, std::size_t N, class Type, std::size_t Size,
          typename std::enable_if_t<(((Start + N) <= Size), N > 0)>* = nullptr>
Array<Type, N> Sub(const Array<Type, Size>& a) {
  Array<Type, N> ret{};
  for (std::size_t i = 0; i < N; ++i) ret[i] = a[Start + i];
  return ret;
}
template <std::size_t Start, class Type, std::size_t Size>
Array<Type, Size - Start> Sub(const Array<Type, Size>& a) {
  return Sub<Start, Size - Start, Type, Size>(a);
}

// Merge case: lhs convertible to rhs || not rhs convertible to lhs
template <class Type, std::size_t Size1, std::size_t Size2,
          std::size_t... Sizes1, std::size_t... Sizes2,
          typename std::enable_if_t<
              (std::is_convertible_v<Tensor<Type, Sizes1...>,
                                     Tensor<Type, Sizes2...>> ||
               !std::is_convertible_v<Tensor<Type, Sizes2...>,
                                      Tensor<Type, Sizes1...>>)>* = nullptr>
Tensor<Type, Size1 + Size2, Sizes2...> Merge(
    const Tensor<Type, Size1, Sizes1...>& lhs,
    const Tensor<Type, Size2, Sizes2...>& rhs) {
  Tensor<Type, Size1 + Size2, Sizes2...> ret{};
  ret = lhs;
  for (std::size_t i = 0; i < Size2; ++i) ret[Size1 + i] = rhs[i];
  return ret;
}
// Merge case: not lhs convertible to rhs && rhs convertible to lhs
template <class Type, std::size_t Size1, std::size_t Size2,
          std::size_t... Sizes1, std::size_t... Sizes2,
          typename std::enable_if_t<
              (!std::is_convertible_v<Tensor<Type, Sizes1...>,
                                      Tensor<Type, Sizes2...>> &&
               std::is_convertible_v<Tensor<Type, Sizes2...>,
                                     Tensor<Type, Sizes1...>>)>* = nullptr>
Tensor<Type, Size1 + Size2, Sizes1...> Merge(
    const Tensor<Type, Size1, Sizes1...>& lhs,
    const Tensor<Type, Size2, Sizes2...>& rhs) {
  Tensor<Type, Size1 + Size2, Sizes1...> ret{};
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
Tensor<Type, Size - 1> Minor(const Tensor<Type, Size>& a,
                             const std::size_t& index) {
  Tensor<Type, Size - 1> ret{};
  for (std::size_t i = 0; i < index; ++i) ret[i] = a[i];
  for (std::size_t i = index + 1; i < Size; ++i) ret[i - 1] = a[i];
  return ret;
}
template <class... Indices, class Type, std::size_t Size, std::size_t... Sizes,
          typename std::enable_if_t<
              (std::conjunction_v<std::is_same<std::size_t, Indices>...> &&
               sizeof...(Indices) == sizeof...(Sizes))>* = nullptr>
Tensor<Type, Size - 1, (Sizes - 1)...> Minor(
    const Tensor<Type, Size, Sizes...>& a, const std::size_t& index,
    Indices... indices) {
  Tensor<Type, Size - 1, (Sizes - 1)...> ret{};
  for (std::size_t i = 0; i < index; ++i) ret[i] = Minor(a[i], indices...);
  for (std::size_t i = index + 1; i < Size; ++i)
    ret[i - 1] = Minor(a[i], indices...);
  return ret;
}
template <class Type, std::size_t N>
Tensor<Type, N> Transpose(const Tensor<Type, N>& a) {
  return a;
}
template <class Type, std::size_t N, std::size_t M>
Tensor<Type, M, N> Transpose(const Tensor<Type, N, M>& a) {
  Tensor<Type, M, N> ret;
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < M; ++j) {
      ret[j][i] = a[i][j];
    }
  }
  return ret;
}
template <class Type, std::size_t N>
Type Trace(const Tensor<Type, N, N>& matrix) {
  Type ret{0};
  for (std::size_t i = 0; i < N; ++i) {
    ret += matrix[i][i];
  }
  return ret;
}
template <class Type, std::size_t N>
Type RevTrace(const Tensor<Type, N, N>& matrix) {
  Type ret{0};
  for (std::size_t i = 0; i < N; ++i) {
    ret += matrix[i][N - i - 1];
  }
  return ret;
}
template <class Type, std::size_t N, std::size_t M, std::size_t U>
Matrix<Type, N, U> operator*(const Matrix<Type, N, M>& lhs,
                             const Matrix<Type, M, U>& rhs) {
  Matrix<Type, N, U> ret{};
  Matrix<Type, U, M> trhs{Transpose(rhs)};
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
Matrix<Type, N, M> TensorProduct(const Vector<Type, N>& lhs,
                                 const Vector<Type, M> rhs) {
  return Transpose(lhs) * Matrix<Type, 1, M>(rhs);
}
// Exterior/Wedge Product - generalized Cross Product
template <class Type, std::size_t N, std::size_t M>
Matrix<Type, N, M> ExteriorProduct(const Vector<Type, N>& lhs,
                                   const Vector<Type, M> rhs) {
  return TensorProduct(lhs, rhs) - TensorProduct(rhs, lhs);
}
// Inner Product - generalized Dot Product
template <class Type, std::size_t N>
Type InnerProduct(const Vector<Type, N>& lhs, const Vector<Type, N> rhs) {
  return Trace(TensorProduct(lhs, rhs));
}
template <class Type, std::size_t N>
Type Determinant(const Matrix<Type, N, N>& matrix) {
  Type ret{};
  std::array<Matrix<Type, N - 1, N - 1>, N> minors{Minors(matrix, 0)};
  for (std::size_t i = 0; i < N; ++i)
    ret += Determinant(minors[i]) * matrix[0][i] * (i % 2 == 1 ? -1 : 1);
  return ret;
}
template <class Type>
Type Determinant(const Matrix<Type, 2, 2>& matrix) {
  return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
}
template <class Type>
Type Determinant(const Matrix<Type, 1, 1>& matrix) {
  return matrix[0][0];
}

template <class Type, std::size_t N, std::size_t M>
std::array<std::array<Matrix<Type, N - 1, M - 1>, M>, N> Minors(
    const Matrix<Type, N, M>& matrix) {
  std::array<std::array<Matrix<Type, N - 1, M - 1>, M>, N> ret;
  for (std::size_t i = 0; i < N; ++i) ret[i] = Minors(matrix, i);
  return ret;
}
template <class Type, std::size_t N, std::size_t M>
std::array<Matrix<Type, N - 1, M - 1>, M> Minors(
    const Matrix<Type, N, M>& matrix, const std::size_t& index) {
  std::array<Matrix<Type, N - 1, M - 1>, M> ret;
  for (std::size_t i = 0; i < M; ++i) ret[i] = Minor(matrix, index, i);
  return ret;
}
template <class Type, std::size_t N>
Vector<Type, N> OrthogonalVector(const Matrix<Type, N - 1, N>& matrix) {
  std::array<Matrix<Type, N - 1, N - 1>, N> arr{
      Minors(static_cast<Matrix<Type, N, N>>(matrix), N - 1)};
  Vector<Type, N> ret;
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
Type DotProduct(const Vector<Type, N>& lhs, const Vector<Type, N>& rhs) {
  Type ret{0};
  for (std::size_t i = 0; i < N; ++i) {
    ret += lhs[i] * rhs[i];
  }
  return ret;
}
template <class Type, std::size_t N>
Matrix<Type, N, 1> Transpose(const Vector<Type, N>& a) {
  return Transpose(Matrix<Type, 1, N>(a));
}
template <class Type, std::size_t N>
Vector<Type, N> Normalize(const Vector<Type, N>& a) {
  Type sum{0};
  for (std::size_t i = 0; i < N; ++i) {
    sum += std::pow(a[i], 2);
  }
  return a / std::sqrt(sum);
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
Vector<Type, N> CrossProduct(const Vector<Type, N>& lhs,
                             const Vector<Type, N>& rhs) {
  Matrix<Type, N, N> m{ExteriorProduct(lhs, rhs)};
  return Vector<Type, N>{m[1][2], m[2][0], m[0][1]};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
Matrix<Type, N + 1, N + 1> TranslationMatrix(Type Tx, Type Ty, Type Tz) {
  return Matrix<Type, N + 1, N + 1>{
      {1, 0, 0, Tx}, {0, 1, 0, Ty}, {0, 0, 1, Tz}, {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
Matrix<Type, N + 1, N + 1> RotationMatrix(Type Rx, Type Ry, Type Rz) {
  return Matrix<Type, N + 1, N + 1>{{1, 0, 0, 0},
                                    {0, std::cos(Rx), -std::sin(Rx), 0},
                                    {0, std::sin(Rx), std::cos(Rx), 0},
                                    {0, 0, 0, 1}} *
         Matrix<Type, N + 1, N + 1>{{std::cos(Ry), 0, std::sin(Ry), 0},
                                    {0, 1, 0, 0},
                                    {-std::sin(Ry), 0, std::cos(Ry), 0},
                                    {0, 0, 0, 1}} *
         Matrix<Type, N + 1, N + 1>{{std::cos(Rz), -std::sin(Rz), 0, 0},
                                    {std::sin(Rz), std::cos(Rz), 0, 0},
                                    {0, 0, 1, 0},
                                    {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
Matrix<Type, N + 1, N + 1> ScaleMatrix(Type Sx, Type Sy, Type Sz) {
  return Matrix<Type, N + 1, N + 1>{
      {Sx, 0, 0, 0}, {0, Sy, 0, 0}, {0, 0, Sz, 0}, {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
Matrix<Type, N + 1, N + 1> ModelMatrix(Type Tx, Type Ty, Type Tz, Type Rx,
                                       Type Ry, Type Rz, Type Sx, Type Sy,
                                       Type Sz) {
  return TranslationMatrix<N>() * RotationMatrix<N>() * ScaleMatrix<N>();
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
Matrix<Type, N + 1, N + 1> ViewMatrix(Vector<Type, N> position,
                                      Vector<Type, N> target,
                                      Vector<Type, N> up) {
  Vector<Type, 3> view_direction{Normalize(position - target)};
  Vector<Type, 3> view_right{Normalize(CrossProduct(up, view_direction))};
  Vector<Type, 3> view_up{CrossProduct(view_direction, view_right)};

  return Matrix<Type, N + 1, N + 1>{
             {view_right[0], view_right[1], view_right[2], 0},
             {view_up[0], view_up[1], view_up[2], 0},
             {view_direction[0], view_direction[1], view_direction[2], 0},
             {0, 0, 0, 1}} *
         Matrix<Type, N + 1, N + 1>{{1, 0, 0, -position[0]},
                                    {0, 1, 0, -position[1]},
                                    {0, 0, 1, -position[2]},
                                    {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
Matrix<Type, N + 1, N + 1> OrthographicMatrix(Type near_val, Type far_val,
                                              Type aspect, Type fov) {
  Type top = near_val *
             std::tan((boost::math::constants::pi<double>() / 180) * fov / 2);
  Type bottom = -top;
  Type right = top * aspect;
  Type left = -right;

  return Matrix<Type, N + 1, N + 1>{
      {1 / (right - left), 0, 0, -((right + left) / (right - left))},
      {0, 2 / (top - bottom), 0, -((top + bottom) / (top - bottom))},
      {0, 0, -(2 / (far_val - near_val)),
       -((far_val + near_val) / (far_val - near_val))},
      {0, 0, 0, 1}};
}
template <class Type, std::size_t N,
          typename std::enable_if_t<(N == 3)>* = nullptr>
Matrix<Type, N + 1, N + 1> PerspectiveMatrix(Type near_val, Type far_val,
                                             Type aspect, Type fov) {
  Type top = near_val *
             std::tan((boost::math::constants::pi<double>() / 180) * fov / 2);
  Type bottom = -top;
  Type right = top * aspect;
  Type left = -right;

  return Matrix<Type, N + 1, N + 1>{
      {(2 * near_val) / (right - left), 0, (right + left) / (right - left), 0},
      {0, (2 * near_val) / (top - bottom), (top + bottom) / (top - bottom), 0},
      {0, 0, -((far_val + near_val) / (far_val - near_val)),
       -((2 * far_val * near_val) / (far_val - near_val))},
      {0, 0, -1, 0}};
}
template <class Type, std::size_t N>
Matrix<Type, N, N> IdentityMatrix() {
  Matrix<Type, N, N> ret{};
  for (std::size_t i = 0; i < N; ++i) {
    ret[i][i] = 1;
  }
  return ret;
}
}  // namespace Tensor

}  // namespace UD