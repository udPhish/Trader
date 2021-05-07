#pragma once

#include <array>
#include <cmath>
#include <functional>
#include <iostream>

#include <boost/math/constants/constants.hpp>

#include "SFINAE.h"

namespace UD {
namespace Tensor {
struct Incrementor {
  virtual ~Incrementor() {}

  virtual bool operator==(const Incrementor& incrementor) = 0;
  virtual bool operator!=(const Incrementor& incrementor) {
    return !(operator==(incrementor));
  }
  virtual Incrementor& operator++() = 0;
  virtual Incrementor operator++(int) = 0;
};
template <class _IteratorType>
struct IteratorWrapper;
template <class _ValueType>
struct Iterator : public Incrementor {
  using ValueType = _ValueType;

  virtual ~Iterator() {}

  virtual ValueType& operator*() = 0;
  template <class IteratorType, class WrapperType = IteratorWrapper>
  static Iterator Wrap(IteratorType iterator) {
    return WrapperType{iterator};
  }
};
template <class _IteratorType>
struct IteratorWrapper : public Iterator<typename _IteratorType::T> {
  using IteratorType = _IteratorType;

 private:
  using Base = Iterator<typename _IteratorType::T>;

 public:
  using typename Base::ValueType;

  IteratorType _iterator;

  virtual ~IteratorWrapper() {}

  bool operator==(const IteratorWrapper& iterator) override {
    return _iterator == iterator._iterator;
  }
  IteratorWrapper& operator++() override {
    ++_iterator;
    return *this;
  }
  IteratorWrapper operator++(int) override {
    IteratorWrapper it{*this};
    _iterator++;
    return it;
  };
  ValueType& operator*() override { return *_iterator; }
};
// template <class _ContainerType>
// struct ReverseIteratorWrapper
//    : public Iterator<typename _ContainerType::value_type> {
//  using ContainerType = _ContainerType;
//
// private:
//  using Base = Iterator<typename ContainerType::value_type>;
//
// public:
//  using typename Base::ValueType;
//
//  typename ContainerType::reverse_iterator _iterator;
//
//  virtual ~ReverseIteratorWrapper() {}
//
//  ReverseIteratorWrapper& operator++() override {
//    ++_iterator;
//    return *this;
//  }
//  ReverseIteratorWrapper operator++(int) override {
//    _iterator++;
//    return *this;
//  };
//  ValueType& operator*() override { return *_iterator; }
//};
template <class _Type>
struct Range {
  using Type = _Type;
  using Iterator = Iterator<Type>;

  virtual ~Range() {}

 protected:
  virtual Type& _at(const std::size_t& index) = 0;

 public:
  virtual std::size_t size() = 0;
  Type& at(const std::size_t& index) { return this->_at(index); }
  const Type& at(const std::size_t& index) const { return this->_at(index); }
  Type& operator[](const std::size_t& index) { return this->_at(index); }
  const Type& operator[](const std::size_t& index) const {
    return this->_at(index);
  }
  virtual Iterator begin() = 0;
  virtual Iterator end() = 0;
  virtual Iterator rbegin() = 0;
  virtual Iterator rend() = 0;
};
template <class _ContainerType>
struct ContainerWrapper : public Range<typename _ContainerType::value_type> {
  using ContainerType = _ContainerType;

 private:
  using Base = Range<typename ContainerType::value_type>;

 public:
  using typename Base::Iterator;
  using typename Base::Type;

  ContainerType _container;

  virtual ~ContainerWrapper() {}
  ContainerWrapper() {}
  ContainerWrapper(ContainerType container) : _container{container} {}
  ContainerWrapper(std::initializer_list<Type> list) : _container{list} {}
  operator ContainerType() { return _container; }

  ContainerType& data() { return _container; }
  std::size_t size() override { return _container.size(); }

 protected:
  Type& _at(const std::size_t& index) override { return _container.at(index); }

 public:
  Iterator begin() override { return Iterator::Wrap(_container.begin()); }
  Iterator end() override { return Iterator::Wrap(_container.end()); }
  Iterator rbegin() override { return Iterator::Wrap(_container.rbegin()); }
  Iterator rend() override { return Iterator::Wrap(_container.rend()); }
};
template <class Type, std::size_t Size>
struct Array : public ContainerWrapper<std::array<Type, Size>> {
 private:
  using Base = ContainerWrapper<std::array<Type, Size>>;

 public:
  virtual ~Array() {}

  using Base::Base;

  using Base::operator typename Base::ContainerType;

  auto front() { return this->data().front(); }
  auto back() { return this->data().back(); }
  auto empty() const { return this->data().empty(); }
  constexpr void fill(const Type& value) { this->data().fill(value); }
};

template <class Type, std::size_t Size, std::size_t... Sizes>
struct Tensor : public Array<Tensor<Type, Sizes...>, Size> {
  using Element = Tensor<Type, Sizes...>;

 private:
  using Base = Array<Element, Size>;

 public:
  virtual ~Tensor() {}
  using Base::Base;
  using Base::fill;
  constexpr void fill(const Type& value) {
    for (auto& element : *this) element.fill(value);
  }
};
template <class Type, std::size_t Size>
struct Tensor<Type, Size> : public Array<Type, Size> {
  using Element = Type;

 private:
  using Base = Array<Element, Size>;

 public:
  virtual ~Tensor() {}
  using Base::Base;
  template <typename std::enable_if_t<(Size >= 1)>* = nullptr>
  Type& x() {
    return this->at(0);
  }
  template <typename std::enable_if_t<(Size >= 2)>* = nullptr>
  Type& y() {
    return this->at(1);
  }
  template <typename std::enable_if_t<(Size >= 3)>* = nullptr>
  Type& z() {
    return this->at(2);
  }
};
template <class Type, std::size_t Size, std::size_t... TensorSizes>
struct Tensor<Tensor<Type, TensorSizes...>, Size>
    : public Tensor<Type, Size, TensorSizes...> {
  using Tensor<Type, Size, TensorSizes...>;
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
template <class Type, std::size_t Size>
void Scale(Tensor<Type, Size>& tensor, const Type& factor) {
  for (auto& ele : tensor) ele *= factor;
}
template <class Type, std::size_t Size, std::size_t... Sizes>
void Scale(Tensor<Type, Size, Sizes...>& tensor, const Type& factor) {
  for (auto& ele : tensor) Scale(ele);
}
template <class Type, std::size_t ASize, std::size_t... ASizes,
          std::size_t... BSizes>
void Scale(Tensor<Type, ASize, ASizes..., BSizes...>& lhs,
           const Tensor<Type, BSizes...>& rhs) {
  for (std::size_t i = 0; i < ASize; ++i) Scale(lhs[i], rhs[i]);
}
template <class Type, std::size_t Size>
void Translate(Tensor<Type, Size>& tensor, const Type& factor) {
  for (auto& ele : tensor) ele += factor;
}
template <class Type, std::size_t Size, std::size_t... Sizes>
void Translate(Tensor<Type, Size, Sizes...>& tensor, const Type& factor) {
  for (auto& ele : tensor) Translate(ele);
}
template <class Type, std::size_t ASize, std::size_t... ASizes,
          std::size_t... BSizes>
void Translate(Tensor<Type, ASize, ASizes..., BSizes...>& lhs,
               const Tensor<Type, BSizes...>& rhs) {
  for (std::size_t i = 0; i < ASize; ++i) Translate(lhs[i], rhs[i]);
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
template <class Type, std::size_t N>
Type Magnitude(const Vector<Type, N>& vector) {
  Type sum = 0;
  for (auto& v : vector) sum += std::pow(v, 2);
  return std::sqrt(sum);
}
template <class Type, std::size_t N>
Vector<Type, N> UnitVector(const Vector<Type, N>& vector) {
  return vector / Magnitude(vector);
}
template <class Type, std::size_t N>
Type Sum(const Vector<Type, N>& vector) {
  Type sum = 0;
  for (auto& v : vector) sum += v;
}
template <class Type, std::size_t... Sizes>
Tensor<Type, Sizes...> Mutliply(const Tensor<Type, Sizes...>& lhs,
                                const Tensor<Type, Sizes...>& rhs) {
  return Transform(lhs, rhs, std::multiplies);
}
template <class Type, std::size_t Size, std::size_t... Sizes>
Tensor<Type, Size, Sizes...> Transform(
    const Tensor<Type, Size, Sizes...>& lhs,
    const Tensor<Type, Size, Sizes...>& rhs,
    std::function<Type(const Type&, const Type&)> op) {
  Tensor<Type, Size, Sizes...> ret;
  for (std::size_t i = 0; i < Size; ++i) ret[i] = Transform(lhs[i], rhs[i], op);
  return ret;
}
template <class Type, std::size_t Size>
Tensor<Type, Size> Transform(const Tensor<Type, Size>& lhs,
                             const Tensor<Type, Size>& rhs,
                             std::function<Type(const Type&, const Type&)> op) {
  Tensor<Type, Size> ret;
  std::transform(lhs.begin(), lhs.end(), rhs.begin(), ret.begin(), op);
  return ret;
}
}  // namespace Tensor
}  // namespace UD