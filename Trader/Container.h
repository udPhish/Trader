#pragma once
#include <array>
#include <cstddef>
#include <span>
#include <type_traits>

#include "Math.h"
#include "Traits.h"

namespace UD {
namespace Container {
template <class _Element, Math::Type::ULong _Size>
struct Array;
}
namespace Traits {
template <class _Element, Math::Type::ULong _Size, class _Leaf>
struct Register<Container::Array<_Element, _Size>, _Leaf> {
  using Element = _Element;
  static constexpr Math::Type::ULong Size() { return _Size; }
  using UnderlyingContainer = std::array<Element, Size()>;
  using Iterator = typename UnderlyingContainer::iterator;
  using ConstIterator = typename UnderlyingContainer::const_iterator;
  using ReverseIterator = typename UnderlyingContainer::reverse_iterator;
  using ConstReverseIterator =
      typename UnderlyingContainer::const_reverse_iterator;
};
}  // namespace Traits
namespace Container {
namespace Concepts {}

template <class _Element, Math::Type::ULong _Size>
struct Array {
  using Traits = Traits::Traits<Array>;

  typename Traits::UnderlyingContainer _container;

  virtual ~Array() {}
  Array(const Array&)       = default;
  Array(Array&&) noexcept = default;
  auto operator=(const Array&) -> Array& = default;
  auto operator=(Array&&) noexcept -> Array& = default;

  Array() : _container{} {}
  template <class... Ts>
  requires requires(Ts... ts) {
    {typename Traits::UnderlyingContainer{{ts...}}};
  }
  Array(Ts... ts) : _container{{ts...}} {}
  template <Math::Type::ULong S>
  requires(S <= _Size) Array(const Array<_Element, S>& arr) {
    for (Math::Type::ULong i = 0; i < S; ++i) this->at(i) = arr[i];
  }
  // Array(const Array& other) : _container{other._container} {}
  // template <class... Ts>
  // Array(Ts... ts) : _container{ts...} {}
  // Array() {}
  // Array(const typename Traits::UnderlyingContainer& container)
  //    : _container{container} {}
  // Array(std::initializer_list<typename Traits::Element> list)
  //    : _container(list) {}

  typename Traits::UnderlyingContainer& underlying_container() {
    return _container;
  }
  const typename Traits::UnderlyingContainer& underlying_container() const {
    return _container;
  }
  auto data() { return underlying_container().data(); }
  auto begin() { return underlying_container().begin(); }
  auto end() { return underlying_container().end(); }
  auto begin() const { return underlying_container().begin(); }
  auto end() const { return underlying_container().end(); }
  auto rbegin() { return underlying_container().rbegin(); }
  auto rend() { return underlying_container().rend(); }
  auto rbegin() const { return underlying_container().rbegin(); }
  auto rend() const { return underlying_container().rend(); }
  auto cbegin() const { return underlying_container().cbegin(); }
  auto cend() const { return underlying_container().cend(); }
  auto crbegin() const { return underlying_container().crbegin(); }
  auto crend() const { return underlying_container().crend(); }
  auto front() { return underlying_container().front(); }
  auto back() { return underlying_container().back(); }
  auto front() const { return underlying_container().front(); }
  auto back() const { return underlying_container().back(); }
  auto empty() const { return underlying_container().empty(); }
  typename Traits::Element& at(Math::Type::ULong index) {
    return underlying_container().at(index);
  }
  const typename Traits::Element& at(Math::Type::ULong index) const {
    return underlying_container().at(index);
  }
  typename Traits::Element& operator[](std::size_t index) {
    return underlying_container()[index];
  }
  template <class T = typename Traits::Element>
      requires (!std::is_fundamental_v<T>) const typename Traits::Element &
      operator[](std::size_t index) const {
    return underlying_container()[index];
  }
  template <class T = typename Traits::Element>
  requires std::is_fundamental_v<T> typename Traits::Element operator[](
      std::size_t index) const {
    return underlying_container()[index];
  }
  constexpr Math::Type::ULong size() const {
    return underlying_container().size();
  }
  constexpr void fill(const typename Traits::Element& value) {
    underlying_container().fill(value);
  }

  // template <
  //    Math::Type::ULong Start, Math::Type::ULong N, class Type,
  //    Math::Type::ULong Size,
  //    typename std::enable_if_t<(((Start + N) <= Size), N > 0)>* = nullptr>
  // Array<Type, N> Sub() {
  //  Array<Type, N> ret{};
  //  for (Math::Type::ULong i = 0; i < N; ++i) ret[i] = this->at(Start + i);
  //  return ret;
  //}
  // template <Math::Type::ULong Start, class Type, Math::Type::ULong Size>
  // Array<Type, Size - Start> Sub() {
  //  return Sub<Start, Size - Start, Type, Size>();
  //}
};
}  // namespace Container
}  // namespace UD