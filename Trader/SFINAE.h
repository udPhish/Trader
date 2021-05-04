#pragma once
#include <memory>
#include <type_traits>
namespace UD {
namespace SFINAE {
template <class T, class = void>
struct is_defined : std::false_type {};
template <class T>
struct is_defined<
    T, std::enable_if_t<std::is_object<T>::value &&
                        !std::is_pointer<T>::value && (sizeof(T) > 0)>>
    : std::true_type {};
template <typename First, typename... Rest>
struct all_same {
  static_assert(std::conjunction_v<std::is_same<First, Rest>...>);
  using type = First;
};

template <class T>
struct is_shared_ptr : std::false_type {};
template <class T>
struct is_shared_ptr<std::shared_ptr<T>> : std::true_type {};
template <class T>
using is_shared_ptr_t = typename is_shared_ptr<T>::type;
template <class T>
inline constexpr bool is_shared_ptr_v = is_shared_ptr<T>::value;

template <class T>
struct is_weak_ptr : std::false_type {};
template <class T>
struct is_weak_ptr<std::weak_ptr<T>> : std::true_type {};
template <class T>
using is_weak_ptr_t = typename is_weak_ptr<T>::type;
template <class T>
inline constexpr bool is_weak_ptr_v = is_weak_ptr<T>::value;

template <class T>
struct is_unique_ptr : std::false_type {};
template <class T>
struct is_unique_ptr<std::unique_ptr<T>> : std::true_type {};
template <class T>
using is_unique_ptr_t = typename is_unique_ptr<T>::type;
template <class T>
inline constexpr bool is_unique_ptr_v = is_unique_ptr<T>::value;

template <class T>
struct is_smart_ptr : std::false_type {};
template <class T>
struct is_smart_ptr<std::shared_ptr<T>> : std::true_type {};
template <class T>
struct is_smart_ptr<std::weak_ptr<T>> : std::true_type {};
template <class T>
struct is_smart_ptr<std::unique_ptr<T>> : std::true_type {};
template <class T>
using is_smart_ptr_t = typename is_smart_ptr<T>::type;
template <class T>
inline constexpr bool is_smart_ptr_v = is_smart_ptr<T>::value;

template <class T>
struct is_any_ptr : std::false_type {};
template <class T>
struct is_any_ptr<T*> : std::true_type {};
template <class T>
struct is_any_ptr<std::shared_ptr<T>> : std::true_type {};
template <class T>
struct is_any_ptr<std::weak_ptr<T>> : std::true_type {};
template <class T>
struct is_any_ptr<std::unique_ptr<T>> : std::true_type {};
template <class T>
using is_any_ptr_t = typename is_any_ptr<T>::type;
template <class T>
inline constexpr bool is_any_ptr_v = is_any_ptr<T>::value;

template <class T>
struct remove_smart_ptr {
  typedef T type;
};
template <class T>
struct remove_smart_ptr<std::shared_ptr<T>> {
  typedef T type;
};
template <class T>
struct remove_smart_ptr<std::unique_ptr<T>> {
  typedef T type;
};
template <class T>
struct remove_smart_ptr<std::weak_ptr<T>> {
  typedef T type;
};
template <class T>
using remove_smart_ptr_t = typename remove_smart_ptr<T>::type;

template <class T>
struct remove_any_ptr {
  typedef T type;
};
template <class T>
struct remove_any_ptr<T*> {
  typedef T type;
};
template <class T>
struct remove_any_ptr<std::shared_ptr<T>> {
  typedef T type;
};
template <class T>
struct remove_any_ptr<std::unique_ptr<T>> {
  typedef T type;
};
template <class T>
struct remove_any_ptr<std::weak_ptr<T>> {
  typedef T type;
};
template <class T>
using remove_any_ptr_t = typename remove_any_ptr<T>::type;

template <class T, std::size_t = 0>
constexpr std::shared_ptr<T> make_shared_ptr() {
  return std::move(std::make_shared<T>());
}

template <class T, size_t... Indices>
constexpr std::array<std::shared_ptr<T>, sizeof...(Indices)>
make_shared_array_sequence(std::index_sequence<Indices...>) {
  return std::move<std::array<std::shared_ptr<T>, sizeof...(Indices)>>(
      {make_shared_ptr<T, Indices>()...});
}

template <class T, size_t Size>
constexpr std::array<std::shared_ptr<T>, Size> make_shared_array() {
  return std::move(
      make_shared_array_sequence<T>(std::make_index_sequence<Size>{}));
}
template <class T, size_t Size,
          typename std::enable_if_t<!is_any_ptr_v<T>>* = nullptr>
constexpr std::array<T, Size> make_array() {
  return std::move(std::array<T, Size>{});
}
template <class T, size_t Size,
          typename std::enable_if_t<is_shared_ptr_v<T>>* = nullptr>
constexpr std::array<T, Size> make_array() {
  return std::move(make_shared_array<remove_any_ptr_t<T>, Size>());
}
// template<size_t... Ss>
// constexpr std::array<std::size_t, sizeof...(Ss)> make_index_array()
// {
//     return
// }
template <typename F, size_t... Is>
auto indices_impl(F f, std::index_sequence<Is...>) {
  return f(std::integral_constant<size_t, Is>()...);
}
template <size_t N, class F>
auto indices(F f) {
  return indices_impl(f, std::make_index_sequence<N>());
}
template <class T>
constexpr T product(const T& t1, const T& t2) {
  return t1 * t2;
}
template <class T, class... Ts>
constexpr T product(const T& t1, const T& t2, const Ts&... ts) {
  return product(t1 * t2, ts...);
}
}  // namespace SFINAE
}  // namespace UD