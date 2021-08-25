#pragma once

#include <wx/wx.h>

#include <wx/colour.h>

#include "libudph/Math/Tensor.h"

namespace UD
{
namespace Colour
{
template<
    class Type    = unsigned char,
    Type MaxValue = std::numeric_limits<Type>::max(),
    Type MinValue = std::numeric_limits<Type>::min()>
struct Colour : public UD::Tensor::Vector<4, Type>
{
  static_assert(!std::numeric_limits<Type>::is_signed, "Type must not be signed");
  static_assert(std::numeric_limits<Type>::is_exact, "Type must have an exact binary representation");

 private:
  using Base = UD::Tensor::Vector<4, Type>;

 public:
  virtual ~Colour() {}
  Colour() {}
  Colour(Type red, Type green, Type blue)
      : Base{red, green, blue, MaxValue} {}
  Colour(Type red, Type green, Type blue, Type alpha)
      : Base{red, green, blue, alpha} {}
  template<class T>
  Colour(Colour<T> other)
  {
    std::size_t bits       = std::numeric_limits<Type>::digits;
    std::size_t other_bits = std::numeric_limits<T>::digits;
    if (bits > other_bits)
    {
      r() = static_cast<Type>(other.r() << (bits - other_bits));
      g() = static_cast<Type>(other.g() << (bits - other_bits));
      b() = static_cast<Type>(other.b() << (bits - other_bits));
      a() = static_cast<Type>(other.a() << (bits - other_bits));
    }
    else if (other_bits > bits)
    {
      r() = static_cast<Type>(other.r() >> (other_bits - bits));
      g() = static_cast<Type>(other.g() >> (other_bits - bits));
      b() = static_cast<Type>(other.b() >> (other_bits - bits));
      a() = static_cast<Type>(other.a() >> (other_bits - bits));
    }
    else
    {
      r() = other.r();
      g() = other.g();
      b() = other.b();
      a() = other.a();
    }
  }
  Colour(wxColour wxcolour)
      : Base{Colour<unsigned char>{wxcolour.Red(), wxcolour.Green(), wxcolour.Blue(), wxcolour.Alpha()}} {}
  operator wxColour() const
  {
    return wxColour(r(), g(), b(), a());
  }

  static constexpr Type max_value()
  {
    return MaxValue;
  }
  static constexpr Type min_value()
  {
    return MinValue;
  }
  static constexpr Type count_values()
  {
    return max_value() - min_value() + 1;
  }

  Type& red()
  {
    return this->at(0);
  }
  const Type& red() const
  {
    return this->at(0);
  }
  Type& r()
  {
    return red();
  }
  const Type& r() const
  {
    return red();
  }
  Type& green()
  {
    return this->at(1);
  }
  const Type& green() const
  {
    return this->at(1);
  }
  Type& g()
  {
    return green();
  }
  const Type& g() const
  {
    return green();
  }
  Type& blue()
  {
    return this->at(2);
  }
  const Type& blue() const
  {
    return this->at(2);
  }
  Type& b()
  {
    return blue();
  }
  const Type& b() const
  {
    return blue();
  }
  Type& alpha()
  {
    return this->at(3);
  }
  const Type& alpha() const
  {
    return this->at(3);
  }
  Type& a()
  {
    return alpha();
  }
  const Type& a() const
  {
    return alpha();
  }
};

template<
    class Type    = unsigned char,
    Type MaxValue = std::numeric_limits<Type>::max(),
    Type MinValue = std::numeric_limits<Type>::min()>
struct White : public Colour<Type, MaxValue, MinValue>
{
  White()
      : Colour<Type, MaxValue, MinValue>{MaxValue, MaxValue, MaxValue} {}
};
template<
    class Type    = unsigned char,
    Type MaxValue = std::numeric_limits<Type>::max(),
    Type MinValue = std::numeric_limits<Type>::min()>
struct Black : public Colour<Type, MaxValue, MinValue>
{
  Black()
      : Colour<Type, MaxValue, MinValue>{MinValue, MinValue, MinValue} {}
};
template<
    class Type    = unsigned char,
    Type MaxValue = std::numeric_limits<Type>::max(),
    Type MinValue = std::numeric_limits<Type>::min()>
struct Grey : public Colour<Type, MaxValue, MinValue>
{
  Grey()
      : Colour<Type, MaxValue, MinValue>{MaxValue / 2, MaxValue / 2, MaxValue / 2} {}
};
template<
    class Type    = unsigned char,
    Type MaxValue = std::numeric_limits<Type>::max(),
    Type MinValue = std::numeric_limits<Type>::min()>
using Gray = Grey<Type, MaxValue, MinValue>;

template<
    class Type    = unsigned char,
    Type MaxValue = std::numeric_limits<Type>::max(),
    Type MinValue = std::numeric_limits<Type>::min()>
struct Red : public Colour<Type, MaxValue, MinValue>
{
  Red()
      : Colour<Type, MaxValue, MinValue>{MaxValue, MinValue, MinValue} {}
};
template<
    class Type    = unsigned char,
    Type MaxValue = std::numeric_limits<Type>::max(),
    Type MinValue = std::numeric_limits<Type>::min()>
struct Green : public Colour<Type, MaxValue, MinValue>
{
  Green()
      : Colour<Type, MaxValue, MinValue>{MinValue, MaxValue, MinValue} {}
};
template<
    class Type    = unsigned char,
    Type MaxValue = std::numeric_limits<Type>::max(),
    Type MinValue = std::numeric_limits<Type>::min()>
struct Blue : public Colour<Type, MaxValue, MinValue>
{
  Blue()
      : Colour<Type, MaxValue, MinValue>{MinValue, MinValue, MaxValue} {}
};
}  // namespace Colour
}  // namespace UD
