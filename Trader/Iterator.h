#pragma once
#include <concepts>
#include <cstddef>
#include <iterator>
#include <memory>
#include <utility>

namespace UD
{
namespace Container
{
namespace Iterator
{
template<class _ValueType,
         class _PointerType    = _ValueType*,
         class _ReferenceType  = _ValueType&,
         class _DifferenceType = std::ptrdiff_t>
struct ContiguousIteratorBase
{
  using iterator_category = std::contiguous_iterator_tag;
  using difference_type   = _DifferenceType;
  using value_type        = _ValueType;
  using pointer           = _PointerType;
  using reference         = _ReferenceType;

  virtual ~ContiguousIteratorBase()                         = default;
  ContiguousIteratorBase()                                  = default;
  ContiguousIteratorBase(const ContiguousIteratorBase&)     = default;
  ContiguousIteratorBase(ContiguousIteratorBase&&) noexcept = default;
  auto operator                  =(const ContiguousIteratorBase&)
      -> ContiguousIteratorBase& = default;
  auto operator                  =(ContiguousIteratorBase&&) noexcept
      -> ContiguousIteratorBase& = default;

  virtual reference               operator*() const          = 0;
  virtual pointer                 operator->()               = 0;
  virtual ContiguousIteratorBase& operator++()               = 0;
  virtual ContiguousIteratorBase* operator++(int)            = 0;
  virtual ContiguousIteratorBase& operator--()               = 0;
  virtual ContiguousIteratorBase* operator--(int)            = 0;
  virtual ContiguousIteratorBase& operator+=(const int& ind) = 0;
  virtual ContiguousIteratorBase& operator-=(const int& ind) = 0;
  virtual reference               operator[](const int& ind) = 0;

 protected:
  virtual bool Equal(const ContiguousIteratorBase& rhs) const    = 0;
  virtual bool NotEqual(const ContiguousIteratorBase& rhs) const = 0;
  virtual ContiguousIteratorBase* Add(const int& ind) const      = 0;
  virtual ContiguousIteratorBase* Subtract(const int& ind) const = 0;

 public:
  friend bool operator==(const ContiguousIteratorBase& lhs,
                         const ContiguousIteratorBase& rhs)
  {
    return lhs.Equal(rhs);
  };
  friend bool operator!=(const ContiguousIteratorBase& lhs,
                         const ContiguousIteratorBase& rhs)
  {
    return lhs.NotEqual(rhs);
  };
  friend ContiguousIteratorBase* operator+(const ContiguousIteratorBase& it,
                                           const int&                    ind)
  {
    return it.Add(ind);
  }
  friend ContiguousIteratorBase* operator-(const ContiguousIteratorBase& it,
                                           const int&                    ind)
  {
    return it.Subtract(ind);
  }
};
template<class T>
struct Iterator;
template<class T>
  requires std::contiguous_iterator<T>
struct Iterator<T>
    : public ContiguousIteratorBase<typename T::value_type,
                                    typename T::pointer,
                                    typename T::reference,
                                    typename T::difference_type>
{
 private:
  using Base = ContiguousIteratorBase<typename T::value_type,
                                      typename T::pointer,
                                      typename T::reference,
                                      typename T::difference_type>;

  T _iterator = {};

 public:
  Iterator(const T& iterator) : _iterator{iterator} {}
  typename Base::reference operator*() const override
  {
    return *_iterator;
  }
  typename Base::pointer operator->() override
  {
    return &(*_iterator);
  }

  Iterator& operator++() override
  {
    _iterator++;
    return *this;
  }
  Iterator* operator++(int) override
  {
    Iterator* tmp = new Iterator(*this);
    ++(*this);
    return tmp;
  }
  Iterator& operator--() override
  {
    _iterator--;
    return *this;
  }
  Iterator* operator--(int) override
  {
    Iterator* tmp = new Iterator(*this);
    --(*this);
    return tmp;
  }
  Iterator& operator+=(const int& ind) override
  {
    _iterator += ind;
    return *this;
  }
  Iterator& operator-=(const int& ind)
  {
    _iterator -= ind;
    return *this;
  }
  // TODO: This is broken...
  typename Base::reference operator[](const int& ind)
  {
    return *((*this) += ind);
  }

 protected:
  bool Equal(const Base& rhs) const override
  {
    return _iterator == static_cast<const Iterator&>(rhs)._iterator;
  }
  bool NotEqual(const Base& rhs) const override
  {
    return _iterator != static_cast<const Iterator&>(rhs)._iterator;
  }
  Iterator* Add(const int& ind) const override
  {
    Iterator* tmp = new Iterator(*this);
    *tmp += ind;
    return tmp;
  }
  Iterator* Subtract(const int& ind) const override
  {
    Iterator* tmp = new Iterator(*this);
    *tmp -= ind;
    return tmp;
  }
};
template<class _ValueType,
         class _PointerType    = _ValueType*,
         class _ReferenceType  = _ValueType&,
         class _DifferenceType = std::ptrdiff_t>
struct GenericContiguousIterator
{
  using iterator_category = std::contiguous_iterator_tag;
  using difference_type   = _DifferenceType;
  using value_type        = _ValueType;
  using pointer           = _PointerType;
  using reference         = _ReferenceType;

 private:
  using IterType
      = ContiguousIteratorBase<value_type, pointer, reference, difference_type>;
  IterType* _iterator = nullptr;
  GenericContiguousIterator(IterType* iterator)
  {
    if (_iterator)
      delete _iterator;
    _iterator = iterator;
  }

 public:
  GenericContiguousIterator(const GenericContiguousIterator& other)
      : GenericContiguousIterator(
          static_cast<IterType*>(new ContiguousIteratorBase(other._iterator)))
  {
  }
  GenericContiguousIterator(GenericContiguousIterator&& other) noexcept
  {
    std::swap(_iterator, other._iterator);
  }
  GenericContiguousIterator& operator=(GenericContiguousIterator other) noexcept
  {
    std::swap(_iterator, other._iterator);
    return *this;
  }
  template<class T>
    requires requires(T t)
    {
      std::contiguous_iterator<T>;
      std::same_as<typename T::iterator_category, iterator_category>;
      std::same_as<typename T::difference_type, difference_type>;
      std::same_as<typename T::value_type, value_type>;
      std::same_as<typename T::pointer, pointer>;
      std::same_as<typename T::reference, reference>;
      std::derived_from<Iterator<T>, IterType>;
    }
  GenericContiguousIterator(const T& iterator)
      : GenericContiguousIterator{
          static_cast<IterType*>(new Iterator<T>(iterator))}
  {
  }
  ~GenericContiguousIterator()
  {
    delete _iterator;
  }

  reference operator*() const
  {
    return **_iterator;
  }
  pointer operator->()
  {
    return &(**_iterator);
  }

  GenericContiguousIterator& operator++()
  {
    ++(*_iterator);
    return *this;
  }
  GenericContiguousIterator operator++(int)
  {
    auto ret = GenericContiguousIterator(*_iterator);
    (*_iterator)++;
    return ret;
  }
  GenericContiguousIterator& operator--()
  {
    --(*_iterator);
    return *this;
  }
  GenericContiguousIterator operator--(int)
  {
    auto ret = GenericContiguousIterator(*_iterator);
    (*_iterator)--;
    return ret;
  }
  GenericContiguousIterator& operator+=(const int& ind)
  {
    (*_iterator) += ind;
    return *this;
  }
  friend GenericContiguousIterator operator+(
      const GenericContiguousIterator& it,
      const int&                       ind)
  {
    return GenericContiguousIterator((*(it._iterator)) + ind);
  }
  GenericContiguousIterator& operator-=(const int& ind)
  {
    (*_iterator) -= ind;
    return *this;
  }
  friend GenericContiguousIterator operator-(
      const GenericContiguousIterator& it,
      const int&                       ind)
  {
    return GenericContiguousIterator((*(it._iterator)) - ind);
  }

  reference operator[](const int& ind)
  {
    return *this += ind;
  }

  friend bool operator==(const GenericContiguousIterator& a,
                         const GenericContiguousIterator& b)
  {
    return *(a._iterator) == *(b._iterator);
  };
  friend bool operator!=(const GenericContiguousIterator& a,
                         const GenericContiguousIterator& b)
  {
    return *(a._iterator) != *(b._iterator);
  };
};
}  // namespace Iterator
}  // namespace Container
}  // namespace UD