#pragma once
#include <algorithm>
#include <array>
#include <iterator>
#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <wx/wx.h>

#include <wx/glcanvas.h>

#include "libudph/Class/Interface.h"
#include "libudph/Container/Container.h"
#include "libudph/Container/Iterator.h"
#include "libudph/Math/udMath.h"
#include "libudph/Math/Tensor.h"

#include "Colour.h"
#include "Logger.h"

namespace UD::World
{
using DEFAULT_TYPE = double;

template<class Type, Math::Type::ULong Dimension>
struct Point : public Tensor::Vector<Dimension, Type>
{
 private:
  using Base = Tensor::Vector<Dimension, Type>;

 public:
  using Colour = Colour::Colour<>;

  Colour _colour = {};

  [[nodiscard]] auto colour()  //
      -> Colour&
  {
    return _colour;
  }
  [[nodiscard]] auto colour() const  //
      -> const Colour&
  {
    return _colour;
  }
  virtual ~Point()        = default;
  Point()                 = default;
  Point(const Point&)     = default;
  Point(Point&&) noexcept = default;
  auto operator=(const Point&) -> Point& = default;
  auto operator=(Point&&) noexcept -> Point& = default;

  explicit Point(const Base& vector) : Base{vector} {}
  Point(const Base& vector, const Colour& colour)
      : Base{vector}
      , _colour{colour}
  {
  }
};
template<class Type, UD::Math::Type::ULong Dimension>
struct Mesh;
template<class Type, Math::Type::ULong Dimension, class _Leaf>
struct UD::Interface::Traits::InheritDescriptor<Mesh<Type, Dimension>, _Leaf>
    : UD::Interface::Traits::Tag<Mesh<Type, Dimension>,
                                 UD::Interface::Modifier::Tag::Abstract>
{
};
template<class Type, Math::Type::ULong Dimension>
struct Mesh : public UD::Interface::Interface<Mesh<Type, Dimension>>
{
  using Point    = Point<Type, Dimension>;
  using Vector   = Tensor::Vector<Dimension, Type>;
  using Iterator = Container::Iterator::GenericContiguousIterator<Point>;

  virtual ~Mesh()       = default;
  Mesh()                = default;
  Mesh(const Mesh&)     = default;
  Mesh(Mesh&&) noexcept = default;
  auto operator=(const Mesh&) -> Mesh& = default;
  auto operator=(Mesh&&) noexcept -> Mesh& = default;

  [[nodiscard]] virtual auto at(const Math::Type::ULong& index)  //
      -> Point& = 0;
  [[nodiscard]] virtual auto at(const Math::Type::ULong& index) const  //
      -> const Point& = 0;
  [[nodiscard]] virtual auto size() const  //
      -> Math::Type::ULong = 0;
  [[nodiscard]] virtual auto begin()  //
      -> Iterator = 0;
  [[nodiscard]] virtual auto end()  //
      -> Iterator = 0;
  // virtual const Iterator begin() const = 0;
  // virtual const Iterator end() const = 0;
  virtual auto operator[](std::size_t index)  //
      -> Point& = 0;
  virtual auto operator[](std::size_t index) const  //
      -> const Point& = 0;
  auto MinX()  //
      -> Type
  {
    if (size() == 0)
      return 0;
    auto min = at(0).x();
    for (Math::Type::ULong i = 1; i < size(); ++i)
      if (at(i).x() < min)
        min = at(i).x();
    return min;
  }
  auto MinXWithin(double bottom, double top)  //
      -> Type
  {
    if (size() == 0)
      return 0;
    auto min = at(0).x();
    for (Math::Type::ULong i = 1; i < size(); ++i)
      if (at(i).x() < min && at(i).y() >= bottom && at(i).y() <= top)
        min = at(i).x();
    return min;
  }
  auto MinY()  //
      -> Type
  {
    if (size() == 0)
      return 0;
    auto min = at(0).y();
    for (Math::Type::ULong i = 1; i < size(); ++i)
      if (at(i).y() < min)
        min = at(i).y();
    return min;
  }
  auto MinYWithin(double left, double right)  //
      -> Type
  {
    if (size() == 0)
      return 0;
    auto min = at(0).y();
    for (Math::Type::ULong i = 1; i < size(); ++i)
      if (at(i).y() < min && at(i).x() >= left && at(i).x() <= right)
        min = at(i).y();
    return min;
  }
  auto MaxX()  //
      -> Type
  {
    if (size() == 0)
      return 0;
    auto max = at(0).x();
    for (Math::Type::ULong i = 1; i < size(); ++i)
      if (at(i).x() > max)
        max = at(i).x();
    return max;
  }
  auto MaxXWithin(double bottom, double top)  //
      -> Type
  {
    if (size() == 0)
      return 0;
    auto max = at(0).x();
    for (Math::Type::ULong i = 1; i < size(); ++i)
      if (at(i).x() > max && at(i).y() >= bottom && at(i).y() <= top)
        max = at(i).x();
    return max;
  }
  auto MaxY()  //
      -> Type
  {
    if (size() == 0)
      return 0;
    auto max = at(0).y();
    for (Math::Type::ULong i = 1; i < size(); ++i)
      if (at(i).y() > max)
        max = at(i).y();
    return max;
  }
  auto MaxYWithin(double left = MinX(), double right = MaxX())  //
      -> Type
  {
    if (size() == 0)
      return 0;
    auto max = at(0).y();
    for (Math::Type::ULong i = 1; i < size(); ++i)
      if (at(i).y() > max && at(i).x() >= left && at(i).x() <= right)
        max = at(i).y();
    return max;
  }
};
template<class Type, Math::Type::ULong Dimension>
struct VariableMesh;
template<class Type, Math::Type::ULong Dimension, class _Leaf>
struct UD::Interface::Traits::InheritDescriptor<VariableMesh<Type, Dimension>, _Leaf>
    : public Inherit<_Leaf, Mesh<Type, Dimension>>
{
};
template<class Type, Math::Type::ULong Dimension>
struct VariableMesh
    : public UD::Interface::Interface<VariableMesh<Type, Dimension>>
{
 private:
  using MeshBase = Mesh<Type, Dimension>;

 public:
  using Point  = typename MeshBase::Point;
  using Vector = typename MeshBase::Vector;

  std::vector<Point> _points = {};

  ~VariableMesh() override              = default;
  VariableMesh()                        = default;
  VariableMesh(const VariableMesh&)     = default;
  VariableMesh(VariableMesh&&) noexcept = default;
  auto operator=(const VariableMesh&) -> VariableMesh& = default;
  auto operator=(VariableMesh&&) noexcept -> VariableMesh& = default;

  explicit VariableMesh(const std::vector<Point>& points) : _points{points} {}

  auto at(const Math::Type::ULong& index)  //
      -> Point& override
  {
    return _points.at(index);
  }
  auto at(const Math::Type::ULong& index) const  //
      -> const Point& override
  {
    return _points.at(index);
  }
  auto operator[](std::size_t index)  //
      -> typename Point& override
  {
    return _points.operator[](index);
  }
  auto operator[](std::size_t index) const  //
      -> const typename Point& override
  {
    return _points.operator[](index);
  }
  [[nodiscard]] auto size() const  //
      -> Math::Type::ULong override
  {
    return _points.size();
  }
  auto begin()  //
      -> typename MeshBase::Iterator override
  {
    return _points.begin();
  }
  auto end()  //
      -> typename MeshBase::Iterator override
  {
    return _points.end();
  }
  void push_back(const Point& point)
  {
    _points.push_back(point);
  }
  // const typename MeshBase::Iterator begin() const override {
  //  return TensorBase::begin();
  //}
  // const typename MeshBase::Iterator end() const override {
  //  return TensorBase::end();
  //}
};

template<class Type, Math::Type::ULong Dimension, Math::Type::ULong CountPoints>
struct StaticMesh;
template<class Type,
         Math::Type::ULong Dimension,
         Math::Type::ULong CountPoints,
         class _Leaf>
struct UD::Interface::Traits::
    InheritDescriptor<StaticMesh<Type, Dimension, CountPoints>, _Leaf>
    : public Inherit<_Leaf,
                     Mesh<Type, Dimension>,
                     Tensor::Vector<CountPoints,
                                    typename Mesh<Type, Dimension>::Point,
                                    Type>>
{
};
template<class Type, Math::Type::ULong Dimension, Math::Type::ULong CountPoints>
struct StaticMesh
    : public UD::Interface::Interface<StaticMesh<Type, Dimension, CountPoints>>
{
 private:
  using MeshBase = Mesh<Type, Dimension>;

 public:
  using Point      = typename MeshBase::Point;
  using Vector     = typename MeshBase::Vector;
  using TensorBase = Tensor::Vector<CountPoints, Point, Type>;

  ~StaticMesh() override            = default;
  StaticMesh()                      = default;
  StaticMesh(const StaticMesh&)     = default;
  StaticMesh(StaticMesh&&) noexcept = default;
  auto operator=(const StaticMesh&) -> StaticMesh& = default;
  auto operator=(StaticMesh&&) noexcept -> StaticMesh& = default;

  template<class... Ts>
  StaticMesh(Ts... ts) : TensorBase{ts...}
  {
  }

  Point& at(const Math::Type::ULong& index) override
  {
    return TensorBase::at(index);
  }
  const Point& at(const Math::Type::ULong& index) const override
  {
    return TensorBase::at(index);
  }
  typename Point& operator[](std::size_t index) override
  {
    return TensorBase::operator[](index);
  }
  const typename Point& operator[](std::size_t index) const override
  {
    return TensorBase::operator[](index);
  }
  Math::Type::ULong size() const override
  {
    return TensorBase::size();
  }
  typename MeshBase::Iterator begin() override
  {
    return TensorBase::begin();
  }
  typename MeshBase::Iterator end() override
  {
    return TensorBase::end();
  }
  // const typename MeshBase::Iterator begin() const override {
  //  return TensorBase::begin();
  //}
  // const typename MeshBase::Iterator end() const override {
  //  return TensorBase::end();
  //}
  // auto rbegin() { return _points.rbegin(); }
  // auto rend() { return _points.rend(); }
  // auto cbegin() { return _points.cbegin(); }
  // auto cend() { return _points.cend(); }
  // auto crbegin() { return _points.crbegin(); }
  // auto crend() { return _points.crend(); }

  // Mesh& Transform(std::function<Type(const Type&)> op) override {
  //  for (auto& point : _points) {
  //  }
  //}
  // Mesh& Transform(std::function<Type(const Type&, const Type&)> op,
  //                const Type& rhs) override {}
  // Mesh& Transform(std::function<Point(const Point&)> op) override {}
  // Mesh& Transform(std::function<Point(const Point&, const Point&)> op,
  //                const Point& rhs) override {}
};
template<class Type, Math::Type::ULong Dimension>
struct Rectangle;

template<class Type, Math::Type::ULong Dimension, class _Leaf>
struct UD::Interface::Traits::InheritDescriptor<Rectangle<Type, Dimension>, _Leaf>
    : UD::Interface::Traits::Inherit<_Leaf, StaticMesh<Type, Dimension, 4>>
{
};
template<class Type, Math::Type::ULong Dimension>
struct Rectangle : public UD::Interface::Interface<Rectangle<Type, Dimension>>
{
 private:
  using Base = StaticMesh<Type, Dimension, 4>;

 public:
  using Point    = typename Base::Point;
  using Vector   = typename Base::Vector;
  using Vector2D = Tensor::Vector<2, Type>;
  using Colour   = typename Point::Colour;

  // using Base::Base;
  Rectangle(const Vector&                             corner,
            const Vector&                             opposite_corner,
            const Container::Array::Array<Colour, 4>& colours
            = {UD::Colour::Black(),
               UD::Colour::Black(),
               UD::Colour::Black(),
               UD::Colour::Black()},
            const Vector& up = {0.0, 1.0})
  {
    Type   parts   = up.Sum();
    Vector diff    = opposite_corner - corner;
    Vector up_part = (diff * up) / parts;
    if (corner.x() <= opposite_corner.x())
    {
      if (corner.y() <= opposite_corner.y())
      {
        bottom_left()          = Point{corner};
        bottom_left().colour() = colours[0];
        top_right()            = Point{opposite_corner};
        top_right().colour()   = colours[2];

        top_left()              = bottom_left() + up_part;
        top_left().colour()     = colours[1];
        bottom_right()          = top_right() - up_part;
        bottom_right().colour() = colours[3];
      }
      else
      {
        top_left()              = Point{corner};
        top_left().colour()     = colours[0];
        bottom_right()          = Point{opposite_corner};
        bottom_right().colour() = colours[2];

        top_right()            = Point{bottom_right() + up_part};
        top_right().colour()   = colours[1];
        bottom_left()          = Point{top_left() - up_part};
        bottom_left().colour() = colours[3];
      }
    }
    else
    {
      *this = Rectangle{opposite_corner,
                        corner,
                        {colours[2], colours[3], colours[0], colours[1]},
                        -up};
    }
  }

  Point& bottom_left()
  {
    return (*this)[0];
  }
  Point& top_left()
  {
    return (*this)[1];
  }
  Point& top_right()
  {
    return (*this)[2];
  }
  Point& bottom_right()
  {
    return this->operator[](3);
  }
  Type left()
  {
    return bottom_left().x();
  }
  void left(Type value)
  {
    bottom_left().x() = value;
    top_left().x()    = value;
  }
  Type right()
  {
    return bottom_right().x();
  }
  void right(Type value)
  {
    bottom_right().x() = value;
    top_right().x()    = value;
  }
  Type top()
  {
    return top_left().y();
  }
  void top(Type value)
  {
    top_left().y()  = value;
    top_right().y() = value;
  }
  Type bottom()
  {
    return bottom_left().y();
  }
  void bottom(Type value)
  {
    bottom_left().y()  = value;
    bottom_right().y() = value;
  }
  Type width()
  {
    return right() - left();
  }
  void width(Type value)
  {
    right(left() + value);
  }
  Type height()
  {
    return top() - bottom();
  }
  void height(Type value)
  {
    top(bottom() + value);
  }

  Vector center()
  {
    return {left() + width() / 2, bottom() + height() / 2};
  }
  void ScaleWidth(const Type& factor)
  {
    auto d = width();
    auto c = left() + d / 2;
    auto w = (d * factor) / 2;
    left(c - w);
    right(c + w);
  }
  void ScaleHeight(const Type& factor)
  {
    auto d = height();
    auto c = bottom() + d / 2;
    auto w = d * factor;
    bottom(c - w);
    top(c + w);
  }
  void Scale(const Vector& factors)
  {
    ScaleWidth(factors.x());
    ScaleHeight(factors.y());
  }
  void Scale(const Type& factor)
  {
    ScaleWidth(factor);
    ScaleHeight(factor);
  }
};

template<class Type, Math::Type::ULong Dimension>
struct Line : public StaticMesh<Type, Dimension, 3>
{
 private:
  using Base = StaticMesh<Type, Dimension, 3>;

 public:
  using Point = typename Base::Point;
  using Base::Base;
  Line(Point a, Point b)
  {
    this->at(0) = a;
    this->at(1) = (a + b) / 2;
    this->at(2) = b;
  }
};
template<class Type, Math::Type::ULong Dimension>
struct VerticalLine : public Line<Type, Dimension>
{
 private:
  using Base = Line<Type, Dimension>;

 public:
  using Point = typename Base::Point;
  VerticalLine(Point a, Type height)
      : Base{{a - Point{{0.0, height}}}, {a + Point{{0.0, height}}}}
  {
  }
};
template<class Type, Math::Type::ULong Dimension>
struct HorizontalLine : public Line<Type, Dimension>
{
 private:
  using Base = Line<Type, Dimension>;

 public:
  using Point = typename Base::Point;
  HorizontalLine(Point a, Type width)
      : Base{{a - Point{{width, 0.0}}}, {a + Point{{width, 0.0}}}}
  {
  }
};

template<class _Type, Math::Type::ULong Dimension>
struct World
{
  using Type   = _Type;
  using Point  = Point<Type, Dimension>;
  using Mesh   = Mesh<Type, Dimension>;
  using Vector = Tensor::Vector<Dimension, Type>;
  template<Math::Type::ULong Size = Dimension>
  using Matrix = Tensor::Matrix<Size, Dimension, Type>;
  template<Math::Type::ULong Size = Dimension>
  using PointMatrix = Tensor::Tensor<Size, Point>;

  using Face      = PointMatrix<Dimension>;
  using Simplex   = PointMatrix<Dimension + 1>;
  using ID        = Math::Type::ULong;
  using Rectangle = Rectangle<Type, Dimension>;

  using Size = Vector;

  // static std::vector<Vector> ConvexHull_JarvisMarch(
  //    std::vector<Vector> points) {
  //  std::vector<Vector> hull;
  //  Vector first = points.begin();
  //  Vector next = first;
  //  hull.push_back(*points.begin());
  //}
  static constexpr Math::Type::ULong dimension()
  {
    return Dimension;
  }
  static Type Slope(Vector a, Vector b)
  {
    b -= a;
    return b[1] /= b[0];
  }
  static Type InterceptY(Vector vec, Type slope)
  {
    return vec[1] - slope * vec[0];
  }
  static Type InterceptY(Vector a, Vector b)
  {
    return InterceptY(a, Slope(a, b));
  }
  // struct IDrawable {
  //  virtual ~IDrawable() {}
  //  void Draw() = 0;
  //};
  struct Object
  {
    Vector _position = {};

    friend void swap(Object& a, Object& b) noexcept
    {
      using std::swap;
      swap(a._position, b._position);
    }
    virtual ~Object()     = default;
    Object()              = default;
    Object(const Object&) = default;
    Object(Object&&)      = default;
    auto operator=(const Object&) -> Object& = default;
    auto operator=(Object&&) noexcept -> Object& = default;

    Object(Vector position) : _position(position) {}
    Vector& position()
    {
      return _position;
    }
    // const Vector& position() const { return _position; }
  };
  struct Entity : public Object
  {
    std::unique_ptr<Mesh> _mesh = nullptr;

    friend void swap(Entity& a, Entity& b) noexcept
    {
      using std::swap;
      swap(static_cast<Object&>(a), static_cast<Object&>(b));
      swap(a._mesh, b._mesh);
    }
    virtual ~Entity() = default;
    Entity(const Entity& entity) : _mesh{_mesh->Clone()}, Object(entity) {}
    Entity(Entity&& entity) noexcept
    {
      swap(*this, entity);
    }
    auto operator=(Entity entity)
    {
      swap(*this, entity);
      return *this;
    };

    Entity(const Mesh& mesh) : _mesh{mesh.Clone()} {}
    Entity(const Vector& position, const Mesh& mesh)
        : Object{position}
        , _mesh{mesh.Clone()}
    {
    }
    // Entity(Entity&& entity) noexcept : Object{std::move(entity)} {
    //  std::swap(_mesh, entity._mesh);
    //}
    // Entity& operator=(Entity entity) {
    //  Object::operator=(entity);
    //  std::swap(_mesh, entity._mesh);
    //  return *this;
    //}
    Mesh& mesh()
    {
      return *_mesh;
    }
    std::unique_ptr<Mesh> ExchangeMesh(const Mesh& mesh)
    {
      return std::exchange(_mesh, mesh.Clone());
    }
    // const Mesh& mesh() const { return *_mesh; }
    void Draw() /* const*/
    {
      DrawImpl();
    }
    void DrawPoint(const Point& point) /* const*/
    {
      glColor4ub(point.colour().r(),
                 point.colour().g(),
                 point.colour().b(),
                 point.colour().a());
      if constexpr (Dimension == 2)
      {
        glVertex2d(point.x(), point.y());
      }
      else if constexpr (Dimension == 3)
      {
        glVertex3d(point.x(), point.y(), point.z());
      }
    }
    template<Math::Type::ULong Size>
    void DrawPoints(Tensor::Tensor<Size, Point> points)
    {
      GLdouble prev_size;
      glGetDoublev(GL_POINT_SIZE, &prev_size);

      glPointSize(5);
      glBegin(GL_POINTS);
      for (auto& point : points)
        DrawPoint(point + this->position());
      glEnd();
      glPointSize(prev_size);
    }
    void DrawPoints()
    {
      GLdouble prev_size;
      glGetDoublev(GL_POINT_SIZE, &prev_size);

      glPointSize(5);
      glBegin(GL_POINTS);
      for (auto& point : mesh())
        DrawPoint(point + this->position());
      glEnd();
      glPointSize(prev_size);
    }
    template<Math::Type::ULong Size>
    void DrawLineStrip(Tensor::Tensor<Size, Point> points)
    {
      glBegin(GL_LINE_STRIP);
      for (auto& point : points)
        DrawPoint(point + this->position());
      glEnd();
    }
    void DrawLineStrip()
    {
      glBegin(GL_LINE_STRIP);
      for (auto& point : mesh())
        DrawPoint(point + this->position());
      glEnd();
    }
    template<Math::Type::ULong Size>
    void DrawPolygon(Tensor::Tensor<Size, Point> points)
    {
      glBegin(GL_POLYGON);
      for (auto& point : points)
        DrawPoint(point + this->position());
      glEnd();
    }
    void ScaleTo(const Rectangle& rect)
    {
      ScaleXTo(rect.left(), rect.right());
      ScaleYTo(rect.bottom(), rect.top());
    }
    void ScaleXTo(double left,
                  double right,
                  double within_bottom = mesh().MinY(),
                  double within_top    = mesh().MaxY())
    {
      auto width      = right - left;
      auto mesh_width = mesh().MaxXWithin(within_bottom, within_top)
                      - mesh().MinXWithin(within_bottom, within_top);
      for (auto& p : mesh())
        p.y() = width * (p.y() / mesh_width);
    }
    void ScaleYTo(double bottom,
                  double top,
                  double within_left  = mesh().MinX(),
                  double within_right = mesh().MaxX())
    {
      auto height      = top - bottom;
      auto mesh_height = mesh().MaxYWithin(within_left, within_right)
                       - mesh().MinYWithin(within_left, within_right);
      for (auto& p : mesh())
        p.y() = height * (p.y() / mesh_height);
    }
    Vector center()
    {
      Vector center = {0.0, 0.0};
      for (auto& point : mesh())
        center += point;
      center /= mesh().size();
      return center + this->position();
    }
    void ReCenter()
    {
      Vector shift = center() - this->position();
      this->position() += shift;
      auto& m = mesh();
      for (auto& point : mesh())
        point -= shift;
    }
    void ReCenter(const Vector& vec)
    {
      auto   c               = center();
      Vector shift_to_center = c - this->position();
      Vector shift_to_vec    = vec - c;
      Vector shift           = shift_to_center + shift_to_vec;
      this->position() += shift;
      for (auto& point : mesh())
        point -= shift_to_center;
    }
    void ReCenterY(const Type& val)
    {
      ReCenter();
      this->position().y() += val - this->position().y();
    }
    void ReBaseY(Type bottom, Type top)
    {
      ReCenter();
      auto left  = _mesh->MinX() + this->position().x();
      auto right = _mesh->MaxX() + this->position().x();
      ReBaseHelper({{left, bottom}, {right, top}},
                   {{left, _mesh->MinY() + this->position().y()},
                    {right, _mesh->MaxY() + this->position().y()}});
    }
    void ReBaseX(Type left, Type right)
    {
      ReCenter();
      auto top    = _mesh->MaxY() + this->position().y();
      auto bottom = _mesh->MinY() + this->position().y();
      ReBaseHelper({{left, bottom}, {right, top}},
                   {{_mesh->MinX() + this->position().x(), bottom},
                    {_mesh->MaxX() + this->position().x(), top}});
    }
    // TODO: make const&
    void ReBase(Rectangle new_base)
    {
      ReCenter();
      auto left   = _mesh->MinX() + this->position().x();
      auto right  = _mesh->MaxX() + this->position().x();
      auto top    = _mesh->MaxY() + this->position().y();
      auto bottom = _mesh->MinY() + this->position().y();
      ReBaseHelper(new_base, {{left, bottom}, {right, top}});
    }
    void ReBase(Rectangle new_base, Rectangle old_base)
    {
      ReCenter();
      ReBaseHelper(new_base, old_base);
    }
    void ReBaseHelper(Rectangle new_base, Rectangle old_base)
    {
      Vector scale{1.0, 1.0};
      if (old_base.width() > 0)
        scale.x() = new_base.width() / old_base.width();
      if (old_base.height() > 0)
        scale.y() = new_base.height() / old_base.height();
      for (auto& point : *_mesh)
        point *= scale;
      this->position() = (this->position() - old_base.bottom_left()) * scale
                       + new_base.bottom_left();
    }
    bool IsInBoundingBox(const Vector& vec)
    {
      auto left   = _mesh->MinX() + this->position().x();
      auto right  = _mesh->MaxX() + this->position().x();
      auto top    = _mesh->MaxY() + this->position().y();
      auto bottom = _mesh->MinY() + this->position().y();
      return vec.x() > left && vec.x() < right && vec.y() > bottom
          && vec.y() < top;
    }

   protected:
    virtual void DrawImpl() /* const*/
    {
      if (mesh().size() == 1)
        DrawPoints();
      else
        DrawLineStrip();
    }
  };
  // template <class MeshType = Mesh>
  // struct MeshEntity : public Entity {
  //  MeshType _mesh;
  //  virtual ~MeshEntity() {}
  //  MeshEntity() {}
  //  MeshEntity(const MeshType& mesh) : _mesh{mesh} {}
  //  MeshType& mesh() { return _mesh; }

  //  virtual void Draw() { Draw(mesh()); }
  //};
  static bool Contains(const Rectangle& rectangle, const Vector& vector)
  {
    return rectangle.top_left()[0] >= vector[0]
        && rectangle.top_left()[1] >= vector[1]
        && rectangle.bottom_right()[0] <= vector[0]
        && rectangle.bottom_right()[1] <= vector[1];
  }
  template<Math::Type::ULong Size, Math::Type::ULong NumberOfPoints>
  static bool Intersect(const Rectangle& rectangle, const Mesh& mesh)
  {
    for (auto& point : mesh)
      if (Contains(rectangle, point))
        return true;
  }
  template<Math::Type::ULong Size, Math::Type::ULong NumberOfPoints>
  static bool Contains(const Rectangle& rectangle, const Mesh& mesh)
  {
    for (auto& point : mesh)
      if (!Contains(rectangle, point))
        return false;
  }
  struct View : public Entity
  {
   private:
    using Base = Rectangle;

   public:
    World* _world  = nullptr;
    Vector _offset = {};
    // TODO: Make world const& (requires implementing generic const iterators)
    virtual ~View() {}
    View(World& world, Vector position, Rectangle rectangle)
        : Entity{position, rectangle}
        , _world{&world}
    {
    }
    View(World& world, Rectangle rectangle)
        : View{world,
               rectangle.center(),
               Rectangle{rectangle.bottom_left() - rectangle.center(),
                         rectangle.top_left() - rectangle.center()}}
    {
    }

    void DrawPrep(Rectangle rect)
    {
      UpdateViewport(rect);
      ProjectionOrtho();
      wxColour clear_colour = *wxBLACK;
      glClearColor(clear_colour.Red(),
                   clear_colour.Green(),
                   clear_colour.Blue(),
                   clear_colour.Alpha());
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      glPointSize(1);
    }
    void Draw(Rectangle rect, std::unordered_set<ID> ids = {})
    {
      DrawPrep(rect);
      if (_world->entities().empty())
        return;
      if (ids.empty())
      {
        for (auto& entity : Visible())
        {
          _world->Get<Entity>(entity).Draw();
        }
      }
      else
      {
        for (auto& id : ids)
        {
          if (IsVisible(id))
            _world->Get<Entity>(id).Draw();
        }
      }
    }
    void DrawExcluding(Rectangle rect, std::unordered_set<ID> ids)
    {
      DrawPrep(rect);
      for (auto& entity_pair : this->_world->entities())
      {
        auto& entity = *entity_pair.second;
        auto& id     = entity_pair.first;
        if (!ids.contains(id))
        {
          if (IsVisible(entity))
            entity.Draw();
        }
      }
    }
    void UpdateViewport(Rectangle rect)
    {
      glEnable(GL_TEXTURE_2D);
      glEnable(GL_COLOR_MATERIAL);
      glEnable(GL_BLEND);
      glDisable(GL_DEPTH_TEST);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glViewport(rect.left(), rect.bottom(), rect.width(), rect.height());
    }
    void ProjectionOrtho()
    {
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glOrtho(left(), right(), bottom(), top(), 0, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
    }
    Rectangle& rect()
    {
      return dynamic_cast<Rectangle&>(this->mesh());
    }
    Type left()
    {
      return rect().left() + this->position().x();
    }
    void left(const Type& val)
    {
      rect().left(val - this->position().x());
    }
    Type right()
    {
      return rect().right() + this->position().x();
    }
    void right(const Type& val)
    {
      rect().right(val - this->position().x());
    }
    Type bottom()
    {
      return rect().bottom() + this->position().y();
    }
    void bottom(const Type& val)
    {
      rect().bottom(val - this->position().y());
    }
    Type top()
    {
      return rect().top() + this->position().y();
    }
    void top(const Type& val)
    {
      rect().top(val - this->position().y());
    }
    Vector& offset()
    {
      return _offset;
    }
    void Zoom(double amount)
    {
      rect().Scale(amount);
    }
    void ZoomWidth(double amount)
    {
      rect().ScaleWidth(amount);
    }
    void ZoomHeight(double amount)
    {
      rect().ScaleHeight(amount);
    }
    void Translate(Vector amount)
    {
      this->position() += amount;
    }
    void UpdateRect(Rectangle rect)
    {
      this->rect() = rect;
    }
    Vector FramePositionToWorldPosition(Rectangle frame, Vector position)
    {
      auto width_ratio  = this->rect().width() / frame.width();
      auto height_ratio = this->rect().height() / frame.height();
      position -= {frame.left(), frame.bottom()};
      position *= {width_ratio, height_ratio};
      position += {this->left(), this->bottom()};
      return position;
    }
    Vector WorldPositionToFramePosition(Rectangle frame, Vector position)
    {
      auto width_ratio  = frame.width() / this->rect().width();
      auto height_ratio = frame.height() / this->rect().height();
      position -= {this->left(), this->bottom()};
      position *= {width_ratio, height_ratio};
      position += {frame.left(), frame.bottom()};
      return position;
    }
    bool IsVisible(ID entity_id)
    {
      auto& entity = _world->Get<Entity>(entity_id);
      for (auto& point : entity.mesh())
      {
        auto t  = point + entity.position();
        auto l  = left();
        auto l1 = right();
        auto l2 = top();
        auto l3 = bottom();
        if (t.x() >= left() && t.x() <= right() && t.y() >= bottom()
            && t.y() <= top())
        {
          return true;
        }
      }
      return false;
    }
    virtual std::vector<ID> Visible()
    {
      std::vector<ID> entities;
      for (auto& entity : _world->entities())
      {
        if (IsVisible(entity.first))
        {
          entities.push_back(entity.first);
        }
      }
      return entities;
    }
    // TODO: replace with inheritance to redefine Visible criteria
    std::vector<ID> VisibleX(std::unordered_set<ID> excluding_ids = {})
    {
      std::vector<ID> entities;
      auto            left  = this->left();
      auto            right = this->right();
      auto            bot   = this->bottom();
      auto            top   = this->top();
      for (auto& entity : _world->entities())
      {
        if (excluding_ids.contains(entity.first))
          continue;
        for (auto& point : entity.second->mesh())
        {
          auto t = point + entity.second->position();
          if (t.x() >= left && t.x() <= right)
          {
            entities.push_back(entity.first);
            break;
          }
        }
      }
      return entities;
    }
    // TODO: replace with inheritance to redefine Visible criteria
    std::vector<Entity> VisibleY(std::unordered_set<ID> excluding_ids = {})
    {
      std::vector<Entity> entities;
      auto                left  = this->left();
      auto                right = this->right();
      auto                bot   = this->bottom();
      auto                top   = this->top();
      for (auto& entity : _world->entities())
      {
        if (excluding_ids.contains(entity.first))
          continue;
        for (auto& point : entity.second->mesh())
        {
          auto t = point + entity.second->position();
          if (t.y() >= bot && t.y() <= top)
          {
            entities.push_back(*entity.second);
            break;
          }
        }
      }
      return entities;
    }
    void Fit(std::vector<ID> entities, std::unordered_set<ID> ids = {})
    {
      FitX(entities, ids);
      FitY(entities, ids);
    }
    void FitX(std::vector<ID> entities, std::unordered_set<ID> ids = {})
    {
      if (entities.size() == 0)
        return;
      Math::Type::ULong start = 0;
      if (!ids.empty())
      {
        while (!ids.contains(entities[start]))
        {
          start++;
          if (start >= entities.size())
            return;
        }
      }
      auto entity    = &_world->Get<Entity>(entities[start]);
      auto max_width = entity->mesh().MaxX() + entity->position().x();
      auto min_width = entity->mesh().MinX() + entity->position().x();
      for (std::size_t i = start + 1; i < entities.size(); ++i)
      {
        if (!ids.empty() && !ids.contains(entities[i]))
          continue;
        entity = &_world->Get<Entity>(entities[i]);

        auto potential_max_width
            = entity->mesh().MaxX() + entity->position().x();
        auto potential_min_width
            = entity->mesh().MinX() + entity->position().x();
        if (potential_max_width > max_width)
          max_width = potential_max_width;
        if (potential_min_width < min_width)
          min_width = potential_min_width;
      }
      left(min_width);
      right(max_width);
      this->ReCenter();
    }
    void FitY(std::vector<ID> entities, std::unordered_set<ID> ids = {})
    {
      if (entities.size() == 0)
        return;
      Math::Type::ULong start = 0;
      if (!ids.empty())
      {
        while (!ids.contains(entities[start]))
        {
          start++;
          if (start >= entities.size())
            return;
        }
      }
      auto entity     = &_world->Get<Entity>(entities[start]);
      auto max_height = entity->mesh().MaxY() + entity->position().y();
      auto min_height = entity->mesh().MinY() + entity->position().y();
      for (std::size_t i = start + 1; i < entities.size(); ++i)
      {
        if (!ids.empty() && !ids.contains(entities[i]))
          continue;
        entity = &_world->Get<Entity>(entities[i]);
        auto potential_max_height
            = entity->mesh().MaxY() + entity->position().y();
        auto potential_min_height
            = entity->mesh().MinY() + entity->position().y();
        if (potential_max_height > max_height)
          max_height = potential_max_height;
        if (potential_min_height < min_height)
          min_height = potential_min_height;
      }
      bottom(min_height);
      top(max_height);
      this->ReCenter();
    }
    void FitExcluding(std::vector<ID> entities, std::unordered_set<ID> ids)
    {
      FitXExcluding(entities, ids);
      FitYExcluding(entities, ids);
    }
    void FitXExcluding(std::vector<ID> entities, std::unordered_set<ID> ids)
    {
      if (entities.size() == 0)
        return;
      Math::Type::ULong start = 0;
      while (ids.contains(entities[start]))
      {
        start++;
        if (start >= entities.size())
          return;
      }
      auto entity    = &_world->Get<Entity>(entities[start]);
      auto max_width = entity->mesh().MaxX() + entity->position().x();
      auto min_width = entity->mesh().MinX() + entity->position().x();
      for (std::size_t i = start + 1; i < entities.size(); ++i)
      {
        if (ids.contains(entities[i]))
          continue;
        entity = &_world->Get<Entity>(entities[i]);

        auto potential_max_width
            = entity->mesh().MaxX() + entity->position().x();
        auto potential_min_width
            = entity->mesh().MinX() + entity->position().x();
        if (potential_max_width > max_width)
          max_width = potential_max_width;
        if (potential_min_width < min_width)
          min_width = potential_min_width;
      }
      rect().left(min_width);
      rect().right(max_width);
    }
    void FitYExcluding(std::vector<ID> entities, std::unordered_set<ID> ids)
    {
      if (entities.size() == 0)
        return;
      Math::Type::ULong start = 0;
      while (ids.contains(entities[start]))
      {
        start++;
        if (start >= entities.size())
          return;
      }
      auto entity     = &_world->Get<Entity>(entities[start]);
      auto max_height = entity->mesh().MaxY() + entity->position().y();
      auto min_height = entity->mesh().MinY() + entity->position().y();
      for (std::size_t i = start + 1; i < entities.size(); ++i)
      {
        if (ids.contains(entities[i]))
          continue;
        entity = &_world->Get<Entity>(entities[i]);
        auto potential_max_height
            = entity->mesh().MaxY() + entity->position().y();
        auto potential_min_height
            = entity->mesh().MinY() + entity->position().y();
        if (potential_max_height > max_height)
          max_height = potential_max_height;
        if (potential_min_height < min_height)
          min_height = potential_min_height;
      }
      rect().bottom(min_height);
      rect().top(max_height);
    }
  };

  ID previous_id = 1;

  std::unordered_map<ID, std::unique_ptr<View>> views = {};

 private:
  std::unordered_map<ID, std::unique_ptr<Entity>> _entities = {};

 public:
  std::unordered_map<ID, std::unique_ptr<Entity>>& entities()
  {
    return _entities;
  }
  World() = default;

  template<class T, class... Ts>
  ID Add(Ts... ts)
  {
    ID id = NextID();
    entities().insert({id, std::make_unique<T>(ts...)});
    return id;
  }
  template<class T = View, class... Ts>
  requires std::derived_from<T, View> T& CreateView(Ts... ts)
  {
    ID id = NextID();
    views.insert({id, std::make_unique<View>(*this, ts...)});
    return Get<T>(id);
  }

  void Clear()
  {
    views.clear();
    entities().clear();
    previous_id = 1;
  }
  void ClearEntities()
  {
    entities().clear();
  }
  void ClearViews()
  {
    views.clear();
  }
  ID NextID()
  {
    ID start_id = previous_id++;

    while (Contains(previous_id))
    {
      previous_id++;
      if (previous_id == start_id)
        throw;
    }
    return previous_id;
  }

  template<class T = Object>
  bool Contains(ID id);
  template<>
  bool Contains<Object>(ID id)
  {
    return Contains<View>(id) || Contains<Entity>(id);
  }
  template<>
  bool Contains<View>(ID id)
  {
    return views.count(id);
  }
  template<>
  bool Contains<Entity>(ID id)
  {
    return entities().count(id);
  }
  template<class T>
  T& Get(ID id);
  template<>
  View& Get<View>(ID id)
  {
    return *views.at(id);
  }
  template<>
  Entity& Get<Entity>(ID id)
  {
    return *entities().at(id);
  }
  Type MaxX(const std::vector<ID>& ids)
  {
    if (ids.empty())
      return {};
    auto& temp_ent = Get<Entity>(*ids.begin());
    Type  max      = temp_ent.mesh().MaxX() + temp_ent.position().x();
    for (auto it = ids.begin() + 1; it != ids.end(); ++it)
    {
      auto& ent = Get<Entity>(*it);
      Type  val = ent.mesh().MaxX() + ent.position().x();
      if (val > max)
        max = val;
    }
    return max;
  }
  Type MinX(const std::vector<ID>& ids)
  {
    if (ids.empty())
      return {};
    auto& temp_ent = Get<Entity>(*ids.begin());
    Type  min      = temp_ent.mesh().MinX() + temp_ent.position().x();
    for (auto it = ids.begin() + 1; it != ids.end(); ++it)
    {
      auto& ent = Get<Entity>(*it);
      Type  val = ent.mesh().MinX() + ent.position().x();
      if (val < min)
        min = val;
    }
    return min;
  }
  Type MaxY(const std::vector<ID>& ids)
  {
    if (ids.empty())
      return {};
    auto& temp_ent = Get<Entity>(*ids.begin());
    Type  max      = temp_ent.mesh().MaxY() + temp_ent.position().y();
    for (auto it = ids.begin() + 1; it != ids.end(); ++it)
    {
      auto& ent = Get<Entity>(*it);
      Type  val = ent.mesh().MaxY() + ent.position().y();
      if (val > max)
        max = val;
    }
    return max;
  }
  Type MinY(const std::vector<ID>& ids)
  {
    if (ids.empty())
      return {};
    auto& temp_ent = Get<Entity>(*ids.begin());
    Type  min      = temp_ent.mesh().MinY() + temp_ent.position().y();
    for (auto it = ids.begin() + 1; it != ids.end(); ++it)
    {
      auto& ent = Get<Entity>(*it);
      Type  val = ent.mesh().MinY() + ent.position().y();
      if (val < min)
        min = val;
    }
    return min;
  }
  ID Intersect(const Vector& vec)
  {
    for (auto& [key, val] : entities())
      if (val->IsInBoundingBox(vec))
        return key;
    return 0;
  }
};  // namespace World
}  // namespace UD::World