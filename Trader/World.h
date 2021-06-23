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

#include "Colour.h"
#include "Container.h"
#include "Iterator.h"
#include "Logger.h"
#include "Math.h"
#include "Tensor.h"

namespace UD {
namespace World {
using DEFAULT_TYPE = double;

template <class Type, Math::Type::ULong Dimension>
struct Point : public Tensor::Vector<Dimension, Type> {
 private:
  using Base = Tensor::Vector<Dimension, Type>;

 public:
  using Colour = Colour::Colour<>;
  Colour _colour;
  Colour& colour() { return _colour; }
  const Colour& colour() const { return _colour; }
  Point() = default;
  Point(const Point&) = default;
  Point(Point&&) = default;
  Point& operator=(const Point&) = default;
  Point& operator=(Point&&) = default;

  Point(const Base& vector) : Base{vector}, _colour{} {}
  Point(const Base& vector, Colour colour) : Base{vector}, _colour{colour} {}
};
template <class Type, Math::Type::ULong Dimension>
struct Mesh {
  using Point = Point<Type, Dimension>;
  using Vector = Tensor::Vector<Dimension, Type>;
  using Iterator = Container::Iterator::GenericContiguousIterator<Point>;

  virtual Point& at(const Math::Type::ULong& index) = 0;
  virtual const Point& at(const Math::Type::ULong& index) const = 0;
  virtual Math::Type::ULong size() const = 0;
  virtual Iterator begin() = 0;
  virtual Iterator end() = 0;
  // virtual const Iterator begin() const = 0;
  // virtual const Iterator end() const = 0;
  virtual Mesh* Clone() const = 0;
  virtual typename Point& operator[](std::size_t index) = 0;
  virtual const typename Point& operator[](std::size_t index) const = 0;

  Type MinX() {
    if (size() == 0) return 0;
    auto min = at(0).x();
    for (Math::Type::ULong i = 1; i < size(); ++i)
      if (at(i).x() < min) min = at(i).x();
    return min;
  }
  Type MinY() {
    if (size() == 0) return 0;
    auto min = at(0).y();
    for (Math::Type::ULong i = 1; i < size(); ++i)
      if (at(i).y() < min) min = at(i).y();
    return min;
  }
  Type MaxX() {
    if (size() == 0) return 0;
    auto max = at(0).x();
    for (Math::Type::ULong i = 1; i < size(); ++i)
      if (at(i).x() > max) max = at(i).x();
    return max;
  }
  Type MaxY() {
    if (size() == 0) return 0;
    auto max = at(0).y();
    for (Math::Type::ULong i = 1; i < size(); ++i)
      if (at(i).y() > max) max = at(i).y();
    return max;
  }

  // auto rbegin() { return _points.rbegin(); }
  // auto rend() { return _points.rend(); }
  // auto cbegin() { return _points.cbegin(); }
  // auto cend() { return _points.cend(); }
  // auto crbegin() { return _points.crbegin(); }
  // auto crend() { return _points.crend(); }
  // public:
  // virtual Mesh& Transform(std::function<Type(const Type&)> op) = 0;
  // virtual Mesh& Transform(std::function<Type(const Type&, const Type&)> op,
  //                         const Type& rhs) = 0;
  // virtual Mesh& Transform(std::function<Point(const Point&)> op) = 0;
  // virtual Mesh& Transform(std::function<Point(const Point&, const Point&)>
  // op,
  //                         const Point& rhs) = 0;

  // private:
  // static Point TranslateHelper(const Point& lhs, const Point& rhs) {
  //   return lhs + rhs;
  // }
  // static Point ScaleHelper(const Point& lhs, const Point& rhs) {
  //   return lhs * rhs;
  // }
  // static Type ScaleHelper(const Type& lhs, const Type& rhs) {
  //   return lhs * rhs;
  // }

  // public:
  // Mesh& Translate(Vector position) { Transform(&TranslateHelper, position);
  // }; Mesh& Scale(Vector factor) { return Transform(&ScaleHelper, factor); }
  // Mesh& Scale(Type factor) { return Transform(&ScaleHelper, factor); }
};
template <class Type, Math::Type::ULong Dimension>
struct VariableMesh : public Mesh<Type, Dimension> {
 private:
  using MeshBase = Mesh<Type, Dimension>;

 public:
  using Point = typename MeshBase::Point;
  using Vector = typename MeshBase::Vector;

  std::vector<Point> _points;

  virtual ~VariableMesh() {}
  template <class... Points>
  VariableMesh(Points... points) {
    _points.insert(_points.end(), {points...});
  }

  Point& at(const Math::Type::ULong& index) override { return _points.at(index); }
  const Point& at(const Math::Type::ULong& index) const override { return _points.at(index); }
  typename Point& operator[](std::size_t index) override { return _points.operator[](index); }
  const typename Point& operator[](std::size_t index) const override { return _points.operator[](index); }
  Math::Type::ULong size() const override { return _points.size(); }
  typename MeshBase::Iterator begin() override { return _points.begin(); }
  typename MeshBase::Iterator end() override { return _points.end(); }
  void push_back(const Point& point) { _points.push_back(point); }
  // const typename MeshBase::Iterator begin() const override {
  //  return TensorBase::begin();
  //}
  // const typename MeshBase::Iterator end() const override {
  //  return TensorBase::end();
  //}
  VariableMesh* Clone() const override { return new VariableMesh(*this); }
};

template <class Type, Math::Type::ULong Dimension, Math::Type::ULong CountPoints>
struct StaticMesh : public Mesh<Type, Dimension>, public Tensor::Vector<CountPoints, typename Mesh<Type, Dimension>::Point, Type> {
 private:
  using MeshBase = Mesh<Type, Dimension>;

 public:
  using Point = typename MeshBase::Point;
  using Vector = typename MeshBase::Vector;
  using TensorBase = Tensor::Vector<CountPoints, Point, Type>;

  virtual ~StaticMesh() {}
  StaticMesh() {}
  template <class... Ts>
  StaticMesh(Ts... ts) : TensorBase{ts...} {}

  Point& at(const Math::Type::ULong& index) override { return TensorBase::at(index); }
  const Point& at(const Math::Type::ULong& index) const override { return TensorBase::at(index); }
  typename Point& operator[](std::size_t index) override { return TensorBase::operator[](index); }
  const typename Point& operator[](std::size_t index) const override { return TensorBase::operator[](index); }
  Math::Type::ULong size() const override { return TensorBase::size(); }
  typename MeshBase::Iterator begin() override { return TensorBase::begin(); }
  typename MeshBase::Iterator end() override { return TensorBase::end(); }
  // const typename MeshBase::Iterator begin() const override {
  //  return TensorBase::begin();
  //}
  // const typename MeshBase::Iterator end() const override {
  //  return TensorBase::end();
  //}
  StaticMesh* Clone() const override { return new StaticMesh(*this); }
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
template <class Type, Math::Type::ULong Dimension>
struct Rectangle : public StaticMesh<Type, Dimension, 4> {
 private:
  using Base = StaticMesh<Type, Dimension, 4>;

 public:
  using Point = typename Base::Point;
  using Vector = typename Base::Vector;
  using Vector2D = Tensor::Vector<2, Type>;
  using Colour = typename Point::Colour;

  using Base::Base;
  Rectangle(const Vector& corner, const Vector& opposite_corner,
            const Container::Array::Array<Colour, 4>& colours = {UD::Colour::Black(), UD::Colour::Black(), UD::Colour::Black(), UD::Colour::Black()},
            const Vector& up = {0.0, 1.0}) {
    Type parts = up.Sum();
    Vector diff = opposite_corner - corner;
    Vector up_part = (diff * up) / parts;
    if (corner.x() <= opposite_corner.x()) {
      if (corner.y() <= opposite_corner.y()) {
        bottom_left() = corner;
        bottom_left().colour() = colours[0];
        top_right() = opposite_corner;
        top_right().colour() = colours[2];

        top_left() = bottom_left() + up_part;
        top_left().colour() = colours[1];
        bottom_right() = top_right() - up_part;
        bottom_right().colour() = colours[3];
      } else {
        top_left() = corner;
        top_left().colour() = colours[0];
        bottom_right() = opposite_corner;
        bottom_right().colour() = colours[2];

        top_right() = bottom_right() + up_part;
        top_right().colour() = colours[1];
        bottom_left() = top_left() - up_part;
        bottom_left().colour() = colours[3];
      }
    } else {
      *this = Rectangle{opposite_corner, corner, {colours[2], colours[3], colours[0], colours[1]}, -up};
    }
  }

  Point& bottom_left() { return (*this)[0]; }
  Point& top_left() { return (*this)[1]; }
  Point& top_right() { return (*this)[2]; }
  Point& bottom_right() { return this->operator[](3); }
  Type left() { return bottom_left().x(); }
  void left(Type value) {
    bottom_left().x() = value;
    top_left().x() = value;
  }
  Type right() { return bottom_right().x(); }
  void right(Type value) {
    bottom_right().x() = value;
    top_right().x() = value;
  }
  Type top() { return top_left().y(); }
  void top(Type value) {
    top_left().y() = value;
    top_right().y() = value;
  }
  Type bottom() { return bottom_left().y(); }
  void bottom(Type value) {
    bottom_left().y() = value;
    bottom_right().y() = value;
  }
  Type width() { return right() - left(); }
  void width(Type value) { right(left() + value); }
  Type height() { return top() - bottom(); }
  void height(Type value) { top(bottom() + value); }

  Vector center() { return {left() + width() / 2, bottom() + height() / 2}; }
  void ScaleWidth(const Type& factor) {
    auto d = width();
    auto c = left() + d / 2;
    auto w = (d * factor) / 2;
    left(c - w);
    right(c + w);
  }
  void ScaleHeight(const Type& factor) {
    auto d = height();
    auto c = bottom() + d / 2;
    auto w = d * factor;
    bottom(c - w);
    top(c + w);
  }
  void Scale(const Vector& factors) {
    ScaleWidth(factors.x());
    ScaleHeight(factors.y());
  }
  void Scale(const Type& factor) {
    ScaleWidth(factor);
    ScaleHeight(factor);
  }
};

template <class Type, Math::Type::ULong Dimension>
struct Line : public StaticMesh<Type, 2, Dimension> {
 private:
  using Base = StaticMesh<Type, 2, Dimension>;

 public:
  using Base::Base;
};

template <class _Type, Math::Type::ULong Dimension>
struct World {
  using Type = _Type;
  using Point = Point<Type, Dimension>;
  using Mesh = Mesh<Type, Dimension>;
  using Vector = Tensor::Vector<Dimension, Type>;
  template <Math::Type::ULong Size = Dimension>
  using Matrix = Tensor::Matrix<Size, Dimension, Type>;
  template <Math::Type::ULong Size = Dimension>
  using PointMatrix = Tensor::Tensor<Size, Point>;

  using Face = PointMatrix<Dimension>;
  using Simplex = PointMatrix<Dimension + 1>;
  using ID = Math::Type::ULong;
  using Rectangle = Rectangle<Type, Dimension>;

  using Size = Vector;

  // static std::vector<Vector> ConvexHull_JarvisMarch(
  //    std::vector<Vector> points) {
  //  std::vector<Vector> hull;
  //  Vector first = points.begin();
  //  Vector next = first;
  //  hull.push_back(*points.begin());
  //}
  static constexpr Math::Type::ULong dimension() { return Dimension; }
  static Type Slope(Vector a, Vector b) {
    b -= a;
    return b[1] /= b[0];
  }
  static Type InterceptY(Vector vec, Type slope) { return vec[1] - slope * vec[0]; }
  static Type InterceptY(Vector a, Vector b) { return InterceptY(a, Slope(a, b)); }
  // struct IDrawable {
  //  virtual ~IDrawable() {}
  //  void Draw() = 0;
  //};
  struct Object {
    Vector _position;
    Object() : _position{} {}
    Object(Vector position) : _position(position) {}
    virtual ~Object() {}
    Vector& position() { return _position; }
    // const Vector& position() const { return _position; }
  };
  struct Entity : public Object {
    std::unique_ptr<Mesh> _mesh;

    virtual ~Entity() {}
    Entity(const Mesh& mesh) : _mesh{mesh.Clone()} {}
    Entity(const Point& position, const Mesh& mesh) : Object{position}, _mesh{mesh.Clone()} {}
    Entity(const Entity& entity) : _mesh{entity._mesh->Clone()}, Object{entity._position} {}
    Entity(Entity&& entity) noexcept : Object{std::move(entity)}, _mesh{std::move(entity._mesh)} {}
    Entity& operator=(const Entity& entity) { return *this = Entity(entity); }
    Entity& operator=(Entity&& entity) noexcept {
      _mesh = std::move(entity._mesh);
      return *this;
    }
    Mesh& mesh() { return *_mesh; }
    std::unique_ptr<Mesh> ExchangeMesh(const Mesh& mesh) {
      std::unique_ptr<Mesh> old_mesh{_mesh.release()};
      _mesh.reset(mesh.Clone());
      return std::move(old_mesh);
    }
    // const Mesh& mesh() const { return *_mesh; }
    void Draw() /* const*/ { DrawImpl(); }
    void DrawPoint(const Point& point) /* const*/ {
      glColor4ub(point.colour().r(), point.colour().g(), point.colour().b(), point.colour().a());
      if constexpr (Dimension == 2) {
        glVertex2d(point.x(), point.y());
      } else if constexpr (Dimension == 3) {
        glVertex3d(point.x(), point.y(), point.z());
      }
    }
    template <Math::Type::ULong Size>
    void DrawPoints(Tensor::Tensor<Size, Point> points) {
      glBegin(GL_POINTS);
      for (auto& point : points) DrawPoint(point + this->position());
      glEnd();
    }
    template <Math::Type::ULong Size>
    void DrawLineStrip(Tensor::Tensor<Size, Point> points) {
      glBegin(GL_LINE_STRIP);
      for (auto& point : points) DrawPoint(point + this->position());
      glEnd();
    }
    void DrawLineStrip() {
      glBegin(GL_LINE_STRIP);
      for (auto& point : mesh()) DrawPoint(point + this->position());
      glEnd();
    }
    template <Math::Type::ULong Size>
    void DrawPolygon(Tensor::Tensor<Size, Point> points) {
      glBegin(GL_POLYGON);
      for (auto& point : points) DrawPoint(point + this->position());
      glEnd();
    }
    void ScaleTo(const Rectangle& rect) {
      ScaleXTo(rect.left(), rect.right());
      ScaleYTo(rect.bottom(), rect.top());
    }
    void ScaleXTo(double left, double right) {
      auto width = right - left;
      auto mesh_width = mesh().MaxX() - mesh().MinX();
      auto rel_pos = mesh().begin()->x() - this->position().x();
      this->position().x() = mesh().begin()->x() + width * (rel_pos / mesh_width);
      for (auto& p : mesh()) {
        p.x() = width * (p.x() / mesh_width);
      }
    }
    void ScaleYTo(double bottom, double top) {
      auto height = top - bottom;
      auto mesh_height = mesh().MaxY() - mesh().MinY();
      // Logger::Log(wxString("Height: ") << height);
      // Logger::Log(wxString("MeshHeight: ") << mesh_height);
      // Logger::Log(wxString("MeshBegin: ") << mesh().begin()->y());
      // Logger::Log(wxString("PosBefore: ") << this->position().y());
      // Logger::Log(wxString("PosAfter: ") << this->position().y());
      // auto rel_pos = mesh().begin()->y() - this->position().y();
      for (auto& p : mesh()) {
        p.y() = height * (p.y() / mesh_height);
      }
      // this->position().y() =
      //    mesh().begin()->y() - height * (rel_pos / mesh_height) / 2;
    }
    Vector center() {
      Vector center = {0.0, 0.0};
      for (auto& point : mesh()) center += point + this->position();
      center /= mesh().size();
      return center;
    }
    void ReCenter() {
      Vector shift = center() - this->position();
      this->position() += shift;
      for (auto& point : mesh()) point -= shift;
    }
    void ReCenter(const Vector& vec) {
      auto c = center();
      Vector shift_to_center = c - this->position();
      Vector shift_to_vec = vec - c;
      Vector shift = shift_to_center + shift_to_vec;
      this->position() += shift;
      for (auto& point : mesh()) point -= shift_to_center;
    }
    void ReCenterY(const Type& val) {
      ReCenter();
      this->position().y() += val - this->position().y();
    }

   protected:
    virtual void DrawImpl() /* const*/ { DrawLineStrip(); }
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
  static bool Contains(const Rectangle& rectangle, const Vector& vector) {
    return rectangle.top_left()[0] >= vector[0] && rectangle.top_left()[1] >= vector[1] && rectangle.bottom_right()[0] <= vector[0] &&
           rectangle.bottom_right()[1] <= vector[1];
  }
  template <Math::Type::ULong Size, Math::Type::ULong NumberOfPoints>
  static bool Intersect(const Rectangle& rectangle, const Mesh& mesh) {
    for (auto& point : mesh)
      if (Contains(rectangle, point)) return true;
  }
  template <Math::Type::ULong Size, Math::Type::ULong NumberOfPoints>
  static bool Contains(const Rectangle& rectangle, const Mesh& mesh) {
    for (auto& point : mesh)
      if (!Contains(rectangle, point)) return false;
  }
  struct View : public Entity {
   private:
    using Base = Rectangle;

   public:
    World* _world;
    Vector _offset;
    // TODO: Make world const& (requires implementing generic const iterators)
    virtual ~View() {}
    View(World& world, Vector position, Rectangle rectangle) : Entity{position, rectangle}, _world{&world} {}
    View(World& world, Rectangle rectangle)
        : View{world, rectangle.center(), Rectangle{rectangle.bottom_left() - rectangle.center(), rectangle.top_left() - rectangle.center()}} {}

    void DrawPrep(Rectangle rect) {
      UpdateViewport(rect);
      ProjectionOrtho();
      wxColour clear_colour = *wxBLACK;
      glClearColor(clear_colour.Red(), clear_colour.Green(), clear_colour.Blue(), clear_colour.Alpha());
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      glPointSize(1);
    }
    void Draw(Rectangle rect, std::unordered_set<ID> ids = {}) {
      DrawPrep(rect);
      if (_world->entities.empty()) return;
      if (ids.empty()) {
        for (auto& entity : Visible()) {
          _world->Get<Entity>(entity).Draw();
        }
      } else {
        for (auto& id : ids) {
          if (IsVisible(id)) _world->Get<Entity>(id).Draw();
        }
      }
    }
    void DrawExcluding(Rectangle rect, std::unordered_set<ID> ids) {
      DrawPrep(rect);
      for (auto& entity_pair : this->_world->entities) {
        auto& entity = *entity_pair.second;
        auto& id = entity_pair.first;
        if (!ids.contains(id)) {
          if (IsVisible(entity)) entity.Draw();
        }
      }
    }
    void UpdateViewport(Rectangle rect) {
      glEnable(GL_TEXTURE_2D);
      glEnable(GL_COLOR_MATERIAL);
      glEnable(GL_BLEND);
      glDisable(GL_DEPTH_TEST);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glViewport(rect.left(), rect.bottom(), rect.width(), rect.height());
    }
    void ProjectionOrtho() {
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glOrtho(left(), right(), bottom(), top(), 0, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
    }
    Rectangle& rect() { return static_cast<Rectangle&>(this->mesh()); }
    Type left() { return rect().left() + this->position().x(); }
    Type right() { return rect().right() + this->position().x(); }
    Type bottom() { return rect().bottom() + this->position().y(); }
    Type top() { return rect().top() + this->position().y(); }
    Vector& offset() { return _offset; }
    void Zoom(double amount) { rect().Scale(amount); }
    void ZoomWidth(double amount) { rect().ScaleWidth(amount); }
    void ZoomHeight(double amount) { rect().ScaleHeight(amount); }
    void Translate(Vector amount) { this->position() += amount; }
    void UpdateRect(Rectangle rect) { this->rect() = rect; }

    bool IsVisible(ID entity_id) {
      auto& entity = _world->Get<Entity>(entity_id);
      for (auto& point : entity.mesh()) {
        auto t = point + entity.position();
        auto l = left();
        auto l1 = right();
        auto l2 = top();
        auto l3 = bottom();
        if (t.x() >= left() && t.x() <= right() && t.y() >= bottom() && t.y() <= top()) {
          return true;
        }
      }
      return false;
    }
    virtual std::vector<ID> Visible() {
      std::vector<ID> entities;
      for (auto& entity : _world->entities) {
        if (IsVisible(entity.first)) {
          entities.push_back(entity.first);
        }
      }
      return entities;
    }
    // TODO: replace with inheritance to redefine Visible criteria
    std::vector<ID> VisibleX() {
      std::vector<ID> entities;
      auto left = this->left();
      auto right = this->right();
      auto bot = this->bottom();
      auto top = this->top();
      for (auto& entity : _world->entities) {
        for (auto& point : entity.second->mesh()) {
          auto t = point + entity.second->position();
          if (t.x() >= left && t.x() <= right) {
            entities.push_back(entity.first);
            break;
          }
        }
      }
      return entities;
    }
    // TODO: replace with inheritance to redefine Visible criteria
    std::vector<Entity> VisibleY() {
      std::vector<Entity> entities;
      auto left = this->left();
      auto right = this->right();
      auto bot = this->bottom();
      auto top = this->top();
      for (auto& entity : _world->entities) {
        for (auto& point : entity.second->mesh()) {
          auto t = point + entity.second->position();
          if (t.y() >= bot && t.y() <= top) {
            entities.push_back(*entity.second);
            break;
          }
        }
      }
      return entities;
    }
    void Fit(std::vector<ID> entities, std::unordered_set<ID> ids = {}) {
      FitX(entities, ids);
      FitY(entities, ids);
    }
    void FitX(std::vector<ID> entities, std::unordered_set<ID> ids = {}) {
      if (entities.size() == 0) return;
      Math::Type::ULong start = 0;
      if (!ids.empty()) {
        while (!ids.contains(entities[start])) {
          start++;
          if (start >= entities.size()) return;
        }
      }
      auto entity = &_world->Get<Entity>(entities[start]);
      auto max_width = entity->mesh().MaxX() + entity->position().x();
      auto min_width = entity->mesh().MinX() + entity->position().x();
      for (std::size_t i = start + 1; i < entities.size(); ++i) {
        if (!ids.empty() && !ids.contains(entities[i])) continue;
        entity = &_world->Get<Entity>(entities[i]);

        auto potential_max_width = entity->mesh().MaxX() + entity->position().x();
        auto potential_min_width = entity->mesh().MinX() + entity->position().x();
        if (potential_max_width > max_width) max_width = potential_max_width;
        if (potential_min_width < min_width) min_width = potential_min_width;
      }
      rect().left(min_width);
      rect().right(max_width);
    }
    void FitY(std::vector<ID> entities, std::unordered_set<ID> ids = {}) {
      if (entities.size() == 0) return;
      Math::Type::ULong start = 0;
      if (!ids.empty()) {
        while (!ids.contains(entities[start])) {
          start++;
          if (start >= entities.size()) return;
        }
      }
      auto entity = &_world->Get<Entity>(entities[start]);
      auto max_height = entity->mesh().MaxY() + entity->position().y();
      auto min_height = entity->mesh().MinY() + entity->position().y();
      for (std::size_t i = start + 1; i < entities.size(); ++i) {
        if (!ids.empty() && !ids.contains(entities[i])) continue;
        entity = &_world->Get<Entity>(entities[i]);
        auto potential_max_height = entity->mesh().MaxY() + entity->position().y();
        auto potential_min_height = entity->mesh().MinY() + entity->position().y();
        if (potential_max_height > max_height) max_height = potential_max_height;
        if (potential_min_height < min_height) min_height = potential_min_height;
      }
      rect().bottom(min_height);
      rect().top(max_height);
    }
    void FitExcluding(std::vector<ID> entities, std::unordered_set<ID> ids) {
      FitXExcluding(entities, ids);
      FitYExcluding(entities, ids);
    }
    void FitXExcluding(std::vector<ID> entities, std::unordered_set<ID> ids) {
      if (entities.size() == 0) return;
      Math::Type::ULong start = 0;
      while (ids.contains(entities[start])) {
        start++;
        if (start >= entities.size()) return;
      }
      auto entity = &_world->Get<Entity>(entities[start]);
      auto max_width = entity->mesh().MaxX() + entity->position().x();
      auto min_width = entity->mesh().MinX() + entity->position().x();
      for (std::size_t i = start + 1; i < entities.size(); ++i) {
        if (ids.contains(entities[i])) continue;
        entity = &_world->Get<Entity>(entities[i]);

        auto potential_max_width = entity->mesh().MaxX() + entity->position().x();
        auto potential_min_width = entity->mesh().MinX() + entity->position().x();
        if (potential_max_width > max_width) max_width = potential_max_width;
        if (potential_min_width < min_width) min_width = potential_min_width;
      }
      rect().left(min_width);
      rect().right(max_width);
    }
    void FitYExcluding(std::vector<ID> entities, std::unordered_set<ID> ids) {
      if (entities.size() == 0) return;
      Math::Type::ULong start = 0;
      while (ids.contains(entities[start])) {
        start++;
        if (start >= entities.size()) return;
      }
      auto entity = &_world->Get<Entity>(entities[start]);
      auto max_height = entity->mesh().MaxY() + entity->position().y();
      auto min_height = entity->mesh().MinY() + entity->position().y();
      for (std::size_t i = start + 1; i < entities.size(); ++i) {
        if (ids.contains(entities[i])) continue;
        entity = &_world->Get<Entity>(entities[i]);
        auto potential_max_height = entity->mesh().MaxY() + entity->position().y();
        auto potential_min_height = entity->mesh().MinY() + entity->position().y();
        if (potential_max_height > max_height) max_height = potential_max_height;
        if (potential_min_height < min_height) min_height = potential_min_height;
      }
      rect().bottom(min_height);
      rect().top(max_height);
    }
  };

  ID previous_id;

  std::unordered_map<ID, std::unique_ptr<View>> views;
  std::unordered_map<ID, std::unique_ptr<Entity>> entities;

  World() : previous_id{0}, entities{} {}

  template <class T, class... Ts>
  ID Add(Ts... ts) {
    ID id = NextID();
    entities.insert({id, std::make_unique<T>(ts...)});
    return id;
  }
  template <class T = View, class... Ts>
  requires std::derived_from<T, View> T& CreateView(Ts... ts) {
    ID id = NextID();
    views.insert({id, std::make_unique<View>(*this, ts...)});
    return Get<T>(id);
  }

  void Clear() {
    views.clear();
    entities.clear();
    previous_id = 0;
  }
  void ClearEntities() { entities.clear(); }
  void ClearViews() { views.clear(); }
  ID NextID() {
    ID start_id = previous_id++;

    while (Contains(previous_id)) {
      previous_id++;
      if (previous_id == start_id) throw;
    }
    return previous_id;
  }

  template <class T = Object>
  bool Contains(ID id);
  template <>
  bool Contains<Object>(ID id) {
    return Contains<View>(id) || Contains<Entity>(id);
  }
  template <>
  bool Contains<View>(ID id) {
    return views.count(id);
  }
  template <>
  bool Contains<Entity>(ID id) {
    return entities.count(id);
  }
  template <class T>
  T& Get(ID id);
  template <>
  View& Get<View>(ID id) {
    return *views.at(id);
  }
  template <>
  Entity& Get<Entity>(ID id) {
    return *entities.at(id);
  }
};  // namespace World
}  // namespace World
}  // namespace UD