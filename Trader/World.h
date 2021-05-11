#pragma once
#include <array>
#include <limits>
#include <memory>
#include <unordered_map>
#include <vector>

#include <wx/wx.h>

#include <wx/glcanvas.h>

#include "Colour.h"
#include "Logger.h"
#include "Tensor.h"

namespace UD {
namespace World {
using DEFAULT_TYPE = double;

template <class Type, std::size_t Dimension>
struct Point : public UD::Tensor::Vector<Type, Dimension> {
 private:
  using Base = UD::Tensor::Vector<Type, Dimension>;

 public:
  using Colour = UD::Colour::Colour<>;
  Colour _colour;
  Colour& colour() { return _colour; }
  Point() : _colour{} {}
  Point(Base vector) : Base{vector}, _colour{} {}
  Point(Base vector, Colour colour) : Base{vector}, _colour{colour} {}
};

template <class Type, std::size_t Dimension>
struct Mesh : public UD::Tensor::Range<Point<Type, Dimension>> {
  using Point = Point<Type, Dimension>;
  using Vector = UD::Tensor::Vector<Type, Dimension>;

 private:
  using Base = UD::Tensor::Range<Point>;

 protected:
  virtual Mesh& _Scale(Vector factor) = 0;

 public:
  virtual Point& _at(const std::size_t& index) override{};
  virtual UD::Tensor::Iterator<Point> begin() override{};
  virtual UD::Tensor::Iterator<Point> end() override{};
  virtual UD::Tensor::Iterator<Point> rbegin() override{};
  virtual UD::Tensor::Iterator<Point> rend() override{};
  virtual Mesh& Translate(Vector position){};
  Mesh& Scale(Vector factor) { return _Scale(std::move(factor)); }
  Mesh& Scale(Type factor) {
    Vector vec;
    vec.fill(factor);
    return _Scale(vec);
  }
};

template <class Type, std::size_t Dimension, std::size_t CountPoints>
struct StaticMesh : public Mesh<Type, Dimension> {
 private:
  using Base = Mesh<Type, Dimension>;

 public:
  using Point = typename Base::Point;
  using Vector = typename Base::Vector;
  using PointMatrix = UD::Tensor::Vector<Point, CountPoints>;

  PointMatrix _points;

  virtual ~StaticMesh() {}
  StaticMesh(PointMatrix points) : _points{points} {}

  virtual Point& _at(const std::size_t& index) override {
    return _points[index];
  }
  virtual std::size_t size() const override { return _points.size(); }
  virtual UD::Tensor::Iterator<Point> begin() override {
    return UD::Tensor::IteratorWrapper<Point>{_points.begin()};
  }
  virtual UD::Tensor::Iterator<Point> end() override {
    return UD::Tensor::IteratorWrapper<Point>{_points.end()};
  }
  virtual UD::Tensor::Iterator<Point> rbegin() override {
    return UD::Tensor::IteratorWrapper<Point>{_points.rbegin()};
  }
  virtual UD::Tensor::Iterator<Point> rend() override {
    return UD::Tensor::IteratorWrapper<Point>{_points.rend()};
  }
  PointMatrix& points() { return _points; }

  virtual StaticMesh& Translate(Vector position) override {
    UD::Tensor::Translate(points(), position);
  }
  virtual StaticMesh& _Scale(Vector factor) override {
    UD::Tensor::Scale(points(), factor);
  }
};
template <class Type, std::size_t Dimension>
struct Rectangle : public StaticMesh<Type, 4, Dimension> {
 private:
  using Base = StaticMesh<Type, 4, Dimension>;

 public:
  using Point = typename Base::Point;
  using Vector = typename Base::Vector;
  using Vector2D = UD::Tensor::Vector<Type, 2>;
  using Colour = typename Point::Colour;

  using Base::Base;
  Rectangle(const Vector& corner, const Vector& opposite_corner,
            const UD::Tensor::Array<Colour, 4>& colours = {UD::Colour::Black(),
                                                           UD::Colour::Black(),
                                                           UD::Colour::Black(),
                                                           UD::Colour::Black()},
            const Vector& up = {0, 1}) {
    Type parts = UD::Tensor::Sum(up);
    Vector diff = opposite_corner - corner;
    Vector up_part = UD::Tensor::Mutliply(diff, up) / parts;
    if (corner.x() < opposite_corner.x()) {
      if (corner.y() < opposite_corner.y()) {
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
      *this = Rectangle{opposite_corner,
                        corner,
                        {colours[2], colours[3], colours[0], colours[1]},
                        -up};
    }
  }

  Point bottom_left() { return (*this)[0]; }
  Point top_left() { return (*this)[1]; }
  Point top_right() { return (*this)[2]; }
  Point bottom_right() { return (*this)[3]; }
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
  void Width(Type value) { right(left() + value); }
  Type height() { return top() - bottom(); }
  void height(Type value) { top(bottom() + value); }
};

template <class Type, std::size_t Dimension>
struct Line : public StaticMesh<Type, 2, Dimension> {
 private:
  using Base = StaticMesh<Type, 2, Dimension>;

 public:
  using Base::Base;
};

static void Draw(Mesh<double, 2>& mesh) {
  glBegin(GL_POINTS);
  for (auto& point : mesh) {
  }
  glEnd();
}
template <class _Type, std::size_t Dimension>
struct World {
  using Type = _Type;
  using Point = Point<Type, Dimension>;
  using Mesh = Mesh<Type, Dimension>;
  using Vector = UD::Tensor::Vector<Type, Dimension>;
  template <std::size_t Size = Dimension>
  using Matrix = UD::Tensor::Matrix<Type, Size, Dimension>;
  template <std::size_t Size = Dimension>
  using PointMatrix = UD::Tensor::Tensor<Point, Size>;

  using Face = PointMatrix<Dimension>;
  using Simplex = PointMatrix<Dimension + 1>;
  using ID = std::size_t;
  using Rectangle = Rectangle<Type, Dimension>;

  using Size = Vector;

  // static std::vector<Vector> ConvexHull_JarvisMarch(
  //    std::vector<Vector> points) {
  //  std::vector<Vector> hull;
  //  Vector first = points.begin();
  //  Vector next = first;
  //  hull.push_back(*points.begin());
  //}
  static constexpr std::size_t dimension() { return Dimension; }
  static Type Slope(Vector a, Vector b) {
    b -= a;
    return b[1] /= b[0];
  }
  static Type InterceptY(Vector vec, Type slope) {
    return vec[1] - slope * vec[0];
  }
  static Type InterceptY(Vector a, Vector b) {
    return InterceptY(a, Slope(a, b));
  }
  struct IDrawable {
    virtual ~IDrawable() {}
    void Draw() = 0;
  };
  struct Object {
    Vector _position;
    Object() : _position{} {}
    Object(Vector position) : _position(position) {}
    virtual ~Object() {}
    Vector& position() { return _position; }
  };
  struct Entity : public Object {
    virtual ~Entity() {}
    virtual void Draw() = 0;
  };
  template <class MeshType = Mesh>
  struct MeshEntity : public Entity {
    MeshType _mesh;
    virtual ~MeshEntity() {}
    MeshEntity() {}
    MeshEntity(const MeshType& mesh) : _mesh{mesh} {}
    MeshType& mesh() { return _mesh; }

    virtual void Draw() { Draw(mesh()); }
  };
  static bool Contains(const Rectangle& rectangle, const Vector& vector) {
    return rectangle.top_left()[0] >= vector[0] &&
           rectangle.top_left()[1] >= vector[1] &&
           rectangle.bottom_right()[0] <= vector[0] &&
           rectangle.bottom_right()[1] <= vector[1];
  }
  template <std::size_t Size, std::size_t NumberOfPoints>
  static bool Intersect(const Rectangle& rectangle, const Mesh& mesh) {
    for (auto& point : mesh)
      if (Contains(rectangle, point)) return true;
  }
  template <std::size_t Size, std::size_t NumberOfPoints>
  static bool Contains(const Rectangle& rectangle, const Mesh& mesh) {
    for (auto& point : mesh)
      if (!Contains(rectangle, point)) return false;
  }
  struct View : public Rectangle {
   private:
    using Base = Rectangle;

   public:
    const World& _world;

    View(const World& world, Rectangle rectangle)
        : Rectangle{rectangle}, _world{world} {}
    View(const World& world, Rectangle size, Vector position)
        : Rectangle{position, size}, _world{world} {}

    void Draw(Rectangle rect) {
      UpdateViewport(rect);
      ProjectionOrtho();
      wxColour clear_colour = *wxBLACK;
      glClearColor(clear_colour.Red(), clear_colour.Green(),
                   clear_colour.Blue(), clear_colour.Alpha());
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      std::size_t iter = 0;
      std::size_t print_iter = 500;
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glPointSize(1);
      for (auto& entity : _world.entities) {
        entity.second.Draw();
      }
      iter++;
    }
    void UpdateViewport(Rectangle rect) {
      Logger::Log(wxString("Viewport{x:")
                  << rect.x() << ",y:" << rect.y() << ",width:" << rect.width()
                  << ",height:" << rect.height() << "}");
      glEnable(GL_TEXTURE_2D);
      glEnable(GL_COLOR_MATERIAL);
      glEnable(GL_BLEND);
      glDisable(GL_DEPTH_TEST);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glViewport(rect.x(), rect.y(), rect.width(), rect.height());
    }
    void ProjectionOrtho() {
      Logger::Log(wxString("Projecting{left:")
                  << this->left() << ",right:" << this->right() << ",bottom:"
                  << this->bottom() << ",top:" << this->top() << "}");
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glOrtho(this->left(), this->right(), this->bottom(), this->top(), 0, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
    }
    void Zoom(double amount) {
      this->width(this->width() * amount);
      this->height(this->height() * amount);
    }
    void Translate(Vector amount) { UD::Tensor::Translate(*this, amount); }
  };

  ID previous_id;

  std::unordered_map<ID, std::unique_ptr<View>> views;
  std::unordered_map<ID, std::unique_ptr<IDrawable>> entities;

  World() : previous_id{0}, entities{} {}

  template <class T, class... Ts>
  ID Add(Ts... ts) {
    ID id = NextID();
    entities.insert({id, std::make_unique<T>(ts...)});
    return id;
  }
  View& CreateView(Rectangle rect) {
    ID id = NextID();
    views.insert({id, std::make_unique<View>(*this, rect)});
    return Get<View>(id);
  }

  void Clear() {
    views.clear();
    entities.clear();
    previous_id = 0;
  }
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
};
}  // namespace World
}  // namespace UD