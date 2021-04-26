#pragma once
#include <array>
#include <memory>
#include <unordered_map>
#include <vector>

#include <wx/wx.h>

#include <wx/colour.h>
#include <wx/glcanvas.h>

#include "Logger.h"

template <std::size_t Dimension>
struct World {
  using Unit = double;
  using Vector = std::array<Unit, Dimension>;
  template<std::size_t Size>
  using Matrix = std::array<Vector, Size>;
  using Face = std::array<Vector, Dimension>;
  using Simplex = std::array<Face, Dimension + 1>;
  using ID = std::size_t;

  static Vector Difference(Vector lhs, Vector rhs) {
    Vector ret{};
    for (std::size_t i = 0; i < lhs.size(); ++i) {
      ret[i] = lhs[i] - rhs[i];
    }
    return ret;
  }
  static Vector Sum(Vector lhs, Vector rhs) {
    Vector ret{};
    for (std::size_t i = 0; i < lhs.size(); ++i) {
      ret[i] = lhs[i] + rhs[i];
    }
    return ret;
  }
  static std::vector<Vector> ConvexHull_JarvisMarch(
      std::vector<Vector> points) {
    std::vector<Vector> hull;
    Vector first = points.begin();
    Vector next = first;
    hull.push_back(*points.begin());
  }
  static Vector DirectionalVector(Vector a, Vector b) {
    return Difference(b, a);
  }
  static Unit Slope(Vector a, Vector b) {
    return (b[1] - a[1]) / (b[0] - a[0]);
  }
  static Unit InterceptY(Vector vec, Unit slope) {
    return vec[1] - slope * vec[0];
  }
  static Unit InterceptY(Vector a, Vector b) {
    return InterceptY(a, Slope(a, b));
  }
  struct Rectangle {
    Vector bottom_left;
    Vector size;
    Unit x() { return bottom_left[0]; }
    Unit y() { return bottom_left[1]; }
    Unit width() { return size[0]; }
    Unit height() { return size[1]; }

    Unit left() { return x(); }
    Unit right() { return x() + width(); }
    Unit bottom() { return y(); }
    Unit top() { return y() + height(); }

    Rectangle Centered(Vector center) {
      Unit half_width = width() / 2;
      Unit half_height = height() / 2;
      return Rectangle{{center[0] - half_width, center[1] - half_height},
                       {width(), height()}};
    }
  };

  struct Mesh {
    std::vector<Vector> points;
    std::vector<wxColour> colours;
    std::vector<std::size_t> indices;

    Mesh(std::vector<Vector> points) : points{points} {
      colours.push_back(*wxWHITE);
      for (std::size_t i = 0; i < points.size(); ++i) {
        indices.push_back(i);
      }
    }
    Mesh(std::vector<Vector> points, std::vector<std::size_t> indices)
        : points{points}, indices{indices} {
      colours.push_back(*wxWHITE);
    }
    Mesh(std::vector<Vector> points, std::vector<wxColour> colours)
        : points{points}, colours{colours} {
      for (std::size_t i = 0; i < points.size(); ++i) {
        indices.push_back(i);
      }
    }
    Mesh(std::vector<Vector> points, std::vector<std::size_t> indices,
         std::vector<wxColour> colours)
        : points{points}, indices{indices}, colours{colours} {}

    Vector& point(std::size_t index) { return points[indices[index]]; }
    wxColour& colour(std::size_t index) {
      return colours[indices[index] % colours.size()];
    }
    std::size_t size() { return indices.size(); }
  };
  struct Object {
    virtual Vector position() = 0;
    virtual ~Object() {}
  };
  struct Entity : public Object {
    virtual Mesh mesh() = 0;
    virtual ~Entity() {}
  };
  struct View : public Object {
    const World& world;
    Rectangle rect;
    Vector m_position;

    View(const World& world, Rectangle rect)
        : world{world}, rect{rect}, m_position{0, 0} {}
    View(const World& world, Rectangle rect, Vector position)
        : world{world}, rect{rect}, m_position{position} {}

    Vector position() override { return m_position; }
    void Draw(Rectangle rect) {
      UpdateViewport(rect);
      ProjectionOrtho();

      // glEnable(GL_TEXTURE_2D);
      // glEnable(GL_COLOR_MATERIAL);
      // glEnable(GL_BLEND);
      // glEnable(GL_SCISSOR_TEST);
      // glDisable(GL_DEPTH_TEST);
      // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      // glDepthFunc(GL_LEQUAL);

      // glEnable(GL_DEPTH_TEST);
      // glDepthFunc(GL_LEQUAL);
      // glEnable(GL_BLEND);
      // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      wxColour clear_colour = *wxBLACK;
      glClearColor(clear_colour.Red(), clear_colour.Green(),
                   clear_colour.Blue(), clear_colour.Alpha());
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      std::size_t iter = 0;
      std::size_t print_iter = 500;
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glPointSize(1);
      for (auto& entity : world.entities) {
        if (entity.second->mesh().size() > 1) {
          glBegin(GL_LINE_STRIP);
          Vector prev = entity.second->mesh().point(0);
          for (std::size_t i = 0; i < entity.second->mesh().size(); ++i) {
            // glColor4f(colour.Red(), colour.Green(), colour.Blue(),
            //          colour.Alpha());
            // glVertex2f(prev[0], prev[1]);
            prev =
                Sum(entity.second->mesh().point(i), entity.second->position());
            wxColour colour = entity.second->mesh().colour(i);
            // glColor4f(colour.Red(), colour.Green(), colour.Blue(),
            //          colour.Alpha());
            if (!(iter % print_iter)) {
              Logger::Log(wxString("Vec[")
                          << i << "]{x:" << prev[0] << ",y:" << prev[1] << "}");
            }
            glVertex2i(prev[0], prev[1]);
            glColor4f(colour.Red(), colour.Green(), colour.Blue(), 1);
          }
          glEnd();
        }
        iter++;
      }
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
      Rectangle rect = VisibleRectangle();
      Logger::Log(wxString("Projecting{left:")
                  << rect.left() << ",right:" << rect.right() << ",bottom:"
                  << rect.bottom() << ",top:" << rect.top() << "}");
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glOrtho(rect.left(), rect.right(), rect.bottom(), rect.top(), 0, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
    }
    Rectangle VisibleRectangle() { return rect.Centered(position()); }
    void Zoom(double amount) {
      rect.size[0] *= amount;
      rect.size[1] *= amount;
    }
    void Translate(Vector amount) {
      m_position[0] += amount[0];
      m_position[1] += amount[1];
    }
  };

  static bool Intersects(const Entity& lhs, const Entity& rhs) {}
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
  View& CreateView(Rectangle rect, Vector position) {
    ID id = NextID();
    views.insert({id, std::make_unique<View>(*this, rect, position)});
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
