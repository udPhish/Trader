#pragma once
#include <algorithm>

#include "Candle.h"
#include "Colour.h"
#include "Strategies.h"
#include "World.h"
struct EntityCandle : public UD::World::World<double, 2>::Entity {
 private:
  using World = UD::World::World<double, 2>;
  using Mesh = UD::World::StaticMesh<double, 2, 8>;
  using Base = World::Entity;
  Position _position;

 public:
  EntityCandle(const Candle& candle, Position strat_position) : Entity{Mesh{}} {
    double top = std::max(candle.open, candle.close);
    double bottom = std::min(candle.open, candle.close);
    double left = static_cast<double>(candle.open_time.count());
    double right = static_cast<double>(candle.close_time.count());
    double mid = left + (right - left) / 2;

    position() = World::Vector{left, bottom};
    UD::Colour::Colour<> col = strat_position == Position::Short
                                   ? (UD::Colour::Colour<>)UD::Colour::Red()
                                   : (UD::Colour::Colour<>)UD::Colour::Green();

    top_left() = World::Point{{left, top}, col} - position();
    top_right() = World::Point{{right, top}, col} - position();
    top_middle() = World::Point{{mid, top}, col} - position();
    top_wick() = World::Point{{mid, candle.high}, col} - position();
    bottom_left() = World::Point{{left, bottom}, col} - position();
    bottom_right() = World::Point{{right, bottom}, col} - position();
    bottom_middle() = World::Point{{mid, bottom}, col} - position();
    bottom_wick() = World::Point{{mid, candle.low}, col} - position();
  }

  World::Point& top_left() { return this->mesh().at(0); }
  World::Point& top_right() { return this->mesh().at(1); }
  World::Point& top_middle() { return this->mesh().at(2); }
  World::Point& top_wick() { return this->mesh().at(3); }
  World::Point& bottom_left() { return this->mesh().at(4); }
  World::Point& bottom_right() { return this->mesh().at(5); }
  World::Point& bottom_middle() { return this->mesh().at(6); }
  World::Point& bottom_wick() { return this->mesh().at(7); }
  // const World::Point& top_left() const { return this->mesh().at(0); }
  // const World::Point& top_right() const { return this->mesh().at(1); }
  // const World::Point& top_middle() const { return this->mesh().at(2); }
  // const World::Point& top_wick() const { return this->mesh().at(3); }
  // const World::Point& bottom_left() const { return this->mesh().at(4); }
  // const World::Point& bottom_right() const { return this->mesh().at(5); }
  // const World::Point& bottom_middle() const { return this->mesh().at(6); }
  // const World::Point& bottom_wick() const { return this->mesh().at(7); }
  void SetColour(UD::Colour::Colour<> colour) {
    top_left().colour() = colour;
    top_right().colour() = colour;
    top_middle().colour() = colour;
    top_wick().colour() = colour;
    bottom_left().colour() = colour;
    bottom_right().colour() = colour;
    bottom_middle().colour() = colour;
    bottom_wick().colour() = colour;
  }

 protected:
   void DrawImpl() /* const*/ override {
    this->DrawPolygon<4>({top_left(), top_right(),
                          bottom_right(),
                          bottom_left()});
    this->DrawLineStrip<2>(
        {top_wick(), top_middle() });
    this->DrawLineStrip<2>(
        {bottom_wick() , bottom_middle() });
  }
};
