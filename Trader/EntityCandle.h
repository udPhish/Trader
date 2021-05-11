#pragma once
#include <algorithm>

#include "Candle.h"
#include "World.h"
struct EntityCandle : public UD::World::World<double, 2>::MeshEntity<
                          UD::World::StaticMesh<double, 2, 6>> {
 private:
  using World = UD::World::World<double, 2>;
  using Mesh = UD::World::StaticMesh<double, 2, 6>;
  using Base = World::MeshEntity<Mesh>;

 public:
  EntityCandle(const Candle& candle) {
    double top = std::max(candle.open, candle.close);
    double bottom = std::min(candle.open, candle.close);
    double left = static_cast<double>(candle.open_time.count());
    double right = static_cast<double>(candle.close_time.count());
    double mid = left + (right - left) / 2;

    position() = World::Vector{left, bottom};
    top_left() = World::Point{{left, top}, UD::Colour::White()} - position();
    top_right() = World::Point{{right, top}, UD::Colour::White()} - position();
    top_middle() = World::Point{{mid, top}, UD::Colour::White()} - position();
    top_wick() =
        World::Point{{mid, candle.high}, UD::Colour::Green()} - position();
    bottom_left() =
        World::Point{{left, bottom}, UD::Colour::White()} - position();
    bottom_right() =
        World::Point{{right, bottom}, UD::Colour::White()} - position();
    bottom_middle() =
        World::Point{{mid, bottom}, UD::Colour::White()} - position();
    bottom_wick() =
        World::Point{{mid, candle.low}, UD::Colour::Red()} - position();
  }

 public:
  World::Point& top_left() { return this->mesh().at(0); }
  World::Point& top_right() { return this->mesh().at(1); }
  World::Point& top_middle() { return this->mesh().at(2); }
  World::Point& top_wick() { return this->mesh().at(3); }
  World::Point& bottom_left() { return this->mesh().at(4); }
  World::Point& bottom_right() { return this->mesh().at(5); }
  World::Point& bottom_middle() { return this->mesh().at(6); }
  World::Point& bottom_wick() { return this->mesh().at(7); }
};
