#pragma once
#include <algorithm>

#include "Candle.h"
#include "World.h"
struct EntityCandle : public World<2>::Entity {
  Candle candle;

  EntityCandle(Candle candle) : candle{candle} {}
  virtual World<2>::Vector position() override {
    double bottom = std::min(candle.open, candle.close);
    return {static_cast<double>(candle.open_time.count()), bottom};
  }
  virtual World<2>::Mesh mesh() override {
    // wxColour colour = *wxWHITE;
    World<2>::Mesh ret(
        {top_left(), top_middle(), top_wick(), top_right(), bottom_right(),
         bottom_middle(), bottom_wick(), bottom_left()},
        {0, 1, 2, 1, 3, 4, 5, 6, 5, 7, 0},
        {*wxWHITE, *wxWHITE, *wxGREEN, *wxWHITE, *wxWHITE, *wxWHITE,
         *wxRED, *wxWHITE});
    return ret;
  }
  World<2>::Vector top_left() {
    double top = std::max(candle.open, candle.close);
    return World<2>::Difference(
        position(), {static_cast<double>(candle.open_time.count()), top});
  }
  World<2>::Vector top_middle() {
    double top = std::max(candle.open, candle.close);
    return World<2>::Difference(
        position(), {static_cast<double>(candle.open_time.count()) +
                         (static_cast<double>(candle.close_time.count()) -
                          static_cast<double>(candle.open_time.count())) /
                             2,
                     top});
  }
  World<2>::Vector top_right() {
    double top = std::max(candle.open, candle.close);
    return World<2>::Difference(
        position(), {static_cast<double>(candle.close_time.count()), top});
  }
  World<2>::Vector top_wick() {
    return World<2>::Difference(
        position(), {static_cast<double>(candle.open_time.count()) +
                         (static_cast<double>(candle.close_time.count()) -
                          static_cast<double>(candle.open_time.count())) /
                             2,
         candle.high});
  }
  World<2>::Vector bottom_left() {
    double bottom = std::min(candle.open, candle.close);
    return World<2>::Difference(
        position(), {static_cast<double>(candle.open_time.count()), bottom});
  }
  World<2>::Vector bottom_middle() {
    double bottom = std::min(candle.open, candle.close);
    return World<2>::Difference(
        position(), {static_cast<double>(candle.open_time.count()) +
                         (static_cast<double>(candle.close_time.count()) -
                          static_cast<double>(candle.open_time.count())) /
                             2,
         bottom});
  }
  World<2>::Vector bottom_right() {
    double bottom = std::min(candle.open, candle.close);
    return World<2>::Difference(
        position(), {static_cast<double>(candle.close_time.count()), bottom});
  }
  World<2>::Vector bottom_wick() {
    return World<2>::Difference(
        position(), {static_cast<double>(candle.open_time.count()) +
                         (static_cast<double>(candle.close_time.count()) -
                          static_cast<double>(candle.open_time.count())) /
                             2,
         candle.low});
  }
};
