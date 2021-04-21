#pragma once
#include <algorithm>

#include "Candle.h"
#include "World.h"
struct EntityCandle : public World<2>::Entity {
   Candle candle;

  EntityCandle( Candle candle) : candle{candle} {}
  virtual World<2>::Vector position() override {
    double bottom = std::min(candle.open, candle.close);
    return {candle.open_time.count(), static_cast<int>(bottom * 10000)};
  }
  virtual World<2>::Mesh mesh() override {
    wxColour colour = *wxWHITE;
    World<2>::Mesh ret(
        {top_left(), top_middle(), top_wick(), top_right(), bottom_right(),
         bottom_middle(), bottom_wick(), bottom_left()},
        {0, 1, 2, 1, 3, 4, 5, 6, 5, 7, 0}, {colour});
    return ret;
  }
  World<2>::Vector top_left() {
    double top = std::max(candle.open, candle.close);
    return World<2>::Difference(
        position(), {candle.open_time.count(), static_cast<int>(top * 10000)});
  }
  World<2>::Vector top_middle() {
    double top = std::max(candle.open, candle.close);
    return World<2>::Difference(
        position(),
        {candle.open_time.count() +
             (candle.close_time.count() - candle.open_time.count()) / 2,
         static_cast<int>(top * 10000)});
  }
  World<2>::Vector top_right() {
    double top = std::max(candle.open, candle.close);
    return World<2>::Difference(
        position(), {candle.close_time.count(), static_cast<int>(top * 10000)});
  }
  World<2>::Vector top_wick() {
    return World<2>::Difference(
        position(),
        {candle.open_time.count() +
             (candle.close_time.count() - candle.open_time.count()) / 2,
         static_cast<int>(candle.high * 10000)});
  }
  World<2>::Vector bottom_left() {
    double bottom = std::min(candle.open, candle.close);
    return World<2>::Difference(
        position(), {candle.open_time.count(),
                                             static_cast<int>(bottom * 10000)});
  }
  World<2>::Vector bottom_middle() {
    double bottom = std::min(candle.open, candle.close);
    return World<2>::Difference(
        position(),
        {candle.open_time.count() +
             (candle.close_time.count() - candle.open_time.count()) / 2,
         static_cast<int>(bottom * 10000)});
  }
  World<2>::Vector bottom_right() {
    double bottom = std::min(candle.open, candle.close);
    return World<2>::Difference(
        position(), {candle.close_time.count(),
                                             static_cast<int>(bottom * 10000)});
  }
  World<2>::Vector bottom_wick() {
    return World<2>::Difference(
        position(),
        {candle.open_time.count() +
             (candle.close_time.count() - candle.open_time.count()) / 2,
         static_cast<int>(candle.low * 10000)});
  }
};
