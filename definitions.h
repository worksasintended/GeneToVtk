#pragma once
#include <vector>

struct Point {
  float phi;
  float theta;
  float x;
  float y;
  float z;
  std::vector<float> time_values;
};

struct Slice {
  std::vector<Point> points;
};


struct Surface {
  std::vector<Slice> slices;
};

struct SmallPoint{
  float x;
  float y;
  float z;
  float value;
};

struct SmallSlice{
  std::vector<SmallPoint> points;
};

struct SmallSurface{
  std::vector<SmallSlice> slices;
};
