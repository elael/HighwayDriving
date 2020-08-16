#pragma once
#include <vector>
#include "Eigen/Dense"
#include "Eigen/StdVector"
#include "spline.h"

struct FrenetFrame
{
  double s, d;
};

struct Map
{
  std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> map_xy;
  std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> map_dxdy;
  std::vector<double> maps_s;
  tk::spline inner_sx;
  tk::spline inner_sy;
  tk::spline outer_sx;
  tk::spline outer_sy;
  static constexpr double smoother_s = 20;

  void smoothMap();

  // Transform from Cartesian x,y coordinates to Frenet s,d coordinates
  FrenetFrame getFrenet(const Eigen::Vector2d &xy, double theta) const;

  // Transform from Frenet s,d coordinates to Cartesian x,y
  Eigen::Vector2d getXY(double s, double d) const;

  Eigen::Vector2d smooth_getXY(double s, double d) const;

private:
  // Calculate closest waypoint to current x, y position
  int ClosestWaypoint(const Eigen::Vector2d &xy) const;

  // Returns next waypoint of the closest waypoint
  int NextWaypoint(const Eigen::Vector2d &xy, double theta) const;

  double road_length_;
};
