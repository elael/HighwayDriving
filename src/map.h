#pragma once
#include <vector>
#include "Eigen/Dense"
#include "Eigen/StdVector"


struct FrenetFrame {
  double s, d;
};

struct Map{
  std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d> > map_xy;
  std::vector<double> maps_x, maps_y, maps_s, maps_dx, maps_dy;
  static constexpr double smoother_s = 10;

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
FrenetFrame getFrenet(const Eigen::Vector2d& xy, double theta);

// Transform from Frenet s,d coordinates to Cartesian x,y
Eigen::Vector2d getXY(double s, double d);


Eigen::Vector2d smooth_getXY(double s, double d){
  return (getXY(s + smoother_s, d) + getXY(s - smoother_s, d))/2;
}

private:

// Calculate closest waypoint to current x, y position
int ClosestWaypoint(const Eigen::Vector2d& xy);

// Returns next waypoint of the closest waypoint
int NextWaypoint(const Eigen::Vector2d& xy, double theta);

};
