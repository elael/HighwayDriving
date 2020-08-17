#include "map.h"

#include <vector>

#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "helpers.h"

using Eigen::Vector2d;

// Calculate closest waypoint to current x, y position
int Map::ClosestWaypoint(const Vector2d& xy) const {
  double closestLen = std::numeric_limits<double>::infinity();  //large number
  int closestWaypoint = 0;

  for (int i = 0; i < map_xy.size(); ++i) {
    double dist = (xy - map_xy[i]).squaredNorm();
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }
  }

  return closestWaypoint;
}

// Returns next waypoint of the closest waypoint
int Map::NextWaypoint(const Vector2d& xy, double theta) const {
  int closestWaypoint = ClosestWaypoint(xy);

  double heading = angle(map_xy[closestWaypoint] - xy);

  double angle = fabs(theta - heading);
  angle = std::min(2 * pi() - angle, angle);

  if (angle > pi() / 2 && ++closestWaypoint == map_xy.size())
    return 0;

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
FrenetFrame Map::getFrenet(const Vector2d& xy, double theta) const {
  int next_wp = NextWaypoint(xy, theta);

  int prev_wp;
  prev_wp = next_wp - 1;
  if (next_wp == 0) {
    prev_wp = map_xy.size() - 1;
  }

  const Vector2d& p_prev = map_xy[prev_wp];

  Vector2d p_n = map_xy[next_wp] - p_prev;
  Vector2d p_x = xy - p_prev;

  // find the projection of x onto n
  Vector2d proj = p_n * (p_x.dot(p_n) / p_n.dot(p_n));

  double frenet_d = (p_x - proj).norm();

  //see if d value is positive or negative by comparing it to a center point
  static const Vector2d center_ref = []() {Vector2d temp; temp << 1000, 2000; return temp; }();
  Vector2d center = center_ref - p_prev;

  double centerToPos = (center - p_x).squaredNorm();
  double centerToRef = (center - proj).squaredNorm();

  if (centerToPos <= centerToRef) frenet_d *= -1;

  // calculate s value
  double frenet_s = maps_s[prev_wp] + proj.norm();

  return {frenet_s, frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
Vector2d Map::getXY(double s, double d) const {
  int prev_wp = -1;

  while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1)))
    ++prev_wp;

  int wp2 = (prev_wp + 1) % map_xy.size();

  double ds = (maps_s[wp2] - maps_s[prev_wp]);

  Vector2d ending = map_xy[wp2] + d * map_dxdy[wp2];
  double ending_weight = (s - maps_s[prev_wp]) / ds;

  Vector2d begining = map_xy[prev_wp] + d * map_dxdy[prev_wp];
  double begining_weight = (maps_s[wp2] - s) / ds;

  return begining_weight * begining + ending_weight * ending;
}

// Transform from Frenet s,d coordinates to Cartesian x,y
Vector2d Map::smooth_getXY(double s, double d) const {
  s = fmod(s, road_length_);

  Vector2d inner_pos(inner_sx(s), inner_sy(s));
  Vector2d outer_pos(outer_sx(s), outer_sy(s));

  return inner_pos + d * (outer_pos - inner_pos);
}

void Map::smoothMap() {
  std::vector<double> s, inner_sxx, inner_syy, outer_sxx, outer_syy;
  for (size_t i = 0; i < maps_s.size(); ++i) {
    s.emplace_back(maps_s[i]);
    inner_sxx.emplace_back(map_xy[i][0]);
    inner_syy.emplace_back(map_xy[i][1]);
    Vector2d outer_xy = map_xy[i] + map_dxdy[i];
    outer_sxx.emplace_back(outer_xy[0]);
    outer_syy.emplace_back(outer_xy[1]);
  }

  Vector2d tg = ccPerp(map_dxdy[0]);

  inner_sx.set_boundary(tk::spline::first_deriv, tg[0], tk::spline::first_deriv, tg[0]);
  inner_sx.set_points(s, std::move(inner_sxx));

  inner_sy.set_boundary(tk::spline::first_deriv, tg[1], tk::spline::first_deriv, tg[1]);
  inner_sy.set_points(s, std::move(inner_syy));

  outer_sx.set_boundary(tk::spline::first_deriv, tg[0], tk::spline::first_deriv, tg[0]);
  outer_sx.set_points(s, std::move(outer_sxx));

  outer_sy.set_boundary(tk::spline::first_deriv, tg[1], tk::spline::first_deriv, tg[1]);
  outer_sy.set_points(std::move(s), std::move(outer_syy));

  road_length_ = maps_s.back() + (map_xy.back() - map_xy.front()).norm();
}
