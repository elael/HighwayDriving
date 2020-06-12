#include <vector>
#include "Eigen/Dense"
#include "map.h"
#include "helpers.h"

using Eigen::Vector2d;

// Calculate closest waypoint to current x, y position
int Map::ClosestWaypoint(const Vector2d& xy) {
  double closestLen = std::numeric_limits<double>::infinity(); //large number
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
int Map::NextWaypoint(const Vector2d& xy, double theta) {
  int closestWaypoint = ClosestWaypoint(xy);

  double heading = angle(map_xy[closestWaypoint] - xy);

  double angle = fabs(theta-heading);
  angle = std::min(2*pi() - angle, angle);

  if (angle > pi()/2 && ++closestWaypoint == map_xy.size())
    return 0;

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
FrenetFrame Map::getFrenet(const Vector2d& xy, double theta) {
  int next_wp = NextWaypoint(xy, theta);

  int prev_wp;
  prev_wp = next_wp-1;
  if (next_wp == 0) {
    prev_wp  = map_xy.size()-1;
  }

  const Vector2d& p_prev = map_xy[prev_wp];

  Vector2d p_n = map_xy[next_wp] - p_prev;
  Vector2d p_x = xy - p_prev;

  // find the projection of x onto n
  Vector2d proj = p_n * (p_x.dot(p_n) / p_n.dot(p_n));

  double frenet_d = (p_x - proj).norm();

  //see if d value is positive or negative by comparing it to a center point
  static const Vector2d center_ref = [](){Vector2d temp; temp << 1000, 2000; return temp;}();
  Vector2d center = center_ref - p_prev;

  double centerToPos = (center - p_x).squaredNorm();
  double centerToRef = (center - proj).squaredNorm();
  
  if (centerToPos <= centerToRef) frenet_d *= -1;
  
  // calculate s value
  double frenet_s = 0;
  for (int i = 0; i < prev_wp; ++i) {
    frenet_s += (map_xy[i] - map_xy[i+1]).norm();
  }

  frenet_s += proj.norm();

  return {frenet_s,frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
Vector2d Map::getXY(double s, double d) {

  int prev_wp = -1;

  while (s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1))) 
    ++prev_wp;
  

  int wp2 = (prev_wp+1)%map_xy.size();

  double heading = angle(map_xy[wp2]-map_xy[prev_wp]);
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  Vector2d seg = map_xy[prev_wp] + seg_s*Vector2d{cos(heading), sin(heading)};

  double perp_heading = heading-pi()/2;

  return seg + d*Vector2d{cos(perp_heading), sin(perp_heading)};
}
