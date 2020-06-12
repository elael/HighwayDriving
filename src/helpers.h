#ifndef HELPERS_H
#define HELPERS_H

#include <math.h>
#include <string>
#include <array>
#include <vector>
#include "Eigen/Dense"

// for convenience
using std::string;
using std::vector;
using Eigen::Vector2d;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
//   else the empty string "" will be returned.
static string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

//
// Helper functions related to waypoints and converting from XY to Frenet
//   or vice versa
//

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
constexpr double deg2rad(const double& x) { return x * pi() / 180; }
constexpr double rad2deg(const double& x) { return x * 180 / pi(); }

// Calculate distance between two points
constexpr double distance(const double& x1, const double& y1, const double& x2, const double& y2) {
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

// Calculate distance between two points
constexpr double sq_distance(const double& x1, const double& y1, const double& x2, const double& y2) {
  return (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);
}

// Calculate closest waypoint to current x, y position
static int ClosestWaypoint(double x, double y, const vector<double> &maps_x, 
                    const vector<double> &maps_y) {
  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for (int i = 0; i < maps_x.size(); ++i) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }
  }

  return closestWaypoint;
}

// Returns next waypoint of the closest waypoint
static int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, 
                 const vector<double> &maps_y) {
  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y-y),(map_x-x));

  double angle = fabs(theta-heading);
  angle = std::min(2*pi() - angle, angle);

  if (angle > pi()/2) {
    ++closestWaypoint;
    if (closestWaypoint == maps_x.size()) {
      closestWaypoint = 0;
    }
  }

  return closestWaypoint;
}

struct FrenetFrame {
  double s, d;
};

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
inline FrenetFrame getFrenet(double x, double y, double theta, 
                         const vector<double> &maps_x, 
                         const vector<double> &maps_y) {
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if (next_wp == 0) {
    prev_wp  = maps_x.size()-1;
  }

  Vector2d p_prev; p_prev << maps_x[prev_wp], maps_y[prev_wp];

  Vector2d p_n; p_n << maps_x[next_wp], maps_y[next_wp];
  Vector2d p_x; p_x << x, y;
  p_n -= p_prev;
  p_x -= p_prev;

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
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += proj.norm();

  return {frenet_s,frenet_d};
}

struct XYFrame {
  double x, y;
};
// Transform from Frenet s,d coordinates to Cartesian x,y
inline Vector2d getXY(double s, double d, const vector<double> &maps_s, 
                     const vector<double> &maps_x, 
                     const vector<double> &maps_y) {
  int prev_wp = -1;

  while (s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1))) {
    ++prev_wp;
  }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),
                         (maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  Eigen::Vector2d xy;
  xy << x, y;
  return xy;
}


inline Vector2d smooth_getXY(double s, double d, const vector<double> &maps_s, 
                     const vector<double> &maps_x, 
                     const vector<double> &maps_y){
  static const double smoother_s = 10;
  return (getXY(s + smoother_s, d, maps_s, maps_x, maps_y) + getXY(s - smoother_s, d, maps_s, maps_x, maps_y))/2;
}

#endif  // HELPERS_H
