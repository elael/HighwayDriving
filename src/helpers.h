#ifndef HELPERS_H
#define HELPERS_H

#include <math.h>
#include <string>
#include <vector>
#include "Eigen/Dense"

// for convenience
using std::string;
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
constexpr double deg2rad(double x) { return x * pi() / 180; }
constexpr double rad2deg(double x) { return x * 180 / pi(); }
constexpr double mph2ms(double x) { return x * 0.44704; }

// Calculate distance between two points
constexpr double distance(const double& x1, const double& y1, const double& x2, const double& y2) {
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

// Calculate distance between two points
constexpr double sq_distance(const double& x1, const double& y1, const double& x2, const double& y2) {
  return (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);
}

inline double angle(const Vector2d& xy){
  return atan2(xy[1],xy[0]);
}

inline Vector2d unitVector(double theta){
  return Vector2d(cos(theta),sin(theta));
}

// counter clockwise perpendicular vector
inline Vector2d ccPerp(const Vector2d& original){
    Vector2d perp; perp << -original[1], original[0];
    return perp;
}
#endif  // HELPERS_H
