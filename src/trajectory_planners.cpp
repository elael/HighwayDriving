#include "trajectory_planners.h"
#include "helpers.h"
#include <limits>
#include <cmath>
#include <tuple>
#include "Eigen/Dense"

/**
 *  define a path made up of (x,y) points that the car will visit
 *   sequentially every .02 seconds
 */
using std::array;
using std::vector;

constexpr double kHalfLane = 2;
constexpr double kUpdatePeriod = 0.02;
constexpr uint kNsamples = 50;
constexpr double kTimeFrame = kNsamples * kUpdatePeriod;

constexpr double kMaxJerk = 9; // coefficients
constexpr double kMaxAcc = 9;
constexpr double kMaxVel = 20;

/* Other cars
 * [0] = x
 * [1] = y
 * [2] = vx
 * [3] = vy
 * [4] = s
 * [5] = d 
 */


void TrajectoryPlanners::keep_lane(const Car& car,  
                                   array<vector<double>,2>& previous_path, 
                                   const vector<vector<double>>& other_cars)
{
  const double current_time = kUpdatePeriod * previous_path[0].size();
  const double time_frame = kTimeFrame - current_time;
  if (time_frame <= 0) return;
  
  const double theta = deg2rad(car.yaw);

  // Select two last positions on FrenetFrame
  FrenetFrame p0, p_old;
  if(previous_path[0].size() < 4)
  for (double t = current_time; t < kTimeFrame; t += kUpdatePeriod) {
    previous_path[0].emplace_back(car.x + t*car.speed * 0.44704*cos(theta));
    previous_path[1].emplace_back(car.y + t*car.speed * 0.44704*sin(theta));
  }
  
  
  
  // if(previous_path.size()>=4)
  // {
  //   p0 = getFrenet(previous_path[0].back(), previous_path[1].back(), theta, map_x, map_y);
  //   p_old = getFrenet(*(previous_path[0].rbegin()+1), *(previous_path[1].rbegin()+1), theta, map_x, map_y);
  // }
  // else{
  //   p0 = {car.s, car.d};
  //   p_old = getFrenet(car.x - car.speed*cos(theta)*kUpdatePeriod, car.y - car.speed*sin(theta)*kUpdatePeriod, theta, map_x, map_y); 
  // }

  // Calculate distante to the next car
  double next_car_distance = std::numeric_limits<double>::infinity();
  for(const auto& other_car: other_cars){
    // predict constant speed for other cars
    const auto [other_s, other_d] = getFrenet(other_car[0] + current_time*other_car[2], other_car[1] + current_time*other_car[3], theta, 
                         map_x, 
                         map_y);
    if (double distance = (other_s - p0.s);
        0 < distance && distance < next_car_distance && fabs(other_d-p0.d) < kHalfLane ) // in front
      next_car_distance = distance;
  }

  static const Eigen::FullPivLU<Eigen::Matrix3d> max_jerk_matrix = [] {
    double t = -kUpdatePeriod;
    Eigen::Matrix3d tmp;
    for (size_t n = 0; n < 3; ++n, t -= kUpdatePeriod)
      tmp.row(n) << t, t*t, t*t*t;
      
    return tmp.fullPivLu();
  }();

  // If no car in front, use max jerk
  if(next_car_distance == std::numeric_limits<double>::infinity())
  { 
            
    auto prev_x = previous_path[0].rbegin();
    const auto& x = previous_path[0].back();
    Eigen::Vector3d xb;    
    xb << (*(prev_x+1) - x), (*(prev_x+2) - x), (*(prev_x+3) - x);
    auto x_coef = max_jerk_matrix.solve(xb);
    double x_jerk_coef = std::min({2*((kMaxVel-x_coef[0])/time_frame - 2*x_coef[1])/time_frame, 
                                   (kMaxAcc-2*x_coef[1])/time_frame, 
                                   kMaxJerk}) / 6;

    
    const auto& y = previous_path[1].back();
    const auto prev_y = previous_path[1].rbegin();
    Eigen::Vector3d yb;
    yb << (*(prev_y+1) - y), (*(prev_y+2) - y), (*(prev_y+3) - y);
    auto y_coef = max_jerk_matrix.solve(yb);
    double y_jerk_coef = std::min({2*((kMaxVel-y_coef[0])/time_frame - 2*y_coef[1])/time_frame, 
                                   (kMaxAcc-2*y_coef[1])/time_frame, 
                                   kMaxJerk}) / 6;

    for (double t = kUpdatePeriod; t < time_frame; t += kUpdatePeriod) {
      previous_path[0].emplace_back(x + t*(x_coef[0] +t*(x_coef[1] + t*x_jerk_coef)));
      previous_path[1].emplace_back(y + t*(y_coef[0] +t*(y_coef[1] + t*y_jerk_coef)));
    }
  }
  
  // static Eigen::Matrix4d keep_lane_matrix = [] {
  //   const auto& T = kTimeFrame;
  //   Eigen::Matrix4d tmp;
  //   tmp.row(0) << T*T, T*T*T, T*T*T*T;
  //   tmp << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
  //   return tmp;
  // }();

}

