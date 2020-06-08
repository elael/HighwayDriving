#include "trajectory_planners.h"
#include "helpers.h"
#include <limits>
#include <cmath>
#include <tuple>
#include "Eigen/Dense"
#include "spline.h"

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

constexpr double kMaxJerk = 10; // coefficients
constexpr double kMaxAcc = 10;
constexpr double kMaxVel = 20;
constexpr double k1 = -kMaxJerk/kMaxAcc, k2 = -k1*k1/4;

/* Other cars
 * [0] = unique_id
 * [1] = x
 * [2] = y
 * [3] = vx
 * [4] = vy
 * [5] = s
 * [6] = d 
 */


void TrajectoryPlanners::keep_lane(const Car& car,  
                                   array<vector<double>,2>& previous_path, 
                                   const vector<vector<double>>& other_cars)
{
  const double current_time = kUpdatePeriod * previous_path[0].size();
  const double time_frame = kTimeFrame - current_time;
  if (time_frame <= 0) return;

  // Select two last positions on FrenetFrame
  double theta;
  FrenetFrame p0;
  if(previous_path[0].size()>=2){
    theta = std::atan2(previous_path[1].back() - *(previous_path[1].rbegin()+1),previous_path[0].back() - *(previous_path[0].rbegin()+1));
    p0 = getFrenet(previous_path[0].back(), previous_path[1].back(), theta, map_x, map_y);
  }
  else
  {
    theta = deg2rad(car.yaw);
    p0 = {car.s, car.d};
    previous_path[0] = {car.x, car.x};
    previous_path[1] = {car.y, car.y};
  }


  // Calculate distante to the next car
  double next_car_distance = 30;
  const std::vector<double>* closest_car_pt;
  for(const auto& other_car: other_cars){
    // predict constant speed for other cars
    const auto [other_s, other_d] = getFrenet(other_car[1] + current_time*other_car[3], other_car[2] + current_time*other_car[4], theta, 
                         map_x, 
                         map_y);
    if (double distance = (other_s - p0.s);
        0 < distance && distance < next_car_distance && fabs(other_d-p0.d) < kHalfLane ){ // in front
          next_car_distance = distance;
          closest_car_pt = &other_car;
        }
        
        
  }

  constexpr double final_ds = 30;
  const auto [last_s, last_d] = getFrenet(previous_path[0].back(), previous_path[1].back(), theta, map_x, map_y);
  const double center_d = 4*std::floor(last_d/4)+2;
  const auto [mid_x,mid_y] = getXY(last_s + final_ds/2, center_d, map_s, map_x, map_y);
  const auto [final_x,final_y] = getXY(last_s + final_ds, center_d, map_s, map_x, map_y);
  const auto [ahead_x,ahead_y] = getXY(last_s + final_ds + 0.01, center_d, map_s, map_x, map_y);

  tk::spline sx;
  sx.set_boundary(tk::spline::first_deriv, previous_path[0].back() - *(previous_path[0].rbegin()+1),
                  tk::spline::first_deriv, ahead_x - final_x);
  sx.set_points({0, final_ds/2, final_ds},{previous_path[0].back(), mid_x, final_x});

  tk::spline sy;
  sy.set_boundary(tk::spline::first_deriv, previous_path[1].back() - *(previous_path[1].rbegin()+1),
                  tk::spline::first_deriv, ahead_y - final_y);
  sy.set_points({0, final_ds/2, final_ds},{previous_path[1].back(), mid_y, final_y});

  
  double x = sx(0);
  double y = sy(0);



  // If no car in front, use max jerk
  double target_speed;
  if(next_car_distance < 30){
    target_speed = sqrt((*closest_car_pt)[3]*(*closest_car_pt)[3] + (*closest_car_pt)[4]*(*closest_car_pt)[4]);
    target_speed /= 30/next_car_distance;
    if(target_speed < 0) target_speed = 0;
    fprintf(stderr, "next_car_distance = %f, target_speed = %f\n", next_car_distance,  target_speed);
    }
  else
    target_speed = kMaxVel;
  

  for (double t = kUpdatePeriod,  n = 0; t < time_frame; t += kUpdatePeriod) {
    double sq_target_distance = kUpdatePeriod * (v_ + kUpdatePeriod/2 * (a_ + kUpdatePeriod/3* j_));
    sq_target_distance *= sq_target_distance;

    while (sq_distance(x,y,sx(n),sy(n)) < sq_target_distance) n += final_ds/10000;
    x = sx(n);
    y = sy(n);
    previous_path[0].emplace_back(x);
    previous_path[1].emplace_back(y);

    j_ = k1*a_ + k2*(v_ - target_speed);
    v_ += kUpdatePeriod * (a_ + kUpdatePeriod/2 * j_ );
    a_ += kUpdatePeriod * j_;
  }

}

