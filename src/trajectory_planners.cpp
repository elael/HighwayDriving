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
using Eigen::Vector2d;
using Vector6d = Eigen::Matrix<double, 6, 1>;

constexpr double kHalfLane = 2;
constexpr double kUpdatePeriod = 0.02;
constexpr uint kNsamples = 50;
constexpr double kTimeFrame = kNsamples * kUpdatePeriod;

constexpr double kMaxJerk = 10; // coefficients
constexpr double kMaxAcc = 10;
constexpr double k = 1;
constexpr double kMaxVel = 20;
constexpr double k1 = -kMaxJerk/kMaxAcc, k2 = -k1*k1/4;

struct MinPath
{
  Vector6d coefs;
  double operator()(const double& s){
    return coefs[0] + s*(coefs[1] + s*(coefs[2] + s*(coefs[3] + s*(coefs[4] + s* coefs[5]))));
  }
};



MinPath minimizer_path(const array<double,3>& starting_state, const array<double,3>& ending_state, double t){
  typedef Eigen::Matrix<double, 6, 6> Matrix6d;
  const auto change_lane_matrix = [=](){
    Matrix6d tmp;
    //start 
    tmp.row(0) << 1, 0, 0, 0,   0,   0; //pos
    tmp.row(1) << 0, 1, 0, 0,   0,  0; //vel
    tmp.row(2) << 0, 0, 2, 0,  0, 0; //acc
    //end
    double t2 = t*t, t3 = t2*t, t4 = t3*t;
    tmp.row(3) << 1, t, t2, t3,  t4,    t4*t; //pos
    tmp.row(4) << 0, 1, 2*t, 3*t2,   4*t3,  5*t4; //vel
    tmp.row(5) << 0, 0, 2, 6*t, 12*t2, 20*t3; //acc
    return tmp.fullPivLu();
  }();

  Vector6d constrains;
  constrains << starting_state[0], starting_state[1], starting_state[2], ending_state[0], ending_state[1], ending_state[2];
  
  return {change_lane_matrix.solve(constrains)};
}

/* Other cars
 * [0] = unique_id
 * [1] = x
 * [2] = y
 * [3] = vx
 * [4] = vy
 * [5] = s
 * [6] = d 
 */

void TrajectoryPlanners::goto_lane(const Car& car, array<vector<double>,2>& previous_path, const vector<vector<double>>& other_cars)
{
  const double current_time = kUpdatePeriod * previous_path[0].size();
  const double time_frame = kTimeFrame - current_time;
  if (time_frame <= 0) return;

  // Select two last positions on FrenetFrame
  double theta;
  if(previous_path[0].size()>=3){
    theta = std::atan2(previous_path[1].back() - *(previous_path[1].rbegin()+1),previous_path[0].back() - *(previous_path[0].rbegin()+1));
  }
  else
  {
    theta = deg2rad(car.yaw);
    v_ = 1;
    const auto initial = smooth_getXY(car.s + v_ * kUpdatePeriod, car.d, map_s, map_x, map_y);
    const auto initial2 = smooth_getXY(car.s + 2 * v_ * kUpdatePeriod, car.d, map_s, map_x, map_y);
    const auto initial3 = smooth_getXY(car.s + 3 * v_ * kUpdatePeriod, car.d, map_s, map_x, map_y);
    previous_path[0] = {initial[0], initial2[0], initial3[0]};
    previous_path[1] = {initial[1], initial2[1], initial3[1]};    
  }

  FrenetFrame p0 = getFrenet(previous_path[0].back(), previous_path[1].back(), theta, map_x, map_y);


  // Calculate distante to the next car
  double next_car_distance = 60;
  const std::vector<double>* closest_car_pt;
  for(const auto& other_car: other_cars){
    // predict constant speed for other cars
    const auto [other_s, other_d] = getFrenet(other_car[1] + current_time*other_car[3], other_car[2] + current_time*other_car[4], theta, 
                         map_x, 
                         map_y);
    if (double distance = (other_s - p0.s);
        0 < distance && distance < next_car_distance && fabs(other_d-lane_center) < kHalfLane ){ // in front
          next_car_distance = distance;
          closest_car_pt = &other_car;
        }
  }

  Vector2d initial; initial << previous_path[0].back(), previous_path[1].back();
  Vector2d behind_initial; behind_initial << *(previous_path[0].rbegin()+1), *(previous_path[1].rbegin()+1);
  Vector2d behind2_initial; behind2_initial << *(previous_path[0].rbegin()+2), *(previous_path[1].rbegin()+2);
  Vector2d v_initial = v_ * (initial-behind_initial).normalized();
  Vector2d v_behind_initial = v_ * (behind_initial-behind2_initial).normalized();
  Vector2d a_initial = (v_initial - v_behind_initial)/kUpdatePeriod/2;

  const auto [last_s, last_d] = getFrenet(previous_path[0].back(), previous_path[1].back(), theta, map_x, map_y);
  double final_ds = v_ * (1 + fabs(lane_center-last_d)/4.0);//kMaxVel * final_dt;
  const Vector2d ahead = smooth_getXY(last_s + final_ds + 0.1, lane_center, map_s, map_x, map_y);
  const Vector2d final = smooth_getXY(last_s + final_ds, lane_center, map_s, map_x, map_y);
  
  Vector2d v_final = v_ * (ahead-final).normalized();
  
  MinPath px = minimizer_path({initial[0], v_initial[0], a_initial[0]},{final[0], v_final[0],  0}, final_ds/v_);
  MinPath py = minimizer_path({initial[1], v_initial[1], a_initial[1]},{final[1], v_final[1],  0}, final_ds/v_);

  // fprintf(stderr, "py = {%f,%f,%f,%f,%f,%f}\n", py.coefs[0],py.coefs[1],py.coefs[2],py.coefs[3],py.coefs[4],py.coefs[5]);
  
  double x = px(0);
  double y = py(0);


  // If no car in front, use max jerk
  double target_speed;
  if(next_car_distance < 30){
    if (lane_center == 2) lane_center = 6;
    else if (lane_center == 6 && car.d < 6) lane_center = 10;
    else if (lane_center == 6 && car.d > 6) lane_center = 2;
    else if (lane_center == 10) lane_center = 6;
    
    target_speed = sqrt((*closest_car_pt)[3]*(*closest_car_pt)[3] + (*closest_car_pt)[4]*(*closest_car_pt)[4]);
    target_speed /= 30/next_car_distance;
    fprintf(stderr, "next_car_distance = %f, target_speed = %f\n", next_car_distance,  target_speed);
    }
  else
    target_speed = kMaxVel;
  

  for (double t = kUpdatePeriod,  n = 0; t < time_frame; t += kUpdatePeriod) {
    double sq_target_distance = kUpdatePeriod * (v_ + kUpdatePeriod/2 * (a_ + kUpdatePeriod/3* j_));
    sq_target_distance *= sq_target_distance;
    while (sq_distance(x,y,px(n),py(n)) < sq_target_distance) n += 1.0/10000;
    if (n > final_ds/v_) return;
    
    
    x = px(n);
    y = py(n);
    previous_path[0].emplace_back(x);
    previous_path[1].emplace_back(y);

    j_ = k1*a_ + k2*(v_ - target_speed);
    v_ += kUpdatePeriod * (a_ + kUpdatePeriod/2 * j_ );
    a_ += kUpdatePeriod * j_;
  }

}

