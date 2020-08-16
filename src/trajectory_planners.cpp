#include "trajectory_planners.h"
#include "helpers.h"
#include "map.h"
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

struct MinPath
{
  Vector6d coefs;
  double operator()(double t) const{
    return coefs[0] + t*(coefs[1] + t*(coefs[2] + t*(coefs[3] + t*(coefs[4] + t* coefs[5]))));
  }
  double speed(double t) const{
    return coefs[1] + t*(2*coefs[2] + t*(3*coefs[3] + t*(4*coefs[4] + t*5*coefs[5])));
  }
  double acc(double t) const{
    return 2*coefs[2] + t*(2*3*coefs[3] + t*(3*4*coefs[4] + t*4*5*coefs[5]));
  }
  double jerk(double t) const{
    return 2*3*coefs[3] + t*(2*3*4*coefs[4] + t*3*4*5*coefs[5]);
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
    double t2 = t*t;
    double t3 = t2*t;
    double t4 = t3*t;
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

inline uint lane2ind(double lane){
  lane = lane > 0? lane : 0;
  lane = lane < 12? lane : 8;
  return static_cast<uint>(lane) / 4; 
}

void TrajectoryPlanners::goto_lane(const Car& car, array<vector<double>,2>& previous_path, const vector<vector<double>>& other_cars)
{
  double theta = deg2rad(car.yaw);

  if(previous_path[0].size() == 0){
    fwd.position = map.getFrenet(Vector2d(car.x, car.y), theta).s;
  }
  const double current_time = kUpdatePeriod * previous_path[0].size();
  const double time_frame = kTimeFrame - current_time;
  if (time_frame <= kUpdatePeriod) return;


  // // Select two last positions on FrenetFrame
  // double speed = mph2ms(car.speed);
  // Vector2d velocity = speed * unitVector(theta);
  // Vector2d initial(car.x, car.y);


  // Calculate distante to the next car
  double next_car_distance = 30;
  array<double, 3> lane_speed = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
  const std::vector<double>* closest_car_pt;
  for(const auto& other_car: other_cars){
    // predict constant speed for other cars
    const auto [other_s, other_d] = map.getFrenet(Vector2d(other_car[1],other_car[2])  + current_time*Vector2d(other_car[3], other_car[4]), theta);
    if (double distance = (other_s - fwd.position);
        -10 < distance && distance < next_car_distance && fabs(other_d-lane_center) < 1.2*kHalfLane ){ // in front
          next_car_distance = distance;
          closest_car_pt = &other_car;
        }
  }

  double final_dt = 3*fwd.velocity; //2*fabs(lane_center - d_initial_)/fwd.velocity;
  
  MinPath pd = minimizer_path({d_initial_, dd_initial_, ddd_initial_},{lane_center, 0,  0}, final_dt);


  // If no car in front, use max jerk
  if(next_car_distance < 30){
    if (lane_center == 2) lane_center = 6;
    else if (lane_center == 6 && car.d < 6) lane_center = 10;
    else if (lane_center == 6 && car.d > 6) lane_center = 2;
    else if (lane_center == 10) lane_center = 6;
    
    double other_speed = Vector2d((*closest_car_pt)[3], (*closest_car_pt)[4]).norm();
    fwd.target_speed = other_speed;
    fwd.target_speed /= 30/next_car_distance;
    fprintf(stderr, "next_car_distance = %f, other_speed = %f, target_speed = %f\n", next_car_distance, other_speed, fwd.target_speed);
    }
  else
    fwd.target_speed = kMaxVel;

  double t;
  const double initial_s = fwd.position;
  for (t = kUpdatePeriod; t < time_frame; t += kUpdatePeriod) {
    fwd.step(kUpdatePeriod);
    double s = fwd.position - initial_s;

    Vector2d xy = map.smooth_getXY(fwd.position, pd(s));//6+4*sin(2*M_PI/50*s));
    
    fprintf(stderr, "At %f: s=%f, d=%f, velocity=%f, accelaration=%f, jerk=%f\n", 
    t, fwd.position, pd(s), fwd.velocity, fwd.accelaration, fwd.jerk);
    
    previous_path[0].emplace_back(xy[0]);
    previous_path[1].emplace_back(xy[1]);
  }
  t -= kUpdatePeriod;
  double s = fwd.position - initial_s;

  fprintf(stderr, "pd(t)=%f, pd.speed(t)=%f, pd.acc(t)=%f\n\n", pd(s), pd.speed(s), pd.acc(s));
  d_initial_ = pd(s);
  dd_initial_ = pd.speed(s);
  ddd_initial_ = pd.acc(s);
}
