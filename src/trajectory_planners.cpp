#include "trajectory_planners.h"

#include <cmath>
#include <limits>
#include <tuple>

#include "Eigen/Dense"
#include "helpers.h"
#include "map.h"
#include "spline.h"

/**
 *  define a path made up of (x,y) points that the car will visit
 *   sequentially every .02 seconds
 */
using Eigen::Vector2d;
using std::array;
using std::vector;
using Vector6d = Eigen::Matrix<double, 6, 1>;

MinPath minimizer_path(const array<double, 3> &starting_state, const array<double, 3> &ending_state, double t) {
  typedef Eigen::Matrix<double, 6, 6> Matrix6d;
  const auto change_lane_matrix = [=]() {
    Matrix6d tmp;
    //start
    tmp.row(0) << 1, 0, 0, 0, 0, 0;  //pos
    tmp.row(1) << 0, 1, 0, 0, 0, 0;  //vel
    tmp.row(2) << 0, 0, 2, 0, 0, 0;  //acc
    //end
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t3 * t;
    tmp.row(3) << 1, t, t2, t3, t4, t4 * t;             //pos
    tmp.row(4) << 0, 1, 2 * t, 3 * t2, 4 * t3, 5 * t4;  //vel
    tmp.row(5) << 0, 0, 2, 6 * t, 12 * t2, 20 * t3;     //acc
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

inline uint lane2ind(double lane) {
  lane = lane > 0 ? lane : 0;
  lane = lane < 12 ? lane : 8;
  return static_cast<uint>(lane) / 4;
}

double TrajectoryPlanners::get_cost(const Car &car, double target_lane, const vector<double> &s_path, const MinPath &candidate, double current_time, const std::vector<std::vector<double>> &other_cars) {
  double max_d_jerk = 0;
  const double dt = kUpdatePeriod * s_path.size();
  for (double t = 0; t <= dt; t += 5 * kUpdatePeriod) {
    double d_jerk = fabs(candidate.jerk(t));
    if (d_jerk > kMaxJerk)
      return std::numeric_limits<double>::infinity();
    if (d_jerk > max_d_jerk)
      max_d_jerk = d_jerk;

    if (candidate(t) < 1 || candidate(t) > 11)
      return std::numeric_limits<double>::infinity();
  }

  double cost = 10 * fabs(max_d_jerk) / kMaxJerk;  //lateral jerk

  cost += 10 * exp((target_lane - candidate(dt)) * (target_lane - candidate(dt)));  // center line distance cost
  cost += 10 * exp(candidate.speed(dt) * candidate.speed(dt));                      // center line alignment cost
  cost += 10 * exp(candidate.acc(dt) * candidate.acc(dt));                         // wobbling cost

  double min_dist = std::numeric_limits<double>::infinity();
  for (const auto &other_car : other_cars) {
    const auto [other_s, other_d] = map.getFrenet(Vector2d(other_car[1], other_car[2]) + (current_time + dt) * Vector2d(other_car[3], other_car[4]), car.yaw);
    double distance = fabs(other_s - fwd.state.pos);
    if (distance < min_dist && fabs(other_d - target_lane) < 1.5 * kHalfLane)
      min_dist = distance;
  }

  const double s_tolerance = t_tolerance * s_path.back();
  double safety_cost = s_tolerance / min_dist;
  safety_cost *= safety_cost;
  cost += 30 * safety_cost;  // safety cost

  return cost;
}

MinPath TrajectoryPlanners::select_candidate(const Car &car, double target_lane, const vector<double> &s_path, vector<MinPath> candidates, double current_time, const std::vector<std::vector<double>> &other_cars) {
  size_t min_i = 0;
  double min_cost = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < candidates.size(); ++i) {
    double cost = get_cost(car, target_lane, s_path, candidates[i], current_time, other_cars);

    if (cost < min_cost) {
      min_cost = cost;
      min_i = i;
    }
  }

  return candidates[min_i];
}

void TrajectoryPlanners::goto_lane(const Car &car, double target_lane, xy_path &previous_path, const std::vector<std::vector<double>> &other_cars) {
  constexpr double dt_step = kUpdatePeriod;
  constexpr uint dt_samples = 100;

  const double current_time = kUpdatePeriod * previous_path[0].size();
  const double time_frame = kTimeFrame - current_time;
  const double initial_dt = 1.5;
  const double final_dt = dt_samples * dt_step + initial_dt;

  // Lateral movement computation
  vector<MinPath> candidates;
  candidates.reserve(dt_samples);
  for (double dt = initial_dt; dt < final_dt; dt += dt_step) {
    candidates.emplace_back(minimizer_path({lat.pos, lat.vel, lat.acc}, {target_lane, 0, 0}, dt));
  }

  // Forward movement computation
  vector<double> s_path;
  s_path.reserve(kNsamples - previous_path[0].size());
  for (size_t n = previous_path[0].size(); n < kNsamples; ++n) {
    fwd.step(kUpdatePeriod);
    s_path.emplace_back(fwd.state.pos);
  }

  // Compound movement selection
  const MinPath &final_path = select_candidate(car, target_lane, s_path, std::move(candidates), current_time, other_cars);

  // Update path with new steps
  double t;
  size_t n;
  for (t = kUpdatePeriod, n = 0; n < s_path.size(); t += kUpdatePeriod, ++n) {
    Vector2d xy = map.smooth_getXY(s_path[n], final_path(t));
    previous_path[0].emplace_back(xy[0]);
    previous_path[1].emplace_back(xy[1]);
  }
  t -= kUpdatePeriod;

  // Update selected lateral movement
  lat = {final_path(t), final_path.speed(t), final_path.acc(t), final_path.jerk(t)};
}

void TrajectoryPlanners::select_state(const Car &car, xy_path &previous_path, const std::vector<std::vector<double>> &other_cars) {
  if (previous_path[0].size() == 0) {
    const auto initial = map.getFrenet(Vector2d(car.x, car.y), car.yaw);
    fwd.state.pos = initial.s;
    lat.pos = initial.d;
  }

  const double target_vel = kMaxVel * 0.95;

  const double current_time = kUpdatePeriod * previous_path[0].size();
  const double time_frame = kTimeFrame - current_time;

  // Lanes speed and distance computation for all lanes
  array lane_speed = {target_vel, target_vel, target_vel};
  array lane_position = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
  for (const auto &other_car : other_cars) {
    // predict constant speed for other cars
    const auto [other_s, other_d] = map.getFrenet(Vector2d(other_car[1], other_car[2]) + kTimeFrame * Vector2d(other_car[3], other_car[4]), car.yaw);
    double speed = Vector2d(other_car[3], other_car[4]).norm();
    const double horizon = t_tolerance * fwd.state.vel;

    if (double distance = (other_s - fwd.state.pos - (fwd.state.vel - speed) * time_frame); - 10 < distance && distance < horizon)
      if (uint index = lane2ind(other_d); speed < lane_speed[index]) {
        lane_speed[index] = speed;
        lane_position[index] = other_s - horizon / 2;
      }
  }

  uint target_lane;
  // Keep lane if it is the fastest
  if (lane_speed[ergo_lane] >= lane_speed[(ergo_lane + 1) % 3] && lane_speed[ergo_lane] >= lane_speed[(ergo_lane + 2) % 3])
    target_lane = ergo_lane;
  // "State Machine" for lane selection
  else switch (ergo_lane) {
    case 0:
      target_lane = lane_position[1] > t_tolerance * fwd.state.vel + fwd.state.pos ? 1 : 0;
      break;
    case 1:
      target_lane = lane_speed[2] > lane_speed[0] ? 2 : 0;
      target_lane = lane_position[target_lane] > t_tolerance * fwd.state.vel + fwd.state.pos ? target_lane : 1;
      break;
    case 2:
      target_lane = lane_position[1] > t_tolerance * fwd.state.vel + fwd.state.pos ? 1 : 2;
      break;
  }

  ergo_lane = target_lane;
  double ds = lane_position[target_lane] - fwd.state.pos;
  ds = ds > 0 ? ds : 0;
  fwd.target_speed = lane_speed[target_lane] * (1 - std::exp(-ds / 30));
  return goto_lane(car, 2 + 4.0 * target_lane, previous_path, other_cars);
}
