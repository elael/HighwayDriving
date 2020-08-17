#include <vector>
#include <array>
#include "map.h"

/**
 *  define a path made up of (x,y) points that the car will visit
 *   sequentially every .02 seconds
 */
using xy_path = std::array<std::vector<double>,2>;


constexpr double kHalfLane = 2;
constexpr double kUpdatePeriod = 0.02;
constexpr uint kNsamples = 50;
constexpr double kTimeFrame = kNsamples * kUpdatePeriod;

constexpr double kMaxJerk = 10; // coefficients
constexpr double kMaxAcc = 10;
constexpr double k = 1;
constexpr double kMaxVel = 22;
constexpr double k1 = -kMaxJerk/kMaxAcc, k2 = -k1*k1/4;

struct Car{
  double x, y, s, d, yaw, speed;
};

struct KineState{
  double pos=0, vel=0, acc=0, jerk=0;
};

struct FeedbackVelocity{
  double k1=0, k2=0, target_speed=0;
  KineState state;
  void step(double dt){
    state.jerk = k1*state.acc + k2*(state.vel - target_speed);
    state.vel += dt * (state.acc + dt/2 * state.jerk );
    state.acc += dt * state.jerk;
    state.pos += dt * (state.vel + dt/2 * (state.acc + dt/3* state.jerk));
  }
};

struct MinPath
{
  Eigen::Matrix<double, 6, 1> coefs;
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

struct MinFrenetPath
{
  MinPath s;
  MinPath d;
  double dt;
  MinFrenetPath(MinPath&& s, MinPath&& d, double dt): s(std::move(s)), d(std::move(d)), dt(dt){}
};

struct TrajectoryPlanners
{
  Map map;
  FeedbackVelocity fwd{k1, k2};
  KineState lat{6,0,0,0};
  uint ergo_lane = 1;
  double theta_ = 0;
  
  void goto_lane(const Car& car, xy_path& previous_path, const std::vector<std::vector<double>>& other_cars);
  void keep_speed(std::array<std::vector<double>,2>& previous_path);
  void follow_car(Car lead_car, std::array<std::vector<double>,2>& previous_path);
  double get_cost(const Car& car, double target_lane, const std::vector<double>& s_path, const MinPath& candidate, double current_time, const std::vector<std::vector<double>>& other_cars);
  MinPath select_candidate(const Car& car, double target_lane, const std::vector<double>& s_path, std::vector<MinPath> candidates, double current_time, const std::vector<std::vector<double>>& other_cars);
  void change_lane(const Car& car, double target_lane, xy_path& previous_path, const std::vector<std::vector<double>>& other_cars);
  void select_state(const Car& car, xy_path& previous_path, const std::vector<std::vector<double>>& other_cars);
};
