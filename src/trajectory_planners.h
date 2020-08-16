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
constexpr double kMaxVel = 20;
constexpr double k1 = -kMaxJerk/kMaxAcc, k2 = -k1*k1/4;

struct Car{
  double x, y, s, d, yaw, speed;
};

struct FeedbackVelocity{
  double k1=0, k2=0, target_speed=0;
  double position=0, velocity=0, accelaration=0, jerk=0;
  void step(double dt){
    jerk = k1*accelaration + k2*(velocity - target_speed);
    velocity += dt * (accelaration + dt/2 * jerk );
    accelaration += dt * jerk;
    position += dt * (velocity + dt/2 * (accelaration + dt/3* jerk));
  }
};

struct TrajectoryPlanners
{
  Map map;
  FeedbackVelocity fwd{k1, k2};
  double lane_center = 6;
  double d_initial_=6, dd_initial_=0, ddd_initial_=0;
  double theta_ = 0;
  
  void goto_lane(const Car& car, xy_path& previous_path, const std::vector<std::vector<double>>& other_cars);
};
