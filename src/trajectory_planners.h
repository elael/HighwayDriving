#include <vector>
#include <array>
#include "map.h"

/**
 *  define a path made up of (x,y) points that the car will visit
 *   sequentially every .02 seconds
 */
using xy_path = std::array<std::vector<double>,2>;


struct Car{
  double x, y, s, d, yaw, speed;
};

struct TrajectoryPlanners
{
  Map map;
  double v_ = 0, j_ = 0.001, a_ = 0, lane_center = 6;
  
  void goto_lane(const Car& car, xy_path& previous_path, const std::vector<std::vector<double>>& other_cars);
};
