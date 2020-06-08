#include <vector>
#include <array>

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
  std::vector<double> map_x, map_y, map_s, map_dx, map_dy;
  double v_ = 0, j_ = 0.001, a_ = 0;
  
  void keep_lane(const Car& car, xy_path& previous_path, const std::vector<std::vector<double>>& other_cars);

};
