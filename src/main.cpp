#include <uWS/uWS.h>

#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Eigen/Core"
#include "Eigen/QR"
#include "helpers.h"
#include "json.hpp"
#include "map.h"
#include "trajectory_planners.h"

// for convenience
using nlohmann::json;
using std::string;
using std::vector;

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  TrajectoryPlanners planner;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  std::ifstream in_map_(map_file_.c_str(), std::ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    std::istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    planner.map.map_xy.emplace_back(x, y);
    planner.map.maps_s.emplace_back(s);
    planner.map.map_dxdy.emplace_back(d_x, d_y);
  }
  planner.map.smoothMap();

  h.onMessage([&planner](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                         uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {
      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

          // Main car's localization Data
          Car main_car{
              j[1]["x"],
              j[1]["y"],
              j[1]["s"],
              j[1]["d"],
              deg2rad(j[1]["yaw"]),
              mph2ms(j[1]["speed"])};

          // Previous path data given to the Planner
          std::array<vector<double>, 2> previous_path = {j[1]["previous_path_x"], j[1]["previous_path_y"]};

          // Previous path's end s and d values
          std::array<double, 2> end_path_frenet = {j[1]["end_path_s"], j[1]["end_path_d"]};

          // Sensor Fusion Data, a list of all other cars on the same side
          //   of the road.
          vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];

          /**
           * TODO: define a path made up of (x,y) points that the car will visit
           *   sequentially every .02 seconds
           */
          // double dist_inc = 0.4;
          // for (int i = 0; i < 50; ++i) {
          //   next_x_vals.emplace_back(main_car.x+(dist_inc*i)*cos(deg2rad(main_car.yaw)));
          //   next_y_vals.emplace_back(main_car.y+(dist_inc*i)*sin(deg2rad(main_car.yaw)));
          // }
          planner.select_state(main_car, previous_path, sensor_fusion);

          json msgJson;
          msgJson["next_x"] = std::move(previous_path[0]);
          msgJson["next_y"] = std::move(previous_path[1]);

          auto msg = "42[\"control\"," + msgJson.dump() + "]";

          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }  // end "telemetry" if
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }  // end websocket if
  });  // end h.onMessage

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }

  h.run();
}
