#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <sys/time.h>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

const int latency_sec = 0.1; /* 100 milliseconds */

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

double polyeval_gradiant(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 1; i < coeffs.size(); i++) {
    result += i * coeffs[i] * pow(x, i-1);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;
  long double predicted_latency_sec = latency_sec;
  double avg_steering_angle = 0;

  h.onMessage([&mpc, &predicted_latency_sec, &avg_steering_angle](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    struct timeval tp; // Keep track of how long we spend in the function
    gettimeofday(&tp, NULL);
    long long starttime_ms = (long long)(tp.tv_sec * 1000) + tp.tv_usec / 1000;

    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object

          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double angle = j[1]["steering_angle"];
          double throttle = j[1]["throttle"];
          vector<double> ptsx_car;
          vector<double> ptsy_car;
          ptsx_car.reserve(ptsx.size());
          ptsy_car.reserve(ptsy.size());


          // Predict the state of the car in the future by approximately
          // 100ms, to compensate for latency.
          // We make the assumption that our latency now will be the latency
          // that we saw in the previous frame.
          px += v * cos(psi) * predicted_latency_sec;
          py += v * sin(psi) * predicted_latency_sec;
          psi += v / Lf * (-angle) * predicted_latency_sec;
          v += throttle*predicted_latency_sec;

          for (int i = 0; i < ptsx.size(); i++) {
            const double x = (ptsx[i] - px);
            const double y = (ptsy[i] - py);
            // Move the x back to compensate for the amount of latency that we had in the
            // last time we called this function.
            ptsx_car.push_back(x * cos(psi) + y * sin(psi));
            ptsy_car.push_back(-x * sin(psi) + y * cos(psi));
          }

          Eigen::VectorXd ptsx_ = Eigen::VectorXd::Map(ptsx_car.data(), ptsx_car.size());
          Eigen::VectorXd ptsy_ = Eigen::VectorXd::Map(ptsy_car.data(), ptsy_car.size());

          auto coeffs = polyfit(ptsx_, ptsy_, 3);
          double cte = -coeffs[0];
          double epsi = atan(-polyeval_gradiant(coeffs, 0));

          Eigen::VectorXd state(6);
          state << 0, 0, 0, v, cte, epsi;

          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          vector<double> solution = mpc.Solve(state, coeffs, mpc_x_vals, mpc_y_vals);
          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          const double steer_value = -solution[0] / deg2rad(25);
          const double throttle_value = solution[1];

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          const double D = 0.4;
          avg_steering_angle = avg_steering_angle * D + steer_value * (1-D);
          msgJson["steering_angle"] = avg_steering_angle;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          /* Draw the raw data */
          //next_x_vals = ptsx_car;
          //next_y_vals = ptsy_car;

          // Or draw the poly best-fit line
          for(int x = 0; x < 100; x+= 10) {
            const double y = polyeval(coeffs, x);
            next_x_vals.push_back(x);
            next_y_vals.push_back(y);
          }

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;



          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(latency_sec*1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

          gettimeofday(&tp, NULL);
          long long endtime_ms = (long long)(tp.tv_sec * 1000) + tp.tv_usec / 1000;
          const double C = 0.5; // low pass filter
          const double time_in_function_sec = (endtime_ms - starttime_ms)/1000.0;
          predicted_latency_sec = C * predicted_latency_sec + (1-C) * time_in_function_sec;
          //cout << "Time diff is" << time_in_function_sec << endl;
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

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
