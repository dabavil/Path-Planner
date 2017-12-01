#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  const double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
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
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  //Set reference velocity in m/s
  double reference_velocity = 0.0;
  double const MAX_SPEED = 21.5;

  //Set reference lane: 0 for left lane, 1 for middle, 2 for right lane
  int lane = 1;
  

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&lane,&reference_velocity,&MAX_SPEED](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;


    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];


            //SITUATIONAL AWARENESS - get clarity on what's what: which lanes are safe to drive? What is the 'speed' that would be possible in each lane?

            



            bool car_left = false;
            bool car_mid = false;
            bool car_right = false;

            double left_maxspeed = MAX_SPEED;
            double mid_maxspeed = MAX_SPEED;
            double right_maxspeed = MAX_SPEED;


              //Cycle through all cars from the sensor fusion vector, find out lane, distance to us, speed difference
              
              //TODO - delete these old variables
              bool front_warning = false;
              bool change_lane;

              for(int i = 0; i < sensor_fusion.size(); i++)
              {

                //other car's d pose
                double other_car_d = sensor_fusion[i][6];
                double d_distance = other_car_d - end_path_d;

                //other car's velocity
                double other_car_vx = sensor_fusion[i][3]; 
                double other_car_vy = sensor_fusion[i][4];
                double other_car_speed = sqrt(other_car_vx*other_car_vx + other_car_vy*other_car_vy);

                //other car's s pose
                double other_car_s = sensor_fusion[i][5]; 
                double other_car_s_future = other_car_s + previous_path_x.size() * other_car_speed * 0.02;
                double s_distance = other_car_s - end_path_s;


                //determine other car's lane
                int other_car_lane = -1;

                //TODO: find more elegant math, e.g., divide by 4 and round to get lane number 
                if (other_car_d > 0 && other_car_d < 4) 
                  {other_car_lane = 0;}
                if (other_car_d > 4 && other_car_d < 8) 
                  {other_car_lane = 1;}
                if (other_car_d > 8 && other_car_d < 12) 
                  {other_car_lane = 2;}

                //debug messages
                //std::cout << "Other car " << sensor_fusion[i][0] << " is in lane " << other_car_lane << " at a distance of " << s_distance << std::endl;

                
                // evaluate what cars position means for the lane attractiveness 
                  const int MAX_DISTANCE_AHEAD = 100;
                  const int MIN_DISTANCE_AHEAD = 0;

                  //for left lane
                  if (other_car_lane == 0 && s_distance < MAX_DISTANCE_AHEAD && s_distance > MIN_DISTANCE_AHEAD)
                  {
                   car_left = true;
                   //std::cout << "Left lane current max speed: " << left_maxspeed << " other_car_speed: " << other_car_speed << std::endl;
                   if (other_car_speed < left_maxspeed)
                   {
                    left_maxspeed = other_car_speed;
                   }
                   //left_maxspeed = std::min(other_car_speed, left_maxspeed); //TODO: Investigate why this throws an error
                   std::cout << "Left lane max speed set to " << left_maxspeed << std::endl;
                  } 
                
                  //for mid lane
                  if (other_car_lane == 1 && s_distance < MAX_DISTANCE_AHEAD && s_distance > MIN_DISTANCE_AHEAD)
                  {
                   car_mid = true;
                   if (other_car_speed < mid_maxspeed)
                   {
                    mid_maxspeed = other_car_speed;
                   }
                   std::cout << "Mid lane max speed set to " << mid_maxspeed << std::endl;
                  }

                  //for right lane
                  if (other_car_lane == 2 && s_distance < MAX_DISTANCE_AHEAD && s_distance > MIN_DISTANCE_AHEAD)
                  {
                   car_right = true;
                   if (other_car_speed < right_maxspeed)
                   {
                    right_maxspeed = other_car_speed;
                   }
                   std::cout << "Right lane max speed set to " << right_maxspeed << std::endl;
                  } 

                //check if car is in our path and in same lane
                if (s_distance > 1 && s_distance < 10.0 && abs(d_distance) < 1.0)
                {
                  
                  front_warning = true;
                  change_lane = true;
                }
              }

            //END of checking if there is a car in front routine

      
    

            // Speed adjustment routine - do we need to speed up or slow down?

            if (front_warning || change_lane)
            {
              
              if(reference_velocity > 14.0)
              {
                reference_velocity -= .1;
              }

              if(reference_velocity < 17.0)
              {

                //here goes the lane changing logic

                bool ok_to_go_left = false;
                bool ok_to_go_right = false;


                //check lane to the left - is it OK to drive


                if(lane > 0){ //first of all check if that we're not in the leftmost lane already

                  //now check all the cars in that lane
                  int collisions_left = 0;           
                  for(int i = 0; i < sensor_fusion.size(); i++)
                  {

                      double other_car_s = sensor_fusion[i][5]; 
                      double s_distance = other_car_s - car_s;

                      double other_car_d = sensor_fusion[i][6];
                      double d_distance = other_car_d - (car_d - 4.0);
                    if (s_distance > -20.0 && s_distance < 25.0 && abs(d_distance) < 1.0)
                    {
                      collisions_left += 1;
                    }
                  }
                  if(collisions_left == 0)
                  {
                    lane -= 1;
                    change_lane = false;
                  }
                }

                else 
                if(lane < 2){ //first of all check if that we're not in the rightmost lane already

                  //now check all the cars in that lane
                  int collisions_right = 0;           
                  for(int i = 0; i < sensor_fusion.size(); i++)
                  {

                      double other_car_s = sensor_fusion[i][5]; 
                      double s_distance = other_car_s - car_s;

                      double other_car_d = sensor_fusion[i][6];
                      double d_distance = other_car_d - (car_d + 4.0);
                    if (s_distance > -20.0 && s_distance < 25.0 && abs(d_distance) < .5)
                    {
                      collisions_right += 1;
                    }
                  }
                  if(collisions_right == 0)
                  {
                    lane += 1;
                    change_lane = false;
                  }
                }


                
                // END of lane - chaning logic

              }
            }
            else {

              if(reference_velocity < 21.5 && !change_lane)
              {
                reference_velocity += .1;
              }
            }


            



            //SPLINE: Using previous path and spline library to generate a smooth, continuous path

              //Lenght of previous path
              int previous_path_len = previous_path_x.size();

              //Define some seed / anchor points for the spline
              vector<double> spline_seed_x;
              vector<double> spline_seed_y;

              //Set reference x, y, and yaw

              double ref_x = car_x;
              double ref_y = car_y;
              double ref_yaw = deg2rad(car_yaw);


              if(previous_path_len > 1)
              {
                ref_x = previous_path_x[previous_path_len-1];
                ref_y = previous_path_y[previous_path_len-1];

                double ref_x_prev = previous_path_x[previous_path_len-2];
                double ref_y_prev = previous_path_y[previous_path_len-2];

                ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

                spline_seed_x.push_back(ref_x_prev);
                spline_seed_y.push_back(ref_y_prev);
              }

              //First reference points are the current location of the car
              

              spline_seed_x.push_back(ref_x);
              spline_seed_y.push_back(ref_y);

              //Next generate some waypoints far up ahead using the Frenet coordinates
              //Spacing used here: 25 meters apart

              vector<double> next_waypoint0 = getXY(car_s+35, (lane*4+2), map_waypoints_s, map_waypoints_x, map_waypoints_y);
              vector<double> next_waypoint1 = getXY(car_s+70, (lane*4+2), map_waypoints_s, map_waypoints_x, map_waypoints_y);
              vector<double> next_waypoint2 = getXY(car_s+105, (lane*4+2), map_waypoints_s, map_waypoints_x, map_waypoints_y);

              spline_seed_x.push_back(next_waypoint0[0]);
              spline_seed_x.push_back(next_waypoint1[0]);
              spline_seed_x.push_back(next_waypoint2[0]);

              spline_seed_y.push_back(next_waypoint0[1]);
              spline_seed_y.push_back(next_waypoint1[1]);
              spline_seed_y.push_back(next_waypoint2[1]);

              //Shift the points to car's frame of reference

              for(int i = 0; i < spline_seed_x.size(); i++)
              {
                double shift_x = spline_seed_x[i] - ref_x;
                double shift_y = spline_seed_y[i] - ref_y;

                spline_seed_x[i] = (shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw));
                spline_seed_y[i] = (shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw));

              }


              //initiate the spline
              tk::spline s;

              //Use the seed points to generate the spline
              s.set_points(spline_seed_x,spline_seed_y);

            //END OF SPLINE


          	json msgJson;

            //Define the vectors that will hold the path x and y values
          	vector<double> next_x_vals;
          	vector<double> next_y_vals;


            //First reuse any values left from previous path
            for(int i = 0; i < previous_path_len; i++)
            {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds


            // Figuring out how to space points from the spline
            double anchor_point_x = 40.0;
            double anchor_point_y = s(anchor_point_x);

            double distance_to_anchor  = distance(0,0,anchor_point_x,anchor_point_y);

            double x_add_on = 0;


            for(int i = 0; i <= (40 - previous_path_len); i++)
            {
              double N = distance_to_anchor / (.02 * reference_velocity);

              double x_pt = x_add_on + anchor_point_x / N;
              double y_pt = s(x_pt);

              x_add_on = x_pt;

              double x_ref = x_pt;
              double y_ref = y_pt;

              //change reference to global again

              x_pt = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
              y_pt = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw)); 

              x_pt += ref_x;
              y_pt += ref_y;

              next_x_vals.push_back(x_pt);
              next_y_vals.push_back(y_pt);

            }
            
            //END
            msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
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
