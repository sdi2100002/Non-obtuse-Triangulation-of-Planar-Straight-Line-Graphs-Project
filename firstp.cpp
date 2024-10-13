#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <boost/json/src.hpp>
#include <QApplication>  // Include QApplication
#include "graphics/graphics.h"

using namespace boost::json;

int main(int argc, char *argv[]) { // Accept command-line arguments
    // Create the QApplication instance
    QApplication app(argc, argv); // This must be the first GUI-related object created

    // Open the JSON file
    std::ifstream file("../input.json");
    if(!file.is_open()){
        std::cerr << "Could not open the file!" << std::endl;
        return 1;
    }

    // Read the entire file into a string
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string input = buffer.str();

    // Parse the JSON
    value parsed_json = parse(input);

    // Get the object from the parsed JSON
    object obj = parsed_json.as_object();

    // Retrieve fields from the object
    std::string instance_uid = obj["instance_uid"].as_string().c_str(); // Use std::string
    int num_points = obj["num_points"].as_int64();

    // Retrieve points_x and points_y lists
    std::vector<int> points_x, points_y; // Use std::vector
    for (const auto& x : obj["points_x"].as_array()) {
        points_x.push_back(x.as_int64());
    }
    for (const auto& y : obj["points_y"].as_array()) {
        points_y.push_back(y.as_int64());
    }

    // Retrieve the region_boundary
    std::vector<int> region_boundary; // Use std::vector
    for (const auto& boundary_point : obj["region_boundary"].as_array()) {
        region_boundary.push_back(boundary_point.as_int64());
    }

    // Retrieve num_constraints
    int num_constraints = obj["num_constraints"].as_int64();

    // Retrieve additional_constraints
    std::vector<std::pair<int, int>> additional_constraints; // Use std::vector and std::pair
    for (const auto& constraint : obj["additional_constraints"].as_array()) {
        auto constraint_pair = constraint.as_array();
        additional_constraints.emplace_back(constraint_pair[0].as_int64(), constraint_pair[1].as_int64());
    }

    // Display the data
    std::cout << "Instance UID: " << instance_uid << std::endl; // Use std::cout
    std::cout << "Number of points: " << num_points << std::endl;

    std::cout << "Points X: ";
    for (const auto& x : points_x) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    std::cout << "Points Y: ";
    for (const auto& y : points_y) {
        std::cout << y << " "; 
    }
    std::cout << std::endl;

    std::cout << "Region Boundary: ";
    for (const auto& boundary_point : region_boundary) {
        std::cout << boundary_point << " ";
    }
    std::cout << std::endl;

    std::cout << "Number of constraints: " << num_constraints << std::endl;

    std::cout << "Additional Constraints: ";
    for (const auto& constraint : additional_constraints) {
        std::cout << "[" << constraint.first << ", " << constraint.second << "] ";
    }
    std::cout << std::endl;

    // Create a vector of pairs for the points
    std::vector<std::pair<double, double>> points;
    for (size_t i = 0; i < points_x.size(); ++i) {
        points.emplace_back(static_cast<double>(points_x[i]), static_cast<double>(points_y[i]));
    }

    // Call the visualization function
    visualizePoints(points);

    return 0; // Ensure the app is executed at the end
}
