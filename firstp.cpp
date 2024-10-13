#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <boost/json/src.hpp>
#include <QApplication>
#include "graphics/graphics.h"

using namespace boost::json;

// Function to format the JSON object into the desired output format
void printJsonFormatted(const object& obj) {
    std::cout << "{\n";
    std::cout << "  \"instance_uid\": \"" << obj.at("instance_uid").as_string() << "\",\n";
    std::cout << "  \"num_points\": " << obj.at("num_points").as_int64() << ",\n";
    
    // Print points_x
    std::cout << "  \"points_x\": [";
    const auto& points_x = obj.at("points_x").as_array();
    for (size_t i = 0; i < points_x.size(); ++i) {
        std::cout << points_x[i].as_int64();
        if (i < points_x.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "],\n";
    
    // Print points_y
    std::cout << "  \"points_y\": [";
    const auto& points_y = obj.at("points_y").as_array();
    for (size_t i = 0; i < points_y.size(); ++i) {
        std::cout << points_y[i].as_int64();
        if (i < points_y.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "],\n";

    // Print region_boundary
    std::cout << "  \"region_boundary\": [";
    const auto& region_boundary = obj.at("region_boundary").as_array();
    for (size_t i = 0; i < region_boundary.size(); ++i) {
        std::cout << region_boundary[i].as_int64();
        if (i < region_boundary.size() - 1) {
            std::cout << ",";
        }
    }
    std::cout << "],\n";

    // Print num_constraints
    std::cout << "  \"num_constraints\": " << obj.at("num_constraints").as_int64() << ",\n";

    // Print additional_constraints
    std::cout << "  \"additional_constraints\": [\n";
    const auto& additional_constraints = obj.at("additional_constraints").as_array();
    for (size_t i = 0; i < additional_constraints.size(); ++i) {
        const auto& constraint = additional_constraints[i].as_array();
        std::cout << "    [" << constraint[0].as_int64() << ", " << constraint[1].as_int64() << "]";
        if (i < additional_constraints.size() - 1) {
            std::cout << ",";
        }
        std::cout << "\n";
    }
    std::cout << "  ]\n";
    std::cout << "}\n";
}

int main(int argc, char *argv[]) {
    // Create the QApplication instance
    QApplication app(argc, argv); // This must be the first GUI-related object created

    // Open the JSON file
    std::ifstream file("../input.json");
    if (!file.is_open()) {
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

    // Print the formatted JSON
    printJsonFormatted(obj);

    // Retrieve fields from the object (as before)
    std::string instance_uid = obj["instance_uid"].as_string().c_str();
    int num_points = obj["num_points"].as_int64();

    // Retrieve points_x and points_y lists
    std::vector<int> points_x, points_y;
    for (const auto& x : obj["points_x"].as_array()) {
        points_x.push_back(x.as_int64());
    }
    for (const auto& y : obj["points_y"].as_array()) {
        points_y.push_back(y.as_int64());
    }

    // Retrieve num_constraints
    int num_constraints = obj["num_constraints"].as_int64();
    
    // Retrieve additional_constraints
    std::vector<std::pair<int, int>> additional_constraints;
    for (const auto& constraint : obj["additional_constraints"].as_array()) {
        auto constraint_pair = constraint.as_array();
        additional_constraints.emplace_back(constraint_pair[0].as_int64(), constraint_pair[1].as_int64());
    }

    // Create a vector of pairs for the points
    std::vector<std::pair<double, double>> points;
    for (size_t i = 0; i < points_x.size(); ++i) {
        points.emplace_back(static_cast<double>(points_x[i]), static_cast<double>(points_y[i]));
    }

    // Call the visualization function for points and constraints
    visualizePoints(points, additional_constraints);

    return 0; // Ensure the app is executed at the end
}
