#include <fstream>
#include <sstream>
#include <iostream>
#include "utils.h"
#include "../triangulation/triangulation.h"


// Function to format the JSON object into the desired output format
void printJsonFormatted(const object& obj) {
    std::cout << "{\n";
    std::cout << "  \"instance_uid\": " << obj.at("instance_uid").as_string() << ",\n";
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

// Load JSON file
bool loadJsonFile(const std::string& filename, std::string& output) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
        return false;
    }

    // Read the entire file into a string
    std::stringstream buffer;
    buffer << file.rdbuf();
    output = buffer.str();
    return true;
}

// Parse JSON input
object parseJson(const std::string& input) {
    value parsed_json = parse(input);
    return parsed_json.as_object();
}

// Retrieve fields from the parsed JSON object
void retrieveFields(const object& obj, std::string& instance_uid, int& num_points, std::vector<int>& points_x, std::vector<int>& points_y,int& num_constraints, std::vector<std::pair<int, int>>& additional_constraints,std::vector<int>& region_boundary) {
    instance_uid = obj.at("instance_uid").as_string().c_str(); // Use at() here
    num_points = obj.at("num_points").as_int64(); // Use at() here

    // Retrieve points_x and points_y lists
    for (const auto& x : obj.at("points_x").as_array()) { // Use at() here
        points_x.push_back(x.as_int64());
    }
    for (const auto& y : obj.at("points_y").as_array()) { // Use at() here
        points_y.push_back(y.as_int64());
    }

    // Retrieve num_constraints
    num_constraints = obj.at("num_constraints").as_int64(); // Use at() here
    
    // Retrieve additional_constraints
    for (const auto& constraint : obj.at("additional_constraints").as_array()) { // Use at() here
        auto constraint_pair = constraint.as_array();
        additional_constraints.emplace_back(constraint_pair[0].as_int64(), constraint_pair[1].as_int64());
    }

    // Retrieve region_boundary
    for (const auto& boundary : obj.at("region_boundary").as_array()) { // Use at() here
        region_boundary.push_back(boundary.as_int64());
    }
}


// Create a vector of points from points_x and points_y
std::vector<std::pair<double, double>> createPointsVector(const std::vector<int>& points_x, const std::vector<int>& points_y) {
    std::vector<std::pair<double, double>> points;
    for (size_t i = 0; i < points_x.size(); ++i) {
        points.emplace_back(static_cast<double>(points_x[i]), static_cast<double>(points_y[i]));
    }
    return points;
}

// Process triangulation
void CallprocessTriangulation(const std::vector<std::pair<double, double>>& points, const std::vector<std::pair<int, int>>& additional_constraints, const std::vector<int>& region_boundary, const std::string& instance_uid) {
    // Pass region_boundary into the CDTProcessor
    Triangulation::CDTProcessor cdtProcessor(points, additional_constraints, region_boundary, instance_uid);
    cdtProcessor.processTriangulation();
}

