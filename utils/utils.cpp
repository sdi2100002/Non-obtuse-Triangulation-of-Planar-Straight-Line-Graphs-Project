#include <fstream>
#include <sstream>
#include <iostream>
#include "utils.h"
#include "../triangulation/triangulation.h"


// This Function used to format the JSON object into the desired output format
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

// This function used to load a JSON file and read its content into a string
bool loadJsonFile(const std::string& filename, std::string& output) {
    std::ifstream file(filename);// Open the specified file

    if (!file.is_open()) { // Checking if the file was opened successfully
        std::cerr << "Could not open the file!" << std::endl;
        return false;
    }

    // Read the entire file into a string
    std::stringstream buffer;
    buffer << file.rdbuf();
    output = buffer.str();

    return true; 
}

// This function is used to parse a JSON string input and return it as an object
object parseJson(const std::string& input) {
    value parsed_json = parse(input);
    return parsed_json.as_object();
}

// This function is used to retrieve various fields from the parsed JSON object
void retrieveFields(const object& obj, std::string& instance_uid, int& num_points, std::vector<int>& points_x, std::vector<int>& points_y,int& num_constraints, std::vector<std::pair<int, int>>& additional_constraints,std::vector<int>& region_boundary) {
    // Retrieve instance uid and num_points as a string
    instance_uid = obj.at("instance_uid").as_string().c_str(); 
    num_points = obj.at("num_points").as_int64();

    // Retrieve points_x and points_y lists
    for (const auto& x : obj.at("points_x").as_array()) { //iterate through x coordinates
        points_x.push_back(x.as_int64());  //Add each x to the vector
    }
    for (const auto& y : obj.at("points_y").as_array()) { //iterate through y coordinates
        points_y.push_back(y.as_int64());//Add each y to the vector
    }

    // Retrieve num_constraints
    num_constraints = obj.at("num_constraints").as_int64(); 
    
    // Retrieve additional_constraints as we did with points 
    for (const auto& constraint : obj.at("additional_constraints").as_array()) { // Use at() here
        auto constraint_pair = constraint.as_array();
        additional_constraints.emplace_back(constraint_pair[0].as_int64(), constraint_pair[1].as_int64());
    }

    // Retrieve region_boundary as we did with points 
    for (const auto& boundary : obj.at("region_boundary").as_array()) { // Use at() here
        region_boundary.push_back(boundary.as_int64());
    }
}


// This function is used to create a vector of points from points_x and points_y arrays
std::vector<std::pair<double, double>> createPointsVector(const std::vector<int>& points_x, const std::vector<int>& points_y) {
    std::vector<std::pair<double, double>> points; //in this vector we store the result points
    for (size_t i = 0; i < points_x.size(); ++i) {
        //Create pairs of (x,y) and add them to point vector
        points.emplace_back(static_cast<double>(points_x[i]), static_cast<double>(points_y[i]));
    }
    return points;//return the vector
}

// This function is used to process triangulation using the specified points and constraints
void CallprocessTriangulation(const std::vector<std::pair<double, double>>& points, const std::vector<std::pair<int, int>>& additional_constraints, const std::vector<int>& region_boundary, const std::string& instance_uid) {
    // Create a CDTProcessor instance with the given points, constraints, boundary, and instance UID
    Triangulation::CDTProcessor cdtProcessor(points, additional_constraints, region_boundary, instance_uid);
    // Call the main process which we are using to perform triangulation
    cdtProcessor.processTriangulation();
}

