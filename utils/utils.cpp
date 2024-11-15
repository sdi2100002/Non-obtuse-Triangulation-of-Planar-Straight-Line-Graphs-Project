#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "utils.h"
#include "../triangulation/triangulation.h"


// This Function used to format the JSON object into the desired output format
void printJsonFormatted(const boost::json::object &obj) {
    std::cout << "{\n";

    // Print instance_uid
    std::cout << "  \"instance_uid\": \"" << obj.at("instance_uid").as_string() << "\",\n";

    // Print num_points
    std::cout << "  \"num_points\": " << obj.at("num_points").as_int64() << ",\n";

    // Print points_x
    std::cout << "  \"points_x\": [";
    const auto &points_x = obj.at("points_x").as_array();
    for (size_t i = 0; i < points_x.size(); ++i) {
        std::cout << points_x[i].as_int64();
        if (i < points_x.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "],\n";

    // Print points_y
    std::cout << "  \"points_y\": [";
    const auto &points_y = obj.at("points_y").as_array();
    for (size_t i = 0; i < points_y.size(); ++i) {
        std::cout << points_y[i].as_int64();
        if (i < points_y.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "],\n";

    // Print region_boundary
    std::cout << "  \"region_boundary\": [";
    const auto &region_boundary = obj.at("region_boundary").as_array();
    for (size_t i = 0; i < region_boundary.size(); ++i) {
        std::cout << region_boundary[i].as_int64();
        if (i < region_boundary.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "],\n";

    // Print num_constraints
    std::cout << "  \"num_constraints\": " << obj.at("num_constraints").as_int64() << ",\n";

    // Print additional_constraints
    std::cout << "  \"additional_constraints\": [\n";
    const auto &additional_constraints = obj.at("additional_constraints").as_array();
    for (size_t i = 0; i < additional_constraints.size(); ++i) {
        const auto &constraint = additional_constraints[i].as_array();
        std::cout << "    [" << constraint[0].as_int64() << ", " << constraint[1].as_int64() << "]";
        if (i < additional_constraints.size() - 1) {
            std::cout << ",";
        }
        std::cout << "\n";
    }
    std::cout << "  ],\n";

    // Print method
    std::cout << "  \"method\": " << obj.at("method").as_string() << ",\n";

    // Print parameters
    std::cout << "  \"parameters\": {\n";
    const auto &parameters = obj.at("parameters").as_object();
    for (auto it = parameters.begin(); it != parameters.end(); ++it) {
        std::cout << "    \"" << it->key() << "\": ";
        if (it->value().is_double()) {
            std::cout << std::fixed << std::setprecision(1) << it->value().as_double(); // Print double with 1 decimal place
        } else if (it->value().is_int64()) {
            std::cout << it->value().as_int64();
        } else if (it->value().is_string()) {
            std::cout << "\"" << it->value().as_string() << "\"";
        }
        if (std::next(it) != parameters.end()) {
            std::cout << ",";
    }
    std::cout << "\n";
}
    std::cout << "  },\n";

    // Print delaunay
    std::cout << "  \"delaunay\": " << (obj.at("delaunay").as_bool() ? "true" : "false") << "\n";

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
void retrieveFields(
    const object &obj,
    std::string &instance_uid,
    int &num_points,
    std::vector<int> &points_x,
    std::vector<int> &points_y,
    int &num_constraints,
    std::vector<std::pair<int, int>> &additional_constraints,
    std::vector<int> &region_boundary,
    std::string &method,
    std::map<std::string,double> &parameters,
    bool &delaunay
) {
    // Ανακτά το instance_uid
    instance_uid = obj.at("instance_uid").as_string().c_str();

    // Ανακτά το num_points
    num_points = obj.at("num_points").as_int64();

    // Ανακτά τις συντεταγμένες x
    for (const auto &x : obj.at("points_x").as_array()) {
        points_x.push_back(x.as_int64());
    }

    // Ανακτά τις συντεταγμένες y
    for (const auto &y : obj.at("points_y").as_array()) {
        points_y.push_back(y.as_int64());
    }

    // Ανακτά τα όρια της περιοχής
    for (const auto &rb : obj.at("region_boundary").as_array()) {
        region_boundary.push_back(rb.as_int64());
    }

    // Ανακτά τον αριθμό περιορισμών
    num_constraints = obj.at("num_constraints").as_int64();

    // Ανακτά τους πρόσθετους περιορισμούς
    for (const auto &constraint : obj.at("additional_constraints").as_array()) {
        auto pair = constraint.as_array();
        additional_constraints.emplace_back(pair[0].as_int64(), pair[1].as_int64());
    }

    // Ανακτά τη μέθοδο
    method = obj.at("method").as_string().c_str();

    // Ανακτά τις παραμέτρους
    const auto& params_obj = obj.at("parameters").as_object();
    for (const auto& param : params_obj) {
        // Use std::string(param.key()) to convert key to std::string
        std::string key = std::string(param.key());

        // Ensure the value is numeric (double or int)
        if (param.value().is_double()) {
            parameters[key] = param.value().as_double();
        } else if (param.value().is_int64()) {
            parameters[key] = static_cast<double>(param.value().as_int64());
        }
    }
    

    // Ανακτά την τιμή του delaunay
    delaunay = obj.at("delaunay").as_bool();
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
void CallprocessTriangulation(const std::vector<std::pair<double, double>>& points, const std::vector<std::pair<int, int>>& additional_constraints, const std::vector<int>& region_boundary, const std::string& instance_uid,const std::string& method,const std::map<std::string, double>& parameters , bool delaunay) {
    // Create a CDTProcessor instance with the given points, constraints, boundary, and instance UID
    Triangulation::CDTProcessor cdtProcessor(points, additional_constraints, region_boundary, instance_uid, method, parameters);
    // Call the main process which we are using to perform triangulation
    if(delaunay){
        cdtProcessor.processTriangulation();
    }
    else{
        cdtProcessor.selectMethod(method,parameters);
    }
}

