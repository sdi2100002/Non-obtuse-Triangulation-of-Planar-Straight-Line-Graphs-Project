#define BOOST_ALLOW_DEPRECATED_HEADERS
#ifndef UTILS_H
#define UTILS_H

#include <boost/json.hpp>


using namespace boost::json;

// This Function used to format the JSON object into the desired output format
void printJsonFormatted(const object& obj);

// This function used to load a JSON file and read its content into a string
bool loadJsonFile(const std::string& filename, std::string& output);

// This function is used to parse a JSON string input and return it as an object
object parseJson(const std::string& input);

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
    object &parameters,
    bool &delaunay
);

// This function is used to create a vector of points from points_x and points_y arrays
std::vector<std::pair<double, double>> createPointsVector(const std::vector<int>& points_x, const std::vector<int>& points_y);

// This function is used to process triangulation using the specified points and constraints
void CallprocessTriangulation(const std::vector<std::pair<double, double>>& points, const std::vector<std::pair<int, int>>& additional_constraints, const std::vector<int>& region_boundary, const std::string& instance_uid);

#endif 