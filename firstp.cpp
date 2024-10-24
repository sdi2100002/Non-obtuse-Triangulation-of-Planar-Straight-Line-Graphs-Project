#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <QApplication>
#include <iostream>
#include <boost/json/src.hpp>
#include "graphics/graphics.h"
#include "triangulation/triangulation.h"
#include "utils/utils.h"


using namespace boost::json;


int main(int argc, char *argv[]) {
    // Create the QApplication instance
    QApplication app(argc, argv); 

    // Load the JSON data from the input.json file
    std::string jsonData;
    if (!loadJsonFile("../input.json", jsonData)) {
        return 1;
    }

    // Parse the JSON data into a JSON object
    object obj = parseJson(jsonData);

    // Print the formatted JSON
    printJsonFormatted(obj);

    // Declare variables to hold the retrieved fields from the JSON object
    std::string instance_uid;
    int num_points;
    std::vector<int> points_x, points_y;
    int num_constraints;
    std::vector<std::pair<int, int>> additional_constraints;
    std::vector<int> region_boundary;

    // Retrieve necessary fields from the JSON object into the declared variables
    retrieveFields(obj, instance_uid, num_points, points_x, points_y, num_constraints, additional_constraints, region_boundary);

    // Create a vector of pairs for the points
    auto points = createPointsVector(points_x, points_y);

    // Call the visualization function for points and constraints
    visualizePoints(points, additional_constraints);

    // Process triangulation
    CallprocessTriangulation(points, additional_constraints, region_boundary, instance_uid);

    return 0; 
}
