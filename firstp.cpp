#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <QApplication>
#include <iostream>
#include <boost/json/src.hpp>
#include "triangulation/triangulation.h"  // Use the typedef defined in triangulation.h
#include "utils/utils.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_triangulation_2.h>

using namespace boost::json;


int main(int argc, char *argv[]) {
    // Create the QApplication instance
    // QApplication app(argc, argv);

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

    // Create a CGAL constrained Delaunay triangulation instance
    CDT cdt;

    // Convert points to CGAL points and add to triangulation
    std::vector<Point> cgal_points;
    for (int i = 0; i < num_points; ++i) {
        cgal_points.emplace_back(points_x[i], points_y[i]);
        cdt.insert(cgal_points.back());
    }

    // Insert constraints into the CGAL triangulation
    for (const auto& constraint : additional_constraints) {
        cdt.insert_constraint(cgal_points[constraint.first], cgal_points[constraint.second]);
    }

    // Use CGAL's draw function to visualize the triangulation
    CGAL::draw(cdt);

    // Create points vector from points_x and points_y if CallprocessTriangulation requires it
    auto points = createPointsVector(points_x, points_y);
    CallprocessTriangulation(points, additional_constraints, region_boundary, instance_uid);


    return 0;
}

