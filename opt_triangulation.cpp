#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <QApplication>
#include <iostream>
#include <boost/json/src.hpp>
#include "triangulation/triangulation.h" 
#include "utils/utils.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_triangulation_2.h>

using namespace boost::json;


int main(int argc, char *argv[]) {

    std::cout<<"argc "<<argc<<std::endl;
    if(argc<6){
        std::cerr << "Usage: " << argv[0] << "-i /path/to/input.json -o /path/to/output.json pre_selected_params " << std::endl;
        return 1;
    }

    std::string inputPath,outputPath;
    std::map<std::string,double> parameters = {
        {"alpha",2.0},
        {"beta",5.0},
        {"xi",1.0},
        {"psi",3.0},
        {"lambda",0.5},
        {"kappa", 50},
        {"L", 100}
    };

    for (int i = 1; i < argc; i++) {
        std::string flag(argv[i]);
        if (flag == "-i") {
            inputPath = argv[i + 1];
        } else if (flag == "-o") {
            outputPath = argv[i + 1];
        }else if (parameters.find(flag)!=parameters.end()){
            if (i+1<argc){
                parameters[flag]=std::stod(argv[++i]);
            }else{
                std::cerr << "Error : Missing value for parameter " << flag << std::endl;
                return 1;
            }
        }
    
    }

    // Load the JSON data from the input.json file
    std::string jsonData;
    if (!loadJsonFile(inputPath, jsonData)) {
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
    std::string method;
    bool delaunay;


    // Retrieve necessary fields from the JSON object into the declared variables
    retrieveFields(obj, instance_uid, num_points, points_x, points_y, num_constraints, additional_constraints, region_boundary,method,parameters,delaunay);

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
    CallprocessTriangulation(points, additional_constraints, region_boundary, instance_uid, method, parameters, delaunay,num_constraints);


    return 0;
}

