#include <iostream>
#include <fstream>
#include "json/single_include/nlohmann/json.hpp" // Include the JSON library
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/IO/io.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <cmath>
#include <vector>
#include <cassert>


// #include <CGAL/draw_polygon_2.h>  // Include CGAL's drawing functions
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Polygon_2.h>
// #include <CGAL/Constrained_Delaunay_triangulation_2.h>

using json = nlohmann::json;

struct CDT {
    std::string instance_uid;
    int num_points;
    std::vector<int> points_x;
    std::vector<int> points_y;
    std::vector<int> region_boundary;
    int num_constraints;
    std::vector<std::vector<int>> additional_constraints;
};

// Function to load the JSON data into a CDT struct
void load_json_to_cdt(const std::string& file_name, CDT& cdt) {
    std::ifstream input_file(file_name);
    
    if (!input_file.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
        return;
    }
    
    // Parse the JSON file
    json j;
    input_file >> j;
    
    // Extract data into CDT struct
    cdt.instance_uid = j["instance_uid"];
    cdt.num_points = j["num_points"];
    cdt.points_x = j["points_x"].get<std::vector<int>>();
    cdt.points_y = j["points_y"].get<std::vector<int>>();
    cdt.region_boundary = j["region_boundary"].get<std::vector<int>>();
    cdt.num_constraints = j["num_constraints"];
    cdt.additional_constraints = j["additional_constraints"].get<std::vector<std::vector<int>>>();
}

int main() {
    CDT cdt;

    // Load the JSON file into the CDT struct
    load_json_to_cdt("input.json", cdt);

    // Output some of the loaded data to check
    std::cout << "Instance UID: " << cdt.instance_uid << std::endl;
    std::cout << "Number of Points: " << cdt.num_points << std::endl;
    std::cout << "First X Point: " << cdt.points_x[0] << std::endl;
    std::cout << "First Y Point: " << cdt.points_y[0] << std::endl;

    CGAL::draw(cdt);

    return 0;
}






// Define the kernel and the types used in the triangulation
// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// typedef CGAL::Polygon_2<K> Polygon;
// typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;

// int main() {
//     CDT cdt;

//     // Load the JSON file into the CDT struct (previous function)
//     load_json_to_cdt("input.json", cdt);

//     // Here, you would typically add vertices and constraints to the CDT
//     for (size_t i = 0; i < cdt.num_points; ++i) {
//         cdt.insert(Point_2(cdt.points_x[i], cdt.points_y[i]));
//     }

//     // Add constraints (edges) based on your region_boundary and additional_constraints
//     for (const auto& constraint : cdt.additional_constraints) {
//         cdt.insert_constraint(cdt.points_x[constraint[0]], cdt.points_y[constraint[0]],
//                               cdt.points_x[constraint[1]], cdt.points_y[constraint[1]]);
//     }

//     // Draw the CDT
//     CGAL::draw(cdt);  // Call CGAL's draw function

//     return 0;
// }
