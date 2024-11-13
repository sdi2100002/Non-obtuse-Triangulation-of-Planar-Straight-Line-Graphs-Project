#ifndef TRIANGULATION_H
#define TRIANGULATION_H
#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <boost/json.hpp>

namespace json = boost::json;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<std::size_t, K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CDT::Point Point;

namespace Triangulation {

    // We are using this class to handle the CDT processing
    class CDTProcessor {
    public:
        //Constructor to initialize a CDTProcessor Object given the points, constraints , region boundary and a unique instance ID (the input file name) 
        CDTProcessor(const std::vector<std::pair<double, double>>& points, const std::vector<std::pair<int, int>>& constraints,const std::vector<int>& region_boundary,const std::string& instance_uid);
        
        //This is the Main function to process the triangulation 
        void processTriangulation();

    private:
        //Variables to store the points , constraints , region boundary and instance uid
        std::vector<std::pair<double, double>> points_;
        std::vector<std::pair<int, int>> constraints_;
        std::vector<int> region_boundary_; 
        std::string instance_uid_; 

        // This function is used to visulize the given triangulation
        void visualizeTriangulation(const CDT& cdt);

        // This function checks if the point pt is on the line segment formed by points p1 and p2
        bool isPointOnSegment(const std::pair<double, double>& pt, const std::pair<double, double>& p1, const std::pair<double, double>& p2);

        //This function checks if a given point lies inside the region boundary 
        bool isPointInsideBoundary(const std::pair<double, double>& point,const std::vector<int>& region_boundary,const std::vector<std::pair<double, double>>& points) const ;

        // This function checks if a point is inside a polygon
        bool isPointInPolygon(const std::pair<double, double>& pt, const std::vector<std::pair<double, double>>& polygon);

        //This function is checking if a triangle is within the given boundary defined by region_boundary
        bool isTriangleWithinBoundary(const CDT::Face_handle& face, const std::vector<int>& boundary);
        
        // This function calculates the squared distance between two points
        double squaredDistance(const Point& p1, const Point& p2);

        // This function checks if a triangle is obtuse
        bool isObtuseTriangle(const Point& p1, const Point& p2, const Point& p3);
        
        //This function counts the obtuse triangles in the triangulation
        int countObtuseTriangles(const CDT& cdt);
        
        // This function attempts to flip an edge returns true if a flip done otherwise false (RETURNS FALSE ALWAYS!!)
        bool tryEdgeFlipping(CDT& cdt, CDT::Face_handle face);
        
        // This function is to simulate the effect of a Steiner point without adding it to the triangulation
        int simulateSteinerEffect(CDT& cdt, const Point& steiner_point);

        // This function is used to get the midpoint of the longest edge
        Point getMidpointOfLongestEdge(const Point& p1, const Point& p2, const Point& p3);

        // This function calculates the incenter of a triangle
        Point calculate_incenter(const Point& a, const Point& b, const Point& c);
        
        //This function finds a point that lies on the perpendicular bisector of the line segment connecting two points a and b
        Point calculate_perpendicular_bisector_point(const Point& a, const Point& b, const Point& c);
        
        // This function computes the orthogonal projection of point p onto the triangle
        Point projectPointOntoTriangle(const Point& p, const Point& p1, const Point& p2, const Point& p3);

        // This function creates a JSON object for output
        json::object createOutputJson(const std::string& instance_uid,const std::vector<std::pair<double, double>>& steiner_points,const std::vector<std::pair<int, int>>& edges);
        
        // This function prints the output JSON to the console
        void printOutputJson(const json::object& output_json);
    };
}

#endif // TRIANGULATION_H
