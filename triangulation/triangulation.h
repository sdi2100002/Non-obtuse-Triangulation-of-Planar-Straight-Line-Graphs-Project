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


struct AntSolution {
    CDT::Face_handle face;      // Το τρίγωνο στο οποίο έγινε η βελτίωση
    Point steiner_point;        // Το επιλεγμένο Steiner σημείο
    double improvement_metric;  // Μετρική βελτίωσης (π.χ., μείωση σε αριθμό οξέων τριγώνων)
};

namespace Triangulation {

    // We are using this class to handle the CDT processing
    class CDTProcessor {
    public:
        //Constructor to initialize a CDTProcessor Object given the points, constraints , region boundary and a unique instance ID (the input file name) 
        CDTProcessor(const std::vector<std::pair<double, double>>& points, const std::vector<std::pair<int, int>>& constraints,const std::vector<int>& region_boundary,const std::string& instance_uid,const std::string& method,const std::map<std::string,double>& parameters,const bool delaunay,int numOfConstraints);
        
        //This is the Main function to process the triangulation 
        void processTriangulation();

        void selectMethod(CDT &cdt, const std::string& method,const std::map<std::string, double>& parameters);

        CDT createCDT();

    private:
        //Variables to store the points , constraints , region boundary and instance uid
        std::vector<std::pair<double, double>> points_;
        std::vector<std::pair<int, int>> constraints_;
        std::vector<int> region_boundary_; 
        std::string instance_uid_; 
        std::string method_;
        std::map<std::string,double> parameters_;
        int num_constraints_;
        bool delaunay_;

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

        //This function processes all convex polygons in the CDT 
        void processConvexPolygon(CDT& cdt);

        //This function calculates the centroid of a set of points
        Point calculateCentroid(const std::vector<Point>& points);

        //This function determines whether a given face in the CDT contains any constrained edges
        bool hasConstrainedEdge(CDT::Face_handle face, const CDT& cdt);

        //This function checks if a given vertex is part of any constrained edge in the CDT
        bool isVertexOnConstrainedEdge(CDT::Vertex_handle vertex, const CDT& cdt);

        //This function ensures that a given set of points forms a convex polygon by computing the convex hull
        std::vector<Point> ensureConvexPolygon(const std::vector<Point>& points);

        //This function is used to compute the mean point of all vertices in a given face
        Point getMeanAdjacentPoint(CDT::Face_handle face, CDT& cdt);

        //This function implements the local Search method
        void localSearch(CDT& cdt, double L);

        //This function implements the simulated Annealing method
        void simulatedAnnealing(CDT& cdt,double alpha,double beta,int L);

        //This function inserts a Steiner point into the CDT at the specified coordinates.
        void insertSteinerPoint(CDT& cdt, const std::pair<double, double>& point);
        
        //This function generates a random number
        double getRandomUniform();

        //This function is used to calculate the energy
        double calculateEnergy(const CDT& cdt, double alpha, double beta,int numberOfSteinerPoints);

        //This function implements the ant Colony method
        void antColonyOptimization(CDT& cdt, double alpha, double beta, double xi, double psi, int lambda, int num_ants, int num_cycles);

        //This function calculates a heuristic value for a triangle formed by three points, based on the chosen strategy
        double calculateHeuristic(const Point& p1,const Point& p2,const Point&p3,const std::string& strategy);
        
        //This function selects a Steiner point from a list of options using their weights and the total weight
        Point selectSteinerPoint(const std::vector<std::pair<Point, double>>& options, double totalWeight);

        //This function generates a list of candidate Steiner points for insertion based on the geometry of a given face
        std::vector<Point> generateSteinerOptions(const CDT::Face_handle& face);
        
        //This function resolves the conflicts between solutions generated by different ants
        std::vector<AntSolution> resolve_conflicts(const std::vector<AntSolution>& solutions);
        
        //This function applies a single ants solution to the cdt
        void applyAntSolution(const AntSolution& solution, CDT& cdt);

        //This function is used to update the pheromones
        void UpdatePheromones(std::map<CDT::Face_handle, std::map<Point, double>>& pheromone,const std::vector<AntSolution>& solutions, double lambda, double alpha, double beta);
        
        //This function is used to copy the triangulation from a source CDT to a destination CDT
        void copyTriangulation(const CDT& source, CDT& destination);

        //This function selects a strategy based on  probabilities of each  method
        std::string selectStrategy(const std::vector<std::pair<std::string, double>>& method_probabilities, double total_probability);
        
        //This function converts a decimal value to its fract as a string.
        std::string toFraction(double value);
        
        //This function is used to write the final output to JSON file
        void writeOutputToFile(const std::string& filepath,const std::vector<std::string>& steiner_points_x,const std::vector<std::string>& steiner_points_y,const std::vector<std::pair<int, int>>& edges,int obtuse_count,bool randomizationUsed);
        
        //This function returns true if the boundary is non convex and not axis aligned otherwise returns false
        bool isIrregularNonConvex();

        //This function returns true if all edges are axis aligned otherwise returns false
        bool isAxisAlignedNonConvex();
        
        //This function returns true if all edges are present in the graph create a closed polygon otherwise returns false
        bool hasClosedPolygonConstraints();

        //This function returns true if there is at least one open constraint otherwise false
        bool hasOpenConstraints();

        //This function returns true if the boundary is convex otherwise returns false
        bool isConvexBoundary();

        //This function returns an int depends on the category the cdt belongs to , 1 for category A , 2 for categor B etc...
        int detectCategory(CDT &cdt);

        //This function generates a random double value within the specified range
        double randomDouble(double min, double max);
    };
}

#endif // TRIANGULATION_H
