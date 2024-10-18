#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// CGAL typedefs
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<std::size_t, K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CDT::Point Point;

namespace Triangulation {

    // Κλάση για την τριγωνοποίηση
    class CDTProcessor {
    public:
        CDTProcessor(const std::vector<std::pair<double, double>>& points, 
                     const std::vector<std::pair<int, int>>& constraints,
                     const std::vector<int>& region_boundary);
        
        // Συνάρτηση για την επεξεργασία της τριγωνοποίησης
        void processTriangulation();

    private:
        std::vector<std::pair<double, double>> points_;
        std::vector<std::pair<int, int>> constraints_;
        std::vector<int> region_boundary_;

        // Συνάρτηση για την οπτικοποίηση της τριγωνοποίησης
        void visualizeTriangulation(const CDT& cdt);

        bool isPointOnSegment(const std::pair<double, double>& pt, const std::pair<double, double>& p1, const std::pair<double, double>& p2);

        bool isPointInPolygon(const std::pair<double, double>& pt, const std::vector<std::pair<double, double>>& polygon);
       

        bool isTriangleWithinBoundary(const CDT::Face_handle& face, const std::vector<int>& boundary);

        double squaredDistance(const Point& p1, const Point& p2);

        bool isObtuseTriangle(const Point& p1, const Point& p2, const Point& p3);
        
        int countObtuseTriangles(const CDT& cdt);

        bool isValidFlip(const Point& p1, const Point& p2, const Point& p3, const Point& p4);
        
        bool tryEdgeFlipping(CDT& cdt, CDT::Face_handle face);
        
        void addSteinerPoint(CDT& cdt, const Point& p1, const Point& p2, const Point& p3);

        void addSteinerAtCircumcenter(CDT& cdt, const Point& p1,const Point& p2, const Point& p3);
        
        void addSteinerInConvexPolygon(CDT& cdt, const Point& p1,const Point& p2, const Point& p3);
        
        int simulateSteinerEffect(CDT& cdt, const Point& steiner_point);

        Point getMidpointOfLongestEdge(const Point& p1, const Point& p2, const Point& p3);

        Point calculate_incenter(const Point& a, const Point& b, const Point& c);
        Point calculate_perpendicular_bisector_point(const Point& a, const Point& b, const Point& c);
    
    };
}

#endif // TRIANGULATION_H
