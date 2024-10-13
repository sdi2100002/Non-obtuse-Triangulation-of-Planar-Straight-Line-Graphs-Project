#include "triangulation.h"
#include "../graphics/graphics.h"  // Assuming the visualizePoints function is in the graphics folder
#include <cmath>
#include <iostream>

namespace Triangulation {

    // Constructor
    CDTProcessor::CDTProcessor(const std::vector<std::pair<double, double>>& points, 
                               const std::vector<std::pair<int, int>>& constraints) 
        : points_(points), constraints_(constraints) {}

    // Function to calculate the squared length of a segment
    double CDTProcessor::squaredDistance(const Point& p1, const Point& p2) {
        return CGAL::square(p1.x() - p2.x()) + CGAL::square(p1.y() - p2.y());
    }

    // Function to check if a triangle is obtuse
    bool CDTProcessor::isObtuseTriangle(const Point& p1, const Point& p2, const Point& p3) {
        double a2 = squaredDistance(p2, p3);
        double b2 = squaredDistance(p1, p3);
        double c2 = squaredDistance(p1, p2);

        return (a2 + b2 < c2) || (a2 + c2 < b2) || (b2 + c2 < a2);
    }

    // Συνάρτηση για την επεξεργασία της τριγωνοποίησης
    void CDTProcessor::processTriangulation() {
        // Δημιουργία CDT αντικειμένου
        CDT cdt;

        // Εισαγωγή σημείων στην τριγωνοποίηση με το index ως πληροφορία
        std::vector<CDT::Vertex_handle> vertices;
        for (size_t i = 0; i < points_.size(); ++i) {
            vertices.push_back(cdt.insert(Point(points_[i].first, points_[i].second)));
            vertices.back()->info() = i;
        }

        // Προσθήκη περιορισμών
        for (const auto& constraint : constraints_) {
            cdt.insert_constraint(vertices[constraint.first], vertices[constraint.second]);
        }
        // Check for obtuse triangles
        bool hasObtuse = false;
        for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
            auto p1 = fit->vertex(0)->point();
            auto p2 = fit->vertex(1)->point();
            auto p3 = fit->vertex(2)->point();

            if (isObtuseTriangle(p1, p2, p3)) {
                hasObtuse = true;
                break;
            }
        }

        if (hasObtuse) {
            std::cout << "The triangulation contains at least one obtuse triangle." << std::endl;
        } else {
            std::cout << "All triangles in the triangulation are acute or right-angled." << std::endl;
        }

        // Οπτικοποίηση τριγωνοποίησης
        visualizeTriangulation(cdt);
    }

    // Συνάρτηση για την οπτικοποίηση της τριγωνοποίησης
    void CDTProcessor::visualizeTriangulation(const CDT& cdt) {
        std::vector<std::pair<double, double>> triangulated_points;
        for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
            auto p1 = fit->vertex(0)->point();
            auto p2 = fit->vertex(1)->point();
            auto p3 = fit->vertex(2)->point();
            triangulated_points.push_back({p1.x(), p1.y()});
            triangulated_points.push_back({p2.x(), p2.y()});
            triangulated_points.push_back({p3.x(), p3.y()});
        }
        // Εμφάνιση της τριγωνοποίησης (μέσω της υλοποίησης visualizePoints)
        visualizePoints(triangulated_points, constraints_);
    }

}
