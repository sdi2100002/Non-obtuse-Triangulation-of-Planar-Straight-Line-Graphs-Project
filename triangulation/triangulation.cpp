#include "triangulation.h"
#include "../graphics/graphics.h"  // Assuming the visualizePoints function is in the graphics folder
#include <cmath>
#include <iostream>
#include <CGAL/algorithm.h>

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

        // bool hasObtuse;
        // int iteration = 0;
        // int maxIterations = 1000;  // Μέγιστος αριθμός επαναλήψεων
        // do {
        //     hasObtuse = false;
        //     for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        //         auto p1 = fit->vertex(0)->point();
        //         auto p2 = fit->vertex(1)->point();
        //         auto p3 = fit->vertex(2)->point();

        //         if (isObtuseTriangle(p1, p2, p3)) {
        //             hasObtuse = true;
        //             addSteinerPoint(cdt, p1, p2, p3);
        //             break;
        //         }
        //     }
        //     iteration++;
        //     if (iteration >= maxIterations) {
        //         std::cerr << "Reached maximum iterations, stopping to avoid infinite loop." << std::endl;
        //         break;  // Σταματάμε για να αποφύγουμε άπειρο βρόχο
        //     }
        // } while (hasObtuse);  // Επανάληψη μέχρι να μην υπάρχουν αμβλυγώνια τρίγωνα

        bool hasObtuse;
        bool flipped;
        int max_iter=1050000;
        int iterationss=0;
        do {
            hasObtuse = false;
            flipped = tryEdgeFlipping(cdt);  // Πρώτη προσπάθεια με flip των ακμών

            for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                auto p1 = fit->vertex(0)->point();
                auto p2 = fit->vertex(1)->point();
                auto p3 = fit->vertex(2)->point();

                if (isObtuseTriangle(p1, p2, p3)) {
                    hasObtuse = true;

                    if (!flipped) {
                        // Αν δεν έγινε περιστροφή, προσθέτουμε σημείο Steiner
                        addSteinerPoint(cdt, p1, p2, p3);
                    }
                    break;
                }
            }
            iterationss++;
            if(iterationss>=max_iter){
                std::cout<<"max iter"<<std::endl;
                break;
            }
        } while (hasObtuse);  // Επανάληψη μέχρι να μην υπάρχουν αμβλυγώνια τρίγωνα

        // Οπτικοποίηση τριγωνοποίησης
        visualizeTriangulation(cdt);
    }

    bool CDTProcessor::tryEdgeFlipping(CDT& cdt) {
        bool flipped = false;
        for (auto eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit) {
            auto face = eit->first;
            int index = eit->second;

            // Παίρνουμε το αντίθετο τρίγωνο για την ακμή
            CDT::Face_handle opposite_face = face->neighbor(index);

            // Ελέγχουμε αν η ακμή είναι περιορισμένη (constrained)
            if (cdt.is_constrained(*eit)) {
                continue; // Αγνοούμε τις περιορισμένες ακμές
            }

            // Ελέγχουμε αν η ακμή είναι εσωτερική (αν υπάρχει γειτονικό τρίγωνο)
            if (cdt.is_infinite(face) || cdt.is_infinite(opposite_face)) continue;

            // Παίρνουμε τις κορυφές των δύο τριγώνων
            auto p1 = face->vertex(cdt.ccw(index))->point();
            auto p2 = face->vertex(cdt.cw(index))->point();
            auto p3 = face->vertex(index)->point();
            auto p4 = opposite_face->vertex(opposite_face->index(face))->point();

            // Ελέγχουμε αν τα τρίγωνα είναι έγκυρα και αν μπορεί να γίνει flip
            if (face != opposite_face && isObtuseTriangle(p1, p2, p3)) {
                // Περιστροφή (flip) της ακμής
                try {
                    cdt.flip(face, index);
                    flipped = true;
                } catch (const CGAL::Assertion_exception& e) {
                    std::cerr << "Edge flip failed: " << e.what() << std::endl;
                }
            }
        }
        return flipped;
    }

    // Συνάρτηση για να προσθέσεις σημείο Steiner στη μέση της μεγαλύτερης ακμής
    void CDTProcessor::addSteinerPoint(CDT& cdt, const Point& p1, const Point& p2, const Point& p3) {
        // Βρίσκουμε τη μεγαλύτερη ακμή
        double a2 = squaredDistance(p2, p3);
        double b2 = squaredDistance(p1, p3);
        double c2 = squaredDistance(p1, p2);

        Point midpoint;
        if (a2 >= b2 && a2 >= c2) {
            // Μεσαίο σημείο της ακμής p2-p3
            midpoint = CGAL::midpoint(p2, p3);
        } else if (b2 >= a2 && b2 >= c2) {
            // Μεσαίο σημείο της ακμής p1-p3
            midpoint = CGAL::midpoint(p1, p3);
        } else {
            // Μεσαίο σημείο της ακμής p1-p2
            midpoint = CGAL::midpoint(p1, p2);
        }

        // Προσθέτουμε το σημείο Steiner
        cdt.insert(midpoint);
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
