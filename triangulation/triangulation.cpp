#include "triangulation.h"
#include "../graphics/graphics.h"  // Assuming the visualizePoints function is in the graphics folder
#include <cmath>
#include <iostream>
#include <CGAL/algorithm.h>

namespace Triangulation {

    // Constructor
    CDTProcessor::CDTProcessor(const std::vector<std::pair<double, double>>& points, 
                           const std::vector<std::pair<int, int>>& constraints, 
                           const std::vector<int>& region_boundary) 
        : points_(points), constraints_(constraints), region_boundary_(region_boundary) {
    
    // Δημιουργία additional boundary constraints
        for (size_t i = 0; i < region_boundary_.size(); ++i) {
            int start = region_boundary_[i];
            int end = region_boundary_[(i + 1) % region_boundary_.size()]; // Wrap around to the first point
            constraints_.emplace_back(start, end);
        }
    }

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

    // Συνάρτηση που μετράει τα αμβλυγώνια τρίγωνα
    int CDTProcessor::countObtuseTriangles(const CDT& cdt) {
        int obtuse_count = 0;

        // Ελέγχουμε κάθε τρίγωνο αν είναι αμβλυγώνιο
        for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
            auto p1 = fit->vertex(0)->point();
            auto p2 = fit->vertex(1)->point();
            auto p3 = fit->vertex(2)->point();

            if (isObtuseTriangle(p1, p2, p3)) {
                ++obtuse_count;
            }
        }

        return obtuse_count;
    }


    void CDTProcessor::processTriangulation() {
        CDT cdt;

        // Εισαγωγή σημείων και ακμών στην τριγωνοποίηση
        std::vector<CDT::Vertex_handle> vertices;
        for (size_t i = 0; i < points_.size(); ++i) {
            vertices.push_back(cdt.insert(Point(points_[i].first, points_[i].second)));
            vertices.back()->info() = i;
        }

        // Εισαγωγή constraints που περιλαμβάνουν και τα boundaries από τον constructor
        for (const auto& constraint : constraints_) {
            cdt.insert_constraint(vertices[constraint.first], vertices[constraint.second]);
        }

        int obtuse_before = countObtuseTriangles(cdt);
        std::cout << "Αμβλυγώνια τρίγωνα πριν την επεξεργασία: " << obtuse_before << std::endl;

        int max_iter = 10000;
        int iterations = 0;
        int strategy = 0;  // 0: μέσο μεγαλύτερης ακμής, 1: περικέντρο, 2: εσωτερικό κυρτού πολυγώνου
        bool hasObtuse = true;

        while (hasObtuse) {
            hasObtuse = false;

            // Έλεγχος για αμβλυγώνια τρίγωνα
            for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                // Αγνοούμε τρίγωνα που βρίσκονται εκτός των ορίων
                if (!isTriangleWithinBoundary(fit, region_boundary_)) {
                    //std::cout << "Triangle is outside boundary." << std::endl;
                    continue;  // Skip this triangle if it is outside the boundary
                }
                std::cout << "Triangle is inside the boundary" << std::endl;

                auto p1 = fit->vertex(0)->point();
                auto p2 = fit->vertex(1)->point();
                auto p3 = fit->vertex(2)->point();

                // Έλεγχος αν το τρίγωνο είναι αμβλυγώνιο
                if (isObtuseTriangle(p1, p2, p3)) {
                    hasObtuse = true;  // Found an obtuse triangle
                    bool flipped = tryEdgeFlipping(cdt, fit);

                    if (!flipped) {
                        Point steiner_point;

                        // Προσομοίωση της προσθήκης Steiner ανάλογα με τη στρατηγική
                        if (strategy == 0) {
                            steiner_point = getMidpointOfLongestEdge(p1, p2, p3);
                        } else if (strategy == 1) {
                            steiner_point = CGAL::circumcenter(p1, p2, p3);  // Περικέντρο
                        } else if (strategy == 2) {
                            steiner_point = CGAL::centroid(p1, p2, p3);  // Εσωτερικό κυρτού πολυγώνου
                        }

                        // Υπολογισμός προσωρινής επίδρασης της προσθήκης Steiner
                        int obtuse_after_sim = simulateSteinerEffect(cdt, steiner_point);

                        if (obtuse_after_sim < obtuse_before) {// Αν βελτιώνει την κατάσταση, προσθέτουμε το σημείο πραγματικά
                            cdt.insert(steiner_point);
                            obtuse_before = obtuse_after_sim;
                            //std::cout << "Προσθήκη Steiner point βελτίωσε την κατάσταση." << std::endl;
                        }
                    }
                    break;  // Exit the loop after processing the first obtuse triangle found
                }
            }

            iterations++;
            if (iterations == max_iter) {
                strategy++;
                iterations = 0;
                if (strategy > 2) {
                    std::cout << "Η διαδικασία ολοκληρώθηκε!" << std::endl;
                    break;
                }
            }
        }

        obtuse_before = countObtuseTriangles(cdt);
        std::cout << "Αμβλυγώνια τρίγωνα μετά την επεξεργασία: " << obtuse_before << std::endl;

        visualizeTriangulation(cdt);
    }

    bool CDTProcessor::isPointInPolygon(const std::pair<double, double>& pt, const std::vector<std::pair<double, double>>& polygon) {
        int n = polygon.size();
        bool inside = false;

        for (int i = 0, j = n - 1; i < n; j = i++) {
            if (((polygon[i].second > pt.second) != (polygon[j].second > pt.second)) &&
                (pt.first < (polygon[j].first - polygon[i].first) * (pt.second - polygon[i].second) / 
                (polygon[j].second - polygon[i].second) + polygon[i].first)) {
                inside = !inside;
            }
        }
        return inside;
    }

    bool CDTProcessor::isTriangleWithinBoundary(const CDT::Face_handle& face, const std::vector<int>& boundary) {
        // Prepare the polygon points from region_boundary
        std::vector<std::pair<double, double>> polygon;
        for (int index : boundary) {
            polygon.emplace_back(points_[index]);
        }

        // Retrieve triangle vertices from the face
        auto p1 = std::make_pair(face->vertex(0)->point().x(), face->vertex(0)->point().y());
        auto p2 = std::make_pair(face->vertex(1)->point().x(), face->vertex(1)->point().y());
        auto p3 = std::make_pair(face->vertex(2)->point().x(), face->vertex(2)->point().y());

        // Check if all triangle vertices are inside the polygon
        return isPointInPolygon(p1, polygon) && isPointInPolygon(p2, polygon) && isPointInPolygon(p3, polygon);
    }

    // Συνάρτηση για υπολογισμό της επίδρασης του Steiner point χωρίς να προστεθεί
    int CDTProcessor::simulateSteinerEffect(CDT& cdt, const Point& steiner_point) {
        // Προσωρινή αντιγραφή του CDT για να ελέγξουμε την επίδραση
        CDT temp_cdt = cdt;

        // Προσθήκη του Steiner point στον προσωρινό πίνακα
        temp_cdt.insert(steiner_point);

        // Υπολογισμός του αριθμού των αμβλυγώνιων τριγώνων μετά την προσθήκη
        return countObtuseTriangles(temp_cdt);
    }

// Συνάρτηση για να πάρουμε το μέσο της μεγαλύτερης ακμής
    Point CDTProcessor::getMidpointOfLongestEdge(const Point& p1, const Point& p2, const Point& p3) {
        double a2 = squaredDistance(p2, p3);
        double b2 = squaredDistance(p1, p3);
        double c2 = squaredDistance(p1, p2);

        if (a2 >= b2 && a2 >= c2) {
            return CGAL::midpoint(p2, p3);  // Μεσο της μεγαλύτερης ακμής p2-p3
        } else if (b2 >= a2 && b2 >= c2) {
            return CGAL::midpoint(p1, p3);  // Μεσο της ακμής p1-p3
        } else {
            return CGAL::midpoint(p1, p2);  // Μεσο της ακμής p1-p2
        }
    }

    bool CDTProcessor::isValidFlip(const Point& p1, const Point& p2, const Point& p3, const Point& p4) {
        // Ελέγχουμε αν τα νέα τρίγωνα είναι obtuse
        if (isObtuseTriangle(p1, p2, p4)) {
            return false; // Το πρώτο τρίγωνο είναι obtuse
        }
        if (isObtuseTriangle(p2, p3, p4)) {
            return false; // Το δεύτερο τρίγωνο είναι obtuse
        }
        return true; // Το flip είναι έγκυρο
    }

    bool CDTProcessor::tryEdgeFlipping(CDT& cdt, CDT::Face_handle face) {
        bool flipped = false;

        // Ελέγχουμε όλες τις ακμές του τριγώνου
        for (int index = 0; index < 3; ++index) {
            CDT::Face_handle opposite_face = face->neighbor(index);

            // Αγνοούμε τις περιορισμένες ακμές
            if (cdt.is_constrained(CDT::Edge(face, index))) {
                continue;
            }

            // Ελέγχουμε αν υπάρχει γειτονικό τρίγωνο και αν η ακμή είναι εσωτερική
            if (cdt.is_infinite(face) || cdt.is_infinite(opposite_face)) continue;

            // Παίρνουμε τις κορυφές του τριγώνου και του γειτονικού
            auto p1 = face->vertex(cdt.ccw(index))->point();
            auto p2 = face->vertex(cdt.cw(index))->point();
            auto p3 = face->vertex(index)->point();

            // Αν το τρίγωνο είναι αμβλυγώνιο, ελέγχουμε για flip
            if (isObtuseTriangle(p1, p2, p3)) {
                // Χρησιμοποιούμε τη συνάρτηση is_flippable της CDT για να ελέγξουμε αν το flip είναι δυνατό
                bool isFlipable=cdt.is_flipable(face, index);
                if (isFlipable) { // Η σωστή κλήση με 2 ορίσματα
                    try {
                        cdt.flip(face, index);
                        flipped = true;
                        std::cout << "Flip επιτυχές!" << std::endl;
                        break;
                    } catch (const CGAL::Assertion_exception& e) {
                        std::cerr << "Flip απέτυχε: " << e.what() << std::endl;
                    }
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

    // Συνάρτηση για να προσθέσουμε σημείο Steiner στο περικεντρο του τριγώνου
    void CDTProcessor::addSteinerAtCircumcenter(CDT& cdt, const Point& p1,const Point& p2, const Point& p3) {
        // Υπολογισμός του περικέντρου
        Point circumcenter = CGAL::circumcenter(p1, p2, p3);
        cdt.insert(circumcenter);
    }

    // Συνάρτηση για να προσθέσουμε σημείο Steiner στο εσωτερικό κυρτού πολυγώνου
    void CDTProcessor::addSteinerInConvexPolygon(CDT& cdt, const Point& p1,const Point& p2, const Point& p3) {
        // Προσδιορισμός του κέντρου βάρους ως απλή λύση για το εσωτερικό του κυρτού πολυγώνου
        Point centroid = CGAL::centroid(p1, p2, p3);
        cdt.insert(centroid);
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
