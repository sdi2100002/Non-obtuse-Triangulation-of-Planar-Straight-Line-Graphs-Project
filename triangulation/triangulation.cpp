#include "triangulation.h"
#include "../graphics/graphics.h"  // Assuming the visualizePoints function is in the graphics folder
#include <cmath>
#include <iostream>
#include <CGAL/algorithm.h>
#include <utility>

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

        int max_iter = 100;
        int iterations = 0;
        bool hasObtuse = true;

        while (hasObtuse && iterations < max_iter) {
            hasObtuse = false;
            int best_obtuse_after_sim = obtuse_before; // Track best obtuse triangle count in this iteration
            Point best_steiner_point;

            // Έλεγχος για αμβλυγώνια τρίγωνα
            for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                //std::cout<<"for loop"<<std::endl;
                // Αγνοούμε τρίγωνα που βρίσκονται εκτός των ορίων
                if (!isTriangleWithinBoundary(fit, region_boundary_)) {
                    //std::cout << "Triangle is outside boundary." << std::endl;
                    continue;  // Skip this triangle if it is outside the boundary
                }

                auto p1 = fit->vertex(0)->point();
                auto p2 = fit->vertex(1)->point();
                auto p3 = fit->vertex(2)->point();
                // Σημείο που θέλεις να προβάλεις
                CGAL::Point_2<CGAL::Epick> p; // Ορισμός ή απόκτηση του σημείου p

                if (!CGAL::is_finite(p1.x()) || !CGAL::is_finite(p1.y()) || !CGAL::is_finite(p2.x()) || !CGAL::is_finite(p2.y()) || !CGAL::is_finite(p3.x()) || !CGAL::is_finite(p3.y())) {
                    std::cerr << "Invalid point encountered in triangle." << std::endl;
                    continue;  // Skip the current triangle
                }

                // Έλεγχος αν το τρίγωνο είναι αμβλυγώνιο
                if (isObtuseTriangle(p1, p2, p3)) {
                    hasObtuse = true;  // Found an obtuse triangle
                    bool flipped = tryEdgeFlipping(cdt, fit);

                    if (!flipped) {
                        // Δοκιμή όλων των στρατηγικών για το τρέχον τρίγωνο
                        for (int strategy = 0; strategy <= 5; ++strategy) {
                            Point steiner_point;

                            if (strategy == 0) {
                                steiner_point = CGAL::circumcenter(p1, p2, p3);  // Περικέντρο
                            } else if (strategy == 1) {
                                steiner_point = calculate_incenter(p1, p2, p3);  // Εσωτερικό τριγώνου
                            } else if (strategy == 2) {
                                steiner_point = getMidpointOfLongestEdge(p1, p2, p3);
                            } else if (strategy == 3) {
                                steiner_point = CGAL::centroid(p1, p2, p3);  // Εσωτερικό κυρτού πολυγώνου
                            } else if (strategy == 4) {
                                steiner_point = calculate_perpendicular_bisector_point(p1, p2, p3);  // Κάθετος διχοτόμος
                            }
                            else if (strategy == 5){
                                steiner_point = projectPointOntoTriangle(p,p1,p2,p3);
                                //std::cout << "Steiner Point: (" << steiner_point.x() << ", " << steiner_point.y() << ")\n";
                            }

                            // Validate steiner_point
                            if (!CGAL::is_finite(steiner_point.x()) || !CGAL::is_finite(steiner_point.y())) {
                                //std::cerr << "Invalid Steiner point calculated: " << steiner_point << std::endl;
                                continue;  // Skip this strategy if the Steiner point is invalid
                            }

                            // Υπολογισμός προσωρινής επίδρασης της προσθήκης Steiner point
                            int obtuse_after_sim = simulateSteinerEffect(cdt, steiner_point);
                    
                            // Αν αυτή η στρατηγική βελτιώνει περισσότερο, κρατάμε αυτή τη στρατηγική
                            if (obtuse_after_sim < best_obtuse_after_sim && isPointInsideBoundary(std::make_pair(steiner_point.x(),steiner_point.y()),region_boundary_,points_)) {
                                //std::cout << "Steiner Point in bounds: (" << steiner_point.x() << ", " << steiner_point.y() << ")\n";
                                best_obtuse_after_sim = obtuse_after_sim;
                                best_steiner_point = steiner_point;
                            }
                            else if (!isPointInsideBoundary(std::make_pair(steiner_point.x(),steiner_point.y()),region_boundary_,points_)){
                                //std::cout << "Steiner Point out of bounds: (" << steiner_point.x() << ", " << steiner_point.y() << ")\n";
                            }                        
                        }
                    }
                }
            }

            // Εισαγωγή του καλύτερου Steiner point
            if (best_obtuse_after_sim < obtuse_before) {
                cdt.insert(best_steiner_point);
                obtuse_before = best_obtuse_after_sim;
                std::cout << "Προσθήκη Steiner point βελτίωσε την κατάσταση. Αμβλυγώνια τρίγωνα: " << obtuse_before << std::endl;
            }
            iterations++;
        }

        obtuse_before = countObtuseTriangles(cdt);
        std::cout << "Αμβλυγώνια τρίγωνα μετά την επεξεργασία: " << obtuse_before << std::endl;

        visualizeTriangulation(cdt);
    }

    bool CDTProcessor::isPointInsideBoundary(const std::pair<double, double>& point,const std::vector<int>& region_boundary,const std::vector<std::pair<double, double>>& points) const {
        int n = region_boundary.size();
        bool inside = false;

        // Helper function to check if a point lies on a segment (p1, p2)
        auto is_on_segment = [](const std::pair<double, double>& p,
                                const std::pair<double, double>& p1,
                                const std::pair<double, double>& p2) -> bool {
            if (std::min(p1.first, p2.first) <= p.first && p.first <= std::max(p1.first, p2.first) &&
                std::min(p1.second, p2.second) <= p.second && p.second <= std::max(p1.second, p2.second)) {
                return (p2.first - p1.first) * (p.second - p1.second) == (p2.second - p1.second) * (p.first - p1.first);
            }
            return false;
        };

        // Ray-casting algorithm: Check how many times a ray from the point intersects the polygon edges
        for (int i = 0, j = n - 1; i < n; j = i++) {
            const auto& p1 = points[region_boundary[i]];  // Current vertex
            const auto& p2 = points[region_boundary[j]];  // Previous vertex

            // Check if the point is on the current edge
            if (is_on_segment(point, p1, p2)) {
                return true;  // If the point is on the boundary, return true
            }

            // Check if the ray crosses the edge (p1, p2)
            if (((p1.second > point.second) != (p2.second > point.second)) &&
                (point.first < (p2.first - p1.first) * (point.second - p1.second) / (p2.second - p1.second) + p1.first)) {
                inside = !inside;
            }
        }

        return inside;
    }

    // Συνάρτηση που ελέγχει αν το σημείο pt είναι πάνω στη γραμμή που ενώνεται από τα σημεία p1 και p2
    bool CDTProcessor::isPointOnSegment(const std::pair<double, double>& pt, const std::pair<double, double>& p1, const std::pair<double, double>& p2) {
        double crossProduct = (pt.second - p1.second) * (p2.first - p1.first) - (pt.first - p1.first) * (p2.second - p1.second);
        if (std::abs(crossProduct) > 1e-10) {
            return false; // Το σημείο δεν είναι στη γραμμή
        }

        // Ελέγχουμε αν το σημείο pt βρίσκεται εντός των ορίων της γραμμής
        if (pt.first < std::min(p1.first, p2.first) || pt.first > std::max(p1.first, p2.first) ||
            pt.second < std::min(p1.second, p2.second) || pt.second > std::max(p1.second, p2.second)) {
            return false; // Το σημείο είναι εκτός των ορίων
        }

        return true; // Το σημείο είναι πάνω στη γραμμή
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
            // Έλεγχος αν το σημείο βρίσκεται πάνω σε μια πλευρά του πολυγώνου
            if (isPointOnSegment(pt, polygon[i], polygon[j])) {
                //std::cout << "Το σημείο (" << pt.first << ", " << pt.second << ") είναι πάνω στο περίγραμμα του πολύγωνου." << std::endl;
                return true; // Επιστρέφουμε true εφόσον το σημείο είναι στο περίγραμμα
            }
        }
        // Αν το σημείο δεν είναι μέσα στο πολύγωνο, εκτύπωσε τις συντεταγμένες του
        if (!inside) {
            std::cout << "Το σημείο (" << pt.first << ", " << pt.second << ") δεν είναι μέσα στο πολύγωνο." << std::endl;
        }
        // std::cout << "Το σημείο (" << pt.first << ", " << pt.second << ")  στο πολύγωνο." << std::endl;
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

    Point CDTProcessor::calculate_incenter(const Point& a, const Point& b, const Point& c) {
        // Υπολογισμός των αποστάσεων A, B, C
        double A = CGAL::squared_distance(b, c);
        double B = CGAL::squared_distance(a, c);
        double C = CGAL::squared_distance(a, b);

        // Υπολογισμός των συντεταγμένων του incenter
        double Ix = (A * a.x() + B * b.x() + C * c.x()) / (A + B + C);
        double Iy = (A * a.y() + B * b.y() + C * c.y()) / (A + B + C);

        return Point(Ix, Iy);
    }

    Point CDTProcessor::calculate_perpendicular_bisector_point(const Point& a, const Point& b, const Point& c) {
        // Υπολογισμός του μέσου της πλευράς a-b
        double mx = (a.x() + b.x()) / 2;
        double my = (a.y() + b.y()) / 2;

        // Υπολογισμός της κλίσης της πλευράς a-b
        double dx = b.x() - a.x();
        double dy = b.y() - a.y();

        // Έλεγχος αν η πλευρά a-b είναι κάθετη (dx == 0)
        if (dx == 0) {
            // Αν η πλευρά είναι κάθετη, η κάθετος διχοτόμος θα είναι οριζόντια
            return Point(mx, my + 1);  // Επιστρέφουμε ένα σημείο πάνω στην οριζόντια κάθετο διχοτόμο
        }

        // Κλίση της πλευράς a-b
        double slope_ab = dy / dx;

        // Κλίση της κάθετης διχοτόμου
        double slope_perpendicular = -1 / slope_ab;

        // Μετατόπιση κατά μια αυθαίρετη απόσταση (π.χ. 1 μονάδα) για να βρούμε σημείο πάνω στη διχοτόμο
        double dx_perpendicular = 1 / sqrt(1 + slope_perpendicular * slope_perpendicular);
        double dy_perpendicular = slope_perpendicular * dx_perpendicular;

        // Το νέο σημείο πάνω στην κάθετο διχοτόμο (μετατοπισμένο από το μέσο της πλευράς a-b)
        return Point(mx + dx_perpendicular, my + dy_perpendicular);
    }   

    // Συνάρτηση για τον υπολογισμό της κάθετης προβολής του σημείου p στο τρίγωνο
    Point CDTProcessor::projectPointOntoTriangle(const Point& p, const Point& p1, const Point& p2, const Point& p3) {
        // Υπολογισμός παραμέτρων για τις πλευρές του τριγώνου
        auto projectOntoLine = [](const Point& p, const Point& a, const Point& b) {
            double ab_x = b.x() - a.x();
            double ab_y = b.y() - a.y();
            double ap_x = p.x() - a.x();
            double ap_y = p.y() - a.y();

            double ab_squared = ab_x * ab_x + ab_y * ab_y;
            if (ab_squared == 0) return a; // a και b είναι το ίδιο σημείο

            // Υπολογισμός του κλάσματος
            double t = (ap_x * ab_x + ap_y * ab_y) / ab_squared;
            t = std::max(0.0, std::min(1.0, t)); // περιορισμός του t στο [0, 1]

            return Point(a.x() + t * ab_x, a.y() + t * ab_y);
        };

        // Υπολογισμός της κάθετης προβολής για κάθε πλευρά του τριγώνου
        Point proj1 = projectOntoLine(p, p1, p2);
        Point proj2 = projectOntoLine(p, p2, p3);
        Point proj3 = projectOntoLine(p, p3, p1);

        double dist1 = std::sqrt((p.x() - proj1.x()) * (p.x() - proj1.x()) + (p.y() - proj1.y()) * (p.y() - proj1.y()));
        double dist2 = std::sqrt((p.x() - proj2.x()) * (p.x() - proj2.x()) + (p.y() - proj2.y()) * (p.y() - proj2.y()));
        double dist3 = std::sqrt((p.x() - proj3.x()) * (p.x() - proj3.x()) + (p.y() - proj3.y()) * (p.y() - proj3.y()));

        if (dist1 <= dist2 && dist1 <= dist3) {
            return proj1;
        } else if (dist2 <= dist1 && dist2 <= dist3) {
            return proj2;
        } else {
            return proj3;
        }
    }
}
