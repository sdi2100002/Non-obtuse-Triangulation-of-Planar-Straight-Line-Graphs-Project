#define BOOST_ALLOW_DEPRECATED_HEADERS
#include "triangulation.h"
#include <cmath>
#include <iostream>
#include <CGAL/algorithm.h>
#include <CGAL/draw_triangulation_2.h>
#include <utility>
#include <numeric> 
#include <set>
#include <random>
#include <gmp.h> 
#include <CGAL/Exact_rational.h>
#include "../CGAL-5.6.1/include/CGAL/convex_hull_2.h"

namespace Triangulation {

    // Constructor 
    CDTProcessor::CDTProcessor(const std::vector<std::pair<double, double>>& points, const std::vector<std::pair<int, int>>& constraints, const std::vector<int>& region_boundary,const std::string& instance_uid,const std::string& method,const std::map<std::string,double>& parameters,const bool delaunay) 
        : points_(points), constraints_(constraints), region_boundary_(region_boundary) , instance_uid_(instance_uid) , method_(method) , parameters_(parameters),delaunay_(delaunay) {
    
        //Adding the boundary constraints (the region boundary is treated as a contraint)
        for (size_t i = 0; i < region_boundary_.size(); ++i) {
            int start = region_boundary_[i];
            int end = region_boundary_[(i + 1) % region_boundary_.size()];
            constraints_.emplace_back(start, end);
        }
    }

    // This function calculates the squared distance between two points
    double CDTProcessor::squaredDistance(const Point& p1, const Point& p2) {
        return CGAL::square(p1.x() - p2.x()) + CGAL::square(p1.y() - p2.y());
    }

    // This function checks if a triangle is obtuse
    bool CDTProcessor::isObtuseTriangle(const Point& p1, const Point& p2, const Point& p3) {
        //Calculates the squared lenghts of the triangles sides
        double a2 = squaredDistance(p2, p3);
        double b2 = squaredDistance(p1, p3);
        double c2 = squaredDistance(p1, p2);

        //Checking the obtuse condition and returns true if the triangle is obtuse false otherwise
        return (a2 + b2 < c2) || (a2 + c2 < b2) || (b2 + c2 < a2);
    }

    //This function counts the obtuse triangles in the triangulation
    int CDTProcessor::countObtuseTriangles(const CDT& cdt) {
        //Init the counter
        int obtuse_count = 0;

        //Loop through all finite faces (triangles) in the cdt
        for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
            //Get the vertices of the current triangle
            auto p1 = fit->vertex(0)->point();
            auto p2 = fit->vertex(1)->point();
            auto p3 = fit->vertex(2)->point();

            //Check if the triangle is obtuse and if it is do the counter++
            if (isObtuseTriangle(p1, p2, p3)) {
                obtuse_count++;
            }
        }
        //return the counter
        return obtuse_count;
    }

    //This function calculates the centroid of a set of points
    Point CDTProcessor::calculateCentroid(const std::vector<Point>& points) {
        double x_sum = 0.0, y_sum = 0.0;
        for (const auto& point : points) {
            x_sum += point.x();
            y_sum += point.y();
        }
        return Point(x_sum / points.size(), y_sum / points.size());
    }

    //This function determines whether a given face in the CDT contains any constrained edges
    bool CDTProcessor::hasConstrainedEdge(CDT::Face_handle face, const CDT& cdt) {
        for (int i = 0; i < 3; ++i) {
            if (cdt.is_constrained(std::make_pair(face, i))) {
                return true;
            }
        }
        return false;
    }

    //This function checks if a given vertex is part of any constrained edge in the CDT
    bool CDTProcessor::isVertexOnConstrainedEdge(CDT::Vertex_handle vertex, const CDT& cdt) {
        CDT::Face_circulator fc_start = cdt.incident_faces(vertex), fc = fc_start;
        do {
            for (int i = 0; i < 3; ++i) {
                if (cdt.is_constrained(std::make_pair(fc, i)) &&
                    (fc->vertex(i) == vertex || fc->vertex((i + 1) % 3) == vertex)) {
                    return true;
                }
            }
        } while (++fc != fc_start);
        return false;
    }

    //This function ensures that a given set of points forms a convex polygon by computing the convex hull
    std::vector<Point> CDTProcessor::ensureConvexPolygon(const std::vector<Point>& points) {
        std::vector<Point> convex_hull_points;
        CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(convex_hull_points));
        return convex_hull_points;
    }

    void CDTProcessor::processConvexPolygon(CDT& cdt) {
        std::vector<CDT::Face_handle> obtuse_faces;

        //Find all obtuse triangles without constrained edges and save them for processing
        for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
            if (isObtuseTriangle(fit->vertex(0)->point(), fit->vertex(1)->point(), fit->vertex(2)->point()) 
                && !hasConstrainedEdge(fit, cdt)) {
                obtuse_faces.push_back(fit);
            }
        }

        //Group obtuse triangles into convex polygons
        for (const auto& face : obtuse_faces) {
            std::vector<Point> polygon_points;
            std::vector<CDT::Vertex_handle> non_constrained_vertices;

            //Find vertices that are not part of constrained edges
            for (int i = 0; i < 3; ++i) {
                auto vertex = face->vertex(i);
                if (!isVertexOnConstrainedEdge(vertex, cdt)) {
                    non_constrained_vertices.push_back(vertex);
                    polygon_points.push_back(vertex->point());
                }
            }

            //Make sure that the collected points form a convex polygon
            polygon_points = ensureConvexPolygon(polygon_points);
            if (polygon_points.size() < 3) continue; //Skip if its not

            //Calculate the centroid of the polygon and use it as a Steiner point
            Point centroid = calculateCentroid(polygon_points);

            //Remove vertices that are not part of constrained edges from the CDT
            for (const auto& vertex : non_constrained_vertices) {
                cdt.remove(vertex);
            }

            //Insert constraint edges around the convex polygon to retain its shape
            std::vector<CDT::Vertex_handle> polygon_vertices;
            for (size_t i = 0; i < polygon_points.size(); ++i) {
                auto v1 = cdt.insert(polygon_points[i]);
                polygon_vertices.push_back(v1);
                auto v2 = cdt.insert(polygon_points[(i + 1) % polygon_points.size()]);
                cdt.insert_constraint(v1, v2);
            }

            //Insert the centroid as a new vertex and connect it with constraint edges to the polygon's vertices
            auto centroid_vertex = cdt.insert(centroid);
            for (const auto& v : polygon_vertices) {
                cdt.insert_constraint(centroid_vertex, v);
            }
        }
    }

    //This is the Main function to process the triangulation 
    void CDTProcessor::processTriangulation() {
        CDT cdt;
        std::vector<std::pair<double, double>> steiner_points;
        std::vector<std::pair<int, int>> edges; 


        // Insert the points into the triangulation
        std::vector<CDT::Vertex_handle> vertices;
        for (size_t i = 0; i < points_.size(); ++i) {
            vertices.push_back(cdt.insert(Point(points_[i].first, points_[i].second)));
            vertices.back()->info() = i; //Vertex index for indentification
        }

        //Insert the constraints (and the boundary constraints)
        for (const auto& constraint : constraints_) {
            cdt.insert_constraint(vertices[constraint.first], vertices[constraint.second]);
        }

        //Count and print the number of obtuse triangles before the process of the triangulation begins
        int obtuse_before = countObtuseTriangles(cdt);
        std::cout << "Αμβλυγώνια τρίγωνα πριν την επεξεργασία: " << obtuse_before << std::endl;

        int max_iter = 100; //var to store max number of loops
        int iterations = 0; //var to store loops done
        bool hasObtuse = true;  //there is at least one more obtuse triangle 

        processConvexPolygon(cdt);

        //Main loop to perform the process that reduces the obtuses triangles 
        //While there is at least one obtuse triangle and we haven't reach the number of max iterations
        while (hasObtuse && iterations < max_iter) {
            hasObtuse = false;
            int best_obtuse_after_sim = obtuse_before;
            Point best_steiner_point;

            //For each triangle check if there is an obtuse 
            for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                
                // Skip the triangles outside the region boundary
                if (!isTriangleWithinBoundary(fit, region_boundary_)) {
                    //std::cout << "Triangle is outside boundary." << std::endl;
                    continue;
                }

                auto p1 = fit->vertex(0)->point();
                auto p2 = fit->vertex(1)->point();
                auto p3 = fit->vertex(2)->point();
                
                CGAL::Point_2<CGAL::Epick> p; //point for project strategy
                
                //Skip triangles with vertices values infinite or NaN coordinates
                if (!CGAL::is_finite(p1.x()) || !CGAL::is_finite(p1.y()) || !CGAL::is_finite(p2.x()) || !CGAL::is_finite(p2.y()) || !CGAL::is_finite(p3.x()) || !CGAL::is_finite(p3.y())) {
                    std::cerr << "Invalid point encountered in triangle." << std::endl;
                    continue;  // Skip the current triangle
                }

                // Checking if the current triangle is obtuse
                if (isObtuseTriangle(p1, p2, p3)) {
                    hasObtuse = true;  // Found at least one obtuse triangle

                    // Try edge-flipping to resolve the obtuse triangle
                    bool flipped = tryEdgeFlipping(cdt, fit);

                    //If edge-flipping fails, proceed to Steiner point insertion strategies
                    //The function tryEdgeFlipping always returns False!! so we will always get in this if statement
                    if (!flipped) {
                        
                        // Try all the strategies we are using to place a steiner point in this triangle
                        for (int strategy = 0; strategy <= 5; ++strategy) {
                            Point steiner_point;

                            if (strategy == 0) {
                                steiner_point = CGAL::circumcenter(p1, p2, p3);  // Circumcenter
                            } else if (strategy == 1) {
                                steiner_point = calculate_incenter(p1, p2, p3);  // Incenter
                            } else if (strategy == 2) {
                                steiner_point = getMidpointOfLongestEdge(p1, p2, p3); // Midpoint of longest edge
                            } else if (strategy == 3) {
                                steiner_point = CGAL::centroid(p1, p2, p3);  // Centroid
                            } else if (strategy == 4) {
                                steiner_point = calculate_perpendicular_bisector_point(p1, p2, p3);  // Perpendicular bisector point
                            }
                            else if (strategy == 5){
                                steiner_point = projectPointOntoTriangle(p,p1,p2,p3); // Projected Steiner point onto triangle
                                //std::cout << "Steiner Point: (" << steiner_point.x() << ", " << steiner_point.y() << ")\n";
                            }

                            // Check if the coordinates of the steiner point are either inf or NaN
                            if (!CGAL::is_finite(steiner_point.x()) || !CGAL::is_finite(steiner_point.y())) {
                                //std::cerr << "Invalid Steiner point calculated: " << steiner_point << std::endl;
                                continue;  // Skip this strategy if the Steiner point is invalid
                            }

                            //Simulate the effect of inserting the Steiner point,count the obtuse triangles after the insertion of this steiner point
                            int obtuse_after_sim = simulateSteinerEffect(cdt, steiner_point);
                    
                            // In case this strategy reduces the obtuse triangles and the steiner point is within the boundary select it
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

            // Insert the best Steiner point found from the  strategies
            if (best_obtuse_after_sim < obtuse_before) {
                auto new_vertex=cdt.insert(best_steiner_point);
                steiner_points.push_back({best_steiner_point.x(),best_steiner_point.y()});
                obtuse_before = best_obtuse_after_sim;
                std::cout << "Προσθήκη Steiner point βελτίωσε την κατάσταση. Αμβλυγώνια τρίγωνα: " << obtuse_before << std::endl;
                
            }
            processConvexPolygon(cdt);
            iterations++;
        }

        //Select the final edges from the triangulation
        for (auto eit=cdt.finite_edges_begin();eit!=cdt.finite_edges_end(); ++eit){
            //Get the vertex index1 and 2
            auto index1 = eit->first->vertex(eit->second)->info();  
            auto index2 = eit->first->vertex(eit->second == 0 ? 2 : eit->second - 1)->info();

            // Check if the indices are valid and within the bounds
            if (index1 >= 0 && index1 < points_.size() && index2 >= 0 && index2 < points_.size()) {
                edges.push_back({index1, index2}); // add the edge to the result list
            }
            
        }

        //Create the output JSON object 
        json::object output_json = createOutputJson(instance_uid_, steiner_points, edges);

        // Print the output JSON
        printOutputJson(output_json);

        //Count and print the number of the obtuse triangles after the process of triangulation
        obtuse_before = countObtuseTriangles(cdt);
        std::cout << "Αμβλυγώνια τρίγωνα μετά την επεξεργασία: " << obtuse_before << std::endl;

        //Visulize the final triangulation
        CGAL::draw(cdt);

        selectMethod(cdt, method_,parameters_);
    }

    //This function checks if a given point lies inside the region boundary 
    bool CDTProcessor::isPointInsideBoundary(const std::pair<double, double>& point,const std::vector<int>& region_boundary,const std::vector<std::pair<double, double>>& points) const {
        int n = region_boundary.size(); // Number of vertices in the boundary
        bool inside = false; //flag init as false

        // Helper function to check if a point lies on a segment (p1, p2)
        auto is_on_segment = [](const std::pair<double, double>& p,const std::pair<double, double>& p1,const std::pair<double, double>& p2) -> bool {
            if (std::min(p1.first, p2.first) <= p.first && p.first <= std::max(p1.first, p2.first) &&
                std::min(p1.second, p2.second) <= p.second && p.second <= std::max(p1.second, p2.second)) {
                //Checks if the point lies exactly on the lie segment  
                return (p2.first - p1.first) * (p.second - p1.second) == (p2.second - p1.second) * (p.first - p1.first);
            }
            return false;
        };

        // Ray-casting algorithm: Check how many times a ray from the point intersects the polygon edges
        for (int i = 0, j = n - 1; i < n; j = i++) {
            //Current and previous vertex
            const auto& p1 = points[region_boundary[i]];  
            const auto& p2 = points[region_boundary[j]];  

            // Check if the point is on the current edge
            if (is_on_segment(point, p1, p2)) {
                return true;  // If the point is on the boundary, return true
            }

            // Check if the ray crosses the edge (p1, p2)
            if (((p1.second > point.second) != (p2.second > point.second)) &&
                (point.first < (p2.first - p1.first) * (point.second - p1.second) / (p2.second - p1.second) + p1.first)) {
                inside = !inside; //change the value of the flag
            }
        }

        return inside;
    }

    // This function checks if the point pt is on the line segment formed by points p1 and p2
    bool CDTProcessor::isPointOnSegment(const std::pair<double, double>& pt, const std::pair<double, double>& p1, const std::pair<double, double>& p2) {
        //Calculate the cross product in order to determine the collinearity
        double crossProduct = (pt.second - p1.second) * (p2.first - p1.first) - (pt.first - p1.first) * (p2.second - p1.second);
        if (std::abs(crossProduct) > 1e-10) {
            return false; // The point is not on the line 
        }

        // Checking if the point pt is whitin the bounds of the line segment
        if (pt.first < std::min(p1.first, p2.first) || pt.first > std::max(p1.first, p2.first) ||
            pt.second < std::min(p1.second, p2.second) || pt.second > std::max(p1.second, p2.second)) {
            return false; // The point is outside
        }

        return true; // The point is on the line segment
    }

    // This function checks if a point is inside a polygon
    bool CDTProcessor::isPointInPolygon(const std::pair<double, double>& pt, const std::vector<std::pair<double, double>>& polygon) {
        int n = polygon.size();// Number of vertices in the boundary
        bool inside = false;//flag init as false

        // Ray-casting algorithm to determine if the point is inside the polygon
        for (int i = 0, j = n - 1; i < n; j = i++) {
            // Check if the ray crosses the edge between polygon[i] and polygon[j]
            if (((polygon[i].second > pt.second) != (polygon[j].second > pt.second)) &&
                (pt.first < (polygon[j].first - polygon[i].first) * (pt.second - polygon[i].second) / 
                (polygon[j].second - polygon[i].second) + polygon[i].first)) {
                inside = !inside; // Change the value of flag
            }
            // Checking if the point is on any edge of the polygon
            if (isPointOnSegment(pt, polygon[i], polygon[j])) {
                //std::cout << "Το σημείο (" << pt.first << ", " << pt.second << ") είναι πάνω στο περίγραμμα του πολύγωνου." << std::endl;
                return true; // Returns true when the point is on the boundary
            }
        }

        // In case the point is not insde the polygon prints its coordinates
        if (!inside) {
            std::cout << "Το σημείο (" << pt.first << ", " << pt.second << ") δεν είναι μέσα στο πολύγωνο." << std::endl;
        }
        
        return inside; 
    }

    //This function is checking if a triangle is within the given boundary defined by region_boundary
    bool CDTProcessor::isTriangleWithinBoundary(const CDT::Face_handle& face, const std::vector<int>& boundary) {
        // Preparing the polygon points from region_boundary
        std::vector<std::pair<double, double>> polygon;
        for (int index : boundary) {
            polygon.emplace_back(points_[index]);
        }

        //Get the vertices of the current triangle
        auto p1 = std::make_pair(face->vertex(0)->point().x(), face->vertex(0)->point().y());
        auto p2 = std::make_pair(face->vertex(1)->point().x(), face->vertex(1)->point().y());
        auto p3 = std::make_pair(face->vertex(2)->point().x(), face->vertex(2)->point().y());

        // Check if all triangle vertices are inside the polygon
        return isPointInPolygon(p1, polygon) && isPointInPolygon(p2, polygon) && isPointInPolygon(p3, polygon);
    }

    // This function is to simulate the effect of a Steiner point without adding it to the triangulation
    int CDTProcessor::simulateSteinerEffect(CDT& cdt, const Point& steiner_point) {
        //Create a copy of the CDT to check the effect of the steiner point
        CDT temp_cdt = cdt;

        //Insertion of the steiner point to the temp CDT
        temp_cdt.insert(steiner_point);

        //Calcuate the number of obtuse triangles after the insertion and return this number
        return countObtuseTriangles(temp_cdt);
    }

    // This function is used to get the midpoint of the longest edge
    Point CDTProcessor::getMidpointOfLongestEdge(const Point& p1, const Point& p2, const Point& p3) {
        // Calculate the squared distances between the vertices
        double a2 = squaredDistance(p2, p3);
        double b2 = squaredDistance(p1, p3);
        double c2 = squaredDistance(p1, p2);

        // Find the biggest edge and return the midpoint of this edge
        if (a2 >= b2 && a2 >= c2) {
            return CGAL::midpoint(p2, p3);  
        } else if (b2 >= a2 && b2 >= c2) {
            return CGAL::midpoint(p1, p3);  
        } else {
            return CGAL::midpoint(p1, p2); 
        }
    }

    // This function attempts to flip an edge
    bool CDTProcessor::tryEdgeFlipping(CDT& cdt, CDT::Face_handle face) {
        bool flipped = false; //flag init to false

        // Checking all the edges of the triangle
        for (int index = 0; index < 3; ++index) {
            CDT::Face_handle opposite_face = face->neighbor(index); // Opposite triangle

            // Ignore the constrained edges
            if (cdt.is_constrained(CDT::Edge(face, index))) {
                continue;
            }

            // Checking if there is a neighboring triangle and if the edge is internal
            if (cdt.is_infinite(face) || cdt.is_infinite(opposite_face)) continue;

            //Get the vertices of the triangle and the neighboring triangle
            auto p1 = face->vertex(cdt.ccw(index))->point();
            auto p2 = face->vertex(cdt.cw(index))->point();
            auto p3 = face->vertex(index)->point();

            //If the triangle is obtuse check for a flip 
            if (isObtuseTriangle(p1, p2, p3)) {
                // Use the function is_flippable from CDT to check if the flip is possible
                //is flipable will always return false!!
                bool isFlipable=cdt.is_flipable(face, index);
                if (isFlipable) { 
                    try {
                        cdt.flip(face, index); //try to flip
                        flipped = true;
                        std::cout << "Flip επιτυχές!" << std::endl;
                        break;
                    } catch (const CGAL::Assertion_exception& e) {
                        std::cerr << "Flip απέτυχε: " << e.what() << std::endl;
                    }
                }
            }
        }
        return flipped; //Return is a flip occured or not 
    }

    // This function calculates the incenter of a triangle
    Point CDTProcessor::calculate_incenter(const Point& a, const Point& b, const Point& c) {
        //Calculating the distances between the points
        double A = CGAL::squared_distance(b, c);
        double B = CGAL::squared_distance(a, c);
        double C = CGAL::squared_distance(a, b);

        // Calculating the coordinates of the incenter
        double Ix = (A * a.x() + B * b.x() + C * c.x()) / (A + B + C);
        double Iy = (A * a.y() + B * b.y() + C * c.y()) / (A + B + C);

        return Point(Ix, Iy); //Return the incenter 
    }

    //This function finds a point that lies on the perpendicular bisector of the line segment connecting two points a and b
    Point CDTProcessor::calculate_perpendicular_bisector_point(const Point& a, const Point& b, const Point& c) {
        // Calculate the midpoint of edge a-b
        double mx = (a.x() + b.x()) / 2;
        double my = (a.y() + b.y()) / 2;

        // Calculate the slope of edge a-b
        double dx = b.x() - a.x();
        double dy = b.y() - a.y();

        // Check if edge a-b is vertical (dx == 0)
        if (dx == 0) {
            // If the edge is vertical, the perpendicular bisector will be horizontal
            return Point(mx, my + 1);  // Return a point above the horizontal perpendicular bisector
        }

        // Slope of edge a-b
        double slope_ab = dy / dx;

        // Slope of the perpendicular bisector
        double slope_perpendicular = -1 / slope_ab; // Inverse slope

        // Offset by an arbitrary distance  to find a point on the bisector
        double dx_perpendicular = 1 / sqrt(1 + slope_perpendicular * slope_perpendicular);
        double dy_perpendicular = slope_perpendicular * dx_perpendicular;

        // The new point on the perpendicular bisector
        return Point(mx + dx_perpendicular, my + dy_perpendicular);
    }   

    // This function computes the orthogonal projection of point p onto the triangle
    Point CDTProcessor::projectPointOntoTriangle(const Point& p, const Point& p1, const Point& p2, const Point& p3) {
        // Calculate parametres for the sides of the triangle
        auto projectOntoLine = [](const Point& p, const Point& a, const Point& b) {
            double ab_x = b.x() - a.x(); //Difference in x 
            double ab_y = b.y() - a.y(); //Differnce in y
            double ap_x = p.x() - a.x(); //Difference in x from point p to a
            double ap_y = p.y() - a.y();//Difference in y from point p to a

            double ab_squared = ab_x * ab_x + ab_y * ab_y;// Squared length of edge ab
            if (ab_squared == 0) return a; // a and b are the same point

            // Calculate the fraction
            double t = (ap_x * ab_x + ap_y * ab_y) / ab_squared;
            t = std::max(0.0, std::min(1.0, t)); // to [0,1]

            return Point(a.x() + t * ab_x, a.y() + t * ab_y);
        };

        // Calculate the orthogonal projection onto each side of the triangle
        Point proj1 = projectOntoLine(p, p1, p2);
        Point proj2 = projectOntoLine(p, p2, p3);
        Point proj3 = projectOntoLine(p, p3, p1);

        // Calculate distances to each projection point
        double dist1 = std::sqrt((p.x() - proj1.x()) * (p.x() - proj1.x()) + (p.y() - proj1.y()) * (p.y() - proj1.y()));
        double dist2 = std::sqrt((p.x() - proj2.x()) * (p.x() - proj2.x()) + (p.y() - proj2.y()) * (p.y() - proj2.y()));
        double dist3 = std::sqrt((p.x() - proj3.x()) * (p.x() - proj3.x()) + (p.y() - proj3.y()) * (p.y() - proj3.y()));

        //Find the closes projection and return that point
        if (dist1 <= dist2 && dist1 <= dist3) {
            return proj1;
        } else if (dist2 <= dist1 && dist2 <= dist3) {
            return proj2;
        } else {
            return proj3;
        }
    }

    // This function creates a JSON object for output
    json::object CDTProcessor::createOutputJson(const std::string& instance_uid,const std::vector<std::pair<double, double>>& steiner_points,const std::vector<std::pair<int, int>>& edges) {
        json::object output; //create the main JSON object

        // Add those fields
        output["content_type"] = "CG_SHOP_2025_Solution";
        output["instance_uid"] = instance_uid;

        // Arrays to store the coordinates of steiner_points
        json::array steiner_points_x;
        json::array steiner_points_y;

        // Loop through each Steiner point and add its coordinates to the JSON arrays
        for (const auto& point : steiner_points) {
            steiner_points_x.push_back(point.first);  
            steiner_points_y.push_back(point.second);
        }

        // Add the coordinates to the output JSON
        output["steiner_points_x"] = steiner_points_x;
        output["steiner_points_y"] = steiner_points_y;

        // Array to store the edges
        json::array edges_array;

        // Loop through each edge and create a pair of indices
        for (const auto& edge : edges) {
            // Create an array for each edge
            json::array edge_array;
            edge_array.push_back(edge.first);
            edge_array.push_back(edge.second);
            edges_array.push_back(edge_array); //Add the edge to the edges array
        }

        // Add the edges to the output JSON
        output["edges"] = edges_array;

        return output;
    }

    // This function prints the output JSON to the console
    void CDTProcessor::printOutputJson(const json::object& output_json) {
        std::cout << "{\n";  

        // Prints
        std::cout << " \"content_type\": " << output_json.at("content_type").as_string() << ",\n";
        std::cout << " \"instance_uid\": " << output_json.at("instance_uid").as_string() << ",\n";

        //Print x coordinates of steiner points
        std::cout << " \"steiner_points_x\": [";
        const auto& steiner_points_x = output_json.at("steiner_points_x").as_array();
        for (size_t i = 0; i < steiner_points_x.size(); ++i) {
            if (steiner_points_x[i].is_string()) {
                std::cout << "\"" << steiner_points_x[i].as_string() << "\"";
            } else {
                std::cout << steiner_points_x[i].as_double();
            }
            if (i < steiner_points_x.size() - 1) {
                std::cout << ", "; 
            }
        }
        std::cout << "],\n";

        //Print y coordinates of steiner points
        std::cout << " \"steiner_points_y\": [";
        const auto& steiner_points_y = output_json.at("steiner_points_y").as_array();
        for (size_t i = 0; i < steiner_points_y.size(); ++i) {
            if (steiner_points_y[i].is_string()) {
                std::cout << "\"" << steiner_points_y[i].as_string() << "\"";
            } else {
                std::cout << steiner_points_y[i].as_double();
            }
            if (i < steiner_points_y.size() - 1) {
                std::cout << ", ";  
            }
        }
        std::cout << "],\n";

        // Print edges
        std::cout << " \"edges\": [\n";
        const auto& edges = output_json.at("edges").as_array();
        for (size_t i = 0; i < edges.size(); ++i) {
            auto edge_array = edges[i].as_array();
            std::cout << "  [";
            std::cout << edge_array[0].as_int64() << ", " << edge_array[1].as_int64() << "]";
            if (i < edges.size() - 1) {
                std::cout << ",\n";  
            }
        }
        std::cout << "\n ]\n"; 

        std::cout << "}\n";
   
    }

    void CDTProcessor::selectMethod(CDT &cdt, const std::string& method,const std::map<std::string, double>& parameters){
        if (method == "local") {
            std::cout << "Selected method: Local Search\n";
            std::cout << "Parameters:\n";
            for (const auto& [key, value] : parameters) {
                std::cout << "  " << key << ": " << value << "\n";
            }
            // Placeholder for local search logic
            std::cout << "Executing local search...\n";
            auto iter = parameters.begin();
            double value = iter->second;
            localSearch(cdt, value);
        } 
        else if (method == "sa") { // Simulated Annealing
            std::cout << "Selected method: Simulated Annealing\n";
            std::cout << "Parameters:\n";
            for (const auto& [key, value] : parameters) {
                std::cout << "  " << key << ": " << value << "\n";
            }
            // Placeholder for simulated annealing logic
            std::cout << "Executing simulated annealing...\n";

            // Extract parameters for Simulated Annealing
            double alpha = (parameters.find("alpha") != parameters.end()) ? parameters.at("alpha") : 2.0;
            double beta = (parameters.find("beta") != parameters.end()) ? parameters.at("beta") : 0.2;
            int L = (parameters.find("L") != parameters.end()) ? static_cast<int>(parameters.at("L")) : 5000;

        
            std::cout << "Using alpha: " << alpha << ", beta: " << beta << ", L: " << L << "\n";

            simulatedAnnealing(cdt,alpha,beta,L);
        } 
        else if (method == "ant") { // Ant Colony Optimization
            std::cout << "Selected method: Ant Colony Optimization\n";
            std::cout << "Parameters:\n";
            for (const auto& [key, value] : parameters) {
                std::cout << "  " << key << ": " << value << "\n";
            }
            double alpha = (parameters.find("alpha") != parameters.end()) ? parameters.at("alpha") : 0.0;
            double beta = (parameters.find("beta") != parameters.end()) ? parameters.at("beta") : 0.0;
            double xi = (parameters.find("xi") != parameters.end()) ? parameters.at("xi") : 1.0;
            double psi = (parameters.find("psi") != parameters.end()) ? parameters.at("psi") : 1.0;
            double lambda = (parameters.find("lambda") != parameters.end()) ? parameters.at("lambda") : 0.1;
            int kappa = (parameters.find("kappa") != parameters.end()) ? static_cast<int>(parameters.at("kappa")) : 5;
            int L = (parameters.find("L") != parameters.end()) ? static_cast<int>(parameters.at("L")) : 10;

            std::cout<<"alpha: " << alpha << std::endl;
            
            // Placeholder for ant colony optimization logic
            std::cout << "Executing ant colony optimization...\n";

            antColonyOptimization(cdt,alpha,beta,xi,psi,lambda,kappa,L);
        } 
        else {
            // Handle invalid method input
            std::cerr << "Error: Unknown method \"" << method << "\".\n";
            std::cerr << "Available methods: local, sa, ant\n";
        }
    }

    void CDTProcessor::localSearch(CDT& cdt, double L) {
        int counter = 0;
        bool hasObtuse = true;

        int obtuseCount = countObtuseTriangles(cdt);
        int newL = static_cast<int>(L);

        // Data for output
        std::vector<std::string> steiner_points_x;
        std::vector<std::string> steiner_points_y;
        std::vector<std::pair<int, int>> steiner_edges;

        std::vector<double> steiner_points_x_double;  // New variable to store x as double
        std::vector<double> steiner_points_y_double;  // New variable to store y as double

        while (hasObtuse && counter < newL) {
            hasObtuse = false;  
            int current_obtuse_count = countObtuseTriangles(cdt);

            std::vector<CDT::Face_handle> obtuse_triangles;

            
            for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                // Skip triangles outside the region boundary
                if (!isTriangleWithinBoundary(fit, region_boundary_)) {
                    continue;
                }

                auto p1 = fit->vertex(0)->point();
                auto p2 = fit->vertex(1)->point();
                auto p3 = fit->vertex(2)->point();

                // Skip triangles with invalid points
                if (!CGAL::is_finite(p1.x()) || !CGAL::is_finite(p1.y()) ||
                    !CGAL::is_finite(p2.x()) || !CGAL::is_finite(p2.y()) ||
                    !CGAL::is_finite(p3.x()) || !CGAL::is_finite(p3.y())) {
                    continue;
                }

                if (isObtuseTriangle(p1, p2, p3)) {
                    hasObtuse = true;  
                    obtuse_triangles.push_back(fit);
                }
            }

            for (const auto& fit : obtuse_triangles){
                auto p1 = fit->vertex(0)->point();
                auto p2 = fit->vertex(1)->point();
                auto p3 = fit->vertex(2)->point();

                CDT best_cdt = cdt;  
                int best_obtuse_count = current_obtuse_count;

                // Evaluate Steiner point strategies
                for (int strategy = 0; strategy <= 6; ++strategy) {
                    Point steiner_point;

                    // Generate Steiner points using different strategies
                    if (strategy == 0) {
                        steiner_point = CGAL::circumcenter(p1, p2, p3);
                    } else if (strategy == 1) {
                        steiner_point = calculate_incenter(p1, p2, p3);
                    } else if (strategy == 2) {
                        steiner_point = getMidpointOfLongestEdge(p1, p2, p3); 
                    } else if (strategy == 3) {
                        steiner_point = CGAL::centroid(p1, p2, p3);
                    } else if (strategy == 4) {
                        steiner_point = calculate_perpendicular_bisector_point(p1, p2, p3); 
                    } else if (strategy == 5) {
                        steiner_point = projectPointOntoTriangle(Point(0, 0), p1, p2, p3);
                    } else if(strategy == 6 && !delaunay_){
                        // TODO if delaunay == true dont call strategy 6
                        steiner_point=getMeanAdjacentPoint(fit,best_cdt); //TODO SEGMETATION AFTER 1K LOOPS
                        if (steiner_point.x()==0 && steiner_point.y()==0){
                            continue;
                        }
                    }

                    // Validate the Steiner point
                    if (!CGAL::is_finite(steiner_point.x()) || !CGAL::is_finite(steiner_point.y())) {
                        continue;
                    }

                    if (!isPointInsideBoundary(
                            std::make_pair(steiner_point.x(), steiner_point.y()),
                            region_boundary_, points_)) {
                        continue;
                    }

                    // Simulate insertion of the Steiner point
                    CDT new_cdt = cdt;  // Copy the current triangulation
                    CDT::Vertex_handle new_vertex;
                    try {
                        new_vertex = new_cdt.insert(steiner_point);
                        processConvexPolygon(new_cdt);  // Ensure local re-triangulation
                    } catch (const CGAL::Precondition_exception& e) {
                        std::cerr << "Failed to insert Steiner point: " << e.what() << "\n";
                        continue;
                    }

                    // Count obtuse triangles in the new triangulation
                    int new_obtuse_count = countObtuseTriangles(new_cdt);

                    // Keep the best triangulation
                    if (new_obtuse_count < best_obtuse_count) {
                        std::cout<<"Placed a steiner : " << steiner_point << std::endl ;
                        best_cdt = new_cdt;
                        best_obtuse_count = new_obtuse_count;
                        // Store Steiner point
                        std::string x_frac = toFraction(steiner_point.x());
                        std::string y_frac = toFraction(steiner_point.y());

                        bool exists = false;
                        for (size_t i = 0; i < steiner_points_x.size(); ++i) {
                            if (steiner_points_x[i] == x_frac && steiner_points_y[i] == y_frac) {
                                exists = true;
                                break;
                            }
                        }
                            
                        if(!exists){
                            steiner_points_x.push_back(x_frac);
                            steiner_points_y.push_back(y_frac);

                            steiner_points_x_double.push_back(steiner_point.x());  // Store x as double
                            steiner_points_y_double.push_back(steiner_point.y());  // Store y as double
                        }
                        break;
                    }
                }

                // Update the CDT to the best found triangulation
                cdt = best_cdt;
                current_obtuse_count = best_obtuse_count;

            }
            counter++;
        }
        std::cout<<"Final obtuse triangles count : " << countObtuseTriangles(cdt) << std::endl;
        CGAL::draw(cdt);

        int index = 0;
        for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
            vit->info() = index++;
        }

        std::cout << "Total vertices in CDT after insertion: " << index << std::endl;

        // Find edges involving Steiner points
        std::set<int> steiner_indices;
        for (size_t i = 0; i < steiner_points_x_double.size(); ++i) {
            double x = steiner_points_x_double[i];
            double y = steiner_points_y_double[i];

            bool found = false;
            for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
                // Use squared distance for comparison
                if (CGAL::squared_distance(vit->point(), Point(x, y)) < 1e-9) {
                    steiner_indices.insert(vit->info());
                    std::cout << "Matched Steiner point (" << x << ", " << y 
                            << ") to vertex index: " << vit->info() << std::endl;
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cerr << "Failed to match Steiner point (" << x << ", " << y << ")" << std::endl;
            }
        }

        steiner_edges.clear();

        // Collect edges involving Steiner points
        for (auto eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit) {
            auto face = eit->first;
            int index = eit->second;

            auto vertex1 = face->vertex((index + 1) % 3);
            auto vertex2 = face->vertex((index + 2) % 3);

            if (vertex1 != nullptr && vertex2 != nullptr) {
                int index1 = vertex1->info();
                int index2 = vertex2->info();

                // Check if either vertex is a Steiner point
                if (steiner_indices.count(index1) > 0 || steiner_indices.count(index2) > 0) {
                    if (index1 > index2) std::swap(index1, index2);
                    steiner_edges.push_back({index1, index2});
                    std::cout << "Edge added: (" << index1 << ", " << index2 << ")" << std::endl;
                }
            }
        }

        std::cout << "Total Steiner edges: " << steiner_edges.size() << std::endl;

        if (steiner_indices.size() > steiner_edges.size()) {
            std::cerr << "Warning: Some Steiner points are not associated with any edge!" << std::endl;
        }


        
        writeOutputToFile("../output/output.json", steiner_points_x, steiner_points_y, steiner_edges, countObtuseTriangles(cdt));
    }


    int CDTProcessor::getPointIndex(const Point& point, const std::vector<std::pair<double, double>>& points) {
        for (size_t i = 0; i < points.size(); ++i) {
            // Compare CGAL Point with std::pair<double, double>
            if (std::abs(points[i].first - point.x()) < 1e-9 &&
                std::abs(points[i].second - point.y()) < 1e-9) {
                return static_cast<int>(i);
            }
        }
        return -1;  // Return -1 if the point is not found
    }


    CDT CDTProcessor::createCDT() {
        CDT cdt;
        std::vector<Point> cgal_points;
        
        for (const auto& point : points_) {
            cgal_points.emplace_back(point.first, point.second);
            cdt.insert(cgal_points.back());
        }

        // Insert constraints into the CGAL triangulation
        for (const auto& constraint : constraints_) {
            cdt.insert_constraint(cgal_points[constraint.first], cgal_points[constraint.second]);
        }
        return cdt;
    }

    void CDTProcessor::simulatedAnnealing(CDT& cdt,double alpha,double beta,int L){
        // Step 1: Compute initial energy of the triangulation
        int countOfSteinerPoints=0;
        double currentEnergy = calculateEnergy(cdt, alpha, beta ,countOfSteinerPoints);
        std::cout << "Initial energy: " << currentEnergy << "\n";

        // Step 2: Initialize temperature
        double temperature = 1.0;

        std::vector<std::string> steiner_points_x;
        std::vector<std::string> steiner_points_y;
        std::vector<std::pair<int, int>> steiner_edges;

        std::vector<double> steiner_points_x_double;  // New variable to store x as double
        std::vector<double> steiner_points_y_double;  // New variable to store y as double
        

        // Step 3: Main SA loop
        while (temperature >= 0) {
            std::cout << "Current temperature: " << temperature << "\n";

            std::vector<CDT::Face_handle> faces;
            for (auto fit = cdt.finite_faces_begin();fit != cdt.finite_faces_end();++fit){
                faces.push_back(fit);
            }

            
            // Step 4: Iterate over all finite faces in the CDT
            for (auto face : faces) {
                
                // Get the vertices of the triangle
                auto p1 = face->vertex(0)->point();
                auto p2 = face->vertex(1)->point();
                auto p3 = face->vertex(2)->point();

                // Compute squared edge lengths
                double a2 = squaredDistance(p2, p3);
                double b2 = squaredDistance(p1, p3);
                double c2 = squaredDistance(p1, p2);

                if (isObtuseTriangle(p1,p2,p3)){

                    // Step 5: Try different Steiner point strategies
                    for (int strategy = 0; strategy <= 5; ++strategy) {
                        Point steinerPoint;

                        // Generate Steiner points using different strategies
                        if (strategy == 0) {
                            steinerPoint = CGAL::circumcenter(p1, p2, p3);
                        } else if (strategy == 1) {
                            steinerPoint = calculate_incenter(p1, p2, p3);
                        } else if (strategy == 2) {
                            steinerPoint = getMidpointOfLongestEdge(p1, p2, p3);
                        } else if (strategy == 3) {
                            steinerPoint = CGAL::centroid(p1, p2, p3);
                        } else if (strategy == 4) {
                            steinerPoint = calculate_perpendicular_bisector_point(p1, p2, p3);
                        } else if (strategy == 5) {
                            steinerPoint = projectPointOntoTriangle(Point(0, 0), p1, p2, p3);
                        }

                        // Validate the Steiner point
                        if (!CGAL::is_finite(steinerPoint.x()) || !CGAL::is_finite(steinerPoint.y())) {
                            continue;
                        }

                        if (!isPointInsideBoundary(std::make_pair(steinerPoint.x(), steinerPoint.y()),region_boundary_, points_)) {
                            continue;
                        }

                        // Step 6: Modify triangulation by adding the Steiner point
                        CDT newCdt = cdt; // Create a copy of the CDT
                        
                        insertSteinerPoint(newCdt, {steinerPoint.x(), steinerPoint.y()});
                        int newCountOfSteinerPoints=countOfSteinerPoints + 1 ;

                        // Step 7: Calculate new energy
                        double newEnergy = calculateEnergy(newCdt, alpha, beta,newCountOfSteinerPoints);
                        double deltaEnergy = newEnergy - currentEnergy;

                        // Step 8: Apply Metropolis criterion
                        if (deltaEnergy < 0 || std::exp(-deltaEnergy / temperature) >= getRandomUniform()) {
                            // Accept the new triangulation
                            
                            cdt = newCdt;
                            countOfSteinerPoints=newCountOfSteinerPoints;
                            currentEnergy = newEnergy;

                            std::string x_frac = toFraction(steinerPoint.x());
                            std::string y_frac = toFraction(steinerPoint.y());
                            bool exists = false;
                            for (size_t i = 0; i < steiner_points_x.size(); ++i) {
                                if (steiner_points_x[i] == x_frac && steiner_points_y[i] == y_frac) {
                                    exists = true;
                                    break;
                                }
                            }
                            
                            if(!exists){
                                steiner_points_x.push_back(x_frac);
                                steiner_points_y.push_back(y_frac);

                                steiner_points_x_double.push_back(steinerPoint.x());  // Store x as double
                                steiner_points_y_double.push_back(steinerPoint.y());  // Store y as double
                                std::cout << "Accepted new configuration. Energy: " << currentEnergy << ", Steiner points: " << countOfSteinerPoints << "\n";  
                            }
                            break;

                        } else {
                            //std::cout << "Rejected configuration. Energy: " << currentEnergy << "\n";
                        }
                    }
                }

            }
            // Step 9: Decrease the temperature
            temperature -= 1.0 / L;
        }
        std::cout << "Final energy: " << currentEnergy << "\n";
        std::cout << "Total Steiner points added: " << countOfSteinerPoints << std::endl;
        
        std::cout<<"Final obtuse triangles count : " << countObtuseTriangles(cdt) << std::endl;
        CGAL::draw(cdt);

        int index = 0;
        for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
            vit->info() = index++;
        }

        std::cout << "Total vertices in CDT after insertion: " << index << std::endl;

        // Find edges involving Steiner points
        std::set<int> steiner_indices;
        for (size_t i = 0; i < steiner_points_x_double.size(); ++i) {
            double x = steiner_points_x_double[i];
            double y = steiner_points_y_double[i];

            bool found = false;
            for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
                // Use squared distance for comparison
                if (CGAL::squared_distance(vit->point(), Point(x, y)) < 1e-9) {
                    steiner_indices.insert(vit->info());
                    std::cout << "Matched Steiner point (" << x << ", " << y 
                            << ") to vertex index: " << vit->info() << std::endl;
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cerr << "Failed to match Steiner point (" << x << ", " << y << ")" << std::endl;
            }
        }

        steiner_edges.clear();

        // Collect edges involving Steiner points
        for (auto eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit) {
            auto face = eit->first;
            int index = eit->second;

            auto vertex1 = face->vertex((index + 1) % 3);
            auto vertex2 = face->vertex((index + 2) % 3);

            if (vertex1 != nullptr && vertex2 != nullptr) {
                int index1 = vertex1->info();
                int index2 = vertex2->info();

                // Check if either vertex is a Steiner point
                if (steiner_indices.count(index1) > 0 || steiner_indices.count(index2) > 0) {
                    if (index1 > index2) std::swap(index1, index2);
                    steiner_edges.push_back({index1, index2});
                    std::cout << "Edge added: (" << index1 << ", " << index2 << ")" << std::endl;
                }
            }
        }

        std::cout << "Total Steiner edges: " << steiner_edges.size() << std::endl;

        if (steiner_indices.size() > steiner_edges.size()) {
            std::cerr << "Warning: Some Steiner points are not associated with any edge!" << std::endl;
        }


        
        writeOutputToFile("../output/output.json", steiner_points_x, steiner_points_y, steiner_edges, countObtuseTriangles(cdt));
    }

    // Helper function to insert a Steiner point and update the CDT
    void CDTProcessor::insertSteinerPoint(CDT& cdt, const std::pair<double, double>& point) {
        //std::cout<<"inserting inside function..." << std::endl; 
        CGAL::Epick::Point_2 steinerPoint(point.first, point.second);
        cdt.insert(steinerPoint);
    }

    // Random number generator for uniform distribution
    double CDTProcessor::getRandomUniform() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> dis(0.0, 1.0);
        return dis(gen);
    }

    double CDTProcessor::calculateEnergy(const CDT& cdt, double alpha, double beta,int numberOfSteinerPoints) {
        int obtuseTriangles = countObtuseTriangles(cdt); 
        int steinerPoints = numberOfSteinerPoints;     
        return alpha * obtuseTriangles * obtuseTriangles  + beta * numberOfSteinerPoints;//TODO maybe change the obtuseTriangles^2 (based on the algorithm its only ^1)
    }

    Point CDTProcessor::getMeanAdjacentPoint(CDT::Face_handle face, CDT& cdt) {
        auto p1 = face->vertex(0)->point();
        auto p2 = face->vertex(1)->point();
        auto p3 = face->vertex(2)->point();

        Point circumcenter = CGAL::circumcenter(p1, p2, p3); // Current face circumcenter

        if (!CGAL::is_finite(circumcenter.x()) || !CGAL::is_finite(circumcenter.y())) {
            return Point(0, 0); // Return an invalid point if circumcenter is invalid
        }
        
        std::vector<Point> adjacentCircumcenters;

        // Check all neighbors for obtuse triangles
        for (int i = 0; i < 3; ++i) { // Each edge of the triangle
            auto neighbor = face->neighbor(i);
            if (cdt.is_infinite(neighbor) || neighbor == face) continue; // Skip infinite faces

            auto q1 = neighbor->vertex(0)->point();
            auto q2 = neighbor->vertex(1)->point();
            auto q3 = neighbor->vertex(2)->point();

            if (isObtuseTriangle(q1, q2, q3)) {
                Point neighborCircumcenter=CGAL::circumcenter(q1,q2,q3);
                if(CGAL::is_finite(neighborCircumcenter.x()) && CGAL::is_finite(neighborCircumcenter.y())){
                    adjacentCircumcenters.push_back(neighborCircumcenter);
                }
            }
        }

        // If fewer than two adjacent obtuse triangles, return an invalid point
        if (adjacentCircumcenters.size() < 2) {
            return Point(0, 0); // Return a "null" point to indicate no mean point is applicable
        }

        // Compute the mean of circumcenters
        Point meanPoint = circumcenter;
        for (const auto& adjCircumcenter : adjacentCircumcenters) {
            meanPoint = Point(
                (meanPoint.x() + adjCircumcenter.x()) / 2,
                (meanPoint.y() + adjCircumcenter.y()) / 2
            );
        }

        if (!isPointInsideBoundary(std::make_pair(meanPoint.x(), meanPoint.y()), region_boundary_, points_) ||
        !CGAL::is_finite(meanPoint.x()) || !CGAL::is_finite(meanPoint.y())) {
            return Point(0, 0); 
        }

        if (!CGAL::is_finite(meanPoint.x()) || !CGAL::is_finite(meanPoint.y()) || 
            !isPointInsideBoundary(std::make_pair(meanPoint.x(), meanPoint.y()), region_boundary_, points_)) {
            return Point(0, 0); 
        }

        return meanPoint;
    }

    void CDTProcessor::UpdatePheromones(std::map<CDT::Face_handle, std::map<Point, double>>& pheromone,const std::vector<AntSolution>& solutions, double lambda, double alpha, double beta) {
        for (const auto& solution : solutions) {
            // Εξάτμιση φερομόνης
            for (auto& pheromone_pair : pheromone[solution.face]) {
                pheromone_pair.second *= (1 - lambda);
            }

            // Ενίσχυση φερομόνης
            pheromone[solution.face][solution.steiner_point] += (1.0 / (1 + alpha * solution.improvement_metric + beta));
        }
    }

    void CDTProcessor::applyAntSolution(const AntSolution& solution, CDT& cdt) {
        cdt.insert(solution.steiner_point);
    }

    std::vector<AntSolution> CDTProcessor::resolve_conflicts(const std::vector<AntSolution>& solutions) {
        std::map<CDT::Face_handle, AntSolution> best_solutions;

        for (const auto& solution : solutions) {
            if (best_solutions.find(solution.face) == best_solutions.end() ||
                solution.improvement_metric > best_solutions[solution.face].improvement_metric) {
                best_solutions[solution.face] = solution;
            }
        }

        std::vector<AntSolution> resolved_solutions;
        for (const auto& pair : best_solutions) {
            resolved_solutions.push_back(pair.second);
        }
        return resolved_solutions;
    }

    std::vector<Point> CDTProcessor::generateSteinerOptions(const CDT::Face_handle& face) {
        std::vector<Point> options;

        Point p1 = face->vertex(0)->point();
        Point p2 = face->vertex(1)->point();
        Point p3 = face->vertex(2)->point();

        // Add options robustly
        try {
            options.push_back(CGAL::circumcenter(p1, p2, p3));
            options.push_back(calculate_incenter(p1, p2, p3));
            options.push_back(getMidpointOfLongestEdge(p1, p2, p3));
            options.push_back(CGAL::centroid(p1, p2, p3));
            options.push_back(calculate_perpendicular_bisector_point(p1, p2, p3));
            options.push_back(projectPointOntoTriangle(Point(0, 0), p1, p2, p3));
        } catch (const std::exception& e) {
            std::cerr << "Exception in generateSteinerOptions: " << e.what() << std::endl;
        }

        if (options.empty()) {
            std::cerr << "generateSteinerOptions: No valid options generated for face." << std::endl;
        }
        
        return options;
    }

    std::string CDTProcessor::selectStrategy(const std::vector<std::pair<std::string, double>>& method_probabilities, double total_probability) {
        if (total_probability <= 0.0) return "";

        double random_value = static_cast<double>(rand()) / RAND_MAX * total_probability;
        double cumulative_probability = 0.0;

        for (const auto& method_prob : method_probabilities) {
            cumulative_probability += method_prob.second;
            if (random_value <= cumulative_probability) {
                return method_prob.first;
            }
        }

        return ""; // Δεν πρέπει να φτάσουμε εδώ
    }

    void CDTProcessor::antColonyOptimization(CDT& cdt, double alpha, double beta, double xi, double psi, int lambda, int num_ants, int num_cycles) {
        std::map<CDT::Face_handle, std::map<Point, double>> pheromone;

        // Αρχικοποίηση φερομονών 
        for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
            auto steiner_options = generateSteinerOptions(face);
            for (const auto& point_option : steiner_options) {
                pheromone[face][point_option] = 1.0; // Αρχική τιμή φερομονών
            }
        }

        CDT best_cdt;
        copyTriangulation(cdt, best_cdt);
        int best_obtuse_count = countObtuseTriangles(best_cdt);
        int totalSteinerInserted=0;

        double best_energy=calculateEnergy(best_cdt,alpha,beta,totalSteinerInserted);

        for (int cycle = 0; cycle < num_cycles; ++cycle) {
            std::vector<AntSolution> all_solutions; // Συνολικές λύσεις μυρμηγκιών

            for (int ant = 0; ant < num_ants; ++ant) {
                std::vector<AntSolution> ant_solutions; // Λύσεις του συγκεκριμένου μυρμηγκιού

                // Επιλογή τριγώνων στην τύχη
                std::vector<CDT::Face_handle> candidate_faces;
                for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
                    if (isObtuseTriangle(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point())) {
                        candidate_faces.push_back(face);
                    }
                }

                if (candidate_faces.empty()) {
                    continue;
                }

                std::shuffle(candidate_faces.begin(), candidate_faces.end(), std::default_random_engine());

                for (auto& face : candidate_faces) {
                    // Υπολογισμός πιθανοτήτων επιλογής Steiner σημείων
                    auto steiner_options = generateSteinerOptions(face);
                    if (steiner_options.empty()) {
                        continue;
                    }

                    double total_probability = 0.0;
                    std::vector<std::pair<Point, double>> probabilities;

                    // Διαθέσιμες μέθοδοι (εξαιρείται η "projection")
                    std::vector<std::string> methods = {"circumcenter", "midpoint", "midpoint"};

                    // Υπολογισμός πιθανότητας για κάθε μέθοδο
                    std::vector<std::pair<std::string, double>> method_probabilities;
                    double total_method_probability = 0.0;

                    for (const auto& method : methods) {
                        double method_tau = 1.0; // Αρχική φερομόνη για τη μέθοδο (αν χρειαστεί, μπορεί να προσαρμοστεί)
                        double method_eta = calculateHeuristic(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point(), method);
                        double method_probability = std::pow(method_tau, xi) * std::pow(method_eta, psi);
                        method_probabilities.emplace_back(method, method_probability);
                        total_method_probability += method_probability;
                    }

                    // Επιλογή μεθόδου με βάση τις πιθανότητες
                    std::string selected_method = selectStrategy(method_probabilities, total_method_probability);

                    if (selected_method.empty()) {
                        //std::cout<<"No method selected fall back to circumcenter\n\n";
                        selected_method = "circumcenter";
                    }

                    // Υπολογισμός πιθανότητας για κάθε Steiner σημείο
                    for (const auto& option : pheromone[face]) {
                        double tau = std::pow(option.second, xi); // Φερομόνη του σημείου
                        double eta = calculateHeuristic(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point(), selected_method);
                        if (tau == 0.0 || eta == 0.0) {
                            continue;
                        }
                        double probability = tau * std::pow(eta, psi);
                        probabilities.emplace_back(option.first, probability);
                        total_probability += probability;
                    }

                    if (total_probability <= 0.0) {
                        continue;
                    }

                    // Επιλογή σημείου με πιθανότητα
                    Point selected_point = selectSteinerPoint(probabilities, total_probability);

                    // Προσομοίωση προσθήκης σημείου
                    CDT temp_cdt;
                    copyTriangulation(cdt, temp_cdt);
                    temp_cdt.insert(selected_point);
                    totalSteinerInserted++;
                    int new_obtuse_count = countObtuseTriangles(temp_cdt);
                    double energy_gain = best_obtuse_count - new_obtuse_count;

                    AntSolution solution = {face, selected_point, energy_gain};
                    ant_solutions.push_back(solution);
                }

                all_solutions.insert(all_solutions.end(), ant_solutions.begin(), ant_solutions.end());
            }

            // Επίλυση συγκρούσεων μεταξύ λύσεων μυρμηγκιών
            std::vector<AntSolution> resolved_solutions = resolve_conflicts(all_solutions);

            // Εφαρμογή λύσεων στο CDT
            CDT merged_cdt;
            copyTriangulation(best_cdt, merged_cdt);
            for (const auto& solution : resolved_solutions) {
                applyAntSolution(solution, merged_cdt);
            }

            // Ενημέρωση του καλύτερου τριγωνισμού
            int current_obtuse_count = countObtuseTriangles(merged_cdt);
            double current_energy=calculateEnergy(merged_cdt,alpha,beta,totalSteinerInserted);

            if(current_energy < best_energy || current_obtuse_count < best_obtuse_count){
                best_energy=current_obtuse_count == 0 ? 0 : current_energy;
                best_obtuse_count=current_obtuse_count;
                copyTriangulation(merged_cdt,best_cdt);
            }

            // Ενημέρωση φερομονών
            UpdatePheromones(pheromone, resolved_solutions, alpha, beta, xi);
        }

        // Επιστροφή του καλύτερου τριγωνισμού
        copyTriangulation(best_cdt, cdt);
        std::cout << "Total obtuse triangles after ant colony: " << countObtuseTriangles(cdt) << std::endl;
        CGAL::draw(cdt);
    }

    void CDTProcessor::copyTriangulation(const CDT& source, CDT& destination) {
        destination.clear(); // Clear the current triangulation

        // Step 1: Copy vertices
        std::map<CDT::Vertex_handle, CDT::Vertex_handle> vertex_map;
        for (auto vertex = source.finite_vertices_begin(); vertex != source.finite_vertices_end(); ++vertex) {
            auto new_vertex = destination.insert(vertex->point());
            vertex_map[vertex] = new_vertex; // Map old vertices to new vertices
        }

        // Step 2: Copy constrained edges only
        for (auto edge = source.constrained_edges_begin(); edge != source.constrained_edges_end(); ++edge) {
            auto face = edge->first; // The face containing this constrained edge
            int index = edge->second; // Index of the edge in the face

            // Get the vertices of the constrained edge
            CDT::Vertex_handle v1 = face->vertex((index + 1) % 3);
            CDT::Vertex_handle v2 = face->vertex((index + 2) % 3);

            // Insert the constrained edge into the destination CDT
            destination.insert_constraint(vertex_map[v1]->point(), vertex_map[v2]->point());
        }
    }

    Point CDTProcessor::selectSteinerPoint(const std::vector<std::pair<Point, double>>& options, double totalWeight) {
        // If options are empty, handle the fallback logic
        if (options.empty()) {
            std::cout<< " return point (0,0)" << std::endl;
            //std::cerr << "No Steiner point options available. Returning a default point." << std::endl;
            return Point(0, 0); // Replace with an appropriate fallback point or behavior
        }

        double randomValue = getRandomUniform() * totalWeight;
        for (const auto& [sp, weight] : options) {
            if ((randomValue -= weight) <= 0) {
                return sp;
            }
        }

        // If no point is selected due to precision issues, fallback to the first option
        std::cout << "SELECT Steiner point: (" << options.front().first.x() << ", " << options.front().first.y() << ")" << std::endl;
        return options.front().first; 
    }

    double CDTProcessor::calculateHeuristic(const Point& p1,const Point& p2,const Point&p3,const std::string& strategy){
        // Calculate circumcenter and circumradius (R)
        Point circumcenter = CGAL::circumcenter(p1, p2, p3);
        double circumradius = std::sqrt(CGAL::squared_distance(circumcenter, p1));  // Circumradius (R)

        // Find the longest edge and calculate the height from that edge
        double a = std::sqrt(CGAL::squared_distance(p2, p3)); // Edge length between p2 and p3
        double b = std::sqrt(CGAL::squared_distance(p1, p3)); // Edge length between p1 and p3
        double c = std::sqrt(CGAL::squared_distance(p1, p2)); // Edge length between p1 and p2
        double longestEdge = std::max({a, b, c});  // The longest side of the triangle

        // Calculate the area of the triangle using the determinant formula (shoelace formula)
        double area = std::abs(CGAL::area(p1, p2, p3));

        // Calculate height from the longest edge
        double height = (2 * area) / longestEdge;

        // Compute the radius-to-height ratio (ρ)
        double rho = circumradius / height;

        // Depending on the strategy, we calculate different heuristics
        if (strategy == "perpendicular") {
            // Vertex projection heuristic, most useful when ρ > 2.0
            if (rho > 2.0) {
                return std::max(0.0, (rho - 1) / rho);
            } else {
                return 0.0;
            }
        } 
        else if (strategy == "circumcenter") {
            // Circumcenter heuristic, probable when ρ ∈ [1.0, 2.0]
            if (rho >= 1.0 && rho <= 2.0) {
                return rho / (2 + rho);  // Prioritize circumcenter for moderately obtuse triangles
            } else {
                return 0.0;
            }
        } 
        else if (strategy == "midpoint") {
            // Midpoint heuristic, favored when ρ < 1.0
            if (rho < 1.0) {
                return std::max(0.0, (3 - 2 * rho) / 3);  // Favor midpoint for low ρ values
            } else {
                return 0.0;
            }
        } 
        else if (strategy == "mean_adjacent") {
            // Mean of adjacent obtuse triangles, if there are at least 2 adjacent obtuse triangles
            // For this case, you would need additional logic to check for adjacent obtuse triangles.
            return 1.0;  // Priority is given to this option in this case
        } 
        else {
            // Default case: no heuristic for unrecognized strategies
            std::cout << "Unkown method " << std::endl;
            return 0.0;
        }
    }



    // Βοηθητική συνάρτηση για τη μετατροπή αριθμών σε κλάσματα
    std::string CDTProcessor::toFraction(double value) {
        constexpr double epsilon = 1e-6;
        int sign = value < 0 ? -1 : 1;  
        value = std::fabs(value);       
        long long numerator = 1, denominator = 1;

        while (std::fabs(value * denominator - std::round(value * denominator)) > epsilon) {
            ++denominator;
        }
        numerator = static_cast<long long>(std::round(value * denominator));

        long long gcd = std::gcd(numerator, denominator);
        return (sign < 0 ? "-" : "") + std::to_string(numerator / gcd) + "/" + std::to_string(denominator / gcd);
    }

  void CDTProcessor::writeOutputToFile(
        const std::string& filepath,
        const std::vector<std::string>& steiner_points_x,
        const std::vector<std::string>& steiner_points_y,
        const std::vector<std::pair<int, int>>& steiner_edges,
        int obtuse_count) 
    {
        namespace json = boost::json;

        // Create JSON object
        json::object output;

        // Add content metadata
        output["content_type"] = "CG_SHOP_2025_Solution";
        output["instance_uid"] = instance_uid_; // Assumes access to `instance_uid_`
        output["method"] = method_;            // Assumes access to `method_`

        // Convert parameters to JSON object
        json::array parameters_json;
        for (const auto& param : parameters_) {
            std::ostringstream param_entry;
            param_entry << param.first << ":" << param.second;
            parameters_json.push_back(boost::json::value(param_entry.str()));
            // parameters_json[param.first] = json::value_from(param.second);
        }
        output["parameters"] = parameters_json;

        // Directly add Steiner points
        json::array steiner_points_x_json(steiner_points_x.begin(), steiner_points_x.end());
        json::array steiner_points_y_json(steiner_points_y.begin(), steiner_points_y.end());

        output["steiner_points_x"] = steiner_points_x_json;
        output["steiner_points_y"] = steiner_points_y_json;

        // Handle edges
        if (steiner_edges.empty()) {
            std::cout << "EMPTY EDGES\n\n";
        }

        json::array edges_json;
        for (const auto& edge : steiner_edges) {
            json::array edge_pair = {edge.first, edge.second};
            edges_json.push_back(edge_pair);
        }
        output["edges"] = edges_json;

        // Add remaining data
        output["obtuse_count"] = obtuse_count;

        // Write to file
        std::ofstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open output file: " + filepath);
        }

        // Write JSON to file with improved formatting
        file << "{\n";
        file << "  \"content_type\": " << json::serialize(output["content_type"]) << "," << std::endl;
        file << "  \"instance_uid\": " << json::serialize(output["instance_uid"]) << "," << std::endl;
        file << "  \"method\": " << json::serialize(output["method"]) << "," << std::endl;

        // Write parameters
        file << "  \"parameters\": " << json::serialize(output["parameters"]) << "," << std::endl;

        // Write Steiner points
        file << "  \"steiner_points_x\": " << json::serialize(output["steiner_points_x"]) << "," << std::endl;
        file << "  \"steiner_points_y\": " << json::serialize(output["steiner_points_y"]) << "," << std::endl;

        // Write edges
        file << "  \"edges\": " << json::serialize(output["edges"]) << "," << std::endl;

        // Write obtuse count
        file << "  \"obtuse_count\": " << output["obtuse_count"] << std::endl;
        file << "}" << std::endl;

        file.close();

        std::cout << "Output written to " << filepath << std::endl;
    }


}
