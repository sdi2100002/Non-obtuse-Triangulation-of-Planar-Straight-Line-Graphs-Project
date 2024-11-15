#define BOOST_ALLOW_DEPRECATED_HEADERS
#include "triangulation.h"
#include <cmath>
#include <iostream>
#include <CGAL/algorithm.h>
#include <CGAL/draw_triangulation_2.h>
#include <utility>
#include "../CGAL-5.6.1/include/CGAL/convex_hull_2.h"

namespace Triangulation {

    // Constructor 
    CDTProcessor::CDTProcessor(const std::vector<std::pair<double, double>>& points, const std::vector<std::pair<int, int>>& constraints, const std::vector<int>& region_boundary,const std::string& instance_uid,const std::string& method,const std::map<std::string,double>& parameters) 
        : points_(points), constraints_(constraints), region_boundary_(region_boundary) , instance_uid_(instance_uid) , method_(method) , parameters_(parameters) {
    
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

        selectMethod(method_,parameters_);
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

    void CDTProcessor::selectMethod(const std::string& method,const std::map<std::string, double>& parameters){
        if (method == "local") {
            std::cout << "Selected method: Local Search\n";
            // std::cout << "Parameters:\n";
            // for (const auto& [key, value] : parameters) {
            //     std::cout << "  " << key << ": " << value << "\n";
            // }
            // // Placeholder for local search logic
            // std::cout << "Executing local search...\n";
        } 
        else if (method == "sa") { // Simulated Annealing
            std::cout << "Selected method: Simulated Annealing\n";
            // std::cout << "Parameters:\n";
            // for (const auto& [key, value] : parameters) {
            //     std::cout << "  " << key << ": " << value << "\n";
            // }
            // // Placeholder for simulated annealing logic
            // std::cout << "Executing simulated annealing...\n";
        } 
        else if (method == "ant") { // Ant Colony Optimization
            std::cout << "Selected method: Ant Colony Optimization\n";
            // std::cout << "Parameters:\n";
            // for (const auto& [key, value] : parameters) {
            //     std::cout << "  " << key << ": " << value << "\n";
            // }
            // // Placeholder for ant colony optimization logic
            // std::cout << "Executing ant colony optimization...\n";
        } 
        else {
            // Handle invalid method input
            std::cerr << "Error: Unknown method \"" << method << "\".\n";
            std::cerr << "Available methods: local, sa, ant\n";
        }
    }












}
