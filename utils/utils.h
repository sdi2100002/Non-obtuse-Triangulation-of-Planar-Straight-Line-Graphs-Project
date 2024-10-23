#ifndef UTILS_H
#define UTILS_H

#include <boost/json.hpp>


using namespace boost::json;

void printJsonFormatted(const object& obj);

bool loadJsonFile(const std::string& filename, std::string& output);

object parseJson(const std::string& input);

void retrieveFields(const object& obj, std::string& instance_uid, int& num_points, std::vector<int>& points_x, std::vector<int>& points_y,int& num_constraints, std::vector<std::pair<int, int>>& additional_constraints,std::vector<int>& region_boundary);

std::vector<std::pair<double, double>> createPointsVector(const std::vector<int>& points_x, const std::vector<int>& points_y);

void CallprocessTriangulation(const std::vector<std::pair<double, double>>& points, const std::vector<std::pair<int, int>>& additional_constraints, const std::vector<int>& region_boundary, const std::string& instance_uid);

#endif 