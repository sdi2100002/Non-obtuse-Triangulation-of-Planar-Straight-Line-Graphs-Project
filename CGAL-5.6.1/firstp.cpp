#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <fstream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef std::istream_iterator< Point_2 > point2_iterator;
typedef std::vector<Point_2> Points;
typedef std::vector<Point_2>::iterator pveciterator;
int main()
{
 Points result;
 std::ifstream in("data.in");
 point2_iterator begin(in);
 point2_iterator end;
 CGAL::convex_hull_2( begin, end, std::back_inserter(result) );
 std::cout << result.size() << " points on the convex hull" << std::endl;
 for (pveciterator iter=result.begin(); iter!=result.end(); ++iter)
 std::cout << *iter << std::endl;
 return 0;
}