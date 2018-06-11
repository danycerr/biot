#ifndef HORIZONHF
#define HORIZONHF
#include <iostream>
#include <fstream>
#include <memory>
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_import.h>

class horizon{
public:
  horizon(const char* file);
  int find_element(std::vector<double>& pt);
  int up_down(std::vector<double>& pt, int el);
private:
  getfem::mesh mesh_;
  void quick_search(std::vector<double>&pt,std::vector<int>& maybe_element, double eps=0.0);
  int fine_search(std::vector<double>&pt,std::vector<int>& maybe_element);
  inline double sign (std::vector<double> p1, std::vector<double> p2, std::vector<double> p3)
{
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]);
}
bool is_in (std::vector<double>  pt, std::vector<double>  v1, std::vector<double>  v2, std::vector<double>  v3);

};





#endif