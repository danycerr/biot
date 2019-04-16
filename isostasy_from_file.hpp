#ifndef ISOSTASY_FROMFILE_HS
#define ISOSTASY_FROMFILE_HS
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "getfem/getfem_config.h"
// // // // // // // // // // // // // // // // // // // // 
class isostasy_from_files {
protected:
  std::vector<double> rotations_;
  std::vector<double> tilts_;
  std::vector<double> disp_z_;
  double dt_; // time steps of the file
private:
public:
  isostasy_from_files(std::string, std::string, std::string, double dt);
  void get_rotation(double time,bgeot::base_matrix& M );
  void get_tilt(double time,std::vector<double> & tilt );
  void get_disp(double time,double & dz );
};


#endif