#ifndef ISOSTASY_HS
#define ISOSTASY_HS
// getfem include
#include "getfem/getfem_mesh.h"

class isostasy{
public:
  isostasy(void):tilt_ {0.,0.}, dz_(0.) {};
  double* get_tilt() {return tilt_;}
  double get_dz() {return dz_;}  
  void move_mesh(getfem::mesh& mesh);
private: 
  double tilt_[2]; /// vector describing tilt_x and tilt_y
  double dz_; /// vector describing tilt_x and tilt_y
};
#endif