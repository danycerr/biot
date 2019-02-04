#ifndef ISOSTASY_DES_HS
#define ISOSTASY_DES_HS

#include "gmm/gmm.h"
//structure describing isostasy
//rotation matrix
struct isostasy_descriptor {
  double dz; //vertical displacement
  bgeot::base_matrix M,Mm1;
  bool active=false;
  std::unique_ptr<bgeot::small_vector<double>> dx,mdx;
  isostasy_descriptor() : M(3,3),Mm1(3,3){}
} ;
#endif