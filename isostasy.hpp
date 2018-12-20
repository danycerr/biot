#ifndef ISOSTASY_HS
#define ISOSTASY_HS
// getfem include
#include "getfem/getfem_mesh.h"
#include "gmm/gmm.h"

class isostasy{
public:
  isostasy(void):tilt_ {3.1415/20,3.1415/10}, dz_(0.), Mx_(N_,N_), Mxm1_(N_,N_)
                 , My_(N_,N_), Mym1_(N_,N_){
      for (int i=0; i<N_;  i++)
	for (int j=0; j<N_; j++)    {Mx_(i,j)=0.;My_(i,j)=0.;}
      {
	Mx_(0,0) = cos(tilt_[0]); Mx_(0,2) = -sin(tilt_[0]);
	Mx_(2,0) = sin(tilt_[0]);Mx_(2,2) = cos(tilt_[0]);
	Mx_(1,1)=1.;
	My_(1,1) = cos(tilt_[1]); My_(1,2) = -sin(tilt_[1]);
	My_(2,1) = sin(tilt_[1]);My_(2,2) = cos(tilt_[1]);
	My_(0,0)=1.;
	invert_transformation();
      }
  };
  double* get_tilt() {return tilt_;}
  double get_dz()    {return dz_;}
  void move_mesh(getfem::mesh *mesh);
  void undo_move_mesh(getfem::mesh *mesh);
  void set_center_of_rotation(getfem::mesh *mesh); // center point of rotation
private:
  double tilt_[2]; /// vector describing tilt_x and tilt_y
  double dz_;      /// vector describing tilt_x and tilt_y
  int N_=3;        /// dimension of the problem
  bgeot::base_matrix Mx_,Mxm1_;
  bgeot::base_matrix My_,Mym1_;
  std::unique_ptr<bgeot::small_vector<double>> dx_,mdx_;
  void invert_transformation(){
    gmm::copy(Mx_, Mxm1_);gmm::lu_inverse(Mxm1_);
    gmm::copy(My_, Mym1_);gmm::lu_inverse(Mym1_);
  }
};
#endif