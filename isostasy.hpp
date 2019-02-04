#ifndef ISOSTASY_HS
#define ISOSTASY_HS
// getfem include
#include "getfem/getfem_mesh.h"
#include "gmm/gmm.h"

// isostasy descriptor
#include "isostasy_descriptor.hpp"

#include "isostasy_from_file.hpp"

class isostasy{
public:
  isostasy(void): dz_(0.), M_(N_,N_), Mm1_(N_,N_) // file isostasy name and timsteps of the file
      {
      for (int i=0; i<N_;  i++)
	for (int j=0; j<N_; j++) M_(i,j)=0.;
      for (int i=0; i<N_;  i++)  M_(i,i)=1.; // initial trivial transvormation
      invert_transformation();
      is_des_.active=true;
      f_isos=std::make_unique<isostasy_from_files> ("mesh/rotation.txt",4.3155e+09);
      };
      
  double get_dz()    {return dz_;}
  void move_mesh(getfem::mesh *mesh);
  void undo_move_mesh(getfem::mesh *mesh);
  void set_center_of_rotation(getfem::mesh *mesh); // center point of rotation
  inline void set_time(double time){time_=time;}
  inline isostasy_descriptor * get_descriptor(){return &is_des_;} // get the current descriptor of isostasy
  void set_transformation();
  std::vector<double>* get_gravity();
private:
  double dz_;      /// vector describing vertical disp
  double time_=0;  /// time
  isostasy_descriptor is_des_;
  std::unique_ptr<isostasy_from_files> f_isos;
  int N_=3;        /// dimension of the problem
  bgeot::base_matrix M_,Mm1_;
  void invert_transformation(){
    gmm::copy(M_, Mm1_);gmm::lu_inverse(Mm1_);
//     gmm::copy(gmm::transposed(Mx_), Mxm1_);
  }
  void dummy_transform(); // dummy transformation
};
#endif