#include "isostasy.hpp"
// ======================================================
void isostasy::move_mesh(getfem::mesh *mesh){
  // displacement vector for displace to orig
  mesh->translation(*dx_);
  mesh->transformation(M_);
  mesh->translation(*mdx_);
}
// ======================================================
// ======================================================
void isostasy::undo_move_mesh(getfem::mesh *mesh){
  // displacement vector for displace to orig
  mesh->translation(*dx_);
  mesh->transformation(Mm1_);
  mesh->translation(*mdx_);
}
// ======================================================

// ======================================================
void isostasy::set_center_of_rotation(getfem::mesh *mesh){
  // A trasformation for the squarred mesh
  bgeot::base_node Pmin, Pmax;
  mesh->bounding_box (Pmin, Pmax);
  std::vector<double> ds_bb; // bounding box of mesh
  for (int idir=0; idir<N_;idir++) ds_bb.push_back((Pmax[idir]-Pmin[idir])/2.);
  dx_ =std::make_unique<bgeot::small_vector<double>> (-(Pmin[0]+ds_bb[0]),-(Pmin[1]+ds_bb[1]),-(Pmin[2]+ds_bb[2])); 
  mdx_=std::make_unique<bgeot::small_vector<double>> ( (Pmin[0]+ds_bb[0]), (Pmin[1]+ds_bb[1]), (Pmin[2]+ds_bb[2])); 
  // displacement vector for displace to orig
}
// ======================================================
// ======================================================
void isostasy::set_transformation(){
  dummy_transform();
//  read_transform)
}
// ======================================================
// ======================================================
std::vector<double>* isostasy::get_gravity(){
std::vector<double> *g=new std::vector<double>{0., 0., -1.};
gmm::mult(Mm1_, *g, *g);
return g;
}
// ======================================================

// ==================dummy transformation ====================================
void isostasy::dummy_transform(){
  double tilt[]={3.1415/4*time_/(80.*1.e+12),3.1415/10,3.1415/50}; /// vector describing a casual rotation
  bgeot::base_matrix Mx_(N_,N_);
  bgeot::base_matrix My_(N_,N_);
  bgeot::base_matrix Mz_(N_,N_);
  for (int i=0; i<N_;  i++)
    for (int j=0; j<N_; j++)    {Mx_(i,j)=0.;My_(i,j)=0.;Mz_(i,j)=0.;M_(i,j)=0.;}
      {
	//=================== xz rotation ===================
	Mx_(0,0) = cos(tilt[0]); Mx_(0,2) = -sin(tilt[0]);
	Mx_(2,0) = sin(tilt[0]);Mx_(2,2) = cos(tilt[0]);
	Mx_(1,1)=1.;
	//=================== yz rotation ===================
	My_(1,1) = cos(tilt[1]); My_(1,2) = -sin(tilt[1]);
	My_(2,1) = sin(tilt[1]);My_(2,2) = cos(tilt[1]);
	My_(0,0)=1.;
	//=================== xy rotation ===================
// 	Mz_(0,0) = cos(tilt[2]); Mz_(0,1) = -sin(tilt[2]);
// 	Mz_(1,0) = sin(tilt[2]);Mz_(1,1) = cos(tilt[2]);
	Mz_(2,2)=1.;
	Mz_(0,0) =1.;Mz_(1,1) =1.;
// 	My_(2,2)=1.;My_(0,0) =1.;My_(1,1) =1.;
	gmm::mult(Mx_, My_, M_); // Mx * My ---> M
	gmm::mult(M_, Mz_, M_); // Mx * My ---> M
        invert_transformation();
        }
}
// ======================================================