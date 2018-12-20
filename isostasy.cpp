#include "isostasy.hpp"
// ======================================================
void isostasy::move_mesh(getfem::mesh *mesh){
  // displacement vector for displace to orig
  mesh->translation(*dx_);
  mesh->transformation(Mx_);
  mesh->transformation(My_);
  mesh->translation(*mdx_);
}
// ======================================================
// ======================================================
void isostasy::undo_move_mesh(getfem::mesh *mesh){
  // displacement vector for displace to orig
  mesh->translation(*dx_);
  mesh->transformation(Mxm1_);
  mesh->transformation(Mym1_);
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