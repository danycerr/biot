#include "isostasy.hpp"
// ======================================================
void isostasy::move_mesh(getfem::mesh *mesh){
  // displacement vector for displace to orig
  mesh->translation(*(is_des_.dx));
  mesh->transformation(M_);
  mesh->translation(*(is_des_.mdx));
}
// ======================================================
// ======================================================
void isostasy::undo_move_mesh(getfem::mesh *mesh){
  // displacement vector for displace to orig
  mesh->translation(*(is_des_.dx));
  mesh->transformation(Mm1_);
  mesh->translation(*(is_des_.mdx));
}
// ======================================================

// ======================================================
void isostasy::set_center_of_rotation(getfem::mesh *mesh){
  // A trasformation for the squarred mesh
  bgeot::base_node Pmin, Pmax;
  mesh->bounding_box (Pmin, Pmax);
  std::vector<double> ds_bb; // bounding box of mesh
  for (int idir=0; idir<N_;idir++) ds_bb.push_back((Pmax[idir]-Pmin[idir])/2.);
//   dx_ =std::make_unique<bgeot::small_vector<double>> (-(Pmin[0]+ds_bb[0]),-(Pmin[1]+ds_bb[1]),-(Pmin[2]+ds_bb[2])); 
//   mdx_=std::make_unique<bgeot::small_vector<double>> ( (Pmin[0]+ds_bb[0]), (Pmin[1]+ds_bb[1]), (Pmin[2]+ds_bb[2])); 
  is_des_.dx =std::make_unique<bgeot::small_vector<double>> (-(Pmin[0]+ds_bb[0]),-(Pmin[1]+ds_bb[1]),-(Pmin[2]+ds_bb[2])); 
  is_des_.mdx=std::make_unique<bgeot::small_vector<double>> ( (Pmin[0]+ds_bb[0]), (Pmin[1]+ds_bb[1]), (Pmin[2]+ds_bb[2])); 
  
  // displacement vector for displace to orig
}
// ======================================================
// ======================================================
void isostasy::set_transformation(){
  dummy_transform();
//  read_transform)
gmm::copy(M_, is_des_.M);
gmm::copy(Mm1_, is_des_.Mm1);
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
  double dt=1.e+12;
  double tilt[]={3.1415/4*(time_ -20*dt)/((70-20)*dt),3.1415/10*(time_ -20*dt)/((70-20)*dt),3.1415/50*(time_ -20*dt)/((70-20)*dt)}; /// vector describing a casual rotation
  if (time_<20*dt) for (int i=0;i<3;i++) tilt[i]=0.;
  else if (time_ > 70*dt){
    tilt[0]=3.1415/4;tilt[1]=3.1415/10;tilt[2]=3.1415/50;
  }
  bgeot::base_matrix Mx(N_,N_);
  bgeot::base_matrix My(N_,N_);
  bgeot::base_matrix Mz(N_,N_);
  for (int i=0; i<N_;  i++)
    for (int j=0; j<N_; j++)    {Mx(i,j)=0.;My(i,j)=0.;Mz(i,j)=0.;M_(i,j)=0.;}
      {
	//=================== xz rotation ===================
	Mx(0,0) = cos(tilt[0]); Mx(0,2) = -sin(tilt[0]);
	Mx(2,0) = sin(tilt[0]);Mx(2,2) = cos(tilt[0]);
	Mx(1,1)=1.;
	//=================== yz rotation ===================
	My(1,1) = cos(tilt[1]); My(1,2) = -sin(tilt[1]);
	My(2,1) = sin(tilt[1]);My(2,2) = cos(tilt[1]);
	My(0,0)=1.;
	//=================== xy rotation ===================
// 	Mz(0,0) = cos(tilt[2]); Mz(0,1) = -sin(tilt[2]);
// 	Mz(1,0) = sin(tilt[2]);Mz(1,1) = cos(tilt[2]);
	Mz(2,2)=1.;
	Mz(0,0) =1.;Mz(1,1) =1.;
// 	My(2,2)=1.;My(0,0) =1.;My(1,1) =1.;
	gmm::mult(Mx, My, M_); // Mx * My ---> M
	gmm::mult(M_, Mz, M_); // Mx * My ---> M
        invert_transformation();
        }
}
// ======================================================