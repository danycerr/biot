
#include "temp_ls_dome.hpp" //class definition
//==============================================================
void temp_ls_dome::import_mesh(){
  std::cout<<"temp_ls_dome::importing for the dome problem"<<std::endl;
  // import mesh
//   getfem::import_mesh("gmsh:mesh/ringmeshes/pinch_2.msh",mesh);// labeled_domain=1; //official lk
  getfem::import_mesh("gmsh:mesh/ringmeshes/layer_cake.msh",mesh);// labeled_domain=1;
  // layer cake mesh from ring need 180 degree rotation
  bgeot::base_matrix M(N_,N_);  M(0,0) = 1.;M(1,1) = 1.;M(2,2) = -1.; // 180degree rotation x
  mesh.transformation(M);
}
//===============================================================
// level set function
base_small_vector temp_ls_dome::ls_function(const base_node P, double time,int num) {
//   std::cout<<"temp_ls_dome::ls_function"<<std::endl;
  scalar_type x = P[0]*p_des.l_ref, y = P[1]*p_des.l_ref, z=0;
  if (N_==3)  z = P[2]*p_des.l_ref;
  y = P[1]*p_des.l_ref;
  double dt=1.e+12;
  // time*=p_des.t_ref;
  base_small_vector res(2);
  // ====== gaussian ls for layer cake from ring
  double s = sqrt(2.5e+5); // dispersion factor
  double pi=3.1415;
  double x0=4954.; double y0=2000.; double z0=1470.; // center of the gaussian
  res[0] =  1.e+3 / (s * sqrt (2*pi)) * exp(-0.5 * pow(((x-x0)/s),2))
                                      * exp(-0.5 * pow(((y-y0)/s),2))
	    - (z-z0+200)/2740.; 
  res[1] = 0.;
  return res;
}

//========= assemby procedures ==========================
//=======================================================
// Assembly solution and iteration matrix for 
// the monolithic problem
//=======================================================
void temp_ls_dome::assembly(double dt,double time) {
  std::cout<< "temp_ls_dome::assembly" <<std::endl;
  // assembly matrix
  // resize and clear vector and matrix
  getfem::size_type nb_dof_p = mf_p.nb_dof();
  getfem::size_type nb_dof_p_x = nb_dof_p+ nb_x_dof_p;
  gmm::resize(B, nb_dof_p_x ); gmm::clear(B);
  gmm::resize(UP,  nb_dof_p_x); //gmm::clear(UP);
  gmm::resize(P, nb_dof_p); gmm::resize(P_old, nb_dof_p);
  gmm::resize(K, nb_dof_p_x, nb_dof_p_x); gmm::clear(K);
  std::cout<<" nb_dof_p "<< nb_dof_p<<std::endl;
  getfem::ga_workspace workspace; // generic workspace for the solution
  configure_workspace(workspace,dt); 
  // getfem::base_vector invdt(1); invdt[0] = 1/dt;
  // workspace.add_fixed_size_constant("invdt", invdt);

  workspace.add_fem_variable("p", mf_p, gmm::sub_interval(0, nb_dof_p), P);
  workspace.add_fem_constant("p_old", mf_p, P_old);
  workspace.add_fem_constant("pres", mf_p, press_);
  // ------------------ expressions --------------------------
  workspace.add_expression( "p.Test_p + tau*Kr*Grad_p.Grad_Test_p"
      , mim_ls_in, UNCUT_REGION);
  workspace.add_expression("penalty*p*Test_p *tau"	, mim_ls_out, UNCUT_REGION);
  workspace.assembly(2);
  gmm::copy(workspace.assembled_matrix(), gmm::sub_matrix( K ,
        gmm::sub_interval(0,  nb_dof_p),
        gmm::sub_interval(0,  nb_dof_p) 
        )
      );
  workspace.clear_expressions();


  //======= RHS =====================
  workspace.add_expression("+1.*Test_p*tau + p_old.Test_p - 1.e-19*C1*Er*Grad_pres.Grad_p*Test_p*tau", mim,UNCUT_REGION_IN);
  // workspace.add_expression("-0.25*Test_p*tau ", mim,LATERAL); //lateral heat sink
  // workspace.add_expression("nls*Test_p*tau ", mim_ls_in,UNCUT_REGION_IN);
  workspace.set_assembled_vector(B);
  workspace.assembly(1);
  workspace.clear_expressions();
  //Boudanry conditions // NITSCHE
  // getfem::base_vector penalty(1); penalty[0] = 2e+18; // 1---10
  // workspace.add_fixed_size_constant("penalty", penalty);
  //Matrix term for nietche bc
 //  workspace.add_expression("penalty/element_size*p*Test_p", mim, LEFT);
 //  workspace.add_expression("penalty/element_size*p*Test_p", mim, RIGHT);// 1 is the region		
 //  workspace.add_expression("-Grad_p.Normal*Test_p - Grad_Test_p.Normal*p ", mim, LEFT); 	
 //  workspace.add_expression("-Grad_p.Normal*Test_p - Grad_Test_p.Normal*p", mim, RIGHT); 	
  workspace.add_expression("2/element_size*p*Test_p*200", mim, TOP);	
  workspace.add_expression("-Grad_p.Normal*Test_p - Grad_Test_p.Normal*p", mim, TOP); 
  if(const_lateral_temp_){
      workspace.add_expression("2/element_size*p*Test_p*200", mim, LATERAL);	
      workspace.add_expression("-Grad_p.Normal*Test_p - Grad_Test_p.Normal*p", mim, LATERAL);
  }
  workspace.assembly(2);
  gmm::add(workspace.assembled_matrix(), gmm::sub_matrix( K ,
        gmm::sub_interval(0,  nb_dof_p),
        gmm::sub_interval(0,  nb_dof_p) 
        )
      );
  workspace.clear_expressions();
  // rhs for niche boudary condition
  // workspace.add_expression(" 0*penalty/element_size*Test_p -Grad_Test_p.Normal*0 ", mim_ls_in, LEFT);
  // workspace.add_expression(" 0*penalty/element_size*Test_p -Grad_Test_p.Normal*0 ", mim_ls_in, RIGHT);
   workspace.add_expression("2/element_size*top_temp*Test_p*200- top_temp*Grad_Test_p.Normal", mim, TOP);
   if(const_lateral_temp_)
       workspace.add_expression("2/element_size*p_old*Test_p*200- p_old*Grad_Test_p.Normal", mim, LATERAL);
   workspace.assembly(1);
   workspace.clear_expressions();
  // end boundary contions
  
  // Kout fotr enriched dof
  sparse_matrix_type K_out(nb_dof_p,nb_dof_p);
  for (int i=0; i< nb_dof_p; i++)  K_out(i,i)=1.e+0;
  std::cout<< "end kout"<< std::endl; 
  //end Kout
  // Kin for enriched dof
  sparse_matrix_type K_in( nb_dof_p,nb_dof_p);
  // NITSCHE for boundary 
//   if(const_lateral_temp_){
//      workspace.add_expression("2/element_size*p*Test_p*tau*2000", mim_ls_bd, CUT_REGION);// 1 is the region		
//      workspace.add_expression("-nlsv.Grad_p*Test_p*tau - nlsv.Grad_Test_p*p*tau ", mim_ls_bd, CUT_REGION);
//   }
  if(const_lateral_temp_) workspace.add_expression("2/element_size*p*Test_p*tau*2000", mim_ls_bd, CUT_REGION);// 1 is the region		
//    workspace.add_expression("-nlsv.Grad_p*Test_p*tau - nlsv.Grad_Test_p*p*tau ", mim_ls_bd, CUT_REGION);
  //NITSCHE
  workspace.add_expression( "+p.Test_p + tau*Kr*Grad_p.Grad_Test_p"
      , mim_ls_in, CUT_REGION);
  workspace.assembly(2);
  gmm::copy(workspace.assembled_matrix(),K_in);
  workspace.clear_expressions();
  workspace.add_expression("+[+1.].Test_p*tau + p_old.Test_p  - 1.e-19*C1*Er*Grad_pres.Grad_Test_p*tau", mim_ls_in,CUT_REGION);
//   if(const_lateral_temp_)
//     workspace.add_expression( "2/element_size*over_p*Test_p*tau*2000- nlsv.Grad_Test_p*over_p*tau", mim_ls_bd, CUT_REGION);
  if(const_lateral_temp_) workspace.add_expression("2/element_size*dome_t*Test_p*2000*tau", mim_ls_bd, CUT_REGION);
  workspace.assembly(1);
  workspace.clear_expressions();
 //  workspace.add_expression("2/element_size*p*Test_p", mim, LEFT);
 //  workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, LEFT); 	
 //  workspace.add_expression("2/element_size*p*Test_p", mim, RIGHT);
 //  workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, RIGHT);
 // workspace.assembly(2);
//  gmm::add(workspace.assembled_matrix(),K_in);
//  workspace.clear_expressions();
  // stabilization term
#ifdef  STAB_P
{
    getfem::mesh_region  inner_faces;
    inner_faces = getfem::inner_faces_of_mesh(mesh, CUT_REGION);

    workspace.add_expression("2*element_size*Grad_p.Normal*Grad_Test_p.Normal*1", mim, inner_faces);// 1 is the region		
    workspace.assembly(2);
    gmm::add(workspace.assembled_matrix(), K_in);
    workspace.clear_expressions();

 }
#endif

   // for (int i=0; i< nb_dof_p; i++)  K_in(i,i)+=1.e+0;
  std::cout<< "end kin"<< std::endl; 


  {// mapping for pressure
    std::cout<<"start mapping pressure"<<std::endl;
    size_type dof_shift = nb_dof_p;
    // On cut elements, the contribution is splitted in two parts,
    // corresponding to the sub-elements In (ls < 0) and Out (ls >= 0).
    // Basis functions on cut elements are copies of standard
    // functions, but are integrated on the In and Out sub-element (since
    // they are extended by zero on the rest of the element).
    // The simple rule to enrich the FE spaces is thus the following:
    // * If a dof index ii is In, its In contribution is 
    // is mapped to the index iIn = ii;
    // * on the other hand, its Out contribution is mapped 
    // in the extended dof range: ii = i + nb_dof_u, where eXt_dof[i] = ii;
    // * If a dof index ii is Out, proceeds analogously.
    // No interface terms here: only In-In and Out-Out dofs.
    for (size_type i = 0; i < eXt_dof.size(); ++i) {
      size_type ii = eXt_dof[i];
      double ls_i = ls_function(mf_p.point_of_basic_dof(ii),time, LS_TYPE)[0];
      // std::cout<< "ls values of dof"<< ii << " is "  << ls_i<< std::endl;
      B[dof_shift + i] =B[ii];
      for (size_type j = 0; j < eXt_dof.size(); ++j) {
        size_type jj = eXt_dof[j];
        double ls_j = ls_function(mf_p.point_of_basic_dof(jj),time, LS_TYPE)[0];
        // std::cout<< "ls_i "<< ls_i << " ls_j "<< ls_j << std::endl;
        if ( (ls_i <= 0) && (ls_j <= 0) ) {
          // i and j are both In
          K(ii , jj) += K_in(ii, jj);
          K(i + dof_shift, j  + dof_shift) += K_out(ii, jj);
        }
        else if ( (ls_i>= 0) && (ls_j >= 0) ) {
          // i and j are both Out
          K(ii, jj) += K_out(ii, jj);
          K(i + dof_shift, j  + dof_shift) += K_in(ii, jj);
        }
        else if ( (ls_i <= 0) && (ls_j>= 0) ) {
          // i is In, j is Out
          K(ii, j  + dof_shift) += K_in(ii, jj);
          K(i + dof_shift, jj) += K_out(ii, jj);
        }
        else {
          // i is Out, j is In
          K(i + dof_shift, jj) += K_in(ii, jj);
          K(ii, j  + dof_shift) += K_out(ii, jj);
        }
      }
    }
  }
  
  std::cout<< "End mapping"<<std::endl;
  //  std::cout<< K<<std::endl;

  // for (int i=0; i<  nb_dof_p; i++)  K(i,i)+=1.;

} // end assembly
//================================end assembly procedures =================0
// ===========================================
// method for generation of bcs zones
// ===========================================
void temp_ls_dome::gen_bc(){
  std::cout << "temp_ls_dome::gen_bc()"<< std::endl;
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);

  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    assert(i.is_face());
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);

//     if (gmm::abs(un[N_-1] - 1.0) < 1.0E-7) { // new Neumann face //for squared domain
    if (un[N_-1]  > 1.0E-1) { // new Neumann face //for squared domain
      mesh.region(TOP).add(i.cv(), i.f());
    } else if (gmm::abs(un[N_-1] + 1.0) < 1.0E-7) {
      mesh.region(BOTTOM).add(i.cv(), i.f());
    } else if (gmm::abs(un[N_-2] + 1.0) < 1.0E-7) {
      mesh.region(LEFT).add(i.cv(), i.f());
      mesh.region(LATERAL).add(i.cv(), i.f());
    } else if (gmm::abs(un[N_-2] - 1.0) < 1.0E-7) {
      mesh.region(RIGHT).add(i.cv(), i.f());
      mesh.region(LATERAL).add(i.cv(), i.f());
    }
    else if(N_=3){
      if (gmm::abs(un[N_-3] + 1.0) < 1.0E-7) {
        mesh.region(LEFTX).add(i.cv(), i.f());
      mesh.region(LATERAL).add(i.cv(), i.f());
      } else if (gmm::abs(un[N_-3] - 1.0) < 1.0E-7) {
        mesh.region(RIGHTX).add(i.cv(), i.f());
        mesh.region(LATERAL).add(i.cv(), i.f());
      }
    }
    else {
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
    }
  }
}
// end bc generations
// =======================================================