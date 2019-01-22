#include "biot_ls_dome.hpp"


//============== import mesh =======================
void biot_ls_dome::import_mesh(){
  std::cout<<"biot_ls_dome::import_mesh()"<<std::endl;
  getfem::import_mesh("gmsh:mesh/ringmeshes/layer_cake.msh",mesh);// labeled_domain=1;
//  bgeot::pgeometric_trans pgt = 
//  bgeot::geometric_trans_descriptor(p_des.MESH_TYPE);
//   std::vector<size_type> nsubdiv(N_);
//   int NX=p_des.nsubdiv;
//   std::fill(nsubdiv.begin(),nsubdiv.end(),NX);
//   getfem::regular_unit_mesh(mesh, nsubdiv, pgt, p_des.noised);
  // layer cake mesh from ring need 180 degree rotation
  bgeot::base_matrix M(N_,N_);  M(0,0) = 1.;M(1,1) = 1.;M(2,2) = -1.; // 180degree rotation x for layer cake
//   bgeot::base_matrix M(N_,N_);  M(0,0) = 4000.;M(1,1) = 4000.;M(2,2) = 4000.; // 180degree rotation x
  mesh.transformation(M);
}

// =============== lev_set function ==================
base_small_vector biot_ls_dome::ls_function(const base_node P, double time,int num) {
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
//   res[0] = -(pow((x-2000),2) + pow((y-2000),2) + pow((z),2) -pow(1000,2)); // dummy dome
//   res[0] = -(P[2]+0.4 -1.2*( exp( -( pow((2*P[0]-1)/0.5 ,2) ))  *exp( -( pow((2*P[1]-1)/0.5 ,2) )) ) );
              
  return res;
}

// ===========================================
// method for generation of bcs zones
// ===========================================
void biot_ls_dome::gen_bc(){
  std::cout << "biot_ls_dome::gen_bc()"<< std::endl;
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);

  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    assert(i.is_face());
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);

   if (un[N_-1]  > 1.0E-1) { // new Neumann face
      mesh.region(TOP).add(i.cv(), i.f());
    } else if (un[N_-1]  < -5.0E-1) {
      mesh.region(BOTTOM).add(i.cv(), i.f());
    } else if (gmm::abs(un[N_-2] + 1.0) < 1.0E-7) {
      mesh.region(LEFT).add(i.cv(), i.f());
    } else if (gmm::abs(un[N_-2] - 1.0) < 1.0E-7) {
      mesh.region(RIGHT).add(i.cv(), i.f());
    }
    else if(N_=3){
      if (gmm::abs(un[N_-3] + 1.0) < 1.0E-7) {
        mesh.region(LEFTX).add(i.cv(), i.f());
      } else if (gmm::abs(un[N_-3] - 1.0) < 1.0E-7) {
        mesh.region(RIGHTX).add(i.cv(), i.f());
      }
    }
    else {
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
    }
  }
}
// end bc generations
// =======================================================

// -------------------------------------------------
// Assembling pressure matrix for fixed stress pressure
// -------------------------------------------------

void biot_ls_dome::assembly_p(double dt, double time){
  gmm::clear(Bp); gmm::clear(Kp);//  gmm::clear(P_old);
  std::cout<<"biot_ls_dome::assembly_p(double dt)" <<std::endl;
  getfem::size_type nb_dof_u = mf_u.nb_dof();
  getfem::size_type nb_dof_p = mf_p.nb_dof();
  getfem::size_type nb_dof_p_x = nb_dof_p + nb_x_dof_p;
  getfem::size_type nb_dof_u_x = nb_dof_u + nb_x_dof_u;
  std::cout<< "total dofs " <<   mf_p.nb_dof() << "normal dof" << mf_p.nb_dof()  <<std::endl;

  gmm::resize(P, nb_dof_p_x); 
  gmm::resize(Kp, nb_dof_p_x, nb_dof_p_x);  gmm::resize(Bp, nb_dof_p_x);


  getfem::ga_workspace workspace;
  configure_workspace(workspace,dt);


  workspace.add_fem_variable("p", mf_p, gmm::sub_interval(0, nb_dof_p), P);
  workspace.add_fem_variable("p_old", mf_p, gmm::sub_interval(0,nb_dof_p), P_old);
  workspace.add_fem_variable("p_iter", mf_p, gmm::sub_interval(0,nb_dof_p), P_iter);
  workspace.add_fem_variable("u_old", mf_u, gmm::sub_interval(0, nb_dof_u), U_old);
  workspace.add_fem_variable("u_iter", mf_u, gmm::sub_interval(0, nb_dof_u), U_iter);

  // Pressure equation
  workspace.add_expression("(1+beta)*p.Test_p + tau*Kr*Grad_p.Grad_Test_p", mim,UNCUT_REGION_IN); // tau
  // workspace.set_assembled_matrix(Kp);
  workspace.assembly(2);
  gmm::copy(workspace.assembled_matrix(), gmm::sub_matrix(Kp,
        gmm::sub_interval(0, nb_dof_p),
        gmm::sub_interval(0, nb_dof_p) 
        ));
  workspace.clear_expressions();

  //======= RHS =====================
  workspace.add_expression("[+0]*Test_p*tau + p_old.Test_p + Test_p*Trace(Sym(Grad_u_old))", mim,UNCUT_REGION_IN); // tau
  workspace.add_expression("(beta)*p_iter.Test_p - Test_p*Trace(Sym(Grad_u_iter))", mim,UNCUT_REGION_IN); // tau
  workspace.set_assembled_vector(Bp);
  workspace.assembly(1);
  workspace.clear_expressions();

// ----- Matrix term --- bcs
  workspace.add_expression("2*penalty/element_size*p*Test_p*200", mim, TOP);
  workspace.add_expression("-Grad_p.Normal*Test_p - Grad_Test_p.Normal*p ", mim, TOP); 	
//   workspace.add_expression("2*penalty/element_size*p*Test_p", mim, LEFT);
//   workspace.add_expression("-Grad_p.Normal*Test_p*tau- Grad_Test_p.Normal*p*tau ", mim, LEFT); 	
//   workspace.add_expression("2*penalty/element_size*p*Test_p", mim, RIGHT);
//   workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, RIGHT); 	
//   if (N_== 3 ){
//     workspace.add_expression("2*penalty/element_size*p*Test_p", mim, LEFTX);
//     workspace.add_expression("-Grad_p.Normal*Test_p*tau- Grad_Test_p.Normal*p*tau ", mim, LEFTX); 	
//     workspace.add_expression("2*penalty/element_size*p*Test_p", mim, RIGHTX);
//     workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, RIGHTX); 	
//   }

  workspace.assembly(2);
  gmm::add(workspace.assembled_matrix(), gmm::sub_matrix(Kp,
        gmm::sub_interval(0, nb_dof_p),
        gmm::sub_interval(0, nb_dof_p) 
        ));
  workspace.clear_expressions();
  //rhs term
  // workspace.add_expression("0*penalty*p*Test_p", mim, TOP); 
  // workspace.assembly(1);
  // workspace.clear_expressions();

  ////dummy part of the matrix
  workspace.add_expression("penalty*p*Test_p ", mim, UNCUT_REGION_OUT);
  //workspace.add_expression("(1/bm + beta)*p.Test_p + tau*permeability*Grad_p.Grad_Test_p"	, mim_ls_out,UNCUT_REGION_OUT);
  //// workspace.add_expression("1.e+18*p*Test_p "	, mim_ls_bd);
  workspace.assembly(2);
  gmm::add(workspace.assembled_matrix(), gmm::sub_matrix(Kp,
        gmm::sub_interval(0, nb_dof_p),
        gmm::sub_interval(0, nb_dof_p) 
        ));
  workspace.clear_expressions();
// =======================================================
  sparse_matrix_type K_out(nb_dof_p,nb_dof_p);
  {
    for (int i=0; i< nb_dof_p; i++)  K_out(i,i)=1.e+0;
    std::cout<< "end kout"<< std::endl; 
  }
//   ======================================================
  // Kin for enriched dof
  std::cout<< "start kin"<< std::endl; 
  sparse_matrix_type K_in(nb_dof_p,nb_dof_p);
  std::vector<scalar_type> Bp_in(nb_dof_p, 0.0);
  {
    // workspace.set_assembled_vector(Bp_in); 
    // NITSCHE
//   workspace.add_expression("2/element_size*p*Test_p*3", mim_ls_bd, CUT_REGION);// 1 is the region // 200		
//   workspace.add_expression("-nlsv.Grad_p*Test_p*tau- nlsv.Grad_Test_p*p*tau ", mim_ls_bd, CUT_REGION); 
    //NITSCHE
//  workspace.add_expression( "permeability*tau*[0,1].Grad_p*Test_p ", mim_ls_bd, CUT_REGION);
    workspace.add_expression( "(1+beta)*p.Test_p + tau*Kr*Grad_p.Grad_Test_p", mim_ls_in, CUT_REGION);
    workspace.assembly(2);
    gmm::copy(workspace.assembled_matrix(),K_in);
    workspace.clear_expressions();
    workspace.add_expression("+[0].Test_p*tau + p_old.Test_p + Test_p*Div_u_old + (beta)*p_iter.Test_p - Test_p*Trace(Sym(Grad_u_iter))", mim_ls_in,CUT_REGION);
//     workspace.add_expression("(beta)*p_iter.Test_p - Test_p*Trace(Sym(Grad_u_iter))", mim_ls_in,CUT_REGION); // tau
    workspace.set_assembled_vector(Bp_in);
    workspace.assembly(1);
    workspace.clear_expressions();
//     workspace.add_expression("2/element_size*p*Test_p", mim, LEFT);
//     workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, LEFT); 	
//     workspace.add_expression("2/element_size*p*Test_p", mim, RIGHT);
//     workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, RIGHT);
//     if (N_==3){
//       workspace.add_expression("2/element_size*p*Test_p", mim, LEFT);
//       workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, LEFT); 	
//       workspace.add_expression("2/element_size*p*Test_p", mim, RIGHT);
//       workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, RIGHT);
//     }
// 
//     workspace.assembly(2);
//     gmm::add(workspace.assembled_matrix(),K_in);
//     workspace.clear_expressions();
    //pstab stabilization term
#ifdef STAB_P 
{
// >>>>>>> ls_temp
    getfem::mesh_region  inner_faces;
    inner_faces = getfem::inner_faces_of_mesh(mesh, CUT_REGION);

    workspace.add_expression("2*element_size*Grad_p.Normal*Grad_Test_p.Normal*0.1", mim, inner_faces);// 1 is the region		
    workspace.assembly(2);
    gmm::add(workspace.assembled_matrix(), K_in);
    workspace.clear_expressions();

}
#endif
    std::cout<< "end kin"<< std::endl; 
    {// mapping for pressure
      std::cout<<"start mapping pressure"<<std::endl;
      size_type dof_shift = nb_dof_p;

      int j_shift=eXt_dof.size()*eXt_dof.size();
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

      double time_1 = gmm::uclock_sec();  
      if (0) for (size_type i = 0; i < eXt_dof.size(); ++i) {
        size_type ii = eXt_dof[i];
        double ls_i = ls_function(mf_p.point_of_basic_dof(ii),time, LS_TYPE)[0];
        // std::cout<< "ls values of dof"<< ii << " is "  << ls_i<< std::endl;
        // Bp[nb_dof_p+i] =Bp[ii];
        if(ls_i<=0) {Bp[ii]+=Bp_in[ii]; }
        else        {Bp[nb_dof_p+i]+= Bp_in[ii];}
        for (size_type j = 0; j < eXt_dof.size(); ++j) {
          size_type jj = eXt_dof[j];
          double ls_j = ls_function(mf_p.point_of_basic_dof(jj),time, LS_TYPE)[0];
          // std::cout<< "ls_i "<< ls_i << " ls_j "<< ls_j << std::endl;
          if ( (ls_i <= 0) && (ls_j <= 0) ) {
            // i and j are both In
            Kp(ii , jj) += K_in(ii, jj);
            Kp(i + dof_shift, j  + dof_shift) += K_out(ii, jj);
            // std::cout<< " i " << ii << " " << pin_index_[i* eXt_dof.size()+j]<<"   ";
            // std::cout<< " j " << jj << " " << pin_index_[i* eXt_dof.size()+j+j_shift]<<std::endl;
          }
          else if ( (ls_i>= 0) && (ls_j >= 0) ) {
            // i and j are both Out
            Kp(ii, jj) += K_out(ii, jj);
            Kp(i + dof_shift, j  + dof_shift) += K_in(ii, jj);
          }
          else if ( (ls_i <= 0) && (ls_j>= 0) ) {
            // i is In, j is Out
            Kp(ii, j  + dof_shift) += K_in(ii, jj);
            Kp(i + dof_shift, jj) += K_out(ii, jj);
          }
          else {
            // i is Out, j is In
            Kp(i + dof_shift, jj) += K_in(ii, jj);
            Kp(ii, j  + dof_shift) += K_out(ii, jj);
          }
        }
      }

      std::cout << "... time to map pressure with if : " << gmm::uclock_sec() - time_1 << " seconds\n";
      double time_2 = gmm::uclock_sec();  
      for (size_type i = 0; i < eXt_dof.size(); ++i) {
        size_type ii = eXt_dof[i]; 
        Bp[pin_index_[i* eXt_dof.size()]]+=Bp_in[ii];
        for (size_type j = 0; j < eXt_dof.size(); ++j) {
          size_type jj = eXt_dof[j];
          Kp(pin_index_[i* eXt_dof.size()+j],pin_index_[i* eXt_dof.size()+j+j_shift])+=K_in(ii, jj);
          Kp(pout_index_[i* eXt_dof.size()+j],pout_index_[i* eXt_dof.size()+j+j_shift])+=K_out(ii, jj);
        }
      }

      std::cout << "... time to map pressure no if : " << gmm::uclock_sec() - time_2 << " seconds\n";
      // std::cin.ignore();
    }
  }
  return;
}
// -------------------------------------------------
// Assembling displacement matrix for fixed stress
// -------------------------------------------------
void biot_ls_dome::assembly_u(double dt,double time){
    std::cout<<"biot_ls_dome::assembly_u(double dt)" <<std::endl;
    getfem::size_type nb_dof_u = mf_u.nb_dof();
    getfem::size_type nb_dof_p = mf_p.nb_dof();
    getfem::size_type nb_dof_p_x = nb_dof_p + nb_x_dof_p;
    getfem::size_type nb_dof_u_x = nb_dof_u + nb_x_dof_u;

    gmm::clear(Bu);  gmm::clear(Ku);
    gmm::resize(U_old, nb_dof_u);
    std::cout<<"nb_dof_u_x "<<nb_dof_u_x<<std::endl;
    gmm::resize(U, nb_dof_u_x); gmm::resize(Ku, nb_dof_u_x, nb_dof_u_x);  gmm::resize(Bu, nb_dof_u_x);
    getfem::ga_workspace workspace; configure_workspace(workspace,dt);
    
    workspace.add_fem_variable("p_iter", mf_p, gmm::sub_interval(0,nb_dof_p), P_iter);
    workspace.add_fem_variable("u", mf_u, gmm::sub_interval(0, nb_dof_u), U);
    workspace.add_fem_variable("u_old", mf_u, gmm::sub_interval(0, nb_dof_u), U_old);
    
    // ----- Boundary condisions displacemenet
    if(N_==2){
      // workspace.add_expression("[0,-2*force].Test_u" , mim, TOP);    //neumann disp
      // workspace.add_expression("[0,+2*force].Test_u" , mim, BOTTOM); //neumann disp
      // workspace.add_expression("penalty*Grad_u(2,1)*Grad_Test_u(2,1)" , mim, TOP); //neumann disp
      // workspace.add_expression("penalty*Grad_u(2,1)*Grad_Test_u(2,1)" , mim, BOTTOM); //neumann disp
    }
    else if(N_==3){
      // workspace.add_expression("[0,0,-2*force].Test_u" , mim, TOP);    //neumann disp
      // workspace.add_expression("[0,0,+2*force].Test_u" , mim, BOTTOM); //neumann disp
    }
    workspace.set_assembled_vector(Bu);
//     workspace.assembly(1);
//     workspace.clear_expressions();
    // ----- momentum equation
    std::cout<<"biot_ls_dome assembling sym grad and bottom fix"<<std::endl;
    workspace.add_expression("2*Er*Sym(Grad_u):Grad_Test_u + C2*Er*Div_u*Div_Test_u", mim , UNCUT_REGION_IN); // stress tensor 
    workspace.add_expression("penalty*u.Test_u" , mim, BOTTOM); //neumann disp
    workspace.assembly(2);   
    gmm::copy(workspace.assembled_matrix(), gmm::sub_matrix(Ku,
          gmm::sub_interval(0, nb_dof_u),
          gmm::sub_interval(0, nb_dof_u) 
          ));
    workspace.clear_expressions();

    std::cout<<"biot_ls_dome assembling gravity term and pressure "<<std::endl;
    if(N_==2) workspace.add_expression("[0,-1].Test_u"  ,  mim_ls_in,  UNCUT_REGION_IN);
    if(N_==3) workspace.add_expression("gravity.Test_u" ,  mim_ls_in,  UNCUT_REGION_IN);
    workspace.add_expression("C1*p_iter*Div_Test_u "    ,  mim_ls_in,  UNCUT_REGION_IN);
    
    if(N_==3) workspace.add_expression("topload*gravity.Test_u" , mim, TOP);    //neumann disp
    workspace.assembly(1);
    workspace.clear_expressions();
    //dummy part of the matrix
    workspace.add_expression("u.Test_u "	, mim_ls_out, UNCUT_REGION_OUT);
    workspace.assembly(2);    
    gmm::add(workspace.assembled_matrix(), gmm::sub_matrix(Ku,
          gmm::sub_interval(0, nb_dof_u),
          gmm::sub_interval(0, nb_dof_u) 
          ));
    workspace.clear_expressions();

    // Enriched dof
    sparse_matrix_type K_out(nb_dof_u ,nb_dof_u );
    for (int i=0; i< nb_dof_u; i++)  K_out(i,i)=1.e+0;
    std::cout<< "end kout"<< std::endl; 
    // Kin for enriched dof
    sparse_matrix_type K_in(nb_dof_u, nb_dof_u );
    // internal for displacement
    std::cout<< "Assembling interal ext dof matrix "<< std::endl; 
    workspace.add_expression( "2*Er*Sym(Grad_u):Grad_Test_u + Er*C2*Div_u*Div_Test_u", mim_ls_in, CUT_REGION);
    workspace.add_expression("penalty*u.Test_u" , mim_ls_bd, CUT_REGION); //neumann disp
    workspace.add_expression("penalty*u.Test_u" , mim_ls_in, BOTTOM); //neumann disp
    workspace.assembly(2);
    gmm::copy(workspace.assembled_matrix(),K_in);
    workspace.clear_expressions();
    std::cout<< "end kin"<< std::endl; 

    std::vector<scalar_type> B_in(nb_dof_u, 0.0);
    workspace.set_assembled_vector(B_in);
    std::cout<< "Assembling interal ext dof rhs "<< std::endl; 
    if(N_==2) workspace.add_expression("[0,-1].Test_u", mim_ls_in,CUT_REGION);
    if(N_==3) workspace.add_expression("gravity.Test_u", mim_ls_in,CUT_REGION);
    if(N_==3) workspace.add_expression("topload*gravity.Test_u" , mim, TOP);    //neumann disp
//     if(N_==3) workspace.add_expression("topload*gravity.Test_u" , mim_ls_in, TOP);    //neumann disp
    workspace.add_expression("C1*p_iter*Div_Test_u ",  mim_ls_in, CUT_REGION);
    workspace.assembly(1);
    workspace.clear_expressions();
    std::cout<< "end Bin"<< std::endl;

    // stabilization term
//     {
//     getfem::mesh_region  inner_faces;
//     inner_faces = getfem::inner_faces_of_mesh(mesh, CUT_REGION);
// 
//     workspace.add_expression("2*element_size*(Sym(Grad_u).Normal).(Sym(Grad_Test_u).Normal)", mim, inner_faces);// 1 is the region		
//     workspace.assembly(2);
//     std::cout<< workspace.assembled_matrix()<<std::endl;
//     gmm::add(workspace.assembled_matrix(), K_in);
//     workspace.clear_expressions();
// 
//   }
// 
    // Mapping
    {// mapping for displacement
      std::cout<<"start mapping displacement"<<std::endl;
      double time_1 = gmm::uclock_sec();  
      size_type dof_shift = nb_dof_u ;
      int j_shift=eXt_dof_u.size()*eXt_dof_u.size();
      if(0)for (size_type i = 0; i < eXt_dof_u.size(); ++i) {
        size_type ii = eXt_dof_u[i];
        double ls_i = ls_function(mf_u.point_of_basic_dof(ii),time, LS_TYPE)[0];
        // B[dof_shift +i ] = B[ii];
        if (ls_i <= 0) Bu[ii]+=B_in[ii];
        else           Bu[dof_shift+i]+=B_in[ii];
        for (size_type j = 0; j < eXt_dof_u.size(); ++j) {
          size_type jj = eXt_dof_u[j];
          double ls_j = ls_function(mf_u.point_of_basic_dof(jj),time, LS_TYPE)[0];
          // std::cout<< "ls_i "<< ls_i << " ls_j "<< ls_j << std::endl;
          if ( (ls_i <= 0) && (ls_j <= 0) ) {
            // i and j are both In
            Ku(ii , jj ) += K_in(ii, jj);
            Ku(i + dof_shift, j  + dof_shift) += K_out(ii, jj);
          }
          else if ( (ls_i>= 0) && (ls_j >= 0) ) {
            // i and j are both Out
            Ku(ii, jj) += K_out(ii, jj);
            Ku(i + dof_shift, j + dof_shift) += K_in(ii, jj);
          }
          else if ( (ls_i <= 0) && (ls_j>= 0) ) {
            // i is In, j is Out
            Ku(ii, j  + dof_shift) += K_in(ii, jj);
            Ku(i + dof_shift, jj) += K_out(ii, jj);
          }
          else {
            // i is Out, j is In
            Ku(i + dof_shift, jj) += K_in(ii, jj);
            Ku(ii, j  + dof_shift) += K_out(ii, jj);
          }
        }
      }


      std::cout << "... time to map disp with if : " << gmm::uclock_sec() - time_1 << " seconds\n";
      double time_2 = gmm::uclock_sec();  
      for (size_type i = 0; i < eXt_dof_u.size(); ++i) {
        size_type ii = eXt_dof_u[i]; 
        Bu[uin_index_[i* eXt_dof_u.size()]]+=B_in[ii];
        for (size_type j = 0; j < eXt_dof_u.size(); ++j) {
          size_type jj = eXt_dof_u[j];
          Ku(uin_index_[i* eXt_dof_u.size()+j],uin_index_[i* eXt_dof_u.size()+j+j_shift])+=K_in(ii, jj);
          Ku(uout_index_[i* eXt_dof_u.size()+j],uout_index_[i* eXt_dof_u.size()+j+j_shift])+=K_out(ii, jj);
        }
      }
      std::cout << "... time to map dipl no if : " << gmm::uclock_sec() - time_2 << " seconds\n";
      // std::cin.ignore();
    }
    // end enriched dof

    //  std::cout<<Ku<<std::endl;std::cin.get();
  } // end assembling momentum matrix for fixed stress approach