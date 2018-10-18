#include "temp_ls.hpp" 
// Preconditioner for dispacement problem
// mPR Bu diag(exdof) 
// PRu diagonal
#define DISP_PRECOND_PARAM PRu
// fix cut height 
#define H_PARAM 2666.67

// #define STAB_P (0)

// #define L2_NORM

void templs_problem::init(void) {

  std::cout<< "templs_problem::init "  << std::endl;

  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(p_des.MESH_TYPE);

  // Mesh generation
  N_ = pgt->dim();
  std::vector<size_type> nsubdiv(N_);
  int NX=p_des.nsubdiv;
  std::fill(nsubdiv.begin(),nsubdiv.end(),NX);
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt, p_des.noised);
  // // A trasformation for the squared mesh
  p_des.l_ref=4000;


  // import mesh
  // getfem::import_mesh("gmsh:mesh/square2k.msh",mesh);
  // getfem::import_mesh("gmsh:mesh/squarepinch.msh",mesh);
  //  getfem::import_mesh("gmsh:mesh/squarepinch_fine.msh",mesh);
  // dal::bit_vector b; b.add(0);
  // mesh.Bank_refine(b);
  // mesh.Bank_refine(mesh.convex_index());
  bgeot::base_matrix M(N_,N_);
  for (size_type i=0; i < N_; ++i) {
    M(i,i) = 1.;// /p_des.l_ref;
  }
  //  if (N>1) { M(0,1) = 0; }
  //
  mesh.transformation(M);
  // // End of mesh generation
  // std::cout << "has 50 "<< mesh.has_region(50)<<std::endl;std::cin.ignore();

  // set  integration methods  
  getfem::pfem pf_p = getfem::fem_descriptor(p_des.FEM_TYPE_P); // for pressure
  getfem::pintegration_method ppi = getfem::int_method_descriptor(p_des.INTEGRATION);
  getfem::pintegration_method simp_ppi = getfem::int_method_descriptor(p_des.SIMPLEX_INTEGRATION);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_p.set_finite_element(mesh.convex_index(), pf_p); //  finite element for pressure
  mf_coef.set_finite_element(mesh.convex_index(),
      getfem::classical_fem(pgt,0));         // p0 for coefficient
  mf_coef_v.set_qdim(N_);          
  mf_coef_v.set_finite_element(mesh.convex_index(),
      getfem::classical_fem(pgt,0));         // p0 for coefficient
  // boundary conditions zones
  gen_bc();

  // initialization of level set
  mls.add_level_set(ls);
  ls.reinit(); 
  std::cout<< "number of dof level set " <<  ls.get_mesh_fem().nb_basic_dof()<<std::endl;
  for (size_type d = 0; d < ls.get_mesh_fem().nb_basic_dof(); ++d) {
    ls.values(0)[d] = ls_function(ls.get_mesh_fem().point_of_basic_dof(d),0, LS_TYPE)[0];
  }
  ls.touch();
  mls.adapt();
  mls.global_cut_mesh(mesh_ls);
  std::cout<<"End of ls evaluation"<<std::endl;
  // creating cut and uncut region
  std::vector<scalar_type> phicell(mesh.convex_index().size(), 0.0);
  dal::bit_vector bv_cv = mesh.convex_index();
  size_type i_cv = 0;
  for (i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv) {
    if (mls.is_convex_cut(i_cv)) {
      mesh.region(CUT_REGION).add(i_cv);
      // std::cout<<"CUT element"<<std::endl;
      getfem::mesh_fem::ind_dof_ct idofs = mf_p.ind_basic_dof_of_element(i_cv);
      phicell[i_cv] = 0.5;
    }
    else {
      mesh.region(UNCUT_REGION).add(i_cv);
      // std::cout<< "creating cut and uncut region ls is"<<
      //         ls_function(ls.get_mesh_fem().point_of_basic_dof(ls.get_mesh_fem().ind_basic_dof_of_element(i_cv)[0]),0, LS_TYPE)[0]<<std::endl;  
      if (ls_function(mf_p.point_of_basic_dof(mf_p.ind_basic_dof_of_element(i_cv)[0]),0, LS_TYPE)[0]<0) 
      { mesh.region(UNCUT_REGION_IN).add(i_cv);phicell[i_cv] = 0.4;}
      else mesh.region(UNCUT_REGION_OUT).add(i_cv);
    }
  }

  { // Just to see what elements are cut by the level set ls:
    std::cout<<"*** *** *** printin cut elements *** *** *** "<<std::endl;
    getfem::vtk_export vtke(p_des.datafilename + ".cut_elements_0.vtk");
    vtke.exporting(ls.get_mesh_fem());
    vtke.write_mesh();
    vtke.write_cell_data(phicell, "CutEl");
  }
  // base_small_vector a(2);

  gmm::resize(normal_ls_v, mf_coef_v.nb_dof()); gmm::clear(normal_ls_v);    // rhs monolithic problem

  //routine for labeling internal materials
  gen_mat();
  //routine for assignment of material propeties
   gen_coefficient();

  {
    dal::bit_vector bv_cv = mesh.convex_index();
    size_type i_cv = 0;
    for (i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv) {
      getfem::mesh_fem::ind_dof_ct idofs = mf_coef_v.ind_basic_dof_of_element(i_cv);
      for (size_type i=0; i < idofs.size(); ++i) {
        if (i==0) normal_ls_v[idofs[i]]=0;
        if (i==1) normal_ls_v[idofs[i]]=1;
        //  std::vector<scalar_type>(mf_coef.nb_dof(), 1.0)
      }
    }

  }

  compute_normal_2_ls();
  // Create vector eXt_dof.  This is the vector handling the mapping
  // from standard i  = 0, ..., nb_dof_u - 1 to extended degrees of freedom
  // ii = nb_dof_u, ..., nb_dof_x - 1 = nb_dof_u + j. The rule is that if element
  // is CUT, than each node corresponding to the standard dof ii is duplicated, and the
  // extended dof ii is created, such that eXt_dof[j] = i (here j = ii - nb_dof_u).
  //
  // for pressure 
  {
    dal::bit_vector ndofs = mf_p.basic_dof_on_region(CUT_REGION);
    for (dal::bv_visitor i(ndofs); !i.finished(); ++i)
      eXt_dof.push_back(i);
    nb_x_dof_p=eXt_dof.size();
  }
  
  std::cout << "# Cut Elements      = " << mesh.region(CUT_REGION).size() << std::endl;
  std::cout << "# dof to extend are  for pressure   = " << eXt_dof.size() << std::endl;


  // integration method on lev set
  std::cout<< "set integration method"<<std::endl;
  mim_ls_all.set_integration_method(mesh.convex_index(), ppi);
  std::cout<< "set integration method simp"<<std::endl;
  mim_ls_all.set_simplex_im(simp_ppi);
  std::cout<< "set integration adapt"<<std::endl;
  // mim_ls_all.adapt();
  std::cout<< "end"<<std::endl;

  mim_ls_in.set_integration_method(mesh.convex_index(), ppi);
  mim_ls_in.set_simplex_im(simp_ppi);
  mim_ls_in.adapt();

  mim_ls_out.set_integration_method(mesh.convex_index(), ppi);
  mim_ls_out.set_simplex_im(simp_ppi);
  mim_ls_out.adapt();

  mim_ls_bd.set_integration_method(mesh.convex_index(), ppi);
  mim_ls_bd.set_simplex_im(simp_ppi);
  mim_ls_bd.adapt();

  mfls_p.adapt();
  // mfls_u_old.adapt();mfls_p_old.adapt();

  mfls_p.set_qdim(bgeot::dim_type(1)); //number of variable
 if(0){
  getfem::mesh mesh_ls;
  mls.global_cut_mesh(mesh_ls);
  std::cout<<"laplace cut fem::save_mesh  "<<std::endl;
  getfem::vtk_export exp_cut("cutmesh.vtk");
  exp_cut.exporting(mesh_ls);
  exp_cut.write_mesh();
 }

  // init vector
  getfem::size_type nb_dof_p = mf_p.nb_dof();

  gmm::resize(B, nb_dof_p); gmm::clear(B);    // rhs monolithic problem
  gmm::resize(UP, nb_dof_p); gmm::clear(UP);  // solution monolithic
  // displacement
  gmm::resize(P, nb_dof_p); gmm::resize(P_old, nb_dof_p); 
  std::fill(P.begin(), P.end(), 0);
  gmm::copy(P,P_old);
  // iteration matrix monolithic
  gmm::resize(K, nb_dof_p,  nb_dof_p); gmm::clear(K);
}
// assembly with ls

// ===========================================
// method for generation of bcs zones
// ===========================================
void templs_problem::gen_bc(){
  std::cout << "templs_problem::gen_bc()"<< std::endl;
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);

  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    assert(i.is_face());
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);

    if (gmm::abs(un[N_-1] - 1.0) < 1.0E-7) { // new Neumann face
      mesh.region(TOP).add(i.cv(), i.f());
    } else if (gmm::abs(un[N_-1] + 1.0) < 1.0E-7) {
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
// ===========================================
// method for generation of bcs zones
// ===========================================
void templs_problem::gen_mat(){
  std::cout << "templs_problem::gen_intmat()"<< std::endl;
  dal::bit_vector bv_cv = mesh.convex_index();
  size_type i_cv = 0;
  for (i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv) {
    double zc=0.;
    for (int inode=0; inode< mesh.structure_of_convex(i_cv)->nb_points(); inode++)
      zc+=(mesh.points_of_convex(i_cv)[inode])[2] / mesh.structure_of_convex(i_cv)->nb_points();
    if (zc>1./3. && zc< 2./3.) {
      mesh.region(MAT_2).add(i_cv);
    }
    else 
      mesh.region(MAT_1).add(i_cv);
  }

  
  std::cout << "Finish templs_problem::gen_intmat()"<< std::endl;
}
// end bc generations
// =======================================================

//=======================================================
//Configure Workspace for all methods
//=======================================================
void templs_problem::configure_workspace(getfem::ga_workspace & workspace,double dt){
  std::cout << "templs_problem::configure_workspace::Configuring workspace " << std::endl;
  p_des.t_ref = p_des.l_ref * p_des.l_ref /(p_des.alpha_temp);
  p_des.p_ref = (p_des.l_ref*p_des.l_ref*p_des.q_rad)/p_des.alpha_temp;
  std::cout<< "***Non dimensional parameters***"<<std::endl;
  std::cout<< "dt "<< dt/p_des.t_ref 
    << " t_ref "<< p_des.t_ref 
    << " p_ref "<< p_des.p_ref 
    <<std::endl;

  tau_[0]= dt/p_des.t_ref ; // dt into  the workspace
  workspace.add_fixed_size_constant("tau", tau_);
  //---------------------------------------------------------
  // alpha_[0] = p_des.alpha;
  // workspace.add_fixed_size_constant("alpha", alpha_);
  //---------------------------------------------------------
  // permeability_[0] = p_des.k;
  // workspace.add_fixed_size_constant("permeability", permeability_);
  //---------------------------------------------------------
  // force_[0] = 1.e+0;
  // workspace.add_fixed_size_constant("force", force_);
  c1_[0]=p_des.alpha*p_des.alpha*p_des.biot_modulus/p_des.mu_s;
  c2_[0]=2*p_des.poisson/(1-2*p_des.poisson);
  workspace.add_fixed_size_constant("C1", c1_);
  workspace.add_fixed_size_constant("C2", c2_);
  std::cout << "end of Configuring workspace " << std::endl;
  //---------------------------------------------------------
  //  getfem::base_vector beta(1); beta[0] = 1*(alpha[0] * alpha[0] ) / (2 * p_des.mu_s + p_des.lambda_l);
  beta_[0] = (p_des.alpha* p_des.alpha ) /( (-2 * p_des.mu_s / 3 + p_des.lambda_l)*p_des.biot_modulus)/p_des.p_ref;
  workspace.add_fixed_size_constant("beta", beta_);

  penalty_[0] = 1.e+10; // 1---10
  workspace.add_fixed_size_constant("penalty", penalty_);

  workspace.add_fem_constant("Kr", mf_coef, Kr_);
  workspace.add_fem_constant("Er", mf_coef, Er_);

  workspace.add_fem_constant("nlsv", mf_coef_v, normal_ls_v);
  over_p_[0]=0.;
  if(time_iter_ > 10 && time_iter_ < 30 ) 
    over_p_[0] = 1000.*9.81*4000.*(time_iter_ - 10. )/(30.-10.);
  else if( time_iter_ < 50 && time_iter_ > 29.5 ) 
    over_p_[0] = 1000.*9.81*4000.*(50-time_iter_ )/(50-30.);
  else
    over_p_[0]=0.;



   over_p_[0]=1.;
   over_p_[0]=over_p_[0]/p_des.p_ref;
   workspace.add_fixed_size_constant("over_p", over_p_);

   over_p_[0]=10./p_des.p_ref;
   workspace.add_fixed_size_constant("top_temp",over_p_);
}
// 


//=======================================================
// Build fix stress preconditioner 
// for the monolithic problem
//=======================================================
// void templs_problem::build_fix_stress_preconditioner(double dt, double time_ls){
//   std::cout<< "templs_problem::build_fix_stress_preconditioner" <<std::endl;
//  assembly_p(dt,time_ls);assembly_u(dt,time_ls);
//  bPR_ = new biot_precond<sparse_matrix_type> (Ku,Kp,nb_x_dof_u,nb_x_dof_p);
// }
//=======================================================
// Assembly solution and iteration matrix for 
// the monolithic problem
//=======================================================
void templs_problem::assembly(double dt,double time) {
  std::cout<< "templs_problem::assembly" <<std::endl;
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
  workspace.add_expression("+1.*Test_p*tau + p_old.Test_p ", mim,UNCUT_REGION_IN);
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
  workspace.add_expression("2/element_size*p*Test_p*20", mim, TOP);	
  workspace.add_expression("-Grad_p.Normal*Test_p - Grad_Test_p.Normal*p", mim, TOP); 
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
  workspace.add_expression("2/element_size*top_temp*Test_p*20", mim, TOP);	
  workspace.assembly(1);
  workspace.clear_expressions();
  /// end boundary contions
 //  // uncut region penalization
 //  workspace.add_expression("penalty*p*Test_p *tau"	, mim_ls_out, UNCUT_REGION);
 //  workspace.assembly(2);
 //  gmm::add(workspace.assembled_matrix(), gmm::sub_matrix( K ,
 //        gmm::sub_interval(0, nb_dof_p),
 //        gmm::sub_interval(0, nb_dof_p) 
 //        )
 //      );
 //  workspace.clear_expressions();
  
  // Kout fotr enriched dof
  sparse_matrix_type K_out(nb_dof_p,nb_dof_p);
   for (int i=0; i< nb_dof_p; i++)  K_out(i,i)=1.e+0;
  std::cout<< "end kout"<< std::endl; 
  // Kin for enriched dof
  sparse_matrix_type K_in( nb_dof_p,nb_dof_p);
  // NITSCHE
  //  workspace.add_expression("2/element_size*p*Test_p*tau*2000", mim_ls_bd, CUT_REGION);// 1 is the region		
  // workspace.add_expression("-nlsv.Grad_p*Test_p*tau - nlsv.Grad_Test_p*p*tau ", mim_ls_bd, CUT_REGION); 
  //NITSCHE
  workspace.add_expression( "+p.Test_p + tau*Kr*Grad_p.Grad_Test_p"
      , mim_ls_in, CUT_REGION);
  workspace.assembly(2);
  gmm::copy(workspace.assembled_matrix(),K_in);
  workspace.clear_expressions();
  workspace.add_expression("+[+1.].Test_p*tau + p_old.Test_p ", mim_ls_in,CUT_REGION);
  //   workspace.add_expression( "2/element_size*over_p*Test_p*tau*2000", mim_ls_bd, CUT_REGION);
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

    workspace.add_expression("2*element_size*Grad_p.Normal*Grad_Test_p.Normal", mim, inner_faces);// 1 is the region		
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

//====================================================
// method for solving the laplace problem
//====================================================
void templs_problem::solve(double time){
  std::cout<<"Start solving monolithic system" <<std::endl;
  // biot_precond<sparse_matrix_type> bPR(Ku,Kp);
  size_type restart = 50;
  scalar_type cond;
  gmm::identity_matrix PM; // no precond
  gmm::iteration iter(1.e-6);  // iteration object with the max residu
  iter.set_noisy(1);               // output of iterations (2: sub-iteration)
  iter.set_maxiter(1000); // maximum number of iterations
  gmm::diagonal_precond<sparse_matrix_type> PR(K);
  {  std::cout<<"!!!!!! Printing Matrix !!!!!"<<std::endl;
    gmm::MatrixMarket_IO::write("K.mm",K);
  }
  
  getfem::size_type nb_dof_p = mf_p.nb_dof();
  sparse_matrix_type Ku_x;                                /// iteration matrix
  gmm::resize(Ku_x, nb_x_dof_p,  nb_x_dof_p);
  gmm::copy(gmm::sub_matrix(K,
  gmm::sub_interval( nb_dof_p,  nb_x_dof_p)), Ku_x); 
		sparse_matrix_type Ku_s;                                /// iteration matrix
  gmm::resize(Ku_s, nb_dof_p, nb_dof_p);
  gmm::copy(gmm::sub_matrix(K,
        gmm::sub_interval(0,nb_dof_p)), Ku_s); 
  // momentum_precond<sparse_matrix_type> mPR (Ku_s,Ku_x); 
  

  // gmm::gmres(K, UP, B, *bPR_, restart, iter);
   gmm::gmres(K, UP, B, PR, restart, iter);

  // gmm::diagonal_precond<sparse_matrix_type> PRu(Ku);
  // gmm::diagonal_precond<sparse_matrix_type> PRp(Kp);

  //   std::cout<< PR.diag << std::endl;

  // gmm::gmres(K, UP, B, *bPR_, restart, iter);
  // gmm::gmres(K, UP, B, mPR, restart, iter);
  // gmm::SuperLU_solve(K, UP , B, cond);
  std::cout << "  Condition number momentum: " << cond << std::endl;
  gmm::copy(gmm::sub_vector(UP,gmm::sub_interval(0 , nb_dof_p)), P);	

  {  
    std::vector<scalar_type> PIn(mf_p.nb_dof(), 0.0);
    for (size_type i = 0; i < eXt_dof.size(); ++i) 
    {
      size_type ii = eXt_dof[i];
      double ls_i = ls_function(mf_p.point_of_basic_dof(ii),time, LS_TYPE)[0];
      if (ls_i >= 0) PIn[ii] = UP[ mf_p.nb_dof()+ i];
    }
    for (size_type i = 0; i < mf_p.nb_dof(); ++i) 
    {
      double ls_i = ls_function(mf_p.point_of_basic_dof(i),time, LS_TYPE)[0];
      if (ls_i < 0)  PIn[i] = UP[i];
    }
    gmm::copy(PIn,P_old);
  } 
  // gmm::copy(P,P_old);  gmm::copy(U,U_old);
}



// method for solving the laplace problem
void templs_problem::print(double time,int istep,double time_ls){
  std::vector<scalar_type> PIn(mf_p.nb_dof(), 0.0);
  std::vector<scalar_type> POut(mf_p.nb_dof(), 0.0);
  std::vector<scalar_type> Pm(mf_p.nb_dof(), 0.0);
  {
    for (size_type i = 0; i < eXt_dof.size(); ++i) {
      size_type ii = eXt_dof[i];
      double ls_i = ls_function(mf_p.point_of_basic_dof(ii),time_ls, LS_TYPE)[0];
      if (ls_i < 0) {
        POut[ii] = +UP[ mf_p.nb_dof() + i];
      }
      else
        PIn[ii] = UP[ mf_p.nb_dof()+ i];
      Pm[ii] = 0.5* (UP[ii] + UP[mf_p.nb_dof()+ i]);
    }
    for (size_type i = 0; i < mf_p.nb_dof(); ++i) {
      double ls_i = ls_function(mf_p.point_of_basic_dof(i),time_ls, LS_TYPE)[0];
      if (ls_i < 0) {
        PIn[i] = UP[i];Pm[i] = UP[i];}
      else {
        POut[i] = +UP[i]; Pm[i] = UP[i];}
    }

  


    getfem::vtk_export vtke(p_des.datafilename +".pressure."+std::to_string(istep)+".vtk");
    vtke.exporting(mf_p);
    vtke.write_mesh();
    vtke.write_point_data(mf_p, PIn, "PIn");
    vtke.write_point_data(mf_p, POut, "POut");
    vtke.write_point_data(mf_p, Pm, "Pm");

  }
  std::cout<<"end print numerical pressure" <<std::endl;
  
  // 3. NUMERICAL SOLUTION
  // Using the "cut mesh" computed by the mesh_level_set object for visualization.

  getfem::mesh mcut;
  mls.global_cut_mesh(mcut);

  // Discontinuous finite element space
  getfem::mesh_fem mfcut(mcut, mf_p.get_qdim());
  mfcut.set_classical_discontinuous_finite_element(3, 0.01);

  // Interpolate the In and Out parts of the solution
  gmm::resize(Pm, mfcut.nb_dof()); gmm::clear(Pm);

  {
    std::vector<scalar_type> Paux(Pm);
    getfem::interpolation(mf_p, mfcut, PIn, Paux);
    getfem::interpolation(mf_p, mfcut, POut, Pm);

    // Check if element cv is "more In" than "Out" (based on dmean)
    // and set values af all local dofs correspondingly.
    for (dal::bv_visitor cv(mcut.convex_index()); !cv.finished(); ++cv) {
      getfem::mesh_fem::ind_dof_ct idofs = mfcut.ind_basic_dof_of_element(cv);
      scalar_type dmean = 0;
      for (size_type i=0; i < idofs.size(); ++i) {
        double ls_val=ls_function(mfcut.point_of_basic_dof( idofs[i] ) ,time_ls, LS_TYPE)[0];
        dmean += ls_val;
      }
      if (dmean < 0) 
        for (size_type i=0; i < idofs.size(); ++i) {
          size_type ii = idofs[i];
          Pm[ii] = Paux[ii];
        }
    }
  }

  std::cout<< "end pressure projection"<<std::endl; 
  bgeot::base_matrix M(N_,N_);
  bgeot::base_matrix Mm1(N_,N_);
  for (size_type i=0; i < N_; ++i) {
    M(i,i) = p_des.l_ref;
    Mm1(i,i) = 1/p_des.l_ref;
  }
  mcut.transformation(M);
  // Export discontinuous solution
  std::string namefile= p_des.datafilename +"d." +  std::to_string(istep) +".vtk";
  getfem::vtk_export vtkd(namefile);
  vtkd.exporting(mfcut);
  vtkd.write_mesh();
  vtkd.write_point_data(mfcut, Pm, "p");
  std::cout<<"end printing"<<std::endl;
  mcut.transformation(Mm1);
}


// level set function
base_small_vector templs_problem::ls_function(const base_node P, double time,int num) {
  scalar_type x = P[0]*p_des.l_ref, y = P[1]*p_des.l_ref, z=0;
  if (N_==3)  z = P[2]*p_des.l_ref;
  y = P[1]*p_des.l_ref;
  // time*=p_des.t_ref;
  base_small_vector res(2);
  switch (num) {
    case 0: {
              res[0] = y-2376.;
              res[1] = -.5 + x;
            } break;
    case 1: {
              res[0] = y - 2400.01 - 200*time/(20*1e+8);
              res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.27;
            } break;
    case 2: {
              res[0] = y - (8.e+2 * time / (1.e+8) * sin(2 * 3.14 *x/4000) + 3300);
              res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
            } break;
    case 3: {
              res[0] = z  - H_PARAM;
//               res[0] = z  -2666.6664;
//               res[0] = z  - 3050 + 500*time/(1e+8);
              res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
            } break;
    case 4: {
              res[0] = z - (8.e+2 * time / (1.e+8) * sin(2 * 3.14 *x/4000)*sin(2 * 3.14 *y/4000) + 3300);
              res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
            } break;
    case 5: {
              res[0] = (z - (0.7+0.25*((y/p_des.l_ref)*(y/p_des.l_ref)-1)
	                             *((x/p_des.l_ref)*(x/p_des.l_ref)-1)
			    )*p_des.l_ref
		       )*(40.-time/(1.e+8))/40.
	              +(z-0.844*p_des.l_ref)
		      * time/((1.e+8)*40);
              res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
            } break;
    case 6: {
              if (time/(1.e+8) > 30)
		res[0] = (z - (0.7+0.25*((y/p_des.l_ref)*(y/p_des.l_ref)-1)
				      *((x/p_des.l_ref)*(x/p_des.l_ref)-1)
			      )*p_des.l_ref
			)*(50.-time/(1.e+8))/(50.-30.)
			+(z-0.844*p_des.l_ref)
			* (time/(1.e+8)-30.)/(50.-30.);
		else
		 res[0] = (z - (0.7+0.25*((y/p_des.l_ref)*(y/p_des.l_ref)-1)
				      *((x/p_des.l_ref)*(x/p_des.l_ref)-1)
			      )*p_des.l_ref);
              res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
            } break;
    case 7: {

               if (P[2] < 0.6)
		res[0] = -((
                            (P[0]-1/2.)*3*(P[0]-1/2.)*3
                           +(P[1]-1/2.)*3*(P[1]-1/2.)*3      )/(0.5*0.5) - (P[2]*3-1)*(P[2]*3-1)/(0.5*0.5)-1);
		else
                res[0] = (-(
                              (P[0]-1/2.)*3*(P[0]-1/2.)*3 
                            + (P[1]-1/2.)*3*(P[1]-1/2.)*3
                            + ((P[2]*3-1)-0.6)*((P[2]*3-1)-0.6)
                            )+0.93);
              res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
            } break;
    default: assert(0);
  }
  return res;
}


void templs_problem::update_ls(double time, int iter){
  std::cout<< "laplacian_problem::update_ls" <<std::endl;
  //  dim_type Q = mfls_u.get_qdim();
  ls.reinit(); 
  for (size_type d = 0; d < ls.get_mesh_fem().nb_basic_dof(); ++d) {
    ls.values(0)[d] = ls_function(ls.get_mesh_fem().point_of_basic_dof(d), time ,LS_TYPE)[0];
  }
  ls.touch();
  mls.adapt();mls.global_cut_mesh(mesh_ls);
  mim_ls_in.adapt(); mim_ls_out.adapt();mim_ls_bd.adapt();  mim_ls_all.adapt();
  // clear regions todo optimize it
  mesh.region(CUT_REGION).clear();
  mesh.region(UNCUT_REGION).clear();
  mesh.region(UNCUT_REGION_IN).clear();
  mesh.region(UNCUT_REGION_OUT).clear();
  // creating cut and uncut region
  std::vector<scalar_type> phicell(mesh.convex_index().size(), 0.0);
  dal::bit_vector bv_cv = mesh.convex_index();
  size_type i_cv = 0;
  for (i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv) {
    if (mls.is_convex_cut(i_cv)) {
      mesh.region(CUT_REGION).add(i_cv);
      getfem::mesh_fem::ind_dof_ct idofs = mf_p.ind_basic_dof_of_element(i_cv);
      phicell[i_cv] = 0.5;
    }
    else {
      mesh.region(UNCUT_REGION).add(i_cv);
      if (ls_function(mf_p.point_of_basic_dof(mf_p.ind_basic_dof_of_element(i_cv)[0]),time, LS_TYPE)[0]<0) 
        {mesh.region(UNCUT_REGION_IN).add(i_cv);phicell[i_cv] = 0.1;}
      else mesh.region(UNCUT_REGION_OUT).add(i_cv);
    }
  }

  { // Just to see what elements are cut by the level set ls:
    std::cout<<"printin cut elements"<<std::endl;
    getfem::vtk_export vtke(p_des.datafilename + ".cut_elements."+std::to_string(iter)+".vtk");
    vtke.exporting(mf_p);
    vtke.write_mesh();
    vtke.write_cell_data(phicell, "CutEl");
  }

  // Create vector eXt_dof.  This is the vector handling the mapping
  // from standard i  = 0, ..., nb_dof_u - 1 to extended degrees of freedom
  // ii = nb_dof_u, ..., nb_dof_x - 1 = nb_dof_u + j. The rule is that if element
  // is CUT, than each node corresponding to the standard dof ii is duplicated, and the
  // extended dof ii is created, such that eXt_dof[j] = i (here j = ii - nb_dof_u).
  //
  // for pressure
  eXt_dof.clear();
  {
    dal::bit_vector ndofs = mf_p.basic_dof_on_region(CUT_REGION);
    for (dal::bv_visitor i(ndofs); !i.finished(); ++i)  eXt_dof.push_back(i);
    nb_x_dof_p=eXt_dof.size();
  }
  // for displacement 

  std::cout << "# Cut Elements      = " << mesh.region(CUT_REGION).size() << std::endl;
  std::cout << "# dof to extend are  for pressure   = " << eXt_dof.size() << std::endl;
  compute_normal_2_ls();
  update_p_index(time);

  // ls.touch();
  // mls.adapt();mls.global_cut_mesh(mesh_ls);
  // mim_ls_in.adapt(); mim_ls_out.adapt();mim_ls_bd.adapt();  mim_ls_all.adapt();

}
// print the solution on active part of level set
void templs_problem::print_crop(double time,int istep,double time_ls){
  std::cout<< "templs_problem::print crop" <<std::endl;


  std::vector<scalar_type> PIn(mf_p.nb_dof(), 0.0);
  std::vector<scalar_type> POut(mf_p.nb_dof(), 0.0);
  std::vector<scalar_type> Pm(mf_p.nb_dof(), 0.0);
  {
    for (size_type i = 0; i < eXt_dof.size(); ++i) {
      size_type ii = eXt_dof[i];
      double ls_i = ls_function(mf_p.point_of_basic_dof(ii),time_ls, LS_TYPE)[0];
      if (ls_i < 0) {
        POut[ii] = -100+UP[mf_p.nb_dof() + i];
      }
      else
        PIn[ii] = UP[ mf_p.nb_dof()+ i];
        Pm[ii] = 0.5* (UP[ii] + UP[mf_p.nb_dof()+ i]);
    }
    for (size_type i = 0; i < mf_p.nb_dof(); ++i) {
      double ls_i = ls_function(mf_p.point_of_basic_dof(i),time_ls, LS_TYPE)[0];
      if (ls_i < 0) {
        PIn[i] = UP[i];Pm[i] = UP[i];}
      else {
        POut[i] = -100+UP[i]; Pm[i] = UP[i];}
    }
  }
  std::cout<<"end print numerical pressure" <<std::endl;
  

  getfem::mesh_slice_cv_dof_data<std::vector<double>> a(ls.get_mesh_fem(), ls.values(0));
  getfem::slicer_isovalues a1(a,0,-1);
  getfem::mesh_slicer slicer(mls);
  getfem::mesh mesh_dim;
  int nrefine = 1;
  getfem::stored_mesh_slice sl;
  getfem::slicer_build_stored_mesh_slice a2(sl);
  getfem::slicer_build_mesh              a3(mesh_dim);
  slicer.push_back_action(a1);
  slicer.push_back_action(a2);
  slicer.push_back_action(a3);
  slicer.exec(nrefine);

  bgeot::base_matrix M(N_,N_);
  bgeot::base_matrix Mm1(N_,N_);
  for (size_type i=0; i < N_; ++i) {
    M(i,i) = p_des.l_ref;
    Mm1(i,i) = 1/p_des.l_ref;
  }
  mesh_dim.transformation(M);
  mesh.transformation(M);

  //sl.build(mesh, , -1),2);
  std::cout<<"end cropping"<< std::endl;
  {  
    //// Export discontinuous solution

    std::vector<scalar_type> PIn_dim(mf_p.nb_dof(), 0.0);gmm::copy(PIn,PIn_dim);gmm::scale(PIn_dim,p_des.p_ref);
    std::string namefile= p_des.datafilename +".crop." +  std::to_string(istep) +".vtk";
    getfem::vtk_export vtkd(namefile);
    vtkd.exporting(mesh_dim);
    vtkd.write_mesh();
    vtkd.write_point_data(mf_p, PIn_dim, "p");
    std::cout<<"end printing"<<std::endl;

  }

  mesh.transformation(Mm1); 
}  

// evaluate birnak to a level set
void templs_problem::compute_normal_2_ls(){
  dal::bit_vector bv_cv = mesh.convex_index();
  size_type i_cv = 0;
  getfem::base_matrix gradU(1, N_);
  getfem::base_matrix interp;
  bgeot::base_vector ls_val;
  for (i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv) {
    ls_val.resize(ls.get_mesh_fem().nb_basic_dof_of_element(i_cv));
    base_tensor t;
    const base_node x= mf_coef_v.point_of_basic_dof(i_cv);
    bgeot::pgeometric_trans pgt = 
      bgeot::geometric_trans_descriptor(p_des.MESH_TYPE);; 
    bgeot::vectors_to_base_matrix
      (interp, ls.get_mesh_fem().linked_mesh().points_of_convex(i_cv));
    getfem::fem_interpolation_context icx(pgt,
        ls.get_mesh_fem().fem_of_element(i_cv),
        x,
        interp ,
        i_cv
        );

    getfem::mesh_fem::ind_dof_ct idofs_ls = ls.get_mesh_fem().ind_basic_dof_of_element(i_cv);
    for (size_type i=0; i < idofs_ls.size(); ++i)ls_val[i] = ls.values(0)[i];
    icx.pf()->interpolation_grad(icx, ls_val, gradU, 1);
    scalar_type norm = gmm::vect_norm2(gmm::mat_row(gradU, 0));
    getfem::mesh_fem::ind_dof_ct idofs = mf_coef_v.ind_basic_dof_of_element(i_cv);

    for (size_type i=0; i < idofs.size(); ++i) {
      normal_ls_v[idofs[i]]=gradU(0, i) / norm;
    }
  }
}

// =================================================
// update p index 
void templs_problem::update_p_index(double time_ls){
  // std::map<int,int[2]> inx_p_map;
  gmm::resize(pin_index_,2*eXt_dof.size()*eXt_dof.size());
  gmm::resize(pout_index_,2*eXt_dof.size()*eXt_dof.size());
  size_type dof_shift = mf_p.nb_dof();
  int j_shift=eXt_dof.size()*eXt_dof.size();

  for (size_type i = 0; i < eXt_dof.size(); ++i) {
    size_type ii = eXt_dof[i];
    double ls_i = ls_function(mf_p.point_of_basic_dof(ii),time_ls, LS_TYPE)[0];
    for (size_type j = 0; j < eXt_dof.size(); ++j) {
      size_type jj = eXt_dof[j];
      double ls_j = ls_function(mf_p.point_of_basic_dof(jj),time_ls, LS_TYPE)[0];
      // std::cout<< "ls_i "<< ls_i << " ls_j "<< ls_j << std::endl;
      if ( (ls_i <= 0) && (ls_j <= 0) ) {
        // i and j are both In
        pin_index_[i*eXt_dof.size()+j]=ii;
        pin_index_[i*eXt_dof.size()+j+j_shift]=jj;
        // Kp(i + dof_shift, j  + dof_shift) += K_out(ii, jj);
        pout_index_[i*eXt_dof.size()+j]=i + dof_shift;
        pout_index_[i*eXt_dof.size()+j+j_shift]=j  + dof_shift;
      }
      else if ( (ls_i>= 0) && (ls_j >= 0) ) {
        // i and j are both Out
        // Kp(ii, jj) += K_out(ii, jj);
        pin_index_[i*eXt_dof.size()+j]=i + dof_shift;
        pin_index_[i*eXt_dof.size()+j+j_shift]=j  + dof_shift;
        // out
        pout_index_[i*eXt_dof.size()+j]=ii;
        pout_index_[i*eXt_dof.size()+j+j_shift]=jj;
      }
      else if ( (ls_i <= 0) && (ls_j>= 0) ) {
        // i is In, j is Out
        pin_index_[i*eXt_dof.size()+j]=ii;
        pin_index_[i*eXt_dof.size()+j+j_shift]=j+ dof_shift;
        // Kp(i + dof_shift, jj) += K_out(ii, jj);
        pout_index_[i*eXt_dof.size()+j]=i + dof_shift;
        pout_index_[i*eXt_dof.size()+j+j_shift]=jj;
      }
      else {
        // i is Out, j is In
        pin_index_[i*eXt_dof.size()+j]=i + dof_shift;
        pin_index_[i*eXt_dof.size()+j+j_shift]=jj;
        // Kp(ii, j  + dof_shift) += K_out(ii, jj);
        pout_index_[i*eXt_dof.size()+j]=ii;
        pout_index_[i*eXt_dof.size()+j+j_shift]=j+dof_shift;
      }
    }
  }

}


//============================================
// routine for the generation of coeffient
//============================================
void templs_problem::gen_coefficient(){ // creating a coefficient
  
  std::vector<scalar_type> Kr_print; // permeability ratio
  gmm::resize(Kr_print, mf_coef.nb_dof()); gmm::fill(Kr_print,1);    // rhs monolithic problem
  gmm::resize(Kr_, mf_coef.nb_dof()); gmm::fill(Kr_,1);    // rhs monolithic problem
  gmm::resize(Er_, mf_coef.nb_dof()); gmm::fill(Er_,1);    // rhs monolithic problem
  std::vector<int> material; material.push_back(MAT_1);material.push_back(MAT_2);

  //std::vector<double> k; k.push_back(1);k.push_back(1.e+2);
  //std::vector<double> E; E.push_back(1);E.push_back(2.e+0);
  std::vector<double> k; k.push_back(1);k.push_back(1.e+0);
  std::vector<double> E; E.push_back(1);E.push_back(1.e+0);
  for (int imat=0; imat< material.size();imat++){
    dal::bit_vector bv_cv = mesh.region(material[imat]).index();
    size_type i_cv = 0;
    for (i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv) {
      getfem::mesh_fem::ind_dof_ct idofs = mf_coef.ind_basic_dof_of_element(i_cv);
      for (size_type i=0; i < idofs.size(); ++i) {
        Kr_[idofs[i]]=k[imat]; Kr_print[(int) i_cv]=k[imat];
        Er_[idofs[i]]=E[imat];
      }
    }
  }
  if(1){ // Just to see what elements are cut by the level set ls:
    std::string namefile= p_des.datafilename +".materials.vtk";
    getfem::vtk_export vtk_data(namefile);
    vtk_data.exporting(mf_coef);
    vtk_data.write_mesh();
    vtk_data.write_cell_data(Kr_print, "K");
  }
}

//============================================
// routine for printing spasity pattern
//============================================
void templs_problem::print_pattern(int iter){
  std::ifstream myfile ("fort." +  std::to_string(iter+1));
  if (myfile.good()){
    std::cout<<"Reading "<< "fort." +  std::to_string(iter+1) <<" for pattern generation"<<std::endl;
    std::string line;
    std::vector<scalar_type> color;
    // 	Read cvs
    if (myfile.is_open())
    {
      while ( getline (myfile,line) )
      {
	std::istringstream iss(line);
	double a;
	iss>>a;color.push_back(a);
      }
      myfile.close();
    }
    
    bgeot::base_matrix M(N_,N_);
    bgeot::base_matrix Mm1(N_,N_);
    for (size_type i=0; i < N_; ++i) {
      M(i,i) = p_des.l_ref;
      Mm1(i,i) = 1/p_des.l_ref;
    }
    //   mesh_dim.transformation(M);
    mesh.transformation(M);
    
    //sl.build(mesh, , -1),2);
    std::cout<<"end cropping"<< std::endl;
    {  
      //// Export discontinuous solution
      std::cout<<"Size of dof"<< color.size()<< std::endl;
      std::string namefile= p_des.datafilename +".pattern." +  std::to_string(iter) +".vtk";
      getfem::vtk_export vtkd(namefile);
      vtkd.exporting(mf_p);
      vtkd.write_mesh();
      vtkd.write_point_data(mf_p, color, "c");
      std::cout<<"end printing pattern in "<< namefile <<std::endl;
      
    }
    
    mesh.transformation(Mm1); 
  }//end if file exist
  else {std::cout<<
    "skipping "<< "fort." +  std::to_string(iter+1)
    <<" for pattern generation"<<std::endl;}
  }

void templs_problem::print_ls(double time,int istep,double time_ls){
  
  std::cout<<"Start printing ls function"<< std::endl;
  std::vector<scalar_type> ls_value_print(mf_p.nb_dof(), 0.0);
  {
    for (size_type i = 0; i < mf_p.nb_dof(); ++i)
    {
      ls_value_print[i]  = ls_function(mf_p.point_of_basic_dof(i),time_ls, LS_TYPE)[0];
    }
  }
  
  bgeot::base_matrix M(N_,N_);
  bgeot::base_matrix Mm1(N_,N_);
  for (size_type i=0; i < N_; ++i) {
    M(i,i) = p_des.l_ref;
    Mm1(i,i) = 1/p_des.l_ref;
  }
//   mesh_dim.transformation(M);
  mesh.transformation(M);
  std::vector<scalar_type> over_p; // permeability ratio
  gmm::resize(over_p, mf_coef.nb_dof()); gmm::fill(over_p,over_p_[0]);    // rhs monolithic problem
  std::string namefile= p_des.datafilename +".ls." +  std::to_string(istep) +".vtk";
  getfem::vtk_export vtkd(namefile);
  vtkd.exporting(mf_p);vtkd.write_mesh();
  vtkd.write_point_data(mf_p, ls_value_print, "ls");
  vtkd.write_cell_data(over_p, "h_ice");
  std::cout<<"end printing ls function"<< std::endl;

  mesh.transformation(Mm1); 
}
