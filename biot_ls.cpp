#include "biot_ls.hpp" 

void biotls_problem::init(void) {

    std::cout<< "biotls_problem::init "  << std::endl;

	bgeot::pgeometric_trans pgt = 
		bgeot::geometric_trans_descriptor(p_des.MESH_TYPE);

	// Mesh generation
	N_ = pgt->dim();
    std::vector<size_type> nsubdiv(N_);
	int NX=p_des.nsubdiv;
	std::fill(nsubdiv.begin(),nsubdiv.end(),NX);
	getfem::regular_unit_mesh(mesh, nsubdiv, pgt, p_des.noised);
    // A trasformation for the squarred mesh
	  bgeot::base_matrix M(N_,N_);
	  for (size_type i=0; i < N_; ++i) {
	    M(i,i) = 1.0;
	  }
      p_des.l_ref=4000;
	//  if (N>1) { M(0,1) = 0; }
	//
	  mesh.transformation(M);
    // // End of mesh generation

   
    // set  integration methods  
    getfem::pfem pf_u = getfem::fem_descriptor(p_des.FEM_TYPE_U); // for displacement 
    getfem::pfem pf_p = getfem::fem_descriptor(p_des.FEM_TYPE_P); // for pressure
	getfem::pintegration_method ppi = getfem::int_method_descriptor(p_des.INTEGRATION);
    getfem::pintegration_method simp_ppi = getfem::int_method_descriptor(p_des.SIMPLEX_INTEGRATION);
    
	mim.set_integration_method(mesh.convex_index(), ppi);
    mf_u.set_qdim(N_);
	mf_u.set_finite_element(mesh.convex_index(), pf_u); // finite element for displacement
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
         for (size_type i=0; i < idofs.size(); ++i) {
             //  std::cout<<"\t dofs "<< i<<" "<< idofs[i]<<std::endl;
              // std::cout<< "x" << mf_p.point_of_basic_dof(idofs[i])[0] 
              // << "y" << mf_p.point_of_basic_dof(idofs[i])[1] 
             //  << "z" << mf_p.point_of_basic_dof(idofs[i])[2] 
              // <<std::endl;
              }
        }
      else {
          mesh.region(UNCUT_REGION).add(i_cv);
          // std::cout<< "creating cut and uncut region ls is"<<
            //         ls_function(ls.get_mesh_fem().point_of_basic_dof(ls.get_mesh_fem().ind_basic_dof_of_element(i_cv)[0]),0, LS_TYPE)[0]<<std::endl;  
          if (ls_function(mf_p.point_of_basic_dof(mf_p.ind_basic_dof_of_element(i_cv)[0]),0, LS_TYPE)[0]<0) 
              mesh.region(UNCUT_REGION_IN).add(i_cv);
          else mesh.region(UNCUT_REGION_OUT).add(i_cv);
      }
    }
   
    { // Just to see what elements are cut by the level set ls:
    getfem::vtk_export vtke("cut_elements.vtk");
    vtke.exporting(ls.get_mesh_fem());
    vtke.write_mesh();
    vtke.write_cell_data(phicell, "CutEl");
    }
    base_small_vector a(2);
    gmm::resize(normal_ls, mf_coef.nb_dof()); gmm::clear(normal_ls);    // rhs monolithic problem
    
    gmm::resize(normal_ls_v, mf_coef_v.nb_dof()); gmm::clear(normal_ls_v);    // rhs monolithic problem
    { // creating a coefficient
    dal::bit_vector bv_cv = mesh.convex_index();
    size_type i_cv = 0;
    for (i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv) {
         getfem::mesh_fem::ind_dof_ct idofs = mf_coef.ind_basic_dof_of_element(i_cv);
         for (size_type i=0; i < idofs.size(); ++i) {
         normal_ls[idofs[i]]=1.e-9;
       //  std::vector<scalar_type>(mf_coef.nb_dof(), 1.0)
      }
   }
   }
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
   
    getfem::vtk_export vtke("normalls.vtk");
    vtke.exporting(mf_coef_v);
    vtke.write_mesh();
    vtke.write_cell_data(normal_ls_v, "nls");
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
    // for displacement 
    {
       dal::bit_vector ndofs = mf_u.basic_dof_on_region(CUT_REGION);
       for (dal::bv_visitor i(ndofs); !i.finished(); ++i)
        eXt_dof_u.push_back(i);
       nb_x_dof_u=eXt_dof_u.size();
     }
    std::cout << "# Cut Elements      = " << mesh.region(CUT_REGION).size() << std::endl;
    std::cout << "# dof to extend are  for pressure   = " << eXt_dof.size() << std::endl;
    std::cout << "# dof to extend are  for displacement   = " << eXt_dof_u.size() << std::endl;
    
    
    // integration method on lev set
    std::cout<< "set integration method"<<std::endl;
    mim_ls_all.set_integration_method(mesh.convex_index(), ppi);
    std::cout<< "set integration method simp"<<std::endl;
    mim_ls_all.set_simplex_im(simp_ppi);
    std::cout<< "set integration adapt"<<std::endl;
    mim_ls_all.adapt();
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
    
    mfls.adapt();mfls_u.adapt();mfls_p.adapt();
    // mfls_u_old.adapt();mfls_p_old.adapt();
    
	mf_u.set_qdim(bgeot::dim_type(N_)); //number of variable
    mfls_u.set_qdim(bgeot::dim_type(N_)); //number of variable
    mfls_p.set_qdim(bgeot::dim_type(1)); //number of variable
    // mfls_u_old.set_qdim(bgeot::dim_type(2)); //number of variable
    // mfls_p_old.set_qdim(bgeot::dim_type(1)); //number of variable
    
    getfem::mesh mesh_ls;
    mls.global_cut_mesh(mesh_ls);
    std::cout<<"laplace cut fem::save_mesh  "<<std::endl;
    getfem::vtk_export exp_cut("cutmesh.vtk");
    exp_cut.exporting(mesh_ls);
    exp_cut.write_mesh();
    
    
    // init vector
    getfem::size_type nb_dof_u = mf_u.nb_dof();
    getfem::size_type nb_dof_p = mf_p.nb_dof();
    
    gmm::resize(B, nb_dof_u + nb_dof_p); gmm::clear(B);    // rhs monolithic problem
    gmm::resize(Bp, nb_dof_p); gmm::resize(Bu, nb_dof_u);  // rhs fixed stress pressure displacement
    gmm::resize(UP, nb_dof_u + nb_dof_p); gmm::clear(UP);  // solution monolithic
    // displacement
    gmm::resize(U, nb_dof_u); gmm::resize(U_old, nb_dof_u); gmm::resize(U_iter, nb_dof_u); 
    gmm::clear(U);gmm::clear(U_old);gmm::clear(U_iter);
    // pressure
    gmm::resize(P, nb_dof_p); gmm::resize(P_old, nb_dof_p); gmm::resize(P_iter, nb_dof_p); 
    std::fill(P.begin(), P.end(), 0);
    gmm::copy(P,P_old);gmm::copy(P,P_iter);
    // iteration matrix monolithic
    gmm::resize(K, nb_dof_u + nb_dof_p, nb_dof_u + nb_dof_p); gmm::clear(K);
    // iteration matrix displacement fixed stress
    gmm::resize(Ku, nb_dof_u , nb_dof_u ); gmm::clear(Ku);
    // iteration matrix pressure fixed stress
    gmm::resize(Kp, nb_dof_p , nb_dof_p ); gmm::clear(Kp);
}
// assembly with ls

// ===========================================
// method for generation of bcs zones
// ===========================================
void biotls_problem::gen_bc(){
	std::cout << "biotls_problem::gen_bc()"<< std::endl;
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

//=======================================================
//Configure Workspace for all methods
//=======================================================
void biotls_problem::configure_workspace(getfem::ga_workspace & workspace,double dt){
    std::cout << "biotls_problem::configure_workspace::Configuring workspace " << std::endl;
    p_des.t_ref = p_des.l_ref/(p_des.k*p_des.biot_modulus);
    p_des.u_ref = (p_des.rho_r-p_des.rho_l)*9.81*p_des.l_ref*p_des.l_ref/p_des.mu_s;
    p_des.p_ref = p_des.alpha*p_des.u_ref*p_des.biot_modulus/p_des.l_ref;
    std::cout<< "***Non dimensional parameters***"<<std::endl;
    std::cout<< "dt "<< dt/p_des.t_ref 
             << " t_ref "<< p_des.t_ref 
             << " u_ref "<< p_des.u_ref 
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
    
    workspace.add_fem_constant("nls", mf_coef, normal_ls);
    
    workspace.add_fem_constant("nlsv", mf_coef_v, normal_ls_v);
    }
// 


//=======================================================
// Build fix stress preconditioner 
// for the monolithic problem
//=======================================================
void biotls_problem::build_fix_stress_preconditioner(){
    std::cout<< "biotls_problem::build_fix_stress_preconditioner" <<std::endl;
    bPR_ = new biot_precond<sparse_matrix_type> (Ku,Kp);
    }
//=======================================================
// Assembly solution and iteration matrix for 
// the monolithic problem
//=======================================================
void biotls_problem::assembly(double dt,double time) {
    std::cout<< "biotls_problem::assembly" <<std::endl;
    // assembly matrix
    // resize and clear vector and matrix
	getfem::size_type nb_dof_u = mf_u.nb_dof();
    getfem::size_type nb_dof_p = mf_p.nb_dof();
    getfem::size_type nb_dof_p_x = nb_dof_p+ nb_x_dof_p;
    getfem::size_type nb_dof_u_x = nb_dof_u+ nb_x_dof_u;
    gmm::resize(B, nb_dof_u_x + nb_dof_p_x ); gmm::clear(B);
    gmm::resize(UP, nb_dof_u_x + nb_dof_p_x); //gmm::clear(UP);
    gmm::resize(U, nb_dof_u); gmm::resize(U_old, nb_dof_u); 
    gmm::resize(P, nb_dof_p); gmm::resize(P_old, nb_dof_p);
    gmm::resize(K, nb_dof_u_x + nb_dof_p_x, nb_dof_u_x + nb_dof_p_x); gmm::clear(K);
 	std::cout<<"nb_dof_u "<< nb_dof_u<<" nb_dof_p "<< nb_dof_p<<std::endl;
	getfem::ga_workspace workspace; // generic workspace for the solution
    configure_workspace(workspace,dt); 
    // getfem::base_vector invdt(1); invdt[0] = 1/dt;
	// workspace.add_fixed_size_constant("invdt", invdt);
    
   	workspace.add_fem_variable("u", mf_u, gmm::sub_interval(0, nb_dof_u), U);
	workspace.add_fem_variable("p", mf_p, gmm::sub_interval(nb_dof_u, nb_dof_p), P);
	workspace.add_fem_constant("u_old", mf_u, U_old);
	workspace.add_fem_constant("p_old", mf_p, P_old);
  	// ------------------ expressions --------------------------
	workspace.add_expression( "2*Sym(Grad_u):Grad_Test_u + C2*Div_u*Div_Test_u"
                              "- C1*p.Div_Test_u ", mim_ls_in, UNCUT_REGION);
    workspace.add_expression( "p.Test_p + tau*Grad_p.Grad_Test_p"
                              "+Test_p*Div_u", mim_ls_in, UNCUT_REGION);
    workspace.assembly(2);
    gmm::copy(workspace.assembled_matrix(), gmm::sub_matrix( K ,
                                                             gmm::sub_interval(0, nb_dof_u+ nb_dof_p),
                                                             gmm::sub_interval(0, nb_dof_u+ nb_dof_p) 
                                                             )
            );
	workspace.clear_expressions();
   
	//======= RHS =====================
    if(N_== 2 )	workspace.add_expression("[0,-1].Test_u", mim_ls_in);
    else        workspace.add_expression("[0,0,-0].Test_u", mim_ls_in);
     workspace.add_expression("+0.0e-6*Test_p*tau + p_old.Test_p + Test_p*Div_u_old", mim_ls_in,UNCUT_REGION_IN);
     // workspace.add_expression("nls*Test_p*tau ", mim_ls_in,UNCUT_REGION_IN);
    workspace.set_assembled_vector(B);
    workspace.assembly(1);
	workspace.clear_expressions();
	//Boudanry conditions // NICHE
	getfem::base_vector penalty(1); penalty[0] = 2e+18; // 1---10
	workspace.add_fixed_size_constant("penalty", penalty);
	//Matrix term for nietche bc
	workspace.add_expression("penalty/element_size*u.Test_u", mim, BOTTOM);
	workspace.add_expression("penalty/element_size*p*Test_p", mim, LEFT);
    workspace.add_expression("penalty/element_size*p*Test_p", mim, RIGHT);// 1 is the region		
	workspace.add_expression("-Grad_p.Normal*Test_p - Grad_Test_p.Normal*p ", mim, LEFT); 	
	workspace.add_expression("-Grad_p.Normal*Test_p - Grad_Test_p.Normal*p", mim, RIGHT); 
	workspace.add_expression("penalty/element_size*p*Test_p", mim, TOP);	
	workspace.add_expression("-Grad_p.Normal*Test_p - Grad_Test_p.Normal*p", mim, TOP); 
	workspace.assembly(2);
    gmm::add(workspace.assembled_matrix(), gmm::sub_matrix( K ,
                                                             gmm::sub_interval(0, nb_dof_u+ nb_dof_p),
                                                             gmm::sub_interval(0, nb_dof_u+ nb_dof_p) 
                                                             )
            );
	workspace.clear_expressions();
    // rhs for niche boudary condition
	workspace.add_expression(" 0*penalty/element_size*Test_p -Grad_Test_p.Normal*0 ", mim_ls_in, LEFT);
	workspace.add_expression(" 0*penalty/element_size*Test_p -Grad_Test_p.Normal*0 ", mim_ls_in, RIGHT);
	workspace.assembly(1);
	workspace.clear_expressions();
    /// end boundary contions
    // uncut region penalization
    workspace.add_expression("penalty*p*Test_p *tau"	, mim_ls_out, UNCUT_REGION);
    workspace.add_expression("penalty*u.Test_u "	, mim_ls_out, UNCUT_REGION);
    workspace.assembly(2);
    gmm::add(workspace.assembled_matrix(), gmm::sub_matrix( K ,
                                                             gmm::sub_interval(0, nb_dof_u+ nb_dof_p),
                                                             gmm::sub_interval(0, nb_dof_u+ nb_dof_p) 
                                                             )
            );
    workspace.clear_expressions();
    workspace.assembly(1); 
    workspace.clear_expressions();
   // Kout fotr enriched dof
   sparse_matrix_type K_out(nb_dof_u + nb_dof_p,nb_dof_u + nb_dof_p);
   workspace.add_expression("penalty*u.Test_u "	,   mim_ls_out, CUT_REGION);
   workspace.add_expression("penalty*u.Test_u "	,   mim_ls_bd, CUT_REGION);
   workspace.add_expression("penalty*p*Test_p*tau ", mim_ls_out, CUT_REGION);
   workspace.assembly(2);
   gmm::copy(workspace.assembled_matrix(),K_out);
   workspace.clear_expressions();
   workspace.assembly(1); 
   workspace.clear_expressions();
   std::cout<< "end kout"<< std::endl; 
   // Kin for enriched dof
   sparse_matrix_type K_in(nb_dof_u + nb_dof_p,nb_dof_u + nb_dof_p);
   // NICHE
     workspace.add_expression("2/element_size*p*Test_p", mim_ls_bd, CUT_REGION);// 1 is the region		
     workspace.add_expression("-nlsv.Grad_p*Test_p*tau - nlsv.Grad_Test_p*p*tau ", mim_ls_bd, CUT_REGION); 
   //NICHE
   // workspace.add_expression( "permeability*tau*[0,1].Grad_p*Test_p ", mim_ls_bd, CUT_REGION);
    workspace.add_expression( "+p.Test_p + tau*Grad_p.Grad_Test_p"
                             "+Test_p*Div_u", mim_ls_in, CUT_REGION);
       // internal for displacement
   workspace.add_expression( "2*Sym(Grad_u):Grad_Test_u + C2*Div_u*Div_Test_u"
                              "- C1*p.Div_Test_u ", mim_ls_in, CUT_REGION);
   workspace.assembly(2);
   gmm::copy(workspace.assembled_matrix(),K_in);
   workspace.clear_expressions();
   workspace.add_expression("+[+0.0e-6].Test_p*tau + p_old.Test_p + Test_p*Div_u_old", mim_ls_in,CUT_REGION);
   workspace.add_expression("[0,-1].Test_u", mim_ls_in,CUT_REGION);
   workspace.assembly(1);
   workspace.clear_expressions();
   workspace.add_expression("penalty*p*Test_p*tau ", mim_ls_in, LEFT);
   workspace.add_expression("penalty*p*Test_p*tau ", mim_ls_in, RIGHT);
   workspace.assembly(2);
   gmm::add(workspace.assembled_matrix(),K_in);
   std::cout<< "end kin"<< std::endl; 
   
   
   {// mapping for pressure
    std::cout<<"start mapping pressure"<<std::endl;
    size_type dof_shift = nb_dof_u+nb_dof_p;
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
      B[dof_shift + i] = B[nb_dof_u+ ii];
      for (size_type j = 0; j < eXt_dof.size(); ++j) {
	    size_type jj = eXt_dof[j];
       double ls_j = ls_function(mf_p.point_of_basic_dof(jj),time, LS_TYPE)[0];
            // std::cout<< "ls_i "<< ls_i << " ls_j "<< ls_j << std::endl;
	        if ( (ls_i <= 0) && (ls_j <= 0) ) {
	         // i and j are both In
	         K(ii + nb_dof_u, jj+ nb_dof_u) += K_in(ii+ nb_dof_u, jj+ nb_dof_u);
	         K(i + dof_shift, j  + dof_shift) += K_out(ii+ nb_dof_u, jj+ nb_dof_u);
            }
            else if ( (ls_i>= 0) && (ls_j >= 0) ) {
	         // i and j are both Out
	         K(ii+ nb_dof_u, jj+ nb_dof_u) += K_out(ii+ nb_dof_u, jj+ nb_dof_u);
	         K(i + dof_shift, j  + dof_shift) += K_in(ii+ nb_dof_u, jj+ nb_dof_u);
	        }
	        else if ( (ls_i <= 0) && (ls_j>= 0) ) {
             // i is In, j is Out
	         K(ii+ nb_dof_u, j  + dof_shift) += K_in(ii+ nb_dof_u, jj+ nb_dof_u);
	         K(i + dof_shift, jj+ nb_dof_u) += K_out(ii+ nb_dof_u, jj+ nb_dof_u);
	        }
	        else {
	        // i is Out, j is In
	        K(i + dof_shift, jj+ nb_dof_u) += K_in(ii+ nb_dof_u, jj+ nb_dof_u);
	        K(ii+ nb_dof_u, j  + dof_shift) += K_out(ii+ nb_dof_u, jj+ nb_dof_u);
	        }
        }
    }
   }
   {// mapping for displacement
    std::cout<<"start mapping displacement"<<std::endl;
    size_type dof_shift = nb_dof_u+ nb_dof_p_x ;
    for (size_type i = 0; i < eXt_dof_u.size(); ++i) {
      size_type ii = eXt_dof_u[i];
      double ls_i = ls_function(mf_u.point_of_basic_dof(ii),time, LS_TYPE)[0];
      // std::cout<< "ls values of dof"<< ii<< " is "  << ls_i<< std::endl;
       B[dof_shift + i] = B[ii];
      for (size_type j = 0; j < eXt_dof_u.size(); ++j) {
	    size_type jj = eXt_dof_u[j];
       double ls_j = ls_function(mf_u.point_of_basic_dof(jj),time, LS_TYPE)[0];
        // std::cout<< "ls_i "<< ls_i << " ls_j "<< ls_j << std::endl;
	        if ( (ls_i <= 0) && (ls_j <= 0) ) {
	         // i and j are both In
	         K(ii , jj ) += K_in(ii, jj);
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
   // std::cout<< K<<std::endl;
   
   // for (int i=0; i< nb_dof_u + nb_dof_p; i++)  K(i,i)=1.e+4;
   
} // end assembly




// -------------------------------------------------
// Assembling pressure matrix for fixed stress
// -------------------------------------------------

void biotls_problem::assembly_p(double dt, double time){
    
    // gmm::clean(P_old, 1E-10);gmm::clean(U_old, 1E-10);
	// gmm::clean(P, 1E-10);gmm::clean(U, 1E-10);
    gmm::clear(Bp); gmm::clear(Kp);//  gmm::clear(P_old);
	std::cout<<"biot_assembler::assembly_p(double dt)" <<std::endl;
   	getfem::size_type nb_dof_u = mf_u.nb_dof();
    getfem::size_type nb_dof_p = mf_p.nb_dof();
    getfem::size_type nb_dof_p_x = nb_dof_p + nb_x_dof_p;
    getfem::size_type nb_dof_u_x = nb_dof_u + nb_x_dof_u;
    std::cout<< "total dofs " <<   mf_p.nb_dof() << "normal dof" << mf_p.nb_dof()  <<std::endl;
    
    gmm::resize(P, nb_dof_p_x); 
    // gmm::resize(P_old, nb_dof_p);gmm::resize(P_iter, nb_dof_p);
    gmm::resize(Kp, nb_dof_p_x, nb_dof_p_x);  gmm::resize(Bp, nb_dof_p_x);
	
    
    getfem::ga_workspace workspace;
    configure_workspace(workspace,dt);


	workspace.add_fem_variable("p", mf_p, gmm::sub_interval(0, nb_dof_p), P);
	workspace.add_fem_variable("p_old", mf_p, gmm::sub_interval(0,nb_dof_p), P_old);
	workspace.add_fem_variable("p_iter", mf_p, gmm::sub_interval(0,nb_dof_p), P_iter);
	// workspace.add_fem_variable("u", mf_u, gmm::sub_interval(0, nb_dof_u), U);
	workspace.add_fem_variable("u_old", mf_u, gmm::sub_interval(0, nb_dof_u), U_old);
	workspace.add_fem_variable("u_iter", mf_u, gmm::sub_interval(0, nb_dof_u), U_iter);
	
	// Pressure equation
	workspace.add_expression("(1+beta)*p.Test_p + tau*Grad_p.Grad_Test_p", mim,UNCUT_REGION_IN); // tau
	// workspace.set_assembled_matrix(Kp);
	workspace.assembly(2);
    gmm::copy(workspace.assembled_matrix(), gmm::sub_matrix(Kp,
                                                             gmm::sub_interval(0, nb_dof_p),
                                                             gmm::sub_interval(0, nb_dof_p) 
                                                             ));
	workspace.clear_expressions();
    
	//======= RHS =====================
	//  workspace.add_expression("[+0.0].Test_p + (1/bm)*invdt*p_old.Test_p + invdt*alpha*Test_p*Trace((Grad_u_old))", mim);// 1/dt
	//  workspace.add_expression("(beta)*invdt*p_iter.Test_p - invdt*alpha*Test_p*Trace((Grad_u_iter))", mim);// 1/dt
	// tau
	workspace.add_expression("[+0]*Test_p*tau + p_old.Test_p + Test_p*Trace(Sym(Grad_u_old))", mim,UNCUT_REGION_IN); // tau
	workspace.add_expression("(beta)*p_iter.Test_p - Test_p*Trace(Sym(Grad_u_iter))", mim,UNCUT_REGION_IN); // tau
	workspace.set_assembled_vector(Bp);
	workspace.assembly(1);
	workspace.clear_expressions();

	//Matrix term
	workspace.add_expression("2*penalty/element_size*p*Test_p", mim, TOP);
	workspace.add_expression("-Grad_p.Normal*Test_p - Grad_Test_p.Normal*p ", mim, TOP); 	
    workspace.add_expression("2*penalty/element_size*p*Test_p", mim, LEFT);
	workspace.add_expression("-Grad_p.Normal*Test_p*tau- Grad_Test_p.Normal*p*tau ", mim, LEFT); 	
   	workspace.add_expression("2*penalty/element_size*p*Test_p", mim, RIGHT);
	workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, RIGHT); 	
    if (N_= 3 ){
    workspace.add_expression("2*penalty/element_size*p*Test_p", mim, LEFTX);
	workspace.add_expression("-Grad_p.Normal*Test_p*tau- Grad_Test_p.Normal*p*tau ", mim, LEFTX); 	
   	workspace.add_expression("2*penalty/element_size*p*Test_p", mim, RIGHTX);
	workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, RIGHTX); 	
        }
    
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
     workspace.add_expression("penalty*penalty*p*Test_p ", mim, UNCUT_REGION_OUT);
    //workspace.add_expression("(1/bm + beta)*p.Test_p + tau*permeability*Grad_p.Grad_Test_p"	, mim_ls_out,UNCUT_REGION_OUT);
    //// workspace.add_expression("1.e+18*p*Test_p "	, mim_ls_bd);
    workspace.assembly(2);
    gmm::add(workspace.assembled_matrix(), gmm::sub_matrix(Kp,
                                                            gmm::sub_interval(0, nb_dof_p),
                                                            gmm::sub_interval(0, nb_dof_p) 
                                                            ));
	workspace.clear_expressions();
    
    
   sparse_matrix_type K_out(nb_dof_p,nb_dof_p);
   {
      workspace.add_expression("penalty*penalty*p*Test_p ", mim_ls_out, CUT_REGION);
      
  //  workspace.add_expression( "(1+beta)*p.Test_p + tau*Grad_p.Grad_Test_p", mim, CUT_REGION);
   // workspace.add_expression("(1/bm + beta)*p.Test_p + tau*permeability*Grad_p.Grad_Test_p"	, mim, CUT_REGION);
    // workspace.add_expression("2/element_size*p*Test_p*tau", mim_ls_bd, CUT_REGION);// 1 is the region		
  //   workspace.add_expression("-permeability*nlsv.Grad_p*Test_p *tau- permeability*nlsv.Grad_Test_p*p*tau ", mim_ls_bd, CUT_REGION); 
   workspace.assembly(2);
   gmm::copy(workspace.assembled_matrix(),K_out);
   workspace.clear_expressions();
   std::cout<< "end kout"<< std::endl; 
   }
   // Kin for enriched dof
   sparse_matrix_type K_in(nb_dof_p,nb_dof_p);
   std::vector<scalar_type> Bp_in(nb_dof_p, 0.0);
   {
   // workspace.set_assembled_vector(Bp_in); 
   // NICHE
    workspace.add_expression("2/element_size*p*Test_p*20", mim_ls_bd, CUT_REGION);// 1 is the region		
    workspace.add_expression("-nlsv.Grad_p*Test_p*tau- nlsv.Grad_Test_p*p*tau ", mim_ls_bd, CUT_REGION); 
   //NICHE
   //  workspace.add_expression( "permeability*tau*[0,1].Grad_p*Test_p ", mim_ls_bd, CUT_REGION);
   workspace.add_expression( "(1+beta)*p.Test_p + tau*Grad_p.Grad_Test_p", mim_ls_in, CUT_REGION);
   workspace.assembly(2);
   gmm::copy(workspace.assembled_matrix(),K_in);
   workspace.clear_expressions();
   workspace.add_expression("+[0].Test_p*tau + p_old.Test_p + Test_p*Div_u_old", mim_ls_in,CUT_REGION);
   workspace.add_expression("(beta)*p_iter.Test_p - Test_p*Trace(Sym(Grad_u_iter))", mim_ls_in,CUT_REGION); // tau
   workspace.set_assembled_vector(Bp_in);
   workspace.assembly(1);
   workspace.clear_expressions();
   workspace.add_expression("2/element_size*p*Test_p", mim, LEFT);
   workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, LEFT); 	
   workspace.add_expression("2/element_size*p*Test_p", mim, RIGHT);
   workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, RIGHT);
   if (N_=3){
   workspace.add_expression("2/element_size*p*Test_p", mim, LEFT);
   workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, LEFT); 	
   workspace.add_expression("2/element_size*p*Test_p", mim, RIGHT);
   workspace.add_expression("-Grad_p.Normal*Test_p*tau - Grad_Test_p.Normal*p*tau ", mim, RIGHT);
     }
    	
   workspace.assembly(2);
   gmm::add(workspace.assembled_matrix(),K_in);
   std::cout<< "end kin"<< std::endl; 
   //gmm::copy(gmm::sub_matrix(Kp,
                                                             //gmm::sub_interval(0, nb_dof_p),
                                                             //gmm::sub_interval(0, nb_dof_p) 
                                                             //),K_in);
                //gmm::copy(K_in,K_out);                                             
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
    // Bp[nb_dof_p+i] =Bp[ii];
     if(ls_i<=0) Bp[ii]+=Bp_in[ii];
     else        Bp[nb_dof_p+i]+= Bp_in[ii];
      for (size_type j = 0; j < eXt_dof.size(); ++j) {
	    size_type jj = eXt_dof[j];
       double ls_j = ls_function(mf_p.point_of_basic_dof(jj),time, LS_TYPE)[0];
            // std::cout<< "ls_i "<< ls_i << " ls_j "<< ls_j << std::endl;
	        if ( (ls_i <= 0) && (ls_j <= 0) ) {
	         // i and j are both In
	         Kp(ii , jj) += K_in(ii, jj);
	         Kp(i + dof_shift, j  + dof_shift) += K_out(ii, jj);
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
   }
}
    
    
    
    
    
    
    
    
    
    
    //dummy part of the matrix
    // workspace.add_expression("-200*1.e+12*Test_p "	, mim_ls_out);
    // workspace.add_expression("-200*1.e+12*Test_p "	, mim_ls_bd);
    //workspace.assembly(1);
	// workspace.clear_expressions();
	return;}



// -------------------------------------------------
// Assembling displacement matrix for fixed stress
// -------------------------------------------------
void biotls_problem::assembly_u(double dt,double time){
	std::cout<<"biot_assembler::assembly_u(double dt)" <<std::endl;
   	getfem::size_type nb_dof_u = mf_u.nb_dof();
    getfem::size_type nb_dof_p = mf_p.nb_dof();
    getfem::size_type nb_dof_p_x = nb_dof_p + nb_x_dof_p;
    getfem::size_type nb_dof_u_x = nb_dof_u + nb_x_dof_u;
    
    gmm::clear(Bu);  gmm::clear(Ku);
    gmm::resize(U_old, nb_dof_u);
    // gmm::resize(P_iter, nb_dof_p);
    gmm::resize(U, nb_dof_u_x); gmm::resize(Ku, nb_dof_u_x, nb_dof_u_x);  gmm::resize(Bu, nb_dof_u_x);
	getfem::ga_workspace workspace; configure_workspace(workspace,dt);
	
    workspace.add_fem_variable("p_iter", mf_p, gmm::sub_interval(0,nb_dof_p), P_iter);
	workspace.add_fem_variable("u", mf_u, gmm::sub_interval(0, nb_dof_u), U);
	workspace.add_fem_variable("u_old", mf_u, gmm::sub_interval(0, nb_dof_u), U_old);
	// workspace.add_fem_variable("u_iter", mf_u, gmm::sub_interval(0, nb_dof_u), U_iter);
    
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
    workspace.assembly(1);
	workspace.clear_expressions();
    // ----- momentum equation
	workspace.add_expression("2*Sym(Grad_u):Grad_Test_u + C2*Div_u*Div_Test_u", mim , UNCUT_REGION_IN); // stress tensor 
    workspace.add_expression("penalty*u.Test_u" , mim, BOTTOM); //neumann disp
	workspace.assembly(2);   
    gmm::copy(workspace.assembled_matrix(), gmm::sub_matrix(Ku,
                                                             gmm::sub_interval(0, nb_dof_u),
                                                             gmm::sub_interval(0, nb_dof_u) 
                                                             ));
	workspace.clear_expressions();
    
	if(N_==2) workspace.add_expression("[0,-1].Test_u"  ,  mim_ls_in,  UNCUT_REGION_IN);
    if(N_==3) workspace.add_expression("[0,0,-1].Test_u",  mim_ls_in,  UNCUT_REGION_IN);
    workspace.add_expression("C1*p_iter*Div_Test_u "    ,  mim_ls_in,  UNCUT_REGION_IN);
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
   workspace.add_expression("u.Test_u "	,mim_ls_out, CUT_REGION);;
   workspace.assembly(2);
   gmm::copy(workspace.assembled_matrix(),K_out);
   workspace.clear_expressions();
   std::cout<< "end kout"<< std::endl; 
   // Kin for enriched dof
   sparse_matrix_type K_in(nb_dof_u, nb_dof_u );
       // internal for displacement
   workspace.add_expression( "2*Sym(Grad_u):Grad_Test_u + C2*Div_u*Div_Test_u", mim_ls_in, CUT_REGION);
   workspace.assembly(2);
   gmm::copy(workspace.assembled_matrix(),K_in);
   workspace.clear_expressions();
   
   std::vector<scalar_type> B_in(nb_dof_u, 0.0);
   workspace.set_assembled_vector(B_in);
   if(N_==2) workspace.add_expression("[0,-1].Test_u", mim_ls_in,CUT_REGION);
   if(N_==3) workspace.add_expression("[0,0,-1].Test_u", mim_ls_in,CUT_REGION);
   workspace.add_expression("C1*p_iter*Div_Test_u ",  mim_ls_in, CUT_REGION);
   workspace.assembly(1);
   //workspace.clear_expressions();
   //workspace.set_assembled_vector(B_in2);
   //workspace.add_expression("[0,-1].Test_u", mim,CUT_REGION);
   //workspace.add_expression("C1*p_iter*Div_Test_u ",  mim, CUT_REGION);
   //workspace.assembly(1);
   workspace.clear_expressions();
   // for(int i=0;i<nb_dof_u; i++)std::cout << B_in[i]<<"<---mimlsin  mim--->"<< B_in2[i]<<std::endl;
   // std::cin.ignore();
   std::cout<< "end kin"<< std::endl; 

   
   // Mapping
   {// mapping for displacement
    std::cout<<"start mapping displacement"<<std::endl;
    size_type dof_shift = nb_dof_u ;
    for (size_type i = 0; i < eXt_dof_u.size(); ++i) {
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
   }
    // end enriched dof
  
   //  std::cout<<Ku<<std::endl;std::cin.get();
} // end assembling momentum matrix for fixed stress approach

//====================================================
// method for solving the laplace problem
//====================================================
void biotls_problem::solve(double time){
    
  // biot_precond<sparse_matrix_type> bPR(Ku,Kp);
  size_type restart = 50;
  scalar_type cond;
  // gmm::identity_matrix PM; // no precond
  gmm::iteration iter(1.e-8);  // iteration object with the max residu
  iter.set_noisy(1);               // output of iterations (2: sub-iteration)
  iter.set_maxiter(1000); // maximum number of iterations
  gmm::diagonal_precond<sparse_matrix_type> PR(K);
  // gmm::gmres(K, UP, B, PR, restart, iter);
  gmm::SuperLU_solve(K, UP , B, cond);
  std::cout << "  Condition number momentum: " << cond << std::endl;
  getfem::size_type nb_dof_u = mf_u.nb_dof();
  getfem::size_type nb_dof_p = mf_p.nb_dof();
  gmm::copy(gmm::sub_vector(UP,gmm::sub_interval(0, nb_dof_u)), U);				
  gmm::copy(gmm::sub_vector(UP,gmm::sub_interval(nb_dof_u , nb_dof_p)), P);	
  
  {  
  std::vector<scalar_type> PIn(mf_p.nb_dof(), 0.0);
  for (size_type i = 0; i < eXt_dof.size(); ++i) 
    {
     size_type ii = eXt_dof[i];
     double ls_i = ls_function(mf_p.point_of_basic_dof(ii),time, LS_TYPE)[0];
     if (ls_i >= 0) PIn[ii] = UP[mf_u.nb_dof() + mf_p.nb_dof()+ i];
    }
  for (size_type i = 0; i < mf_p.nb_dof(); ++i) 
    {
     double ls_i = ls_function(mf_p.point_of_basic_dof(i),time, LS_TYPE)[0];
     if (ls_i < 0)  PIn[i] = UP[mf_u.nb_dof() +i];
    }
  gmm::copy(PIn,P_old);
  } 
  
  {  
   std::vector<scalar_type> UIn(mf_u.nb_dof(), 0.0);
   for (size_type i = 0; i < eXt_dof_u.size(); ++i) 
    {
     size_type ii = eXt_dof_u[i];
     double ls_i = ls_function(mf_u.point_of_basic_dof(ii),time, LS_TYPE)[0];
     if (ls_i >= 0) UIn[ii] = UP[mf_u.nb_dof() + mf_p.nb_dof() + eXt_dof.size() + i];
    }
  for (size_type i = 0; i < mf_u.nb_dof(); ++i) 
    {
     double ls_i = ls_function(mf_u.point_of_basic_dof(i),time, LS_TYPE)[0];
     if (ls_i < 0)  UIn[i] = UP[i];
    }
  gmm::copy(UIn,U_old);
  } 
 // gmm::copy(P,P_old);  gmm::copy(U,U_old);
 }
 
//====================================================
// method for solving the fixed stress  problem
//====================================================
void biotls_problem::solve_fix_stress(double dt, int max_iter,double time){

  
  double epsu=1.e-6; double epsp=1.e-6;
  double rel_unorm=1; double rel_pnorm=1; int fix_count=0;
  double old_unorm=1; double new_unorm=1;
  double old_pnorm=1; double new_pnorm=1;
  getfem::size_type nb_dof_u = mf_u.nb_dof();
  getfem::size_type nb_dof_p = mf_p.nb_dof();
  int min_iter=2;
  while( ( fix_count < max_iter  && ( rel_unorm>epsu ||  rel_pnorm > epsp)) ||  fix_count < min_iter  )
		{
		 fix_count++; 
		 std::cout<<"\033[1;34m***** iteration " << fix_count 
                  << " norm p " <<  rel_pnorm 
                 << " norm u " <<  rel_unorm << std::endl;
		 
                 
         std::cout<< " \033[1;32m Assembling pressure"<<std::endl;
         assembly_p(dt,time);
         std::cout<< " \033[1;32m Start solving pressure"<<std::endl;
         
         // precond
         gmm::identity_matrix PM; // no precond
         gmm::diagonal_precond<sparse_matrix_type> PRp(Kp);
         // gmm::clear(P);
         {
          size_type restart = 50;
          gmm::iteration iter(1.e-8);  // iteration object with the max residu
          iter.set_noisy(1);               // output of iterations (2: sub-iteration)
          iter.set_maxiter(1000); // maximum number of iterations
          
          // std::cout<<P.size()<<std::endl;
          // std::cout<<Bp.size()<<std::endl;std::cin.get();
        gmm::gmres(Kp, P, Bp, PRp, restart, iter);
        scalar_type cond;
		 //  gmm::SuperLU_solve(Kp, P , Bp, cond);
        std::cout << "  Condition number pressure: " << cond << std::endl;
        std::vector<scalar_type> PIn(mf_p.nb_dof(), 0.0);
        std::cout<<"Updating P_iter"<<std::endl;
         int nb_exdof_p=eXt_dof.size();
        {  
         for (size_type i = 0; i < eXt_dof.size(); ++i) 
          {
             size_type ii = eXt_dof[i];
             double ls_i = ls_function(mf_p.point_of_basic_dof(ii),time, LS_TYPE)[0];
             if (ls_i >= 0) PIn[ii] = P[mf_p.nb_dof()+ i];
          }
        for (size_type i = 0; i < mf_p.nb_dof(); ++i) 
         {
          double ls_i = ls_function(mf_p.point_of_basic_dof(i),time, LS_TYPE)[0];
          if (ls_i < 0)  PIn[i] = P[i];
         }
         gmm::copy(PIn,P_iter);
         } 
         // updating p_iter
         }
         //--------------------------------------------------------------
         #ifdef PRINT_MATRIX
         std::cout<<"Pressure print matrix"<< std::endl;
         gmm::MatrixMarket_IO::write(("nKp." +  std::to_string(fix_count) + ".mm").c_str() , Kp);
         std::fstream f(("nBp." +  std::to_string(fix_count) + ".mm").c_str(),std::ios::out);
         std::fstream f1(("nP." +  std::to_string(fix_count) + ".mm").c_str(),std::ios::out);
         for (unsigned i=0; i < gmm::vect_size(Bp); ++i) {f<< Bp[i] << "\n"; f1<< P[i] << "\n";} 
         #endif // PRINTMATRIX
         //----------------------------------------------------------------
         

        { // for displacement
         new_pnorm = gmm::vect_norm2(P);
         rel_pnorm=fabs(new_pnorm - old_pnorm)/ (old_pnorm+1.e-18);
         old_pnorm = new_pnorm;
         // updating u
         //////warrnincg
         // gmm::copy(P,P_iter);
         std::cout<< " \033[1;31m Start solving momentum balance"<<std::endl;
         std::cout<<"assembling"<<std::endl;
         assembly_u(dt,time);
         //solving u
         gmm::diagonal_precond<sparse_matrix_type> PRu(Ku);
         {
          size_type restart = 50;
          gmm::iteration iter(1.e-7);  // iteration object with the max residu
          iter.set_noisy(1);               // output of iterations (2: sub-iteration)
          iter.set_maxiter(1000); // maximum number of iterations
          // gmm::MatrixMarket_load("km",Ku);
          // gmm::clear(U);
          gmm::gmres(Ku, U, Bu, PRu, restart, iter);
          scalar_type cond;
		  // gmm::SuperLU_solve(Ku, U , Bu, cond);
		  std::cout << "  Condition number momentum: " << cond << std::endl;
          {
             std::cout << " Updating U_iter"<<std::endl;
             int nb_exdof_u=eXt_dof_u.size();        
             std::vector<scalar_type> UIn(mf_u.nb_dof(), 0.0);
             for (size_type i = 0; i < eXt_dof_u.size(); ++i) 
             {
                size_type ii = eXt_dof_u[i];
                double ls_i = ls_function(mf_u.point_of_basic_dof(ii),time, LS_TYPE)[0];
                if (ls_i >= 0) UIn[ii] = U[mf_u.nb_dof() + i];
            }
            for (size_type i = 0; i < mf_u.nb_dof(); ++i) 
            {
             double ls_i = ls_function(mf_u.point_of_basic_dof(i),time, LS_TYPE)[0];
             if (ls_i < 0)  UIn[i] = U[i];
            }
           gmm::copy(UIn,U_iter);
           } 
        }
         //--------------------------------------------------------------
         #ifdef PRINT_MATRIX
         std::cout<<"Displacement print matrix"<< std::endl;
         gmm::MatrixMarket_IO::write(("nKu." +  std::to_string(fix_count) + ".mm").c_str() , Ku);
         std::fstream f2(("nBu." +  std::to_string(fix_count) + ".mm").c_str(),std::ios::out);
         std::fstream f3(("nU." +  std::to_string(fix_count) + ".mm").c_str(),std::ios::out);
         for (unsigned i=0; i < gmm::vect_size(Bu); ++i) {f2 << Bu[i] << "\n"; f3<< U[i] << "\n";} 
         #endif // PRINT_MATRIX
         //----------------------------------------------------------------
         new_unorm = gmm::vect_norm2(U);
         rel_unorm=fabs(new_unorm - old_unorm)/ (old_unorm+1e-20);
         old_unorm = new_unorm;
       }
     }
  	 std::cout<<"\033[1;34m***** last iteration " << fix_count 
              << " norm p " <<  rel_pnorm 
              << " norm u " <<  rel_unorm << std::endl;
     
     
     std::cout<<"\033[1;34m updating old vectors "<< std::endl;
      std::vector<scalar_type> PIn(mf_p.nb_dof(), 0.0);
      int nb_exdof_p=eXt_dof.size();
      int nb_exdof_u=eXt_dof_u.size();
    {  
         for (size_type i = 0; i < eXt_dof.size(); ++i) 
          {
             size_type ii = eXt_dof[i];
             double ls_i = ls_function(mf_p.point_of_basic_dof(ii),time, LS_TYPE)[0];
             if (ls_i >= 0) PIn[ii] = P[mf_p.nb_dof()+ i];
          }
        for (size_type i = 0; i < mf_p.nb_dof(); ++i) 
         {
          double ls_i = ls_function(mf_p.point_of_basic_dof(i),time, LS_TYPE)[0];
          if (ls_i < 0)  PIn[i] = P[i];
         }
         gmm::copy(PIn,P_old);
    } 
  
    {  
      std::vector<scalar_type> UIn(mf_u.nb_dof(), 0.0);
      for (size_type i = 0; i < eXt_dof_u.size(); ++i) 
      {
       size_type ii = eXt_dof_u[i];
       double ls_i = ls_function(mf_u.point_of_basic_dof(ii),time, LS_TYPE)[0];
       if (ls_i >= 0) UIn[ii] = U[mf_u.nb_dof() + i];
      }
      for (size_type i = 0; i < mf_u.nb_dof(); ++i) 
      {
       double ls_i = ls_function(mf_u.point_of_basic_dof(i),time, LS_TYPE)[0];
       if (ls_i < 0)  UIn[i] = U[i];
      }
      gmm::copy(UIn,U_old);
    } 
     
     gmm::clear(U_iter); gmm::clear(P_iter);
     gmm::resize(UP, nb_dof_u + nb_exdof_u + nb_dof_p + nb_exdof_p);  
     gmm::copy( gmm::sub_vector(U ,gmm::sub_interval(0, nb_dof_u)),
                gmm::sub_vector(UP,gmm::sub_interval(0, nb_dof_u))
               );				
     gmm::copy(
                gmm::sub_vector(U ,gmm::sub_interval(nb_dof_u, nb_exdof_u)),
                gmm::sub_vector(UP,gmm::sub_interval(nb_dof_u+nb_dof_p + nb_exdof_p, nb_exdof_u))
                );				
     gmm::copy(P, gmm::sub_vector(UP,gmm::sub_interval(nb_dof_u , nb_dof_p + nb_exdof_p)));	
  std::cout<<"end updating old"<<std::endl;
 }
 


// method for solving the laplace problem
void biotls_problem::print(double time,int istep,double time_ls){
    std::vector<scalar_type> PIn(mf_p.nb_dof(), 0.0);
    std::vector<scalar_type> POut(mf_p.nb_dof(), 0.0);
    std::vector<scalar_type> Pm(mf_p.nb_dof(), 0.0);
    {
        for (size_type i = 0; i < eXt_dof.size(); ++i) {
         size_type ii = eXt_dof[i];
          double ls_i = ls_function(mf_p.point_of_basic_dof(ii),time_ls, LS_TYPE)[0];
         if (ls_i < 0) {
          POut[ii] = +UP[mf_u.nb_dof() + mf_p.nb_dof() + i];
      }
         else
          PIn[ii] = UP[mf_u.nb_dof() + mf_p.nb_dof()+ i];
          Pm[ii] = 0.5* (UP[mf_u.nb_dof() +ii] + UP[mf_u.nb_dof() +mf_p.nb_dof()+ i]);
      }
      for (size_type i = 0; i < mf_p.nb_dof(); ++i) {
           double ls_i = ls_function(mf_p.point_of_basic_dof(i),time_ls, LS_TYPE)[0];
        if (ls_i < 0) {
          PIn[i] = UP[mf_u.nb_dof() +i];Pm[i] = UP[mf_u.nb_dof() +i];}
        else {
          POut[i] = +UP[mf_u.nb_dof() +i]; Pm[i] = UP[mf_u.nb_dof() +i];}
      }
        getfem::vtk_export vtke(p_des.datafilename +".pressure."+std::to_string(istep)+".vtk");
        vtke.exporting(mf_p);
        vtke.write_mesh();
        vtke.write_point_data(mf_p, PIn, "NSolution_PIn");
        vtke.write_point_data(mf_p, POut, "NSolution_POut");
        vtke.write_point_data(mf_p, Pm, "NSolution_Pm");
    
    }
    std::cout<<"end print numerical pressure" <<std::endl;
    // dispalcement
    std::vector<scalar_type> UIn(mf_u.nb_dof(), 0.0);
    std::vector<scalar_type> UOut(mf_u.nb_dof(), 0.0);
    std::vector<scalar_type> Um(mf_u.nb_dof(), 0.0);
    int shift_xdof_u = mf_u.nb_dof() + mf_p.nb_dof() + eXt_dof.size(); // dof for extended dofs
    {

      for (size_type i = 0; i < mf_u.nb_dof(); ++i) {
           double ls_i = ls_function(mf_u.point_of_basic_dof(i),time_ls, LS_TYPE)[0];
        if (ls_i < 0) {
          UIn[i] = UP[i];Um[i] = UP[i];}
        else {
          UOut[i] = +UP[i]; Um[i] = UP[i];}
      }
       for (size_type i = 0; i < eXt_dof_u.size(); ++i) {
         size_type ii = eXt_dof_u[i];
          double ls_i = ls_function(mf_u.point_of_basic_dof(ii),time_ls, LS_TYPE)[0];
         if (ls_i < 0) {
          UOut[ii] = +UP[shift_xdof_u + i];
      }
         else
          UIn[ii] = UP[shift_xdof_u + i];
          Um[ii] = 0.5* (UP[ii] + UP[shift_xdof_u + i]);
      }
        getfem::vtk_export vtke(p_des.datafilename +".dislpacement."+std::to_string(istep)+".vtk");
        vtke.exporting(mf_u);
        vtke.write_mesh();
        vtke.write_point_data(mf_u, UIn, "NSolution_UIn");
        vtke.write_point_data(mf_u, UOut, "NSolution_UOut");
        vtke.write_point_data(mf_u, Um, "NSolution_Um");
    
    }
    std::cout<<"end print numerical pressure" <<std::endl;
      // 3. NUMERICAL SOLUTION
  // Using the "cut mesh" computed by the mesh_level_set object for visualization.

  getfem::mesh mcut;
  mls.global_cut_mesh(mcut);

  getfem::vtk_export vtkc("cmesh.vtk");
  vtkc.exporting(mcut);
  vtkc.write_mesh();
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
   getfem::mesh_fem mfcut_u(mcut, mf_u.get_qdim());
   mfcut_u.set_classical_discontinuous_finite_element(3, 0.01);

  // Interpolate the In and Out parts of the solution
  gmm::resize(Um, mfcut_u.nb_dof()); gmm::clear(Um);
  
  {
    
      std::vector<scalar_type> Uaux(Um);
      std::cout<< "start interpolation of velocity"<<std::endl; 
      getfem::interpolation(mf_u, mfcut_u, UIn, Uaux);
      getfem::interpolation(mf_u, mfcut_u, UOut, Um);
      std::cout<< "end interpolation of velocity"<<std::endl; 
    
      // Check if element cv is "more In" than "Out" (based on dmean)
      // and set values af all local dofs correspondingly.
      for (dal::bv_visitor cv(mcut.convex_index()); !cv.finished(); ++cv) {
        getfem::mesh_fem::ind_dof_ct idofs = mfcut_u.ind_basic_dof_of_element(cv);
        scalar_type dmean = 0;
        for (size_type i=0; i < idofs.size(); ++i) 
          dmean += ls_function(mfcut_u.point_of_basic_dof( idofs[i] ) ,time_ls, LS_TYPE)[0];
        if (dmean < 0) 
          for (size_type i=0; i < idofs.size(); ++i) {
        size_type ii = idofs[i];
        Um[ii] = Uaux[ii];
          }
       }
  }
  // Export discontinuous solution
  std::string namefile= p_des.datafilename +"d." +  std::to_string(istep) +".vtk";
  getfem::vtk_export vtkd(namefile);
  vtkd.exporting(mfcut);
  vtkd.write_mesh();
  vtkd.write_point_data(mfcut, Pm, "p");
  vtkd.write_point_data(mfcut_u, Um, "u");
  std::cout<<"end printing"<<std::endl;
 }
 
 
  // level set function
  base_small_vector biotls_problem::ls_function(const base_node P, double time,int num) {
  scalar_type x = P[0]*p_des.l_ref, y = P[1]*p_des.l_ref, z=0;
  if (N_=3)  z = P[2]*p_des.l_ref;
  y = P[1]*p_des.l_ref;
  // time*=p_des.t_ref;
  base_small_vector res(2);
  switch (num) {
    case 0: {
      res[0] = y-3050  ;
      res[1] = -.5 + x;
    } break;
    case 1: {
      res[0] = y - 3050 + 500*time/(1e+8);
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.27;
    } break;
    case 2: {
      res[0] = y - (4.e+2 * time / (1.e+8) * sin(2 * 3.14 *x/4000) + 3300);
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
    } break;
    case 3: {
      res[0] = z  - 3050 + 500*time/(1e+8);
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
    } break;
    case 4: {
      res[0] = z - (4.e+2 * time / (1.e+8) * sin(2 * 3.14 *x/4000)*sin(2 * 3.14 *y/4000) + 3300);
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
    } break;
    default: assert(0);
  }
  return res;
}


void biotls_problem::update_ls(double time, int iter){
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
               mesh.region(UNCUT_REGION_IN).add(i_cv);
          else mesh.region(UNCUT_REGION_OUT).add(i_cv);
      }
    }
   
    { // Just to see what elements are cut by the level set ls:
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
    
    eXt_dof_u.clear();
    {
       dal::bit_vector ndofs = mf_u.basic_dof_on_region(CUT_REGION);
       for (dal::bv_visitor i(ndofs); !i.finished(); ++i)
        eXt_dof_u.push_back(i);
       nb_x_dof_u=eXt_dof_u.size();
     }
    std::cout << "# Cut Elements      = " << mesh.region(CUT_REGION).size() << std::endl;
    std::cout << "# dof to extend are  for pressure   = " << eXt_dof.size() << std::endl;
    std::cout << "# dof to extend are  for displacement   = " << eXt_dof_u.size() << std::endl;
    compute_normal_2_ls();
    
    
    // ls.touch();
    // mls.adapt();mls.global_cut_mesh(mesh_ls);
    // mim_ls_in.adapt(); mim_ls_out.adapt();mim_ls_bd.adapt();  mim_ls_all.adapt();

    
    }
// print the solution on active part of level set
void biotls_problem::print_crop(double time,int istep,double time_ls){
    std::cout<< "biotls_problem::print crop" <<std::endl;
    
    
    std::vector<scalar_type> PIn(mf_p.nb_dof(), 0.0);
    std::vector<scalar_type> POut(mf_p.nb_dof(), 0.0);
    std::vector<scalar_type> Pm(mf_p.nb_dof(), 0.0);
    {
        for (size_type i = 0; i < eXt_dof.size(); ++i) {
         size_type ii = eXt_dof[i];
          double ls_i = ls_function(mf_p.point_of_basic_dof(ii),time_ls, LS_TYPE)[0];
         if (ls_i < 0) {
          POut[ii] = -100+UP[mf_u.nb_dof() + mf_p.nb_dof() + i];
      }
         else
          PIn[ii] = UP[mf_u.nb_dof() + mf_p.nb_dof()+ i];
          Pm[ii] = 0.5* (UP[mf_u.nb_dof() +ii] + UP[mf_u.nb_dof() +mf_p.nb_dof()+ i]);
      }
      for (size_type i = 0; i < mf_p.nb_dof(); ++i) {
           double ls_i = ls_function(mf_p.point_of_basic_dof(i),time_ls, LS_TYPE)[0];
        if (ls_i < 0) {
          PIn[i] = UP[mf_u.nb_dof() +i];Pm[i] = UP[mf_u.nb_dof() +i];}
        else {
          POut[i] = -100+UP[mf_u.nb_dof() +i]; Pm[i] = UP[mf_u.nb_dof() +i];}
      }
    }
    std::cout<<"end print numerical pressure" <<std::endl;
    // dispalcement
    std::vector<scalar_type> UIn(mf_u.nb_dof(), 0.0);
    int shift_xdof_u = mf_u.nb_dof() + mf_p.nb_dof() + eXt_dof.size(); // dof for extended dofs
    {

      for (size_type i = 0; i < mf_u.nb_dof(); ++i) {
           double ls_i = ls_function(mf_u.point_of_basic_dof(i),time_ls, LS_TYPE)[0];
        if (ls_i < 0) UIn[i] = UP[i];
      }
       for (size_type i = 0; i < eXt_dof_u.size(); ++i) {
         size_type ii = eXt_dof_u[i];
          double ls_i = ls_function(mf_u.point_of_basic_dof(ii),time_ls, LS_TYPE)[0];
          if (ls_i >= 0)  UIn[ii] = UP[shift_xdof_u + i];
         
        }
     
    }
   
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
          std::vector<scalar_type> UIn_dim(mf_u.nb_dof(), 0.0);gmm::copy(UIn,UIn_dim);gmm::scale(UIn_dim,p_des.u_ref);
          std::string namefile= p_des.datafilename +".crop." +  std::to_string(istep) +".vtk";
          getfem::vtk_export vtkd(namefile);
          vtkd.exporting(mesh_dim);
          vtkd.write_mesh();
          vtkd.write_point_data(mf_p, PIn_dim, "p");
          vtkd.write_point_data(mf_u, UIn_dim, "u");
          std::cout<<"end printing"<<std::endl;
       
        }
        
    mesh.transformation(Mm1); 
      }  
    
// evaluate birnak to a level set
void biotls_problem::compute_normal_2_ls(){
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
    

