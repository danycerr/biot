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
	getfem::regular_unit_mesh(mesh, nsubdiv, pgt, 1);
    // A trasformation for the squarred mesh
	  bgeot::base_matrix M(N_,N_);
	  for (size_type i=0; i < N_; ++i) {
	    M(i,i) = 4000.0;
	  }
	//  if (N>1) { M(0,1) = 0; }
	//
	  mesh.transformation(M);
    // // End of mesh generation

  
    // set  integration methods  
    getfem::pfem pf_u = getfem::fem_descriptor(p_des.FEM_TYPE_U);
    getfem::pfem pf_p = getfem::fem_descriptor(p_des.FEM_TYPE_P);
	getfem::pintegration_method ppi = getfem::int_method_descriptor(p_des.INTEGRATION);
    getfem::pintegration_method simp_ppi = getfem::int_method_descriptor(p_des.SIMPLEX_INTEGRATION);
    
	mim.set_integration_method(mesh.convex_index(), ppi);
    mf_u.set_qdim(N_);
	mf_u.set_finite_element(mesh.convex_index(), pf_u); // finite element for displacement
    mf_p.set_finite_element(mesh.convex_index(), pf_p); //  finite element for pressure
    // boundary conditions zones
    gen_bc();
    
    // initialization of level set
    mls.add_level_set(ls);
    ls.reinit(); 
    for (size_type d = 0; d < ls.get_mesh_fem().nb_basic_dof(); ++d) {
        ls.values(0)[d] = ls_function(ls.get_mesh_fem().point_of_basic_dof(d),0, 2)[0];
    }
    mls.adapt();
    mls.global_cut_mesh(mesh_ls);
    
    // integration method on lev set
    
    mim_ls_all.set_integration_method(mesh.convex_index(), ppi);
    mim_ls_all.set_simplex_im(simp_ppi);
    mim_ls_all.adapt();
    
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
    mfls_u_old.adapt();mfls_p_old.adapt();
    
	mf_u.set_qdim(bgeot::dim_type(2)); //number of variable
    mfls_u.set_qdim(bgeot::dim_type(2)); //number of variable
    mfls_p.set_qdim(bgeot::dim_type(1)); //number of variable
    mfls_u_old.set_qdim(bgeot::dim_type(2)); //number of variable
    mfls_p_old.set_qdim(bgeot::dim_type(1)); //number of variable
    
    getfem::mesh mesh_ls;
    mls.global_cut_mesh(mesh_ls);
    std::cout<<"laplace cut fem::save_mesh  "<<std::endl;
    getfem::vtk_export exp_cut("cutmesh.vtk");
    exp_cut.exporting(mesh_ls);
    exp_cut.write_mesh();
    
    
    // init vector
    getfem::size_type nb_dof_u = mfls_u.nb_dof();
    getfem::size_type nb_dof_p = mfls_p.nb_dof();
    
    gmm::resize(B, nb_dof_u + nb_dof_p); gmm::clear(B);    // rhs monolithic problem
    gmm::resize(Bp, nb_dof_p); gmm::resize(Bu, nb_dof_u);  // rhs fixed stress pressure displacement
    gmm::resize(UP, nb_dof_u + nb_dof_p); gmm::clear(UP);  // solution monolithic
    // displacement
    gmm::resize(U, nb_dof_u); gmm::resize(U_old, nb_dof_u); gmm::resize(U_iter, nb_dof_u); 
    gmm::clear(U);gmm::clear(U_old);gmm::clear(U_iter);
    // pressure
    gmm::resize(P, nb_dof_p); gmm::resize(P_old, nb_dof_p); gmm::resize(P_iter, nb_dof_p); 
    std::fill(P.begin(), P.end(), 1);
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
		else {
			mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
		}
	  }
      
      
      
      
     dal::bit_vector bv_cv = mesh.convex_index();
     size_type i_cv = 0;
    for (i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv) {
      if (mls.is_convex_cut(i_cv)) {
       mesh.region(CUT_REGION).add(i_cv);
       std::cout<<"CUT element"<<std::endl;
       // Fill phi with the jump of the flux across ls:
       getfem::mesh_fem::ind_dof_ct idofs = mfls_u.ind_basic_dof_of_element(i_cv);
       std::cout<< "Dof for element "<< i_cv << " are " <<   idofs.size()<< std::endl;
       // for (size_type i=0; i < idofs.size(); ++i) 
       // phi[idofs[i]] = phivalue;      
       ////  Cell value
       // phicell[i_cv] = 0.5;
    }
    else mesh.region(UNCUT_REGION).add(i_cv);
    }
      
      
      
      
      
      
      
      
    }
// end bc generations
// =======================================================

//=======================================================
//Configure Workspace for all methods
//=======================================================
void biotls_problem::configure_workspace(getfem::ga_workspace & workspace,double dt){
    std::cout << "biotls_problem::configure_workspace::Configuring workspace " << std::endl;
    tau_[0] = dt; // dt into  the workspace
    workspace.add_fixed_size_constant("tau", tau_);
    
    vmu_[0] = p_des.mu_s;
	workspace.add_fixed_size_constant("mu", vmu_);
    //---------------------------------------------------------
	bm_[0] =  p_des.biot_modulus; //biot modoulus
	workspace.add_fixed_size_constant("bm", bm_);
    //---------------------------------------------------------
	lambda_[0] = p_des.lambda_l ;
	workspace.add_fixed_size_constant("lambda", lambda_);
    //---------------------------------------------------------
	alpha_[0] = p_des.alpha;
	workspace.add_fixed_size_constant("alpha", alpha_);
    //---------------------------------------------------------
	permeability_[0] = p_des.k;
	workspace.add_fixed_size_constant("permeability", permeability_);
    //---------------------------------------------------------
	force_[0] = 1.e+0;
	workspace.add_fixed_size_constant("force", force_);
    std::cout << "end of Configuring workspace " << std::endl;
    //---------------------------------------------------------
   	//  getfem::base_vector beta(1); beta[0] = 1*(alpha[0] * alpha[0] ) / (2 * p_des.mu_s + p_des.lambda_l);
	beta_[0] = 1*(alpha_[0] * alpha_[0]  ) / (-2 * p_des.mu_s / 3 + p_des.lambda_l);
	workspace.add_fixed_size_constant("beta", beta_);
    
    penalty_[0] = 1.e+12; // 1---10
	workspace.add_fixed_size_constant("penalty", penalty_);
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
void biotls_problem::assembly(double dt) {
    std::cout<< "biotls_problem::assembly" <<std::endl;
    // assembly matrix
    // resize and clear vector and matrix
	getfem::size_type nb_dof_u = mfls_u.nb_dof();
    getfem::size_type nb_dof_p = mfls_p.nb_dof();
    gmm::resize(B, nb_dof_u + nb_dof_p); gmm::clear(B);
    gmm::resize(UP, nb_dof_u + nb_dof_p); //gmm::clear(UP);
    gmm::resize(U, nb_dof_u); gmm::resize(U_old, nb_dof_u); 
    gmm::resize(P, nb_dof_p); gmm::resize(P_old, nb_dof_p);
    gmm::resize(K, nb_dof_u + nb_dof_p, nb_dof_u + nb_dof_p); gmm::clear(K);
 	
	getfem::ga_workspace workspace; // generic workspace for the solution
    configure_workspace(workspace,dt); 
    // getfem::base_vector invdt(1); invdt[0] = 1/dt;
	// workspace.add_fixed_size_constant("invdt", invdt);
    // 
   
   	workspace.add_fem_variable("u", mfls_u, gmm::sub_interval(0, nb_dof_u), U);
	workspace.add_fem_variable("p", mfls_p, gmm::sub_interval(nb_dof_u, nb_dof_p), P);
	workspace.add_fem_variable("u_old", mfls_u, gmm::sub_interval(0, nb_dof_u), U_old);
	workspace.add_fem_variable("p_old", mfls_p, gmm::sub_interval(nb_dof_u,nb_dof_p), P_old);
  	// ------------------ expressions --------------------------
	workspace.add_expression( "2*mu*Sym(Grad_u):Grad_Test_u + lambda*Div_u*Div_Test_u"
                              "- alpha*p.Div_Test_u ", mim_ls_in);
    workspace.add_expression( "+(1/bm)*p.Test_p + tau*permeability*Grad_p.Grad_Test_p"
                              "+ alpha*Test_p*Div_u", mim_ls_in);
    workspace.assembly(2);
    gmm::copy(workspace.assembled_matrix(), K);
	workspace.clear_expressions();
   

	//======= RHS =====================
    if(N_== 2 )	workspace.add_expression("[0,-0].Test_u", mim_ls_in);
    else        workspace.add_expression("[0,0,-0].Test_u", mim_ls_in);
    workspace.add_expression("+[+0.0].Test_p*tau + (0/bm)*p_old.Test_p + alpha*Test_p*Div_u_old", mim_ls_in);
    workspace.set_assembled_vector(B);
    workspace.assembly(1);
    // gmm::copy(workspace.assembled_vector(),B);
	workspace.clear_expressions();

	//Boudanry conditions // NICHE
	getfem::base_vector penalty(1); penalty[0] = 2e+4; // 1---10
	workspace.add_fixed_size_constant("penalty", penalty);
	//Matrix term
	workspace.add_expression("penalty/element_size*p*Test_p", mim_ls_in, LEFT);
    workspace.add_expression("penalty/element_size*p*Test_p", mim_ls_in, RIGHT);// 1 is the region	
	workspace.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim_ls_in, LEFT); 	
	workspace.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p", mim_ls_in, RIGHT); 
	workspace.assembly(2);
    gmm::add(workspace.assembled_matrix(), K);
	workspace.clear_expressions();
		
	//rhs term
    if(N_ == 2){
       // workspace.add_expression("0*penalty*u.Test_u" "+ penalty*0*p*Test_p", mim, TOP); // 1 is the region
       workspace.add_expression("[0,-2*force].Test_u" , mim_ls_in, TOP); //neumann disp
	   // workspace.add_expression("0*penalty*u.Test_u" "+ 0*penalty*0*p*Test_p", mim, BOTTOM); // 1 is the region
	   workspace.add_expression("[0,+2*force].Test_u" , mim_ls_in, BOTTOM); //neumann disp
    }
	else if(N_== 3){
       // workspace.add_expression("0*penalty*u.Test_u" "+ penalty*0*p*Test_p", mim, TOP); // 1 is the region
       workspace.add_expression("[0,0,-2*force].Test_u" , mim_ls_in, TOP); //neumann disp
	   // workspace.add_expression("0*penalty*u.Test_u" "+ 0*penalty*0*p*Test_p", mim, BOTTOM); // 1 is the region
	   workspace.add_expression("[0,0,+2*force].Test_u" , mim_ls_in, BOTTOM); //neumann disp
    }
	
	workspace.add_expression(" 0*penalty/element_size*Test_p - permeability*Grad_Test_p.Normal*0 ", mim_ls_in, LEFT);
	workspace.add_expression(" 0*penalty/element_size*Test_p - permeability*Grad_Test_p.Normal*0 ", mim_ls_in, RIGHT);
	workspace.assembly(1);
    // gmm::add(workspace.assembled_vector(),B);
	workspace.clear_expressions();
   /// end boundary contions
   
   workspace.add_expression("1.e+18*p*Test_p "	, mim_ls_out);
   workspace.add_expression("1.e+12*u.Test_u "	, mim_ls_out);
   workspace.assembly(2);
   gmm::add(workspace.assembled_matrix(), K);
   workspace.clear_expressions();
} // end assembly




// -------------------------------------------------
// Assembling pressure matrix for fixed stress
// -------------------------------------------------

void biotls_problem::assembly_p(double dt){
    
    gmm::clean(P_old, 1E-10);gmm::clean(U_old, 1E-10);
	gmm::clean(P, 1E-10);gmm::clean(U, 1E-10);
    gmm::clear(Bp); gmm::clear(Kp);
	std::cout<<"biot_assembler::assembly_p(double dt)" <<std::endl;
   	getfem::size_type nb_dof_u = mfls_u.nb_dof();
    getfem::size_type nb_dof_p = mfls_p.nb_dof();
    
    //  gmm::resize(P, nb_dof_p); gmm::resize(P_old, nb_dof_p);gmm::resize(P_iter, nb_dof_p);
    //  gmm::resize(Kp, nb_dof_p, nb_dof_p);  gmm::resize(Bp, nb_dof_p);
	
    
    getfem::ga_workspace workspace;
    configure_workspace(workspace,dt);


	workspace.add_fem_variable("p", mfls_p, gmm::sub_interval(0, nb_dof_p), P);
	workspace.add_fem_variable("p_old", mfls_p, gmm::sub_interval(0,nb_dof_p), P_old);
	workspace.add_fem_variable("p_iter", mfls_p, gmm::sub_interval(0,nb_dof_p), P_iter);
	workspace.add_fem_variable("u", mfls_u, gmm::sub_interval(0, nb_dof_u), U);
	workspace.add_fem_variable("u_old", mfls_u, gmm::sub_interval(0, nb_dof_u), U_old);
	workspace.add_fem_variable("u_iter", mfls_u, gmm::sub_interval(0, nb_dof_u), U_iter);
	
	// Pressure equation
	workspace.add_expression("(1/bm + beta)*p.Test_p + tau*permeability*Grad_p.Grad_Test_p", mim_ls_in); // tau
	// workspace.set_assembled_matrix(Kp);
	workspace.assembly(2);
    gmm::copy(workspace.assembled_matrix(), Kp);
	workspace.clear_expressions();
    
	//======= RHS =====================
	//  workspace.add_expression("[+0.0].Test_p + (1/bm)*invdt*p_old.Test_p + invdt*alpha*Test_p*Trace((Grad_u_old))", mim);// 1/dt
	//  workspace.add_expression("(beta)*invdt*p_iter.Test_p - invdt*alpha*Test_p*Trace((Grad_u_iter))", mim);// 1/dt
	// tau
	workspace.add_expression("[+0.0]*Test_p + (1/bm)*p_old.Test_p + alpha*Test_p*Trace(Sym(Grad_u_old))", mim_ls_in); // tau
	workspace.add_expression("(beta)*p_iter.Test_p - alpha*Test_p*Trace(Sym(Grad_u_iter))", mim_ls_in); // tau
	workspace.set_assembled_vector(Bp);
	workspace.assembly(1);
	workspace.clear_expressions();

	//Boudanry conditions
	// getfem::base_vector penalty(1); penalty[0] = 1.e+4;
	// workspace.add_fixed_size_constant("penalty", penalty);
	// Matrix term
	// workspace.add_expression("0*penalty*u.Test_u" "+ 0*penalty*p*Test_p", mim, TOP); // 1 is the region
	// workspace.add_expression("penalty*u.Test_u" "+ 0*penalty*p*Test_p", mim, BOTTOM); // 1 is the region
	// workspace.add_expression("0*penalty*u.Test_u" "+ 0*penalty*p*Test_p", mim, BOTTOM); // 1 is the region	
	// good one
    // workspace.add_expression("penalty*p*Test_p", mim, LEFT); workspace.add_expression("penalty*p*Test_p", mim, RIGHT);// 1 is the region	
		//Boudanry conditions // NICHE

	//Matrix term
	workspace.add_expression("penalty/element_size*p*Test_p", mim_ls_in, TOP);
	workspace.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim_ls_in, TOP); 	
    workspace.add_expression("penalty/element_size*p*Test_p", mim_ls_in, LEFT);
	workspace.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim_ls_in, LEFT); 	
   	workspace.add_expression("penalty/element_size*p*Test_p", mim_ls_in, RIGHT);
	workspace.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim_ls_in, RIGHT); 	
    
    workspace.assembly(2);
    gmm::add(workspace.assembled_matrix(), Kp);
	workspace.clear_expressions();
	//rhs term
	workspace.add_expression("0*penalty*p*Test_p", mim_ls_in, TOP); 
	workspace.assembly(1);
	workspace.clear_expressions();
    
    //dummy part of the matrix
    workspace.add_expression("1.e+18*p*Test_p "	, mim_ls_out);
    workspace.add_expression("1.e+18*p*Test_p "	, mim_ls_bd);
    workspace.assembly(2);
    gmm::add(workspace.assembled_matrix(), Kp);
	workspace.clear_expressions();
    //dummy part of the matrix
    // workspace.add_expression("-200*1.e+12*Test_p "	, mim_ls_out);
    // workspace.add_expression("-200*1.e+12*Test_p "	, mim_ls_bd);
    //workspace.assembly(1);
	// workspace.clear_expressions();
     for (int i=0; i< nb_dof_p; i++) if(fabs(Kp(i,i)) < 1.e-19 ) Kp(i,i)=1.e+4;
	return;}



// -------------------------------------------------
// Assembling displacement matrix for fixed stress
// -------------------------------------------------
void biotls_problem::assembly_u(double dt){
	std::cout<<"biot_assembler::assembly_u(double dt)" <<std::endl;
   	getfem::size_type nb_dof_u = mfls_u.nb_dof();
    getfem::size_type nb_dof_p = mfls_p.nb_dof();
    gmm::clear(Bu);  gmm::clear(Ku);
    // gmm::resize(P, nb_dof_p); gmm::resize(P_old, nb_dof_p);gmm::resize(P_iter, nb_dof_p);
    // gmm::resize(Ku, nb_dof_u, nb_dof_u);  gmm::resize(Bu, nb_dof_u);
	getfem::ga_workspace workspace; configure_workspace(workspace,dt);
	
    workspace.add_fem_variable("p", mfls_p, gmm::sub_interval(0, nb_dof_p), P);
	workspace.add_fem_variable("p_old", mfls_p, gmm::sub_interval(0,nb_dof_p), P_old);
	workspace.add_fem_variable("p_iter", mfls_p, gmm::sub_interval(0,nb_dof_p), P_iter);
	workspace.add_fem_variable("u", mfls_u, gmm::sub_interval(0, nb_dof_u), U);
	workspace.add_fem_variable("u_old", mfls_u, gmm::sub_interval(0, nb_dof_u), U_old);
	workspace.add_fem_variable("u_iter", mfls_u, gmm::sub_interval(0, nb_dof_u), U_iter);
    
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
	workspace.add_expression("2*mu*Sym(Grad_u):Grad_Test_u + lambda*Div_u*Div_Test_u", mim_ls_in); // stress tensor 
    workspace.add_expression("penalty*u.Test_u" , mim_ls_in, BOTTOM); //neumann disp
	workspace.assembly(2);
    gmm::copy(workspace.assembled_matrix(), Ku);
	workspace.clear_expressions();
	if(N_==2) workspace.add_expression("(2200*0.8 + 1000*0.2 -1000 )*[0,-1].Test_u", mim_ls_in);
    if(N_==3) workspace.add_expression("[0,0,-0].Test_u", mim_ls_in);
    workspace.add_expression("alpha*p*Div_Test_u ", mim_ls_in);
	workspace.assembly(1);
	workspace.clear_expressions();
    // std::cout<< Bu<< std::endl; 
    //dummy part of the matrix
    workspace.add_expression("1.e+12*u.Test_u "	, mim_ls_out);
    workspace.assembly(2);
    gmm::add(workspace.assembled_matrix(), Ku);
	workspace.clear_expressions();
    
    for (int i=0; i< nb_dof_u; i++) if(fabs(Ku(i,i)) < 1.e-16 ) Ku(i,i)=1.e+4;
  
    
} // end assembling momentum matrix for fixed stress approach

//====================================================
// method for solving the laplace problem
//====================================================
void biotls_problem::solve(){
    
  // biot_precond<sparse_matrix_type> bPR(Ku,Kp);
  size_type restart = 50;

  // gmm::identity_matrix PM; // no precond
  // gmm::diagonal_precond<sparse_matrix_type> PR(K);
  gmm::iteration iter(1.e-8);  // iteration object with the max residu
  iter.set_noisy(1);               // output of iterations (2: sub-iteration)
  iter.set_maxiter(1000); // maximum number of iterations
  gmm::gmres(K, UP, B, *bPR_, restart, iter);
  getfem::size_type nb_dof_u = mf_u.nb_dof();
  getfem::size_type nb_dof_p = mf_p.nb_dof();
  gmm::copy(gmm::sub_vector(UP,gmm::sub_interval(0, nb_dof_u)), U);				
  gmm::copy(gmm::sub_vector(UP,gmm::sub_interval(nb_dof_u , nb_dof_p)), P);	
  gmm::copy(U,U_old);gmm::copy(P,P_old);
 }
 
//====================================================
// method for solving the fixed stress  problem
//====================================================
void biotls_problem::solve_fix_stress(double dt, int max_iter){

  
  double epsu=1.e-6; double epsp=1.e-6;
  double rel_unorm=1; double rel_pnorm=1; int fix_count=0;
  double old_unorm=1; double new_unorm=1;
  double old_pnorm=1; double new_pnorm=1;
  
  while(fix_count < max_iter && ( rel_unorm>epsu ||  rel_pnorm > epsp) )
		{
		 fix_count++; 
		 std::cout<<"\033[1;34m***** iteration " << fix_count 
                  << " norm p " <<  rel_pnorm 
                   << " norm u " <<  rel_unorm << std::endl;
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
          
          // gmm::gmres(Kp, P, Bp, PRp, restart, iter);
          scalar_type cond;
		  gmm::SuperLU_solve(Kp, P , Bp, cond);
		  std::cout << "  Condition number (momentum: " << cond << std::endl;
          
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
         
         getfem::size_type nb_dof_u = mf_u.nb_dof();
         getfem::size_type nb_dof_p = mf_p.nb_dof();
         
         new_pnorm = gmm::vect_norm2(P);
         rel_pnorm=fabs(new_pnorm - old_pnorm)/ (old_pnorm+1.e-18);
         old_pnorm = new_pnorm;
         // updating u
         gmm::copy(P,P_iter);
         assembly_u(dt);
         //solving u
         std::cout<< " \033[1;31m Start solving momentum balance"<<std::endl;
         gmm::diagonal_precond<sparse_matrix_type> PRu(Ku);
         {
          size_type restart = 50;
          gmm::iteration iter(1.e-7);  // iteration object with the max residu
          iter.set_noisy(1);               // output of iterations (2: sub-iteration)
          iter.set_maxiter(1000); // maximum number of iterations
          // gmm::MatrixMarket_load("km",Ku);
          // gmm::clear(U);
          // gmm::gmres(Ku, U, Bu, PRu, restart, iter);
          scalar_type cond;
		  gmm::SuperLU_solve(Ku, U , Bu, cond);
		  std::cout << "  Condition number (momentum: " << cond << std::endl;
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
	     
         
         gmm::copy(U,U_iter);
         new_unorm = gmm::vect_norm2(U);
         rel_unorm=fabs(new_unorm - old_unorm)/ old_unorm;
         old_unorm = new_unorm;
         // updating u
         assembly_p(dt);
     }
  	 std::cout<<"\033[1;34m***** last iteration " << fix_count 
              << " norm p " <<  rel_pnorm 
              << " norm u " <<  rel_unorm << std::endl;
     gmm::copy(U,U_old); gmm::copy(P,P_old);
     gmm::clear(U_iter); gmm::clear(P_iter);
  
 }
 


// method for solving the laplace problem
void biotls_problem::print(int time){
    std::cout<< "biotls_problem::print" <<std::endl;
    // getfem::vtk_export exp(p_des.datafilename + "." +  std::to_string(time) + ".vtk");
    // exp.exporting(mf_u);  	exp.write_point_data(mf_u, U, "u"); 
    //                     	exp.write_point_data(mf_p, P, "p"); 
                            
    dim_type Q = mfls_u.get_qdim();
    getfem::mesh_fem mf(mesh_ls, Q);
    mf.set_classical_discontinuous_finite_element(2, 0.01);
    getfem::model_real_plain_vector Up(mf.nb_dof());
    getfem::interpolation(mfls_u, mf, U, Up);
    getfem::vtk_export exp(p_des.datafilename + "." +  std::to_string(time) + ".vtk");
    exp.exporting(mf);  	exp.write_point_data(mf, Up, "u"); 
    
    dim_type Qp = mfls_p.get_qdim();
    getfem::mesh_fem mfp(mesh_ls, Qp);
    mfp.set_classical_discontinuous_finite_element(2, 0.01);
    getfem::model_real_plain_vector Pp(mfp.nb_dof());
    getfem::interpolation(mfls_p, mfp, P, Pp);
    exp.exporting(mfp);  	exp.write_point_data(mfp, Pp, "p"); 
 }
 
 
  // level set function
  base_small_vector biotls_problem::ls_function(const base_node P, double time,int num) {
  scalar_type x = P[0], y = P[1];
  base_small_vector res(2);
  switch (num) {
    case 0: {
      res[0] = y-0.55  ;
      res[1] = -.5 + x;
    } break;
    case 1: {
      res[0] = y - 3500 + time/(2e+9)*100;
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.27;
    } break;
    case 2: {
      res[0] = y - (4.e+2 * time / (2.e+9) * sin(2 * 3.14 *x/4000) + 3300);
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
    } break;
    default: assert(0);
  }
  return res;
}


void biotls_problem::update_ls(double time){
    std::cout<< "laplacian_problem::update_ls" <<std::endl;
   //  dim_type Q = mfls_u.get_qdim();
    ls.reinit(); 
    for (size_type d = 0; d < ls.get_mesh_fem().nb_basic_dof(); ++d) {
        ls.values(0)[d] = ls_function(ls.get_mesh_fem().point_of_basic_dof(d), time ,2)[0];
         }
    //dim_type Q = mfls_u.get_qdim();
    //getfem::mesh_fem mf(mesh_ls, Q);
    //mf.set_classical_discontinuous_finite_element(2, 0.01);
    //getfem::model_real_plain_vector Up(mf.nb_dof());
    //getfem::interpolation(mfls_u, mf, U, Up);
    
    //dim_type Qp = mfls_p.get_qdim();
    //getfem::mesh_fem mfp(mesh_ls, Qp);
    //mfp.set_classical_discontinuous_finite_element(2, 0.01);
    //getfem::model_real_plain_vector Pp(mfp.nb_dof());
    //getfem::interpolation(mfls_p, mfp, P, Pp);
    mls.adapt();mls.global_cut_mesh(mesh_ls);
    
    mim_ls_in.adapt(); mim_ls_out.adapt();mim_ls_bd.adapt();  mim_ls_all.adapt();
    mfls.adapt();mfls_u.adapt();mfls_p.adapt();
    getfem::size_type nbdofu  = mfls_u.nb_dof();
	getfem::size_type nbdofp = mfls_p.nb_dof();
    // 
    gmm::resize(U_old, nbdofu); gmm::resize(P_old, nbdofp);
    // getfem::interpolation( mf, mfls_u, Up, U_old);
    // getfem::interpolation( mfp, mfls_p, Pp, P_old);
    
    
    
    
    gmm::resize(B, nbdofu + nbdofp);
    
    gmm::resize(P, nbdofp);  gmm::resize(P_iter, nbdofp);
    gmm::resize(Kp, nbdofp , nbdofp ); gmm::resize(Bp, nbdofp);
     
    gmm::resize(U, nbdofu);  gmm::resize(U_iter, nbdofu);
    gmm::resize(Ku, nbdofu , nbdofu ); gmm::resize(Bu, nbdofu);
    // gmm::clear(U_old);
    // gmm::clear(P_old);
    // for(int i=0; i < nbdofp; i++) (P_old[i]>1.e+4? 0:P_old[i]);
    // for(int i=0; i < nbdofu; i++) (U_old[i]>1.e+4? 0:U_old[i]);
    
}

