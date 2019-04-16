#include "biot.hpp" 

// #define LATERAL_INJECTION
void biot_problem::init(void) {

	std::cout<< "biot_problem::init "  << std::endl;

	bgeot::pgeometric_trans pgt = 
		bgeot::geometric_trans_descriptor(p_des.MESH_TYPE);

	// Mesh generation
	N_ = pgt->dim();
	std::vector<size_type> nsubdiv(N_);
	int NX=p_des.nsubdiv;
	std::fill(nsubdiv.begin(),nsubdiv.end(),NX);
	// import labeled mesh
	int labeled_domain=0;
// 	getfem::regular_unit_mesh(mesh, nsubdiv, pgt, 0);
        // getfem::import_mesh("gmsh:mesh/basin.msh",mesh);
        // getfem::import_mesh("gmsh:mesh/3dbasin_mult_dom.msh",mesh);
//         getfem::import_mesh("gmsh:mesh/3dbasin_dom_monomat.msh",mesh);
//         getfem::import_mesh("gmsh:mesh/patch_6_glued.msh",mesh);
// =============================================================
//         getfem::import_mesh("gmsh:mesh/patch_7light.msh",mesh);
//         getfem::import_mesh("gmsh:mesh/layer_cake/bounding.msh",mesh);
       //  getfem::import_mesh("gmsh:mesh/patch_6fp2.msh",mesh);
        /////////////////////////////////////////////
// 	getfem::import_mesh("gmsh:mesh/layer_cake/lk_fs.msh",mesh);labeled_domain=1; //official lk
// 	getfem::import_mesh("gmsh:mesh/ringmeshes/layer_cake.msh",mesh);labeled_domain=1; //official lk

	// 	//////////////////////////////////////////////
//      getfem::import_mesh("gmsh:mesh/pichout/patch_6.msh",mesh);
// 	getfem::import_mesh("gmsh:mesh/pinchout2/patch_6e.msh",mesh);
// 	getfem::import_mesh("gmsh:mesh/pinchout3/patch_7.msh",mesh);
	getfem::import_mesh("gmsh:mesh/ringmeshes/pinch_2.msh",mesh);labeled_domain=1; //official lk
// 	=============================================
// 	if(1){
// //  	 mesh.read_from_file("mesh/pinchout3/labeled_mesh_fp2"); //good pinchout
// //	mesh.read_from_file("mesh/pinchout2/labeled_mesh_fp2");
//  	mesh.read_from_file("mesh/pichout/labeled_mesh_fp2");
// 	mesh.read_from_file("mesh/labeled_mesh_fp2");// layar cake
//          labeled_domain=1;
// 	}
	//refinement
 	{
// 		// dal::bit_vector b; b.add(0);
//  		mesh.Bank_refine(mesh.convex_index());
 	}


	// A trasformation for the squarred mesh
	bgeot::base_matrix M(N_,N_);
	for (size_type i=0; i < N_; ++i) {
// 		M(i,i) = 4000.0;
		M(i,i) = 1.e+4;
	}
	//  if (N>1) { M(0,1) = 0; }
	M(0,0) = 1.;M(1,1) = 1.;M(2,2) = -1.; // 180degree rotation x for pichout ring
// 	mesh.transformation(M);
	//
	g_=new std::vector<double>{0., 0., -1.};
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

        mf_coef.set_finite_element(mesh.convex_index(),
        getfem::classical_fem(pgt,0));         // p0 for coefficient
 

	// boundary conditions zones
// 	if(!labeled_domain)
	  gen_bc();
        gmm::resize(Kr_, mf_coef.nb_dof()); gmm::fill(Kr_,1);    // rhs monolithic problem
        gmm::resize(Er_, mf_coef.nb_dof()); gmm::fill(Er_,1);    // rhs monolithic problem
        if(!labeled_domain) mesh_labeling();
        gen_coefficient();/// Generation of coefficient for different materials
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
	std::cout<<"number of dof "<< nb_dof_u<< " and "<<  nb_dof_p<<std::endl;
}
// assembly with ls

// ===========================================
// method for generation of bcs zones
// ===========================================
void biot_problem::gen_bc(){
	std::cout << "biot_problem::gen_bc()"<< std::endl;
	getfem::mesh_region border_faces;
	getfem::outer_faces_of_mesh(mesh, border_faces);

	for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
		assert(i.is_face());
		base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
		un /= gmm::vect_norm2(un);

		if ((un[N_-1] ) > 1.0E-1  && (mesh.points_of_convex(i.cv())[0])[2]>3500. ) { // new Neumann face
// 			if ((mesh.points_of_convex(i.cv())[0])[0]>5000. ) 
// 			  mesh.region(TOP).add(i.cv(), i.f());
// 			else
// 			  mesh.region(TOP_P).add(i.cv(), i.f());
			mesh.region(TOP_P).add(i.cv(), i.f()); //all top iced
		} else if ((un[N_-1] ) < -9.0E-1) {  //the bottom surface is the most sharp
			mesh.region(BOTTOM).add(i.cv(), i.f());
		} else if ((un[N_-2] ) < -1.0E-1) {
			mesh.region(LEFT).add(i.cv(), i.f());
		} else if ((un[N_-2] ) > 1.0E-1) {
		   mesh.region(RIGHT).add(i.cv(), i.f());
		}
		else if(N_=3){
			if ((un[N_-3] ) < -1.0E-1) {
			  if ((mesh.points_of_convex(i.cv())[0])[2]>2000. 
		               && 
		             (mesh.points_of_convex(i.cv())[0])[2]<2300.)
			mesh.region(RIGHTP).add(i.cv(), i.f());
				mesh.region(LEFTX).add(i.cv(), i.f());
			} else if ((un[N_-3] ) > 1.0E-1) {
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
void biot_problem::configure_workspace(getfem::ga_workspace & workspace,double dt){
	std::cout << "biot_problem::configure_workspace::Configuring workspace " << std::endl;
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
// 	int t1=10*2; int t2=20*2;
// 	int t3=30*2; int t4=35*2;
	
// //         std::vector<scalar_type> ice_force(1);ice_force[0] = 1.e+0;
// 	if(iter_<t1)       force_[0]= 0.;
// 	else if(iter_<t2)  force_[0]= (iter_ -((double) t1) )/(((double) t2)-((double) t1))
// 	                              *1000*9.81*1700;
// 	else if(iter_<t3)  force_[0]= 1000*9.81*1700;
// 	else if(iter_<t4)  force_[0]= (iter_ -((double) t4))/(((double) t3)-((double) t4))
// 	                              *1000*9.81*1700;
// 	else               force_[0]= 0.;
// 	workspace.add_fixed_size_constant("topload",force_);
	workspace.add_fixed_size_constant("gravity",*g_);
// 	
// 	if(iter_<t1)       overpres_[0]= 0.;
// 	else if(iter_<t2)  overpres_[0]= (iter_ -((double) t1))/(((double) t2)-((double) t1))* 1000*9.81*4000;
// 	else if(iter_<t3)  overpres_[0]= 1000*9.81*4000;
// 	else if(iter_<t4)  overpres_[0]= (iter_ -((double) t4))/(((double) t3)-((double) t4))*1000*9.81*4000;
// 	else               overpres_[0]= 0.;
// 	workspace.add_fixed_size_constant("overpres",overpres_);
//         std::cout<<"ssssssss iter is " <<iter_<< "  and force is "<<force_[0]<<std::endl;
	  //   ================ from Paleoisostesy =========================
  int t1=0; int t2=30;
  int t3=60; int t4=90;
//         std::vector<scalar_type> ice_force(1);ice_force[0] = 1.e+0;
  double p_buf=0.;
  if(iter_<t1)       p_buf= 0.;
  else if(iter_<t2)  p_buf= (iter_ -((double) t1) )/(((double) t2)-((double) t1))
                            *1000*9.81*4000;
  else if(iter_<t3)  p_buf= 1000*9.81*4000;
  else if(iter_<t4)  p_buf= (iter_ -((double) t4))/(((double) t3)-((double) t4))
                            *1000*9.81*4000;
  else               p_buf= 0.;
  overpres_[0]= p_buf; 
  std::cout<<"top load is "<<  p_buf/p_des.p_ref <<std::endl;
  workspace.add_fixed_size_constant("topload",overpres_);
  workspace.add_fixed_size_constant("overpres",overpres_);
	
	
	workspace.add_fem_constant("Kr", mf_coef, Kr_);
        workspace.add_fem_constant("Er", mf_coef, Er_);
	// workspace.add_fem_constant("f", mf_data, F);
}
// 


//=======================================================
// Build fix stress preconditioner 
// for the monolithic problem
//=======================================================
void biot_problem::build_fix_stress_preconditioner(double dt){
	std::cout<< "biot_problem::build_fix_stress_preconditioner" <<std::endl;
	assembly_p(dt);assembly_u(dt);
	bPR_ = new biot_precond<sparse_matrix_type> (Ku,Kp);
}
//=======================================================
// Assembly solution and iteration matrix for 
// the monolithic problem
//=======================================================
void biot_problem::assembly(double dt) {
	std::cout<< "biot_problem::assembly" <<std::endl;
	// assembly matrix
	// resize and clear vector and matrix
	getfem::size_type nb_dof_u = mf_u.nb_dof();
	getfem::size_type nb_dof_p = mf_p.nb_dof();
	// gmm::resize(B, nb_dof_u + nb_dof_p); 
	gmm::clear(B);
	// gmm::resize(UP, nb_dof_u + nb_dof_p); //gmm::clear(UP);
	// gmm::resize(U, nb_dof_u); // gmm::resize(U_old, nb_dof_u); 
	// gmm::resize(P, nb_dof_p); // gmm::resize(P_old, nb_dof_p);
	// gmm::resize(K, nb_dof_u + nb_dof_p, nb_dof_u + nb_dof_p); 
	gmm::clear(K);

	getfem::ga_workspace workspace; // generic workspace for the solution
	configure_workspace(workspace,dt); 
	// getfem::base_vector invdt(1); invdt[0] = 1/dt;
	// workspace.add_fixed_size_constant("invdt", invdt);
	// 

	workspace.add_fem_variable("u", mf_u, gmm::sub_interval(0, nb_dof_u), U);
	workspace.add_fem_variable("p", mf_p, gmm::sub_interval(nb_dof_u, nb_dof_p), P);
	workspace.add_fem_variable("u_old", mf_u, gmm::sub_interval(0, nb_dof_u), U_old);
	workspace.add_fem_variable("p_old", mf_p, gmm::sub_interval(nb_dof_u,nb_dof_p), P);
	// ------------------ expressions --------------------------
	workspace.add_expression("2*mu*Er*Sym(Grad_u):Grad_Test_u + lambda*Er*Div_u*Div_Test_u- alpha*p.Div_Test_u ", mim);
	workspace.add_expression( "+(1/bm)*p.Test_p + tau*Kr*permeability*Grad_p.Grad_Test_p+ alpha*Test_p*Div_u", mim);
	workspace.assembly(2);
	gmm::copy(workspace.assembled_matrix(), K);
	workspace.clear_expressions();


	//======= RHS =====================
	if(N_==2) workspace.add_expression("(2200*0.8 + 1000*0.2 -1000 )*[0,-1].Test_u", mim);
	if(N_==3) workspace.add_expression("(2200*0.8 + 1000*0.2 -1000 )*[0,0,-1].Test_u", mim);
	workspace.add_expression("+[+0.0e-6].Test_p + (1/bm)*p_old.Test_p + alpha*Test_p*Div_u_old", mim);
	workspace.set_assembled_vector(B);
	workspace.assembly(1);
	// gmm::copy(workspace.assembled_vector(),B);
	workspace.clear_expressions();

	//Boudanry conditions // NITSCHE
	// getfem::base_vector penalty(1); penalty[0] = 1e+42; // 1---10
	// workspace.add_fixed_size_constant("penalty", penalty);
	//Matrix term
	workspace.add_expression("penalty/element_size*p*Test_p", mim, TOP);
	workspace.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim, TOP); 	
	workspace.add_expression("penalty/element_size*p*Test_p", mim, LEFT);
	workspace.add_expression("penalty/element_size*p*Test_p", mim, RIGHT);// 1 is the region	
	workspace.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim, LEFT); 	
	workspace.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p", mim, RIGHT); 
	
	workspace.add_expression("penalty/element_size*p*Test_p", mim, LEFTX);
	workspace.add_expression("penalty/element_size*p*Test_p", mim, RIGHTX);// 1 is the region	
	if(N_==3){	 	
	workspace.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim, LEFTX); 	
	workspace.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p", mim, RIGHTX);
	}	
	workspace.add_expression("penalty*u.Test_u" , mim, BOTTOM); //neumann disp
	workspace.assembly(2);
	gmm::add(workspace.assembled_matrix(), K);
	workspace.clear_expressions();

	//rhs term
	if(N_ == 2){
		// workspace.add_expression("0*penalty*u.Test_u" "+ penalty*0*p*Test_p", mim, TOP); // 1 is the region
		// workspace.add_expression("[0,-2*force].Test_u" , mim, TOP); //neumann disp
		// workspace.add_expression("0*penalty*u.Test_u" "+ 0*penalty*0*p*Test_p", mim, BOTTOM); // 1 is the region
		// workspace.add_expression("[0,+2*force].Test_u" , mim, BOTTOM); //neumann disp
	}
	else if(N_== 3){
		// workspace.add_expression("0*penalty*u.Test_u" "+ penalty*0*p*Test_p", mim, TOP); // 1 is the region
	// 	workspace.add_expression("[0,0,-2*force].Test_u" , mim, TOP); //neumann disp
		// workspace.add_expression("0*penalty*u.Test_u" "+ 0*penalty*0*p*Test_p", mim, BOTTOM); // 1 is the region
	// 	workspace.add_expression("[0,0,+2*force].Test_u" , mim, BOTTOM); //neumann disp
	}

	workspace.add_expression(" 0*penalty/element_size*Test_p - permeability*Grad_Test_p.Normal*0 ", mim, LEFT);
	workspace.add_expression(" 0*penalty/element_size*Test_p - permeability*Grad_Test_p.Normal*0 ", mim, RIGHT);
	// workspace.assembly(1);
	// gmm::add(workspace.assembled_vector(),B);
	workspace.clear_expressions();
	/// end boundary contions
} // end assembly




// -------------------------------------------------
// Assembling pressure matrix for fixed stress
// -------------------------------------------------

void biot_problem::assembly_p(double dt){

	// gmm::clean(P_old, 1E-10);gmm::clean(U_old, 1E-10);
	gmm::clean(P, 1E-10);gmm::clean(U, 1E-10);
	gmm::clear(Bp); gmm::clear(Kp);
	std::cout<<"biot_assembler::assembly_p(double dt)" <<std::endl;
	getfem::size_type nb_dof_u = mf_u.nb_dof();
	getfem::size_type nb_dof_p = mf_p.nb_dof();

	//  gmm::resize(P, nb_dof_p); gmm::resize(P_old, nb_dof_p);gmm::resize(P_iter, nb_dof_p);
	//  gmm::resize(Kp, nb_dof_p, nb_dof_p);  gmm::resize(Bp, nb_dof_p);


	getfem::ga_workspace workspace;
	configure_workspace(workspace,dt);


	workspace.add_fem_variable("p", mf_p, gmm::sub_interval(0, nb_dof_p), P);
	workspace.add_fem_variable("p_old", mf_p, gmm::sub_interval(0,nb_dof_p), P_old);
	workspace.add_fem_variable("p_iter", mf_p, gmm::sub_interval(0,nb_dof_p), P_iter);
	workspace.add_fem_variable("u", mf_u, gmm::sub_interval(0, nb_dof_u), U);
	workspace.add_fem_variable("u_old", mf_u, gmm::sub_interval(0, nb_dof_u), U_old);
	workspace.add_fem_variable("u_iter", mf_u, gmm::sub_interval(0, nb_dof_u), U_iter);

	// Pressure equation
	workspace.add_expression("(1/bm + beta)*p.Test_p + tau*permeability*Kr*Grad_p.Grad_Test_p", mim); // tau
	// workspace.set_assembled_matrix(Kp);
	workspace.assembly(2);
	gmm::copy(workspace.assembled_matrix(), Kp);
	workspace.clear_expressions();

	//======= RHS =====================
	//  workspace.add_expression("[+0.0].Test_p + (1/bm)*invdt*p_old.Test_p + invdt*alpha*Test_p*Trace((Grad_u_old))", mim);// 1/dt
	//  workspace.add_expression("(beta)*invdt*p_iter.Test_p - invdt*alpha*Test_p*Trace((Grad_u_iter))", mim);// 1/dt
	// tau
	workspace.add_expression("[+0.0e-6]*Test_p + (1/bm)*p_old.Test_p + alpha*Test_p*Trace(Sym(Grad_u_old))", mim); // tau
	workspace.add_expression("(beta)*p_iter.Test_p - alpha*Test_p*Trace(Sym(Grad_u_iter))", mim); // tau
	workspace.set_assembled_vector(Bp);
	workspace.assembly(1);
	workspace.clear_expressions();


	//Matrix term
	workspace.add_expression("penalty/element_size*p*Test_p", mim, TOP);
	workspace.add_expression("-permeability*Kr*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim, TOP); 
	workspace.add_expression("penalty/element_size*p*Test_p", mim, TOP_P);
	workspace.add_expression("-permeability*Kr*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim, TOP_P); 
	
	
#ifdef LATERAL_INJECTION
	
	workspace.add_expression("penalty*p*Test_p", mim, RIGHTP);
	workspace.add_expression("-permeability*Kr*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim, RIGHTP); 	
// 	
#endif
	
	
	
	
	
	
// 	workspace.add_expression("penalty/element_size*p*Test_p", mim, LEFT);
// 	workspace.add_expression("-permeability*Kr*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim, LEFT); 	
// 	workspace.add_expression("penalty/element_size*p*Test_p", mim, RIGHT);
// 	workspace.add_expression("-permeability*Kr*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim, RIGHT);
// 	if(N_==3){	 	
// 		workspace.add_expression("penalty/element_size*p*Test_p", mim, LEFTX);
// 		workspace.add_expression("-permeability*Kr*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim, LEFTX);
// 		workspace.add_expression("penalty/element_size*p*Test_p", mim, RIGHTX);
// 		workspace.add_expression("-permeability*Kr*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim, RIGHTX);    
// 	}

	workspace.assembly(2);
	gmm::add(workspace.assembled_matrix(), Kp);
	workspace.clear_expressions();
	//rhs term
	workspace.add_expression("0*penalty*p*Test_p", mim, TOP); 
#ifdef LATERAL_INJECTION
	workspace.add_expression("0.1*overpres*penalty*Test_p", mim, RIGHTP); 
#endif
// 	workspace.add_expression("0.4*1000*9.81*4000*penalty*Test_p", mim, TOP_P); 
	workspace.assembly(1);
	workspace.clear_expressions();
	return;}



	// -------------------------------------------------
	// Assembling displacement matrix for fixed stress
	// -------------------------------------------------
void biot_problem::assembly_u(double dt){
		std::cout<<"biot_assembler::assembly_u(double dt)" <<std::endl;

		getfem::size_type nb_dof_u = mf_u.nb_dof();
		getfem::size_type nb_dof_p = mf_p.nb_dof();
		gmm::clear(Bu);  gmm::clear(Ku);
		// gmm::resize(P, nb_dof_p); gmm::resize(P_old, nb_dof_p);gmm::resize(P_iter, nb_dof_p);
		// gmm::resize(Ku, nb_dof_u, nb_dof_u);  gmm::resize(Bu, nb_dof_u);
		getfem::ga_workspace workspace; configure_workspace(workspace,dt);

		workspace.add_fem_variable("p", mf_p, gmm::sub_interval(0, nb_dof_p), P);
		workspace.add_fem_variable("p_old", mf_p, gmm::sub_interval(0,nb_dof_p), P_old);
		workspace.add_fem_variable("p_iter", mf_p, gmm::sub_interval(0,nb_dof_p), P_iter);
		workspace.add_fem_variable("u", mf_u, gmm::sub_interval(0, nb_dof_u), U);
		workspace.add_fem_variable("u_old", mf_u, gmm::sub_interval(0, nb_dof_u), U_old);
		workspace.add_fem_variable("u_iter", mf_u, gmm::sub_interval(0, nb_dof_u), U_iter);

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
		workspace.add_expression("2*mu*Er*Sym(Grad_u):Grad_Test_u + lambda*Er*Div_u*Div_Test_u", mim); // stress tensor 
		workspace.add_expression("penalty*u.Test_u" , mim, BOTTOM); //neumann disp
		workspace.add_expression("penalty*u(2).Test_u(2)" , mim, LEFT); //neumann disp
		workspace.add_expression("penalty*u(2).Test_u(2)" , mim, RIGHT); //neumann disp
		workspace.add_expression("penalty*u(1).Test_u(1)" , mim, LEFTX); //neumann disp
		workspace.add_expression("penalty*u(1).Test_u(1)" , mim, RIGHTX); //neumann disp
// 		workspace.add_expression("penalty*u.Test_u" , mim, TOP_P); //Height og the ice disp
		workspace.assembly(2);
		gmm::copy(workspace.assembled_matrix(), Ku);
		workspace.clear_expressions();
		if(N_==2) workspace.add_expression("(2200*0.8 + 1000*0.2 -1000 )*[0,-1].Test_u", mim);
		if(N_==3) workspace.add_expression("(2200*0.8 + 1000*0.2 -1000 )*gravity.Test_u", mim);
		workspace.add_expression("alpha*p*Div_Test_u ", mim);
		if(N_==3)workspace.add_expression("topload*gravity.Test_u" , mim, TOP_P); //Height og the ice disp
#ifdef LATERAL_INJECTION
		if(N_==3)workspace.add_expression("overpres*[1,0,0].Test_u" , mim, RIGHTP); //Height og the ice disp
#endif
		workspace.assembly(1);
		workspace.clear_expressions();
		// std::cout<< Bu<< std::endl; 


	} // end assembling momentum matrix for fixed stress approach

//====================================================
// method for solving the laplace problem
//====================================================
void biot_problem::solve(){

	// biot_precond<sparse_matrix_type> bPR(Ku,Kp);
	size_type restart = 50;

	gmm::identity_matrix PM; // no precond
	gmm::diagonal_precond<sparse_matrix_type> PR(K);
	// gmm::mr_approx_inverse_precond<sparse_matrix_type> P(K, 10, 10E-17);
	gmm::iteration iter(1.e-8);  // iteration object with the max residu
	iter.set_noisy(1);               // output of iterations (2: sub-iteration)
	iter.set_maxiter(1000); // maximum number of iterations
	gmm::gmres(K, UP, B, PR, restart, iter);
	// gmm::gmres(K, UP, B, PM, restart, iter);
	scalar_type cond;
	// gmm::SuperLU_solve(K, UP , B, cond);
	std::cout<<"condition number "<< cond<< std::endl;
	getfem::size_type nb_dof_u = mf_u.nb_dof();
	getfem::size_type nb_dof_p = mf_p.nb_dof();
	gmm::copy(gmm::sub_vector(UP,gmm::sub_interval(0, nb_dof_u)), U);				
	gmm::copy(gmm::sub_vector(UP,gmm::sub_interval(nb_dof_u , nb_dof_p)), P);	
	gmm::copy(U,U_old);gmm::copy(P,P_old);
}

//====================================================
// method for solving the fixed stress  problem
//====================================================
void biot_problem::solve_fix_stress(double dt, int max_iter){


	double epsu=1.e-6; double epsp=1.e-6;
	double rel_unorm=1; double rel_pnorm=1; int fix_count=0;
	double old_unorm=1; double new_unorm=1;
	double old_pnorm=1; double new_pnorm=1;
        assembly_u(dt);
        assembly_p(dt);
	AMG amg_("momentum");
	gmm::csr_matrix<scalar_type> ku_csr;
//      gmm::clean(ku_csr, 1E-12);
        gmm::copy(Ku, ku_csr);
	amg_.convert_matrix(ku_csr);
	std::cout<<"End samg conversion"<<std::endl;
	
	AMG amg_p_("Pressure");
	gmm::csr_matrix<scalar_type> kp_csr;
//      gmm::clean(ku_csr, 1E-12);
        gmm::copy(Kp, kp_csr);
	amg_p_.convert_matrix(kp_csr);
	std::cout<<"End samg conversion"<<std::endl;
	while(fix_count < max_iter && ( rel_unorm>epsu ||  rel_pnorm > epsp) )
	{
		fix_count++; 
		std::cout<<"\033[1;34m***** iteration " << fix_count 
			<< " norm p " <<  rel_pnorm 
			<< " norm u " <<  rel_unorm << std::endl;

		assembly_p(dt);
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
			gmm::gmres(Kp, P, Bp, PRp, restart, iter);
// 			amg_p_.solve(kp_csr, P , Bp , 1);
// 			gmm::copy(amg_p_.getsol(), P);

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
	// 	gmm::ilu_precond<sparse_matrix_type> PRu(Ku); 
		{
			size_type restart = 50;
			gmm::iteration iter(1.e-8);  // iteration object with the max residu
			iter.set_noisy(1);               // output of iterations (2: sub-iteration)
			iter.set_maxiter(1000); // maximum number of iterations
			// gmm::MatrixMarket_load("km",Ku);
			// gmm::clear(U);
			gmm::gmres(Ku, U, Bu, PRu, restart, iter);
		        
	                scalar_type cond;
// 			amg_.solve(ku_csr, U , Bu , 1);
// 			gmm::copy(amg_.getsol(), U);
// 			gmm::SuperLU_solve(Ku, U , Bu, cond);
	std::cout<<"condition number "<< cond<< std::endl;
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
	}
	std::cout<<"\033[1;34m***** last iteration " << fix_count 
		<< " norm p " <<  rel_pnorm 
		<< " norm u " <<  rel_unorm << std::endl;
	gmm::copy(U,U_old); gmm::copy(P,P_old);
	gmm::clear(U_iter); gmm::clear(P_iter);

}



// method for solving the laplace problem
void biot_problem::print(int time){
	std::cout<< "biot_problem::print  "<< 
		p_des.datafilename + "." +  std::to_string(time) + ".vtk"<<std::endl;
	getfem::vtk_export exp(p_des.datafilename + "." +  std::to_string(time) + ".vtk");
	exp.exporting(mf_u);  	exp.write_point_data(mf_u, U, "u"); 
	exp.write_point_data(mf_p, P, "p"); 
// 	getfem::vtk_export exp_data(p_des.datafilename + ".data." +  std::to_string(time) + ".vtk");
// 	exp_data.exporting(mf_coef);  	exp_data.write_cell_data(mf_coef, Er_, "E"); 

{ // Just to see what elements are cut by the level set ls:
    getfem::vtk_export vtk_data("data.vtk");
    vtk_data.exporting(mf_u);
    vtk_data.write_mesh();
    vtk_data.write_cell_data(Kr_, "K");
  }


}



//============================================
// routine for the generation of coeffient
//============================================
void biot_problem::gen_coefficient(){ // creating a coefficient
  std::cout<<"Generating materials coefficients"<<std::endl;
  std::cout<<"Dof num  "<< mf_coef.nb_dof()<<std::endl;
	
  gmm::resize(Kr_, mf_coef.nb_dof()); gmm::fill(Kr_,1);    // rhs monolithic problem
  gmm::resize(Er_, mf_coef.nb_dof()); gmm::fill(Er_,1);    // rhs monolithic problem
  gmm::resize(Kr_print_, mf_coef.nb_dof()); gmm::fill(Kr_print_,1);    // rhs monolithic problem
  gmm::resize(Er_print_, mf_coef.nb_dof()); gmm::fill(Er_print_,1);    // rhs monolithic problem
  std::vector<int> material; 
  // material.push_back(1);material.push_back(5);material.push_back(3); // for layerckae
  // material.push_back(1);material.push_back(2);material.push_back(3); // for ring_pinch2
   material.push_back(0);material.push_back(1);material.push_back(2); // for ring_pinch2
  // std::vector<double> k; k.push_back(1.e-0);k.push_back(1.e-4);k.push_back(1.e-0);
  // std::vector<double> E; E.push_back(1.e+0);E.push_back(1.e+0);E.push_back(1.e+0);
   /////////////////////////////////////////////////////////////////////////////////////////
   std::vector<double> k; k.push_back(1.e-0);k.push_back(1.e+3);k.push_back(1.e-2); // pinch trimat
   std::vector<double> E; E.push_back(1.e+0);E.push_back(2.e+0);E.push_back(1.e+1);
  //////////////////ring mesh layer cake 
//    std::vector<double> k;k.push_back(1.e+3);k.push_back(1.e-2); k.push_back(1.e-0); // 
//    std::vector<double> E;E.push_back(2.e+0);E.push_back(1.e+1); E.push_back(1.e+0);
   ///////////////////////////////////////////////////
     //////////////////mesh cgal cake
//   std::vector<double> k;k.push_back(1.e-2);k.push_back(1.e-0);k.push_back(1.e+3);  // 
//    std::vector<double> E;E.push_back(1.e+1);E.push_back(1.e+0);E.push_back(2.e+0); 
   ///////////////////////////////////////////////////
//   std::vector<double> k; k.push_back(1);k.push_back(1.e+0);
//   std::vector<double> E; E.push_back(1);E.push_back(1.e+0);
/////////////////layercacke from cgal
  for (int imat=0; imat< material.size();imat++)
//   int imat=2;
  {
   
    dal::bit_vector bv_cv = mesh.region(material[imat]).index();
    size_type i_cv = 0;
    for (i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv) {
	    getfem::mesh_fem::ind_dof_ct idofs = mf_coef.ind_basic_dof_of_element(i_cv);
 	 //    std::cout<<"Material  "<<material[imat] << " with "<< idofs.size()<<std::endl;
      for (size_type i=0; i < idofs.size(); ++i) {
        Kr_[idofs[i]]=k[imat]; Kr_print_[(int) i_cv]=k[imat];
        Er_[idofs[i]]=E[imat];// Er_print_[i_cv]=E[imat];
      }
    }
  }
  if(1){ // Just to see what elements are cut by the level set ls:
    // getfem::vtk_export vtk_data("ring_data_gen_3mat_pinch.vtk");
    getfem::vtk_export vtk_data("pinch_cgal_3mat.vtk");
    vtk_data.exporting(mf_coef);
    vtk_data.write_mesh();
    vtk_data.write_cell_data(Kr_print_, "K");
  }
}



//============================================
// routine for the labeling mesh
//============================================
void biot_problem::mesh_labeling(){
//    horizon h  ("mesh/layer_cake/mesh_horizon_2.msh");
//   horizon h3 ("mesh/layer_cake/mesh_horizon_3.msh");
//=======================================================
//   horizon h  ("mesh/pichout/pinched_horizon_2.msh");
//   horizon h3 ("mesh/pichout/pinched_horizon_3.msh");
//=======================================================
//   horizon h  ("mesh/pinchout2/mesh_pinched_horizon_1.msh");
//=======================================================
//=============pinchout====================
  horizon h3  ("mesh/pinchout3/pinched_horizon_0.msh");
  horizon h ("mesh/pinchout3/pinched_crop_large_horizon_1.msh");
//=======================================================
   dal::bit_vector bv_cv = mesh.convex_index();
   size_type i_cv = 0;
     std::cout<<"Looping internal element "<< i_cv <<std::endl;
   for (i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv) {
     int dof_el=mesh.structure_of_convex(i_cv)->nb_points();
//      std::cout<< mf_coef.ind_basic_dof_of_element(i_cv)[0]<<" <-> "<< i_cv<<std::endl;
       const std::vector<long unsigned int> indicies(
// 	 mesh.ind_points_of_convex(mf_coef.ind_basic_dof_of_element(i_cv)[0]) // wrong inside correct visualization
	 mesh.ind_points_of_convex(i_cv)                                      // correct inside wrong visualization
       );
	  std::vector<double> bc={0.,0.,0.};
          for(int idof=0;idof < dof_el; idof++ ){
	    for(int idir=0; idir<3; idir++)
	       bc[idir] += mesh.points()[indicies[idof]][idir]/dof_el;
    }
	    int el_on_h=h.find_element(bc);
// 	    std::cout<< "barcenter at " << bc[0]<< " "<<bc[1]<< " "<<bc[2]<< 
// 	    "and corresponding element is "<< el_on_h << std::endl;
	    int group =  h.up_down(bc,el_on_h); //1 under 2 above
	    std::cout<< "Group h1 is "<<  group<<std::endl;
	    
	    int el_on_h3=h3.find_element(bc);
	    int group_h3=  h3.up_down(bc,el_on_h3);
	    std::cout<< "Group h2 is "<<  group_h3<<std::endl;
	    //////////////////////////////////
// 	    if (group ==1 ) mesh.region(1).add(i_cv);
// 	    if (group ==2 ) mesh.region(2).add(i_cv);
// 	    int group_h3=2;
// 	    mesh.region(group).add(i_cv);
	    /////////////////////////////////
            if (group_h3==2 && group ==2) mesh.region(1).add(i_cv);
	    if (group_h3==1 && group ==2) mesh.region(2).add(i_cv);
	    if (group_h3==1 && group ==1) mesh.region(3).add(i_cv);
// // 	    if (mesh.points()[indicies[0]][2]<2000) mesh.region(1).add(i_cv);
//        std::cout<< "zed is " <<mesh.points()[mesh.ind_points_of_convex(i_cv)[0]][2]<<std::endl;
// 	    if (mesh.points()[indicies[0]][2]>2000) mesh.region(14).add(i_cv);
   }
//    mesh.write_to_file("mesh/pichout/labeled_mesh_fp2");
   mesh.write_to_file("mesh/pinchout3/labeled_mesh_fp2");
}

// Method for printing auxiliar data
void biot_problem::print_aux_data(int istep){
 std::cout<<"Start printing auxiliar data function"<< std::endl;
 bgeot::pgeometric_trans pgt = 
 bgeot::geometric_trans_descriptor("GT_PK(2,1)");
 getfem::mesh mesh_htop;
 // getfem::import_mesh("gmsh:mesh/layer_cake/horizons/horizon_4.msh",mesh_htop);// layer_cake
 getfem::mesh_region &rg = mesh.region(TOP_P);
 for (getfem::mr_visitor i(rg); !i.finished(); ++i) {
   std::vector<bgeot::base_node> nodel;
if (i.is_face()) {
  for(int ip=0;ip<mesh.structure_of_convex(i.cv())->ind_points_of_face(i.f()).size(); ip++){
  nodel.push_back(mesh.points_of_convex(i.cv())
       [mesh.structure_of_convex(i.cv())->ind_points_of_face(i.f())[ip]]);
  }
}
 mesh_htop.add_convex_by_points(bgeot::simplex_geotrans(2,1),nodel.begin());
 }
  getfem::mesh_fem mf_top(mesh_htop);    /// the main mesh_fem, for the pressure solution
  mf_top.set_finite_element(mesh_htop.convex_index(),
      getfem::classical_fem(pgt,1)); 
  
  
  
  {
  std::vector<scalar_type> over_p; // permeability ratio
  gmm::resize(over_p, mf_p.nb_dof()); gmm::fill(over_p,overpres_[0]);    // rhs monolithic problem
  std::string namefile= p_des.datafilename +".aux." +  std::to_string(istep) +".vtk";
  getfem::vtk_export vtkd(namefile);
  vtkd.exporting(mf_p);vtkd.write_mesh();
  vtkd.write_point_data(mf_p, over_p, "h_ice");
  vtkd.write_cell_data(Kr_print_, "kr");
  }  
  {
  std::vector<scalar_type> over_p; // permeability ratio
  gmm::resize(over_p, mf_top.nb_dof()); gmm::fill(over_p,overpres_[0]);
  
  std::string namefile= p_des.datafilename +".topaux." +  std::to_string(istep) +".vtk";
  getfem::vtk_export vtkd(namefile);
  vtkd.exporting(mf_top);vtkd.write_mesh();
  vtkd.write_point_data(mf_top, over_p, "h_ice");
  }
  std::cout<<"end printing auxiliar function function"<< std::endl;
}

