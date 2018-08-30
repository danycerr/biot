#include "temp.hpp" 

void temperature_problem::init(void) {

	std::cout<< "temperature_problem::init "  << std::endl;

	bgeot::pgeometric_trans pgt = 
		bgeot::geometric_trans_descriptor(p_des.MESH_TYPE);

	// Mesh generation
	N_ = pgt->dim();
	std::vector<size_type> nsubdiv(N_);
	int NX=p_des.nsubdiv;
	std::fill(nsubdiv.begin(),nsubdiv.end(),NX);
	// import labeled mesh
	int labeled_domain=0;
	getfem::regular_unit_mesh(mesh, nsubdiv, pgt, 0);
        // getfem::import_mesh("gmsh:mesh/basin.msh",mesh);
        // getfem::import_mesh("gmsh:mesh/3dbasin_mult_dom.msh",mesh);
//         getfem::import_mesh("gmsh:mesh/3dbasin_dom_monomat.msh",mesh);
//         getfem::import_mesh("gmsh:mesh/patch_6_glued.msh",mesh);
// =============================================================
//         getfem::import_mesh("gmsh:mesh/patch_7light.msh",mesh);
//         getfem::import_mesh("gmsh:mesh/layer_cake/bounding.msh",mesh);
//         getfem::import_mesh("gmsh:mesh/patch_6fp.msh",mesh);
// 	 getfem::import_mesh("gmsh:mesh/pichout/patch_6.msh",mesh);
// 	=============================================
	if(1){
// 	mesh.read_from_file("mesh/pinchout3/labeled_mesh_fp2");
// 	mesh.read_from_file("mesh/labeled_mesh_fp2");// layar cake
// 	mesh.read_from_file("mesh/pichout/labeled_mesh_fp2");
// 	mesh.read_from_file("mesh/pinchout2/labeled_mesh_fp2");
	getfem::import_mesh("gmsh:mesh/layer_cake/lk_fs.msh",mesh);labeled_domain=1;
        labeled_domain=1;
	}
	//refinement
// 	{
// 		// dal::bit_vector b; b.add(0);
// 		mesh.Bank_refin(mesh.convex_index());
// 	}


	// A trasformation for the squarred mesh
	bgeot::base_matrix M(N_,N_);
	for (size_type i=0; i < N_; ++i) {
// 		M(i,i) = 4000.0;
		M(i,i) = 1.e+4;
	}
	//  if (N>1) { M(0,1) = 0; }
	//
// 	mesh.transformation(M);
	// // End of mesh generation


	// set  integration methods  
	getfem::pfem pf_u = getfem::fem_descriptor(p_des.FEM_TYPE_U);
	getfem::pintegration_method ppi = getfem::int_method_descriptor(p_des.INTEGRATION);
	getfem::pintegration_method simp_ppi = getfem::int_method_descriptor(p_des.SIMPLEX_INTEGRATION);

	mim.set_integration_method(mesh.convex_index(), ppi);
	mf_u.set_finite_element(mesh.convex_index(), pf_u); // finite element for displacement

        mf_coef.set_finite_element(mesh.convex_index(),
        getfem::classical_fem(pgt,0));         // p0 for coefficient
 

	// boundary conditions zones
// 	if(!labeled_domain)
	  gen_bc();
        gmm::resize(Kr_, mf_coef.nb_dof()); gmm::fill(Kr_,1);    // rhs monolithic problem
        gmm::resize(Er_, mf_coef.nb_dof()); gmm::fill(Er_,1);    // rhs monolithic problem
//         if(!labeled_domain) 	  mesh_labeling();
        gen_coefficient();/// Generation of coefficient for different materials
	// init vector
	getfem::size_type nb_dof_u = mf_u.nb_dof();

	gmm::resize(B, nb_dof_u); gmm::clear(B);    // rhs monolithic problem
	// displacement
	gmm::resize(U, nb_dof_u); gmm::resize(U_old, nb_dof_u); gmm::resize(U_iter, nb_dof_u); 
	gmm::clear(U);gmm::clear(U_old);gmm::clear(U_iter);
	std::fill(U.begin(), U.end(), 0.);
	std::fill(U_old.begin(), U_old.end(), 10.);
	// iteration matrix monolithic
	gmm::resize(K, nb_dof_u , nb_dof_u ); gmm::clear(K);
	std::cout<<"number of dof "<< nb_dof_u<<std::endl;
}
// assembly with ls

// ===========================================
// method for generation of bcs zones
// ===========================================
void temperature_problem::gen_bc(){
	std::cout << "temperature_problem::gen_bc()"<< std::endl;
	getfem::mesh_region border_faces;
	getfem::outer_faces_of_mesh(mesh, border_faces);

	for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
		assert(i.is_face());
		base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
		un /= gmm::vect_norm2(un);

		if ((un[N_-1] ) > 4.0E-1 && (mesh.points_of_convex(i.cv())[0])[2]>3500) { // new Neumann face
			if ((mesh.points_of_convex(i.cv())[0])[0]>5000. )
			  mesh.region(TOP).add(i.cv(), i.f());
			else
			  mesh.region(TOP_P).add(i.cv(), i.f());
		} 
// 			else if ((un[N_-1] ) < -9.0E-1) {  //the bottom surface is the most sharp
// 			mesh.region(BOTTOM).add(i.cv(), i.f());
// 		} else if ((un[N_-2] ) < -1.0E-1) {
// 			mesh.region(LEFT).add(i.cv(), i.f());
// 		} else if ((un[N_-2] ) > 1.0E-1) {
// 			mesh.region(RIGHT).add(i.cv(), i.f());
// 		}
// 		else if(N_=3){
// 			if ((un[N_-3] ) < -1.0E-1) {
// 				mesh.region(LEFTX).add(i.cv(), i.f());
// 			} else if ((un[N_-3] ) > 1.0E-1) {
// 				mesh.region(RIGHTX).add(i.cv(), i.f());
// 			}
// 		}
// 		else {
// 			mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
// 		}
	}
}
// end bc generations
// =======================================================

//=======================================================
//Configure Workspace for all methods
//=======================================================
void temperature_problem::configure_workspace(getfem::ga_workspace & workspace,double dt){
	std::cout << "temperature_problem::configure_workspace::Configuring workspace " << std::endl;
	tau_[0] = dt; // dt into  the workspace
	workspace.add_fixed_size_constant("tau", tau_);

	vmu_[0] = p_des.diff;
	workspace.add_fixed_size_constant("diffusivity", vmu_);
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

	penalty_[0] = 1.e+5; // 1---10
	workspace.add_fixed_size_constant("penalty", penalty_);
	
	if(iter_<10)       temp_bc_[0]= 10.;
	else if(iter_<20)  temp_bc_[0]= (iter_ -20.)/(10.-20.)*10;
	else if(iter_<30)  temp_bc_[0]= 0;
	else if(iter_<35)  temp_bc_[0]= (iter_ -30.)/(35.-30.)*10;
	else               temp_bc_[0]= 10.;
	std::cout << "Temperature is " <<temp_bc_[0] << std::endl;
	workspace.add_fixed_size_constant("temp_bc",temp_bc_);

  workspace.add_fem_constant("Kr", mf_coef, Kr_);
  workspace.add_fem_constant("lr", mf_coef, Er_);
	// workspace.add_fem_constant("f", mf_data, F);
}
// 


//=======================================================
// Assembly solution and iteration matrix for 
// the monolithic problem
//=======================================================
void temperature_problem::assembly(double dt, getfem::mesh_fem& mf_pressure, std::vector<scalar_type>& p) {
	std::cout<< "temperature_problem::assembly" <<std::endl;
	// assembly matrix
	// resize and clear vector and matrix
	getfem::size_type nb_dof_u = mf_u.nb_dof();
	// gmm::resize(B, nb_dof_u + nb_dof_p); 
	gmm::clear(B);
	gmm::clear(K);

	getfem::ga_workspace workspace; // generic workspace for the solution
	configure_workspace(workspace,dt); 

	workspace.add_fem_variable("T", mf_u, gmm::sub_interval(0, nb_dof_u), U);
	workspace.add_fem_variable("T_old", mf_u, gmm::sub_interval(0, nb_dof_u), U_old);
	workspace.add_fem_variable("p", mf_u, gmm::sub_interval(0, nb_dof_u), p);
	// ------------------ expressions --------------------------
// 	workspace.add_expression("tau*Kr*permeability*Grad_T.Grad_Test_T", mim);
	workspace.add_expression("tau*lr*diffusivity*Grad_T.Grad_Test_T", mim);
	workspace.add_expression( "T.Test_T - tau*Kr*permeability*Grad_p.Grad_Test_T", mim);
	workspace.add_expression("tau*penalty*T*Test_T", mim, TOP);
	workspace.add_expression("tau*penalty*T*Test_T", mim, TOP_P);
	workspace.assembly(2);
	gmm::copy(workspace.assembled_matrix(), K);
	workspace.clear_expressions();


	//======= RHS =====================
	
	

	//Boudanry conditions // NICHE
	// getfem::base_vector penalty(1); penalty[0] = 1e+42; // 1---10
	// workspace.add_fixed_size_constant("penalty", penalty);
	//Matrix term
	
// 	workspace.add_expression("-permeability*Grad_T.Normal*Test_T - permeability*Grad_Test_T.Normal*T ", mim, TOP); 	
// 	workspace.add_expression("penalty*u.Test_u" , mim, BOTTOM); //neumann disp


	//rhs term
// 	workspace.add_expression(" 0*penalty/element_size*Test_p - permeability*Grad_Test_p.Normal*0 ", mim, LEFT);
// 	workspace.add_expression(" 0*penalty/element_size*Test_p - permeability*Grad_Test_p.Normal*0 ", mim, RIGHT);
	workspace.set_assembled_vector(B);
	workspace.add_expression(" tau*10*penalty*Test_T", mim, TOP);
	workspace.add_expression(" tau*temp_bc*penalty*Test_T", mim, TOP_P);
	workspace.add_expression("+tau*[+1.0e-5].Test_T + T_old.Test_T", mim);
	workspace.assembly(1);
// // 	gmm::add(workspace.assembled_vector(),B);
	workspace.clear_expressions();
	/// end boundary contions
} // end assembly

//====================================================
// method for solving the laplace problem
//====================================================
void temperature_problem::solve(){

	// biot_precond<sparse_matrix_type> bPR(Ku,Kp);
	size_type restart = 50;

	gmm::identity_matrix PM; // no precond
	gmm::diagonal_precond<sparse_matrix_type> PR(K);
	// gmm::mr_approx_inverse_precond<sparse_matrix_type> P(K, 10, 10E-17);
	gmm::iteration iter(1.e-8);  // iteration object with the max residu
	iter.set_noisy(1);               // output of iterations (2: sub-iteration)
	iter.set_maxiter(1000); // maximum number of iterations
	gmm::gmres(K, U, B, PR, restart, iter);
	// gmm::gmres(K, UP, B, PM, restart, iter);
	scalar_type cond;
	// gmm::SuperLU_solve(K, UP , B, cond);
	std::cout<<"condition number "<< cond<< std::endl;	
	gmm::copy(U,U_old);
}


// method for solving the laplace problem
void temperature_problem::print(int time){
	std::cout<< "temperature_problem::print  "<< 
		p_des.datafilename + "." +  std::to_string(time) + ".vtk"<<std::endl;
	getfem::vtk_export exp(p_des.datafilename + "." +  std::to_string(time) + ".vtk");
	exp.exporting(mf_u);  	exp.write_point_data(mf_u, U, "T"); 
// 	exp.write_point_data(mf_p, P, "p"); 
// 	getfem::vtk_export exp_data(p_des.datafilename + ".data." +  std::to_string(time) + ".vtk");
// 	exp_data.exporting(mf_coef);  	exp_data.write_cell_data(mf_coef, Er_, "E"); 

// { // Just to see what elements are cut by the level set ls:
//     getfem::vtk_export vtk_data("data.vtk");
//     vtk_data.exporting(mf_u);
//     vtk_data.write_mesh();
//     vtk_data.write_cell_data(Kr_, "K");
//   }


}

//============================================
// routine for the generation of coeffient
//============================================
void temperature_problem::gen_coefficient(){ // creating a coefficient
  std::cout<<"Generating materials coefficients"<<std::endl;
  std::cout<<"Dof num  "<< mf_coef.nb_dof()<<std::endl;
	
  gmm::resize(Kr_, mf_coef.nb_dof()); gmm::fill(Kr_,1);    // rhs monolithic problem
  gmm::resize(Er_, mf_coef.nb_dof()); gmm::fill(Er_,1);    // rhs monolithic problem
  gmm::resize(Kr_print_, mf_coef.nb_dof()); gmm::fill(Kr_print_,1);    // rhs monolithic problem
  gmm::resize(Er_print_, mf_coef.nb_dof()); gmm::fill(Er_print_,1);    // rhs monolithic problem
  std::vector<int> material; 
//   material.push_back(1);material.push_back(2);material.push_back(3);
  material.push_back(1);material.push_back(5);material.push_back(3);
  std::vector<double> k; k.push_back(1.e-0);k.push_back(1.e+3);k.push_back(1.e-2);
  std::vector<double> E; E.push_back(1.2e+0);E.push_back(2.e+0);E.push_back(1.1e+0);
  
//   std::vector<double> k; k.push_back(1);k.push_back(1.e+0);
//   std::vector<double> E; E.push_back(1);E.push_back(1.e+0);
  for (int imat=0; imat< material.size();imat++){
    dal::bit_vector bv_cv = mesh.region(material[imat]).index();
    size_type i_cv = 0;
    for (i_cv << bv_cv; i_cv != size_type(-1); i_cv << bv_cv) {
	  //   std::cout<<"Material  "<<material[imat]<<std::endl;
	    getfem::mesh_fem::ind_dof_ct idofs = mf_coef.ind_basic_dof_of_element(i_cv);
      for (size_type i=0; i < idofs.size(); ++i) {
        Kr_[idofs[i]]=k[imat];// Kr_print_[i_cv]=k[imat];
        Er_[idofs[i]]=E[imat];// Er_print_[i_cv]=E[imat];
      }
    }
  }
     if(0){ // Just to see what elements are cut by the level set ls:
    getfem::vtk_export vtk_data("data_gen_3mat.vtk");
    vtk_data.exporting(mf_coef);
//     vtk_data.write_mesh();
    vtk_data.write_cell_data(Kr_print_, "K");
  }
}



//============================================
// routine for the labeling mesh
//============================================
void temperature_problem::mesh_labeling(){
  horizon h  ("mesh/layer_cake/mesh_horizon_2.msh");
  horizon h3 ("mesh/layer_cake/mesh_horizon_3.msh");
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
	    
// 	    mesh.region(group).add(i_cv);
	    
            if (group_h3==2 && group ==2) mesh.region(1).add(i_cv);
	    if (group_h3==1 && group ==2) mesh.region(2).add(i_cv);
	    if (group_h3==1 && group ==1) mesh.region(3).add(i_cv);
// // 	    if (mesh.points()[indicies[0]][2]<2000) mesh.region(1).add(i_cv);
//        std::cout<< indicies[0] <<" "<< indicies[1] <<" "<< indicies[2] <<" "<< indicies[3] <<std::endl;
//        std::cout<< "zed is " <<mesh.points()[mesh.ind_points_of_convex(i_cv)[0]][2]<<std::endl;
// 	    if (mesh.points()[indicies[0]][2]>2000) mesh.region(14).add(i_cv);
   }
   mesh.write_to_file("mesh/labeled_mesh_fp");
   
}
