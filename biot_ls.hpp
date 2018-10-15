#ifndef BIOTLS_H
#define BIOTLS_H
/*===========================================================================

  Copyright (C) 2002-2016 Yves Renard, Julien Pommier.

  This file is a part of GetFEM++

  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
  under  the  terms  of the  GNU  Lesser General Public License as published
  by  the  Free Software Foundation;  either version 3 of the License,  or
  (at your option) any later version along with the GCC Runtime Library
  Exception either version 3.1 or (at your option) any later version.
  This program  is  distributed  in  the  hope  that it will be useful,  but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License and GCC Runtime Library Exception for more details.
  You  should  have received a copy of the GNU Lesser General Public License
  along  with  this program;  if not, write to the Free Software Foundation,
  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

  ===========================================================================*/

/**@file laplacian.cpp
  @brief Laplacian (Poisson) problem with generic assembly.
  */

#include "getfem/getfem_generic_assembly.h"
#include "getfem/getfem_export.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_superlu.h"
#include "getfem/getfem_models.h"
#include "getfem/getfem_model_solvers.h" // for preconditioners
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"

#include "getfem/getfem_mesh_slice.h"
#include "gmm/gmm.h"
#include "getfem/getfem_import.h"

#include "gmm/gmm_inoutput.h"

// Preconditioner
#include "biot_precond.hpp" 
#include "momentum_precond.hpp" 

/* some GetFEM++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node; /* geometrical nodes (derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::dim_type; 
using bgeot::base_tensor;

typedef gmm::rsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef gmm::col_matrix<sparse_vector_type> col_sparse_matrix_type;
typedef std::vector<scalar_type> plain_vector;

#define LS_TYPE 5 // 5 classic
// Right hand side. Allows an interpolation for the source term.
// scalar_type sol_f(const base_node &x) { return 10.; }

struct problem_descriptor_tri{    
	std::string MESH_TYPE =         "GT_PK(2,1)" ; // triangular elements
	std::string FEM_TYPE_U  =         "FEM_PK(2,1)";
	std::string FEM_TYPE_P  =         "FEM_PK(2,1)";
	std::string INTEGRATION =       "IM_TRIANGLE(6)";
	std::string SIMPLEX_INTEGRATION="IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),6)"; 
	std::string datafilename="resu/laplace2"; 
	int noised =0;  // noise on mesh
	int nsubdiv=32; // subdivision of the sqaured mesh
	double E=1.e+10;
	double poisson =0.3;
	double mu_s = E/( 2 * ( 1 + poisson) ) ;
	double lambda_l= E*poisson/ ( ( 1+poisson ) * (1 - 2 * poisson)) ;
	double biot_modulus=1.e+9;
	double k =1.e-15; //permeability
	double alpha=1; // Biot coefficient
	double rho_l=1000; // Biot coefficient
	double rho_r=2200; // Biot coefficient
	// non dimensional quantities
	double l_ref=1;
	double sigma_ref=1;
	double t_ref=1;
	double u_ref=1;
	double p_ref=1;

};

struct problem_descriptor_quad{    
	std::string MESH_TYPE =         "GT_LINEAR_QK(2)" ; // triangular elements
	std::string FEM_TYPE_U  =         "FEM_QK(2,2)";
	std::string FEM_TYPE_P  =         "FEM_QK(2,1)";
	std::string INTEGRATION =       "IM_GAUSS_PARALLELEPIPED(2,6)";
	std::string SIMPLEX_INTEGRATION="IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),6)"; 
	std::string datafilename="laplace"; 
	int nsubdiv=40; // subdivision of the sqaured mesh
	int noised =1;  // noise on mesh
	double E=1.e+4;
	double poisson =0.0;
	double mu_s = E/( 2 * ( 1 + poisson) ) ;
	double lambda_l= E*poisson/ ( ( 1+poisson ) * (1 - 2 * poisson)) ;
	double biot_modulus=1.e+18;
	double k =1.e-4; //permeability
	double alpha=1; // Biot coefficient
};

struct problem_descriptor_quad_3d{    
	std::string MESH_TYPE =         "GT_LINEAR_QK(3)" ; // triangular elements
	std::string FEM_TYPE_U  =         "FEM_QK(3,2)";
	std::string FEM_TYPE_P  =         "FEM_QK(3,1)";
	std::string INTEGRATION =       "IM_GAUSS_PARALLELEPIPED(3,6)";
	std::string SIMPLEX_INTEGRATION="IM_STRUCTURED_COMPOSITE(IM_TETRAHEDRON(6),3)"; 
	std::string datafilename="resu/laplace"; 
	int noised =0;  // noise on mesh
	int nsubdiv=6; // subdivision of the sqaured mesh
	double E=1.e+10;
	double poisson =0.3;
	double mu_s = E/( 2 * ( 1 + poisson) ) ;
	double lambda_l= E*poisson/ ( ( 1+poisson ) * (1 - 2 * poisson)) ;
	double biot_modulus=1.e+9;
	double k =1.e-15; //permeability
	double alpha=1; // Biot coefficient
	double rho_l=1000; // Biot coefficient
	double rho_r=2200; // Biot coefficient
	// non dimensional quantities
	double l_ref=1;
	double sigma_ref=1;
	double t_ref=1;
	double u_ref=1;
	double p_ref=1;
};
struct problem_descriptor_tetra_3d{    
	std::string MESH_TYPE =         "GT_PK(3,1)" ; // triangular elements
	std::string FEM_TYPE_U  =         "FEM_PK(3,1)";
	std::string FEM_TYPE_P  =         "FEM_PK(3,1)";
	std::string INTEGRATION =       "IM_TETRAHEDRON(6)";
	std::string SIMPLEX_INTEGRATION="IM_STRUCTURED_COMPOSITE(IM_TETRAHEDRON(6),3)"; 
	std::string datafilename="resu/ls_lk_p_nostnopr_p1"; 
	int noised =0;  // noise on mesh
	int nsubdiv=6; //9 // subdivision of the sqaured mesh
	double E=1.e+10;
	double poisson =0.3;
	double mu_s = E/( 2 * ( 1 + poisson) ) ;
	double lambda_l= E*poisson/ ( ( 1+poisson ) * (1 - 2 * poisson)) ;
	double biot_modulus=1.e+9;
	double k =1.e-15; //permeability
	double alpha=1; // Biot coefficient
	double rho_l=1000; // Biot coefficient
	double rho_r=2200; // Biot coefficient
	// non dimensional quantities
	double l_ref=1;
	double sigma_ref=1;
	double t_ref=1;
	double u_ref=1;
	double p_ref=1;
};
// Defining unit normal on a level set ------------------------------------

//   structure for the Laplacian problem
class biotls_problem {
	private:
		getfem::mesh mesh;
		getfem::mesh_im mim;      /// the integration methods
		getfem::mesh_fem mf_u;    /// the main mesh_fem, for the displacement solution
		getfem::mesh_fem mf_p;    /// the main mesh_fem, for the pressure solution
		getfem::mesh_fem mf_rhs;  /// the mesh_fem for the right hand side(f(x),..)
		getfem::mesh_fem mf_coef; /* the mesh_fem to represent pde coefficients    */
		getfem::mesh_fem mf_coef_v; /* the mesh_fem to represent pde vectorial coefficients    */
		// Levelset
		getfem::level_set ls; // create a levelset
		getfem::mesh_level_set mls;
		getfem::mesh mesh_ls;
		getfem::mesh_im_level_set mim_ls_all;
		getfem::mesh_im_level_set mim_ls_in;
		getfem::mesh_im_level_set mim_ls_out;
		getfem::mesh_im_level_set mim_ls_bd;
		getfem::mesh_fem_level_set mfls, mfls_u, mfls_p;
		// getfem::mesh_fem_level_set mfls_u_old, mfls_p_old;
		std::vector<size_type> eXt_dof, eXt_dof_u;  // The extended dofs
		std::vector<size_type> pin_index_, pout_index_;  // The extended dofs
		std::vector<size_type> uin_index_, uout_index_;  // The extended dofs
		size_type nb_x_dof_p, nb_x_dof_u;
// 		problem_descriptor_tri p_des;
		problem_descriptor_tetra_3d p_des;
		// problem_descriptor_tetra_3d p_des;
		enum { DIRICHLET_BOUNDARY_NUM = 10, NEUMANN_BOUNDARY_NUM = 11}; // descriptor for bcs flag
		enum { BOTTOM = 2, TOP = 1 , LEFT = 3, RIGHT =4, LEFTX = 5, RIGHTX =6}; // descriptor for zones
		enum { CUT_REGION = 100, UNCUT_REGION = 200, UNCUT_REGION_IN = 201, UNCUT_REGION_OUT = 202, CUT_EDGE=203};
		enum { MAT_1 = 50,MAT_2=60};
		size_type N_;             /// dimension of the problem

		///  workspace configuration parameters---------------------
		std::vector<scalar_type> tau_, vmu_, bm_ ,lambda_, beta_,
			alpha_, permeability_, force_,penalty_, c1_, c2_;
		// ---------------------------------------------------------
		sparse_matrix_type K;                                /// iteration matrix
		std::vector<scalar_type> U, U_old, P,  Px,           /// diplacement, disp old, pressure
			P_old, B, UP;               /// main unknown, and right hand side
		sparse_matrix_type Kp, Ku;                           /// iteration matrix for fixed steres of tpreconditioner
		std::vector<scalar_type> U_iter, P_iter, Bp, Bu;     /// main unknown, and right hand side
		biot_precond<sparse_matrix_type> *bPR_;               /// preconditioner based on fixed stress

		std::vector<scalar_type> Kr_; // permeability ratio
		std::vector<scalar_type> Er_; // young ratio
		std::vector<scalar_type> normal_ls_v;
                int step_=0;
		/// Methods
		void gen_bc(void);                                /// create zones for boundary conditions
		void gen_mat(void);                                /// create zones for internal conditions
		void compute_normal_2_ls(void);                   /// create normal to ls as a cell field
		void gen_coefficient();                         /// generate coefficient p0

		void configure_workspace                          /// configure the workspace add constants
			(getfem::ga_workspace & workspace,                /// the workspace
			 double dt);                                      /// timestep

		base_small_vector ls_function(const base_node P, 
				double time=0,
				int num = 0);

	public:
		void assembly(double dt,double time);                         /// assemble the monolithic iteration matrix for the problem
		void assembly_p(double dt,double time);                       /// assemble the iteration matrix for pressure, can be used as preconditioner
		void assembly_u(double dt,double time);                       /// assemble the iteration matrix for pressure, can be used as preconditioner
		void build_fix_stress_preconditioner(double dt, double time_ls);
		void solve(double time);                                 /// solves the monolithic system 
		void solve_fix_stress(double dt, int max_iter,double time);   /// solves the system with classic fixed stress approach
		void init(void);                                  /// initial configuration for the problem 
		void print(double time=0,int istep=0,double time_ls=0);

		void print_crop(double time=0,int istep=0,double time_ls=0);
                void print_pattern(int istep=0);
		void update_ls(double time=0, int iter=0);
		void update_p_index(double timels=0);
		void update_u_index(double timels=0);
                void set_step(int step){step_=step;}
		biotls_problem(void): mim(mesh), mf_u(mesh), mf_rhs(mesh), mf_p(mesh),mf_coef(mesh),mf_coef_v(mesh)
				      ,tau_(1), vmu_(1), bm_(1), lambda_(1),alpha_(1), permeability_(1), force_(1), beta_(1),penalty_(1),
				      c1_(1),c2_(1)
						   // level set 
						   ,ls(mesh,2),mls(mesh),
						   mim_ls_all(mls, getfem::mesh_im_level_set::INTEGRATE_ALL),
						   mim_ls_in(mls, getfem::mesh_im_level_set::INTEGRATE_INSIDE),
						   mim_ls_out(mls, getfem::mesh_im_level_set::INTEGRATE_OUTSIDE),
						   mim_ls_bd(mls, getfem::mesh_im_level_set::INTEGRATE_BOUNDARY),
						   mfls(mls, mf_u), mfls_u(mls, mf_u), mfls_p(mls, mf_p)
	{}
};




#endif // LAPLACE_GA_H
