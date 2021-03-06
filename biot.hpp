#ifndef BIOT_H
#define BIOT_H
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

#include "gmm/gmm.h"

#include "gmm/gmm_inoutput.h"

// Preconditioner
#include "biot_precond.hpp" 

/* some GetFEM++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node; /* geometrical nodes (derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::dim_type; 

typedef gmm::rsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef gmm::col_matrix<sparse_vector_type> col_sparse_matrix_type;
typedef std::vector<scalar_type> plain_vector;

// Right hand side. Allows an interpolation for the source term.
// scalar_type sol_f(const base_node &x) { return 10.; }

struct problem_descriptor_tri{    
    std::string MESH_TYPE =         "GT_PK(2,1)" ; // triangular elements
    std::string FEM_TYPE_U  =         "FEM_PK(2,2)";
    std::string FEM_TYPE_P  =         "FEM_PK(2,1)";
	std::string INTEGRATION =       "IM_TRIANGLE(6)";
    std::string SIMPLEX_INTEGRATION="IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),6)"; 
    std::string datafilename="laplace"; 
    int nsubdiv=10; // subdivision of the sqaured mesh
    double E=1.e+4;
	double poisson =0.3;
	double mu_s = E/( 2 * ( 1 + poisson) ) ;
	double lambda_l= E*poisson/ ( ( 1+poisson ) * (1 - 2 * poisson)) ;
	double biot_modulus=1.e+18;
	double k =1.e-4; //permeability
	double alpha=1; // Biot coefficient
    };
    
struct problem_descriptor_quad{    
    std::string MESH_TYPE =         "GT_LINEAR_QK(2)" ; // triangular elements
    std::string FEM_TYPE_U  =         "FEM_QK(2,2)";
    std::string FEM_TYPE_P  =         "FEM_QK(2,1)";
	std::string INTEGRATION =       "IM_GAUSS_PARALLELEPIPED(2,6)";
    std::string SIMPLEX_INTEGRATION="IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),6)"; 
    std::string datafilename="laplace"; 
    int nsubdiv=10; // subdivision of the sqaured mesh
    double E=1.e+4;
	double poisson =0.0;
	double mu_s = E/( 2 * ( 1 + poisson) ) ;
	double lambda_l= E*poisson/ ( ( 1+poisson ) * (1 - 2 * poisson)) ;
	double biot_modulus=1.e+18;
	double k =1.e-4; //permeability
	double alpha=1; // Biot coefficient
    };

//   structure for the Laplacian problem
class biot_problem {
    private:
      getfem::mesh mesh;
      getfem::mesh_im mim;      /// the integration methods
      getfem::mesh_fem mf_u;    /// the main mesh_fem, for the displacement solution
      getfem::mesh_fem mf_p;    /// the main mesh_fem, for the pressure solution
      getfem::mesh_fem mf_rhs;  /// the mesh_fem for the right hand side(f(x),..)
      problem_descriptor_tri p_des;
      enum { DIRICHLET_BOUNDARY_NUM = 10, NEUMANN_BOUNDARY_NUM = 11}; // descriptor for bcs flag
      enum { BOTTOM = 2, TOP = 1 , LEFT = 3, RIGHT =4}; // descriptor for zones
      size_type N_;             /// dimension of the problem
      
      ///  workspace configuration parameters---------------------
      std::vector<scalar_type> tau_, vmu_, bm_ ,lambda_, beta_,
                               alpha_, permeability_, force_;
     // ---------------------------------------------------------
     sparse_matrix_type K;                                /// iteration matrix
     std::vector<scalar_type> U, U_old, P,                /// diplacement, disp old, pressure
                              P_old, B, UP;               /// main unknown, and right hand side
     sparse_matrix_type Kp, Ku;                           /// iteration matrix for fixed steres of tpreconditioner
     std::vector<scalar_type> U_iter, P_iter, Bp, Bu;     /// main unknown, and right hand side
     biot_precond<sparse_matrix_type> *bPR_;               /// preconditioner based on fixed stress
     
     /// Methods
      void gen_bc(void);                                /// create zones for boundary conditions
      void configure_workspace                          /// configure the workspace add constants
          (getfem::ga_workspace & workspace,                /// the workspace
           double dt);                                      /// timestep
    public:
      void assembly(double dt);                         /// assemble the monolithic iteration matrix for the problem
      void assembly_p(double dt);                       /// assemble the iteration matrix for pressure, can be used as preconditioner
      void assembly_u(double dt);                       /// assemble the iteration matrix for pressure, can be used as preconditioner
      void build_fix_stress_preconditioner();
      void solve(void);                                 /// solves the monolithic system 
      void solve_fix_stress(double dt, int max_iter);   /// solves the system with classic fixed stress approach
      void init(void);                                  /// initial configuration for the problem 
      void print(int time=0);
      biot_problem(void): mim(mesh), mf_u(mesh), mf_rhs(mesh), mf_p(mesh)
      ,tau_(1), vmu_(1), bm_(1), lambda_(1),alpha_(1), permeability_(1), force_(1), beta_(1)
      {}
};



 
#endif // LAPLACE_GA_H
