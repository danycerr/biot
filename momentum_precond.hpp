#ifndef MOMENTUM_PRECOND_HS
#define MOMENUTM_PRECOND_HS

#include "gmm_fix.hpp"
#include <vector>
#include <getfem/getfem_generic_assembly.h>
#include "gmm/gmm.h"
#include <gmm/gmm_precond_diagonal.h>
#include <gmm/gmm_superlu_interface.h>
// #include "getfem/getfem_mesher.h"
// #define USE_SAMG 1
//#define USE_MP

// preconditioning strategies SOLVE invert by solving Ax=b
// LUMP mass lumping
#define SOLVE_A_MOMENTUM
//   #define SOLVE_SHUR
//     #define LUMP_A
#define LUMP_SHUR


#ifdef USE_SAMG
#include "AMG_Interface.hpp"
#endif
typedef gmm::rsvector<bgeot::scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;

using bgeot::scalar_type; 

template <class MATRIX>
class momentum_precond
{
public:
  // TODO boundary conditions (Dirichlet, Robin) and coefficient (kappa, beta)
  momentum_precond(MATRIX &A, MATRIX &B, int ind=0);
  
  getfem::size_type nrows() const {
    return gmm::mat_nrows(A_) + gmm::mat_nrows(S_);
  }
  
  getfem::size_type ncols() const {
    return gmm::mat_ncols(A_) + gmm::mat_ncols(S_);
  }
  
  template <class L2, class L3>
  void mult(const L2 &src, L3 &dst) const
  {
    
    //         std::cout<<"My mult"<<std::endl;
    const getfem::size_type n1 =nb_dof_, // gmm::mat_ncols(A_),
    n2 =nb_dof_x_; // gmm::mat_ncols(S_);
    #ifdef SOLVE_A_MOMENTUM
    {
      bgeot::size_type restart = 50;
      std::vector<double> x(n1),b(n1);
      gmm::copy(gmm::sub_vector(src, gmm::sub_interval(0, n1)),b);
      gmm::clear(x);                       
      gmm::identity_matrix PM; // no precond
      // gmm::diagonal_precond<sparse_matrix_type> PR(A_);
      gmm::iteration iter(1.e-12);  // iteration object with the max residu
      iter.set_noisy(0);               // output of iterations (2: sub-iteration)
      iter.set_maxiter(1000); // maximum number of iterations
      // gmm::gmres(A_, x, b, PM, restart, iter);
      #ifdef USE_SAMG
      gmm::clear(x);
      amg_.solve(Ku_csr_, x , b , 1);
      // 	    std::cout<<"my mult"<<std::endl;
      gmm::copy(amg_.getsol(),gmm::sub_vector(dst, gmm::sub_interval(0, n1)));
      #else
      slu_.solve(x,b);
      gmm::copy(x,gmm::sub_vector(dst, gmm::sub_interval(0, n1)));
      #endif
      
      
    }
    
    #endif
    #ifdef SOLVE_SHUR_MOMENTUM
    {
      bgeot::size_type restart = 50;
      std::vector<double> x(n2),b(n2);
      // gmm::diagonal_precond<sparse_matrix_type> PR(S_);
      gmm::copy(gmm::sub_vector(src, gmm::sub_interval(n1, n2))
      ,gmm::sub_vector(b, gmm::sub_interval(0, n2))
      );
      gmm::clear(x);                       
      gmm::identity_matrix PM; // no precond
      gmm::iteration iter(1.e-12);  // iteration object with the max residu
      iter.set_noisy(0);               // output of iterations (2: sub-iteration)
      iter.set_maxiter(1000); // maximum number of iterations
      // gmm::gmres(S_, x, b, PM, restart, iter);
      slup_.solve(x,b);
      gmm::copy(gmm::sub_vector(x, gmm::sub_interval(0, n2)),
		gmm::sub_vector(dst, gmm::sub_interval(n1-n1_x, n2))
      );
    }
    #endif
    #ifdef LUMP_A
    
    gmm::mult(pA_, gmm::sub_vector(src, gmm::sub_interval(0, n1)),
	      gmm::sub_vector(dst, gmm::sub_interval(0, n1)));
    
    #endif
    #ifdef LUMP_SHUR
    gmm::mult(pS_, gmm::sub_vector(src, gmm::sub_interval(n1, n2)),
	      gmm::sub_vector(dst, gmm::sub_interval(n1, n2)));
    #endif
    
  }
  
private:
  MATRIX &A_;
  gmm::diagonal_precond<MATRIX> pA_;
  gmm::diagonal_precond<MATRIX> pS_;
  MATRIX &S_;
  bgeot::size_type nb_dof_, nb_dof_x_;
  int idx_;
  gmm::csr_matrix<scalar_type> Ku_csr_;
  #ifdef USE_SAMG
  mutable AMG amg_;
  
  #else
  gmm::SuperLU_factor<double> slu_;
  #endif
  
  gmm::SuperLU_factor<double> slup_;
};


namespace gmm {
  template <class MATRIX>
  struct linalg_traits<::momentum_precond<MATRIX>> {
    using this_type = ::momentum_precond<MATRIX>;
    using sub_orientation = owned_implementation;
    
    static size_type nrows(const this_type &m) { return m.nrows(); }
    static size_type ncols(const this_type &m) { return m.ncols(); }
  };
} // namespace gmm


template <class MATRIX>
momentum_precond<MATRIX>::momentum_precond(MATRIX &A, MATRIX &B, int ind)
// :nb_dof_(nb_dof), nb_dof_x_(nb_dof_x)
: A_(A),S_(B), idx_(ind)
// :A_(sub_matrix( A ,
//         gmm::sub_interval(0, nb_dof),
//         gmm::sub_interval(0, nb_dof) )),
// S_(sub_matrix( A ,
//         gmm::sub_interval(0, nb_dof_x),
//         gmm::sub_interval(0, nb_dof_x) ))
#ifdef USE_SAMG
, amg_("Schur")
#endif
{
  //       gmm::sub_matrix( Ku ,
  //       gmm::sub_interval(nb_dof_u,nb_x_dof_u  ),
  //       gmm::sub_interval(nb_dof_u,nb_x_dof_u ) ) 
  std::cout<<"Building biot preconditioner"<<std::endl;
  ////  std::cout<<"nb_dof "<<nb_dof<<" nb_dof_x "<<nb_dof_x<<std::endl;
  std::cout<<"nb_dof "<<A.ncols()<<std::endl;
  std::cout<<"nb_dof "<<A.nrows()<<std::endl;
  // std::cout<<"**********"<<A<<std::endl;
  nb_dof_=A.ncols();
  nb_dof_x_=B.ncols();
  std::cout<<"nb_dof A_ "<<A_.ncols()<<std::endl;
  std::cout<<"nb_dof S_ "<<S_.ncols()<<std::endl;
  gmm::MatrixMarket_IO::write("Ku_s",A_);
  gmm::MatrixMarket_IO::write("Ku_x",S_);
  #ifdef LUMP_A
  pA_.build_with(A_);
  #endif
  #ifdef LUMP_SHUR
  pS_.build_with(S_);
  #endif
  
  #ifdef USE_SAMG 
  gmm::copy(A_, Ku_csr_);
  amg_.convert_matrix(Ku_csr_);
  amg_.setid(idx_);
  #else
  slu_.build_with(A_);    slup_.build_with(S_);
  #endif
}








#endif // ifndef darcyprecond

