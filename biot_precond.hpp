#ifndef BIOT_PRECOND_HS
#define BIOT_PRECOND_HS

#include "gmm_fix.hpp"
#include <vector>
#include <getfem/getfem_generic_assembly.h>
#include <gmm/gmm_precond_diagonal.h>
#include <gmm/gmm_superlu_interface.h>

//#define USE_SAMG 1
//#define USE_MP

// preconditioning strategies SOLVE invert by solving Ax=b
// LUMP mass lumping
  #define SOLVE_A 
  #define SOLVE_SHUR
//    #define LUMP_A
//   #define LUMP_SHUR


#ifdef USE_SAMG
#include "AMG_Interface.hpp"
#endif
typedef gmm::rsvector<bgeot::scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;

template <class MATRIX>
class biot_precond
{
public:
    // TODO boundary conditions (Dirichlet, Robin) and coefficient (kappa, beta)
    biot_precond(MATRIX &A,MATRIX &C, int nx1=0, int nx2=0, int ind=0);

    getfem::size_type nrows() const {
        return gmm::mat_nrows(A_) + gmm::mat_nrows(S_);
    }

    getfem::size_type ncols() const {
        return gmm::mat_ncols(A_) + gmm::mat_ncols(S_);
    }

    template <class L2, class L3>
    void mult(const L2 &src, L3 &dst) const
    {

        // std::cout<<"My mult"<<std::endl;
        const getfem::size_type n1 = gmm::mat_ncols(A_),
                                n2 = gmm::mat_ncols(S_);
       #ifdef SOLVE_A
        {
            bgeot::size_type restart = 50;
            std::vector<double> x(n1),b(n1);
            if (n1_x==0) gmm::copy(gmm::sub_vector(src, gmm::sub_interval(0, n1)),b);
            else {
                gmm::copy(gmm::sub_vector(src, gmm::sub_interval(0, n1-n1_x)),
                          gmm::sub_vector(b,   gmm::sub_interval(0, n1-n1_x)));
                gmm::copy(gmm::sub_vector(src, gmm::sub_interval(n1-n1_x+n2,n1_x)),
                          gmm::sub_vector(b,   gmm::sub_interval(n1-n1_x,n1_x)));
                }
            gmm::clear(x);                       
            gmm::identity_matrix PM; // no precond
            // gmm::diagonal_precond<sparse_matrix_type> PR(A_);
            gmm::iteration iter(1.e-12);  // iteration object with the max residu
            iter.set_noisy(0);               // output of iterations (2: sub-iteration)
            iter.set_maxiter(1000); // maximum number of iterations
            // gmm::gmres(A_, x, b, PM, restart, iter);
             slu_.solve(x,b);
            if (n1_x==0) gmm::copy(x,gmm::sub_vector(dst, gmm::sub_interval(0, n1)));
            else{
                gmm::copy(gmm::sub_vector(x, gmm::sub_interval(0, n1-n1_x)),
                         gmm::sub_vector(dst, gmm::sub_interval(0, n1-n1_x)));
                         
                gmm::copy(gmm::sub_vector(x, gmm::sub_interval(n1-n1_x,n1_x)),
                         gmm::sub_vector(dst, gmm::sub_interval(n1-n1_x+n2,n1_x)));
                }
        } 
        #endif
         #ifdef SOLVE_SHUR
          {
            bgeot::size_type restart = 50;
            std::vector<double> x(n2),b(n2);
            // gmm::diagonal_precond<sparse_matrix_type> PR(S_);
            gmm::copy(gmm::sub_vector(src, gmm::sub_interval(n1-n1_x, n2))
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
        
        if(n1_x==0) gmm::mult(pA_, gmm::sub_vector(src, gmm::sub_interval(0, n1)),
                                   gmm::sub_vector(dst, gmm::sub_interval(0, n1)));
        else{ 
                    std::vector<double> src1(n1);
                    std::vector<double> dst1(n1);
                    gmm::copy(
                              gmm::sub_vector(src, gmm::sub_interval(0, n1-n1_x)),
                              gmm::sub_vector(src1, gmm::sub_interval(0, n1-n1_x))
                              );
                    gmm::copy(
                              gmm::sub_vector(src, gmm::sub_interval(n1+n2-n1_x, n1_x)),
                              gmm::sub_vector(src1, gmm::sub_interval(n1-n1_x,n1_x))
                              );
                    gmm::mult(pA_, src1,dst1);
                    gmm::copy(
                              gmm::sub_vector(dst1, gmm::sub_interval(0, n1-n1_x)),
                              gmm::sub_vector(dst, gmm::sub_interval(0, n1-n1_x))
                              );
                    gmm::copy(
                              gmm::sub_vector(dst1, gmm::sub_interval(n1-n1_x,n1_x)),
                              gmm::sub_vector(dst, gmm::sub_interval(n1+n2-n1_x, n1_x))
                              );

            }
        #endif
        #ifdef LUMP_SHUR
        gmm::mult(pS_, gmm::sub_vector(src, gmm::sub_interval(n1-n1_x, n2)),
                  gmm::sub_vector(dst, gmm::sub_interval(n1-n1_x, n2)));
        #endif
//#ifdef USE_SAMG
        //std::vector<double> x(n2),b(n2);
        //gmm::copy(gmm::sub_vector(src, gmm::sub_interval(n1, n2)),b);
        //gmm::clear(x);
        //amg_.solve(S_, x , b , 1);
        //gmm::copy(amg_.getsol(),gmm::sub_vector(dst, gmm::sub_interval(n1, n2)));
//#else
        //slu_.solve(gmm::sub_vector(dst, gmm::sub_interval(n1, n2)),
                   //gmm::sub_vector(src, gmm::sub_interval(n1, n2)));
//#endif
        //// gmm::copy(gmm::sub_vector(src, gmm::sub_interval(n1, n2)),gmm::sub_vector(dst, gmm::sub_interval(n1, n2)));
        //gmm::mult(pAv_, gmm::sub_vector(src, gmm::sub_interval(n1+n2, n3)),
        //gmm::sub_vector(dst, gmm::sub_interval(n1+n2, n3)));

        //sluv_.solve(gmm::sub_vector(dst, gmm::sub_interval(n1+n2+n3, n4)),
                   //gmm::sub_vector(src, gmm::sub_interval(n1+n2+n3, n4)));

    // gmm::copy(src,dst);
    }

private:
    MATRIX &A_;
    gmm::diagonal_precond<MATRIX> pA_;
    gmm::diagonal_precond<MATRIX> pS_;
    MATRIX &S_;
    int n1_x, n2_x;
    int idx;
#ifdef USE_SAMG
    mutable AMG amg_;

#else
    gmm::SuperLU_factor<double> slu_;
#endif

gmm::SuperLU_factor<double> slup_;
};


namespace gmm {
    template <class MATRIX>
    struct linalg_traits<::biot_precond<MATRIX>> {
        using this_type = ::biot_precond<MATRIX>;
        using sub_orientation = owned_implementation;

        static size_type nrows(const this_type &m) { return m.nrows(); }
        static size_type ncols(const this_type &m) { return m.ncols(); }
    };
} // namespace gmm


template <class MATRIX>
biot_precond<MATRIX>::biot_precond(MATRIX &A, MATRIX &C, int nx1, int nx2, int ind)
: A_(A),S_(C), n1_x(nx1), n2_x(nx2), idx(ind)
#ifdef LUMP_A
,pA_(A)
#endif
#ifdef LUMP_SHUR
,pS_(C)
#endif
#ifdef USE_SAMG
, amg_("Schur")
#endif
{
    
    std::cout<<"Building biot preconditioner"<<std::endl;
    slu_.build_with(A_);    slup_.build_with(S_);

    //const getfem::size_type nb_dof_p = mf_p.nb_dof();
    //const getfem::mesh &mesh = mf_p.linked_mesh();
    //getfem::mesh_region inner_faces = getfem::inner_faces_of_mesh(mesh);
    //getfem::mesh_region outer_faces;
    //getfem::outer_faces_of_mesh(mesh, outer_faces);
    //getfem::ga_workspace wp;

    //std::vector<double> p(nb_dof_p);
    //wp.add_fem_variable("p", mf_p, gmm::sub_interval(0, nb_dof_p), p);  

    //#ifdef USE_MP
    //wp.add_expression("p*Test_p",mim);
    //#else
    //wp.add_expression("Grad_p.Grad_Test_p", mim);
    //wp.add_expression("-0.5 * (Grad_p + Interpolate(Grad_p, neighbour_elt)).Normal"
                        //" * (Test_p - Interpolate(Test_p, neighbour_elt))"
                      //"-0.5 * (Grad_Test_p + Interpolate(Grad_Test_p, neighbour_elt)).Normal"
                        //" * (p - Interpolate(p, neighbour_elt))"
                      //"+2 / element_size * (p - Interpolate(p, neighbour_elt))"
                        //" * (Test_p - Interpolate(Test_p, neighbour_elt))",
                      //mim, inner_faces);
    //// to remove in case of mix/neumann condition
    //wp.add_expression("1/element_size*p*Test_p",
                      //mim, outer_faces);
     //#endif
    //wp.assembly(2);

    //gmm::copy(wp.assembled_matrix(), S_);


//#ifdef USE_SAMG
    //amg_.convert_matrix(S_);
//#else
    //slu_.build_with(S_);
//#endif
}








#endif // ifndef darcyprecond

