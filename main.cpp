/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/
#include "biot_ls.hpp"


int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {    
    biotls_problem p;
    p.init();


     double dt=0.5e+8;
     p.assembly_p(dt); 
     p.assembly_u(dt);
     p.build_fix_stress_preconditioner();
     int n_step=80;
     for(int istep=1; istep<n_step; istep++)
     {
         
       if (istep<40) p.update_ls(istep*dt);
       p.solve_fix_stress(dt, 2000);
       
      //    p.assembly(dt);
      //   p.solve();
       
       
       p.print(istep);
   }
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
