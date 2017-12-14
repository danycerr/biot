/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/
#include "biot.hpp"


int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {    
    biot_problem p;
     p.init();


     double dt=1.e-3;
     p.assembly_p(dt); 
     p.assembly_u(dt);
     p.build_fix_stress_preconditioner();
     int n_step=200;
     for(int istep=0; istep<n_step; istep++)
     {
       // p.solve_fix_stress(dt, 2000);
       
        p.assembly(dt);
        p.solve();
       
       
       p.print(istep);
   }
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
