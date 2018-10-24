/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

// #define BIOT

#define BIOT_LS
 #define TEMP_LS

#ifdef BIOT_LS
#include "biot_ls.hpp"
#ifdef TEMP_LS
#include "temp_ls.hpp"
#endif
#endif
#ifdef BIOT
#include "biot.hpp"
#endif
//  #include "temp.hpp"


int main(int argc, char *argv[]) {

	GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
	FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

	try {    
#ifdef BIOT_LS
		biotls_problem p;
#ifdef TEMP_LS
		templs_problem t;
#endif
#endif
#ifdef BIOT
		 biot_problem p;
#endif

//                 temperature_problem t;
		p.init();
#ifdef TEMP_LS
 		t.init();
#endif
		// double dt=1e-3;
		double dt=1e+12;
		// p.build_fix_stress_preconditioner(dt,0);
		//    p.assembly_p(dt,0); 
		//    p.assembly_u(dt);
		int n_step=80;int erosion_limit=0; // usually 50 steps
		double time=0*dt;
		double time_ls  =0;
		p.update_ls(time_ls, 0);
#ifdef TEMP_LS
		t.update_ls(time_ls, 0);
#endif
		for(int istep=0; istep<n_step; istep++)
		{
#ifdef BIOT_LS
// // 		  //level set
			time=istep*dt;
			std::cout<< "*** Time step "<< istep << "/"<<n_step<<std::endl;
			{ // level set
				if (istep<erosion_limit && istep > 30){ 
					time_ls  =istep*dt;
					p.update_ls(time_ls, istep);
#ifdef TEMP_LS
					t.update_ls(time_ls, istep);
#endif

					// p.build_fix_stress_preconditioner(dt,time_ls);
				}
                                 p.set_step(istep);
#ifdef TEMP_LS
                                 t.set_step(istep);
#endif
				 p.solve_fix_stress(dt, 2000,time_ls);
#if defined TEMP_LS && defined BIOT_LS
                                 t.set_pressure(p.get_pressure());
#endif
#ifdef TEMP_LS
                                 t.assembly(dt,time_ls);
		  	   	 t.solve(time_ls);
#endif
                         
			 	//  p.build_fix_stress_preconditioner(dt,time_ls);
				// ===============================
// 			 	  p.assembly(dt,time_ls);
				//   p.solve(time_ls);
				// ===========================
				 p.print(istep*dt,istep,time_ls);      
				 p.print_crop(istep*dt,istep,time_ls);
#ifdef TEMP_LS
			 	 t.print_crop(istep*dt,istep,time_ls);
				 t.print(istep*dt,istep,time_ls);   
#endif   
				// p.print_ls(istep*dt,istep,time_ls); //pint the ls slice
// 				p.print_pattern(istep);
				// p.update_time_iter(istep);
			} // endl of lev_set biot
#endif
// //==================================================//
#ifdef BIOT
				{ // classic biot
// 				          t.assembly(dt, p.get_pressure_fem(), p.get_pressure());
// 				          t.solve();
// 					  t.print(istep);
					  if(istep!=0) p.solve_fix_stress(dt, 100);
					  else  p.solve_fix_stress(dt, 5);
					 //  p.assembly(dt);
					 //  p.solve();
					 p.print(istep);
				}
#endif
// //==================================================//
		}
	}
	GMM_STANDARD_CATCH_ERROR;

	return 0; 
}
