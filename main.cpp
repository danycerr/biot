/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

#define BIOT

#define TEMPERATURE
#define ISOSTASY
#include "gmm/gmm_except.h"

// #define BIOT_LS
// #define TEMP_LS

#ifdef BIOT_LS
#include "biot_ls.hpp"
#endif
#ifdef TEMP_LS
#include "temp_ls.hpp"
#endif
#ifdef BIOT
#include "biot.hpp"
#endif
#ifdef TEMPERATURE
#include "temp.hpp"
#endif
#ifdef ISOSTASY
#include "isostasy.hpp"
#endif

int main(int argc, char *argv[]) {

	GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
	FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

	try {    
#ifdef BIOT_LS
		biotls_problem p;
#endif
#ifdef TEMP_LS
		templs_problem t;
#endif
#ifdef BIOT
		 biot_problem p;
#endif
#ifdef TEMPERATURE
                 temperature_problem t;
#endif
#ifdef ISOSTASY
                 isostasy isos;
#endif

//                 temperature_problem t;
// 		t.init();
// 		// double dt=1e-3;
// 		double dt=0.5e+10;
// 		// p.build_fix_stress_preconditioner(dt,0);
// 		//    p.assembly_p(dt,0); 
// 		//    p.assembly_u(dt);
// 		int n_step=80;int erosion_limit=50; // usually 50 steps
// 
// 		double time=0*dt;
// #ifdef BIOT_LS
// 		double time_ls  =0;p.update_ls(time_ls, 0);
// =======
//                 temperature_problem t;
#if defined BIOT || defined BIOT_LS
		p.init();
#endif
#if defined TEMPERATURE ||defined TEMP_LS 
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
#ifdef BIOT_LS
		p.update_ls(time_ls, 0);
#endif
#ifdef TEMP_LS
		t.update_ls(time_ls, 0);
#endif
#ifdef ISOSTASY
// 		set center of rotation of isostasy
		isos.set_center_of_rotation( t.get_mesh()); 
#endif
		for(int istep=0; istep<n_step; istep++)
		{
#if defined TEMP_LS || defined BIOT_LS
// // 		  //level set
			time=istep*dt;
			std::cout<< "*** Time step "<< istep << "/"<<n_step<<std::endl;
			{ // level set
				if (istep<erosion_limit && istep > 30){ 
					time_ls  =istep*dt;
#ifdef BIOT_LS
					p.update_ls(time_ls, istep);
#endif
#ifdef TEMP_LS
					t.update_ls(time_ls, istep);
#endif

					// p.build_fix_stress_preconditioner(dt,time_ls);
				}
#ifdef BIOT_LS
                                 p.set_step(istep);
#endif
#ifdef TEMP_LS
                                 t.set_step(istep);
#endif
#ifdef BIOT_LS
				 p.solve_fix_stress(dt, 2000,time_ls);
#endif
#if defined TEMP_LS && defined BIOT_LS
                                 t.set_pressure(p.get_pressure());
#endif
#ifdef TEMP_LS
                                 if(istep>0) t.set_clt(true); // imposing constant temperature lateral
				 t.assembly(dt,time_ls);
		  	   	 t.solve(time_ls);
#endif
                         
			 	//  p.build_fix_stress_preconditioner(dt,time_ls);
				// ===============================
// 			 	  p.assembly(dt,time_ls);
				//   p.solve(time_ls);
				// ===========================
#ifdef BIOT_LS
				 p.print(istep*dt,istep,time_ls);      
				 p.print_crop(istep*dt,istep,time_ls);
#endif   
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
#if defined BIOT || defined TEMPERATURE
				{ // classic biot
				         
#ifdef TEMPERATURE
                                          t.set_iter(istep);
//                                           t.move_mesh();
#ifdef ISOSTASY
					  isos.set_time(istep*dt);
					  isos.set_transformation();
#ifdef TEMPERATURE
					  isos.move_mesh( t.get_mesh()); 
#endif
					  std::cout<< "Gravity is "<< *(isos.get_gravity())<<std::endl;
#ifdef BIOT
					  p.set_gravity(isos.get_gravity());
#endif
#endif
#if defined BIOT && defined TEMPERATURE
				          t.assembly(dt, p.get_pressure_fem(), p.get_pressure());
#elif defined TEMPERATURE
					  t.assembly(dt);			  
#endif
					  t.solve();
					  t.print(istep);
					  
// 					  isos.undo_move_mesh( t.get_mesh()); 
#endif
#ifdef BIOT
				          p.set_iter(istep);
					  p.assembly_p(dt);//dummy assembling just for ice
					  if(istep!=0) p.solve_fix_stress(dt, 100); // 100
					  else  p.solve_fix_stress(dt, 5); //5
					  // p.assembly(dt);
					  // p.solve();
#ifdef BIOT
					  isos.move_mesh( p.get_mesh()); 
#endif
					  p.print(istep);
//                                           p.print_aux_data(istep);
#endif
#ifdef ISOSTASY
#ifdef TEMPERATURE
					 isos.undo_move_mesh( t.get_mesh()); 
#endif
#ifdef BIOT
					 isos.undo_move_mesh( p.get_mesh()); 
#endif
#endif
				}
#endif
// //==================================================//
		}
	}
	GMM_STANDARD_CATCH_ERROR;

	return 0; 
}
