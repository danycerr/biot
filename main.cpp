/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/
// #include "biot_ls.hpp"
 #include "biot.hpp"
 #include "temp.hpp"


int main(int argc, char *argv[]) {

	GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
	FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

	try {    
	// 	biotls_problem p;
		 biot_problem p;
		p.init();

                temperature_problem t;
		t.init();
		// double dt=1e-3;
		double dt=1e+8;
		// p.build_fix_stress_preconditioner(dt,0);
		//    p.assembly_p(dt,0); 
		//    p.assembly_u(dt);
		int n_step=80;int erosion_limit=20; // usually 80 steps
		double time=0*dt;
		double time_ls  =0;
		for(int istep=0; istep<n_step; istep++)
		{//level set
	// 		time=istep*dt;
	// 		std::cout<< "*** Time step "<< istep << "/"<<n_step<<std::endl;
	// 		{ // level set
	// 			if (istep<erosion_limit){ 
	// 				time_ls  =istep*dt;
	// 				p.update_ls(time_ls, istep);

	// 				// p.build_fix_stress_preconditioner(dt,time_ls);
	// 			}
	// 	  	   	 p.solve_fix_stress(dt, 2000,time_ls);

	// 		 	//  p.build_fix_stress_preconditioner(dt,time_ls);
	// 			// ===============================
	// 		 	//   p.assembly(dt,time_ls);
	// 			//   p.solve(time_ls);
	// 			// ===========================
	// 			p.print(istep*dt,istep,time_ls);      
	// 			p.print_crop(istep*dt,istep,time_ls);
	// 		} // endl of lev_set biot
				{ // classic biot
				          t.assembly(dt, p.get_pressure_fem(), p.get_pressure());
				          t.solve();
					  t.print(istep);
					  if(istep!=0) p.solve_fix_stress(dt, 100);
					  else  p.solve_fix_stress(dt, 5);
					 //  p.assembly(dt);
					 //  p.solve();
					 p.print(istep);
				}
		}
	}
	GMM_STANDARD_CATCH_ERROR;

	return 0; 
}
