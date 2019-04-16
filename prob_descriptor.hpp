
#ifndef PDES_H
#define PDES_H

struct generic_problem_descriptor_tri{    
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
        double k_temp=1;

};

struct generic_problem_descriptor_quad{    
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

struct generic_problem_descriptor_quad_3d{    
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
struct generic_problem_descriptor_tetra_3d{    
	std::string MESH_TYPE =         "GT_PK(3,1)" ; // triangular elements
	std::string FEM_TYPE_U  =         "FEM_PK(3,1)";
	std::string FEM_TYPE_P  =         "FEM_PK(3,1)";
	std::string INTEGRATION =       "IM_TETRAHEDRON(6)";
	std::string SIMPLEX_INTEGRATION="IM_STRUCTURED_COMPOSITE(IM_TETRAHEDRON(6),3)"; 

// 	std::string datafilename="resu/lk_ls8_2mat_temp_coup"; 
	std::string datafilename="resu/ls_ring_pinch_dtemp"; 
	int noised =0;  // noise on mesh
	int nsubdiv=6; // subdivision of the sqaured mesh

	double E=1.e+10;
	double poisson =0.3;
	double mu_s = E/( 2 * ( 1 + poisson) ) ;
	double lambda_l= E*poisson/ ( ( 1+poisson ) * (1 - 2 * poisson)) ;
	double biot_modulus=1.e+9;
	double k =1.e-19; //permeability
	double alpha=1; // Biot coefficient
	double rho_l=1000; // Biot coefficient
	double rho_r=2200; // Biot coefficient
	// non dimensional quantities
	double l_ref=1;
	double sigma_ref=1;
	double t_ref=1;
	double u_ref=1;
	double p_ref=1;
        double alpha_temp=1.2e-6; //thermal conductivity 0.8e-5
        double q_rad=0.1e-6/(2200); // 1.e-6
        // double q_rad=1.e-10;
};

#endif
