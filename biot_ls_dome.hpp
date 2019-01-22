#ifndef BIOTLSDOME_H
#define BIOTLSDOME_H
#include "biot_ls.hpp"

class biot_ls_dome : public biotls_problem{
private:
protected:
  void import_mesh(void); //method for picking the mesh 
  base_small_vector ls_function(const base_node P, 
				double time=0,
				int num = 0);
  void gen_bc(void);
public:
  void assembly_p(double dt,double time);                       /// assemble the iteration matrix for pressure, can be used as preconditioner
  void assembly_u(double dt,double time);                       /// assemble the iteration matrix for pressure, can be used as preconditioner
		
};
#endif