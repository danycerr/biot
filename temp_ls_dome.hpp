#ifndef TEMPLSDOME_H
#define TEMPLSDOME_H

#include "temp_ls.hpp"

class temp_ls_dome: public templs_problem {
private:
  base_small_vector ls_function(const base_node P, 
				double time=0,
				int num = 0);
  void gen_bc(void);
protected:
  void import_mesh(void);
public:
  void assembly(double dt,double time);
};
#endif