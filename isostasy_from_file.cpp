#include "isostasy_from_file.hpp"
isostasy_from_files::isostasy_from_files(std::string filename, double dt):dt_(dt){
  
  // open a file in read mode.
  std::ifstream infile; 
  infile.open(filename.c_str()); 
  std::cout << "Reading isostasy from file "<<filename << std::endl;
  double buf;
  while(infile){
    infile>>buf;
    rotations_.push_back(buf);
  }
  infile.close();
  std::cout<< " size of rotation matrix "<< rotations_.size()<<std::endl;
//   std::cin.ignore();
}

void isostasy_from_files::get_rotation(double time,bgeot::base_matrix& M ){
  int time_idx=(int) (time/dt_);
  double weight, whole;
  weight=std::modf(time/dt_,&whole);
  int shift=9;
  for(int i =0; i<3;i++)
    for(int j =0; j<3; j++)
      M[i,j]=rotations_[(time_idx)*9 + i*3 + j]*(1.-weight) + rotations_[(time_idx+1)*9 + i*3 +j]*(weight);
  
}