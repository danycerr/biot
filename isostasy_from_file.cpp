#include "isostasy_from_file.hpp"
//===========================================================
isostasy_from_files::isostasy_from_files(std::string filename,std::string filename_tilt,std::string filename_vdisp, double dt):dt_(dt){
  
  // open a file in read mode.
  {
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
  }
  {
    std::ifstream infile; 
    infile.open(filename_tilt.c_str()); 
    std::cout << "Reading tilt from file "<<filename_tilt << std::endl;
    double buf;
    while(infile){
      infile>>buf;
      tilts_.push_back(buf);
    }
    infile.close();
    std::cout<< " size of tilt matrix "<< tilts_.size()<<std::endl;
  }
    {
    std::ifstream infile; 
    infile.open(filename_vdisp.c_str()); 
    std::cout << "Reading tilt from file "<<filename_vdisp << std::endl;
    double buf;
    std::string line;
    std::getline(infile,line);
//     infile.getline(filename_vdisp, 100);
    char comma;
    while(std::getline(infile,line)) // vertical dispacement is the third value
    {
        std::istringstream iss(line);
        iss >> buf >> comma >> buf >> comma >> buf ;
        disp_z_.push_back(buf*100);
    }
    infile.close();
    std::cout<< " size of disp vector "<< disp_z_.size()<<std::endl;
  }
  //   std::cin.ignore();
}
//===========================================================
void isostasy_from_files::get_rotation(double time,bgeot::base_matrix& M ){
  int time_idx=(int) (time/dt_);
  double weight, whole;
  weight=std::modf(time/dt_,&whole);
  int shift=9;
  for(int i =0; i<3;i++)
    for(int j =0; j<3; j++){
      M(i,j)=(double)(rotations_[(time_idx)*9 + i*3 + j]*(1.-weight)) + rotations_[(time_idx+1)*9 + i*3 +j]*(weight);
      
    }
    std::cout<<"sent matrix from file is   "<<M<<std::endl;
  std::cout<<"and weight    "<< weight <<"  index  "<<time_idx<<std::endl;
}
//===========================================================
void isostasy_from_files::get_tilt(double time,std::vector<double> &  tilt ){
  int time_idx=(int) (time/dt_);
  double weight, whole;
  weight=std::modf(time/dt_,&whole);
  int shift=3;
  for(int i =0; i<3;i++)
    tilt[i]=(tilts_[(time_idx)*shift + i]*(1.-weight)) + tilts_[(time_idx+1)*shift + i]*(weight);
}//===========================================================
void isostasy_from_files::get_disp(double time,double &  dz ){
  int time_idx=(int) (time/dt_);
  double weight, whole;
  weight=std::modf(time/dt_,&whole);
  dz=(disp_z_[time_idx]*(1.-weight)) + disp_z_[time_idx+1]*(weight);
}