#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>


class time_step{
private:
   int ts_;
   int num_fs_it_;
   double avg_gmres_u_;
   double avg_gmres_p_;
   std::vector<double> norm_p_,norm_u_;
   std::vector<int> iter_p_,iter_u_;
public:
   time_step(int ts)
               :ts_(ts),num_fs_it_(0),
               avg_gmres_u_(0),avg_gmres_p_(0){};
   void add_subiter(double norm_p, double norm_u, int iter_p, int iter_u){
              norm_p_.push_back(norm_p);norm_u_.push_back(norm_u);
              iter_p_.push_back(iter_p);iter_u_.push_back(iter_u);
   }
   void finalize(){
   num_fs_it_=norm_p_.size();
   for(int i=0; i<num_fs_it_; i++){
     avg_gmres_p_+= (double) iter_p_[i] / num_fs_it_;
     avg_gmres_u_+= (double) iter_u_ [i]/ num_fs_it_;
     }
   }; // end of finalize
  
   void summary(){
   std::cout<< "ts " <<ts_<<" num_fs_it " << num_fs_it_ 
            << " avg_gmres_u " << avg_gmres_u_ 
            << " avg_gmres_p " << avg_gmres_p_ <<std::endl;
   };
   

   int get_ts() {return ts_;};
   int get_fs() {return num_fs_it_;};
   double get_avg_u() {return  avg_gmres_u_;};
   double get_avg_p() {return avg_gmres_p_;};

};

int main(int argc, char* argv[]){
   std::ifstream infile;
   std::string line;
   std::vector<time_step> ts;
   double buf;
   std::cout<<"Opening file"<<std::endl;
   infile.open(argv[1]);
   if(!infile) 
      {
      std::cout<<"Error reading"<<std::endl;
      return 13;
      }
   while (std::getline(infile, line))
   {
     std::istringstream iss(line);
     int a, b;
     double buf;
     std::string buf_str;
     iss >> buf_str >> buf_str;
     if (!buf_str.compare("last")){
         int indx = ts.size()-1;
         ts[indx].finalize();
     }
     else{
         int l_ts, f_s_i;
         double norm_p, norm_u;
         int iter_p,iter_u;
         iss >> buf_str >> l_ts >>buf_str >> buf_str >> buf_str >> f_s_i >> buf_str >> buf_str >> 
         norm_p >> buf_str >>buf_str >> norm_u >> buf_str >> buf_str >> iter_p >> buf_str >>
         buf_str >> iter_u; 
         
         std::cout <<"l_ts " << l_ts
                   <<" f_s_i " << f_s_i
                   <<" norm_p " << norm_p
                   <<" norm_u " << norm_u
                   <<" iter_p " << iter_p
                   <<" iter_u " << iter_u
         <<std::endl;
         if (f_s_i == 1) {
                          time_step buf_ts(l_ts);
                          ts.push_back(buf_ts);
                          }
         int indx = ts.size()-1;
         ts[indx].add_subiter(norm_p,norm_u,iter_p,iter_u);
     }
   }
   infile.close();
   std::cout<<"Number of time steps "<<ts.size()<<std::endl;
   std::ofstream outfile;
   outfile.open (argv[2]);
   outfile << "# ts fsi iter_u iter_p"<<std::endl;
    for (int i=0; i < ts.size(); i++) 
      {  
          // ts[i].summary();
       outfile<< ts[i].get_ts()<<" "<< ts[i].get_fs() / 100.<<" " << ts[i].get_avg_u() / 1000.<<" "
       << ts[i].get_avg_p() / 100. <<std::endl;
      }

  outfile.close();
return 1;
}
