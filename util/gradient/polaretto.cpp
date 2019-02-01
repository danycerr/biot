#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry> 

#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <string>      // std::ifstream
#include <sstream>      // std::stringstream

using Eigen::MatrixXd;
void extractRotation(const Eigen::Matrix3d &A, Eigen::Quaterniond &q,
		     const unsigned int maxIter);

int main()
{
  std::ofstream outfile;
  outfile.open ("rotation.txt");
  for (int ifile =0; ifile<190; ifile++ )
  {
    Eigen::Matrix3d m; 
    std::stringstream namefile;
    namefile << "cell_grad/cel_gradi."<<ifile<<".csv";
    std::ifstream file ( namefile.str() );//  declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
    std::string line;
    //     {
    //       getline ( file, value, ',' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
    //       std::cout << std::string( value, 1, value.length()-2 ); // display value removing the first and the last character from it
    //     }
    getline(file,line);
    getline(file,line,',');
    for (int i =0; i< 3 ; i++)
      for (int j =0; j< 3 ; j++){
	getline(file,line,',');
	m(i,j)=std::stof(line)+ ((i==j)? 1:0);
      }//
      getline(file,line,',');
    
    std::cout<<m<<std::endl;
    Eigen::Quaterniond q(m);
    extractRotation(m,q,500);
    for (int i =0; i< 3 ; i++)
      for (int j =0; j< 3 ; j++) outfile<<q.matrix()(i,j)<<" ";
      outfile<<std::endl;
//     std::cout << "=================" << std::endl;
//     std::cout << m << std::endl;
//     std::cout << "========***=========" << std::endl;
//     std::cout << q.matrix() << std::endl;
  }//end for ifile
  
  
}

void extractRotation(const Eigen::Matrix3d &A, Eigen::Quaterniond &q,
		     const unsigned int maxIter)
{
  for (unsigned int iter = 0; iter < maxIter; iter++)
  {
    std::cout << "num iteration "<< iter << std::endl;
    Eigen::Matrix3d R = q.matrix();
    Eigen::Vector3d omega = (R.col(0).cross(A.col(0)) + R.col
    (1).cross(A.col(1)) + R.col(2).cross(A.col(2))
    )
    *
    (1.0 / fabs(R.col(0).dot(A.col(0)) + R.col
    (1).dot(A.col(1)) + R.col(2).dot(A.col(2))) +
    1.0e-9);
    double w = omega.norm();
    if (w < 1.0e-9)
      break;
    q = Eigen::Quaterniond(Eigen::AngleAxisd(w, (1.0 / w)
    *
    omega))
    *
    q;
    q.normalize();
    
  }
}