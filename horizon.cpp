#include "horizon.hpp"

horizon::horizon(const char* file){
std::cout << "storing horizon "<<std::endl;
std::cout<<"reading points clouds"<<std::endl;
std::string namefile=file;
std::string namefile_cmd="gmsh:" + namefile;
std::cout<< "command "<< namefile_cmd<<std::endl;
getfem::import_mesh(namefile_cmd.c_str(),mesh_);
// A trasformation for the squarred mesh
bgeot::base_matrix M(3,3);
for (bgeot::size_type i=0; i < 3; ++i) {
  M(i,i) = 1.e+4;
}
//  if (N>1) { M(0,1) = 0; }
//
mesh_.transformation(M);
}
/////////////////////////////////////////////
int horizon::find_element(std::vector<double>& pt){
//   std::cout<<"Finding elemet"<<std::endl;
  std::vector<int> maybe_element;
  quick_search(pt, maybe_element);
  if(maybe_element.empty()){
    double eps=1.e-2;
    while (maybe_element.empty()){
    std::cout<< "Bounding box not found, using eps "<<eps<<std::endl;
      quick_search(pt, maybe_element,eps);
      eps*=2.;
    }
  }
  int found = fine_search(pt, maybe_element);
  std::cout<< "found return "<<found<<std::endl; 
  if (found == -1) {
    std::cout <<"Corresponding element not found, using "<<maybe_element[0]<<std::endl;
    found=maybe_element[0];

  }
 return found;
}
/////////////////////////////////////////////
void horizon::quick_search(std::vector<double>&pt,std::vector<int>& maybe_element, double eps){
  bgeot::size_type i_cv = 0;
  dal::bit_vector bv_cv = mesh_.convex_index();
  for (i_cv << bv_cv; i_cv !=bgeot::size_type(-1); i_cv << bv_cv) {
    bgeot::size_type i_pt = 0; 
    std::vector<double> b_box;b_box.resize(4);
    b_box[0]=1.e+10;b_box[2]=1.e+10;
    b_box[1]=-1.e+10;b_box[3]=-1.e+10;
    std::vector<bgeot::base_node> pts;
    for (int i = 0; i < mesh_.points_of_convex(i_cv).size(); i++)
      for (int idim=0; idim<2; idim++){
	// pts.push_back(mesh_.points_of_convex(i_cv)[i]);
	if((mesh_.points_of_convex(i_cv)[i])[idim] < b_box[2*idim])b_box[2*idim]= mesh_.points_of_convex(i_cv)[i][idim];     /// min 
	if(mesh_.points_of_convex(i_cv)[i][idim] > b_box[2*idim+1])b_box[2*idim+1]= mesh_.points_of_convex(i_cv)[i][idim]; /// max
      }
    if( pt[0] > b_box[0] - eps && pt[0] < b_box[1] + eps && pt[1] > b_box[0+2] - eps && pt[1] < b_box[1+2] + eps){
//       std::cout<<"Found possible element "<< i_cv<<std::endl;
      maybe_element.push_back(i_cv);
    }
     
  }
}
/////////////////////////////////////////////
int horizon::fine_search(std::vector<double>&pt,std::vector<int>& maybe_element){
  
  int shift_p0=0;int shift_p1=3;int shift_p2=6;
  int element=-1;
  for (int i_el=0; i_el<maybe_element.size();i_el++){
    std::vector<double> p; // p1x p1y p1z...p2x p2y p2z
//     std::cout<<"Testing element  "<<  maybe_element[i_el]<< std::endl;
    for (int ipt=0; ipt< mesh_.points_of_convex(maybe_element[i_el]).size(); ipt++)
      for (int idim=0; idim < mesh_.points_of_convex(maybe_element[i_el])[ipt].size(); idim++) 
	p.push_back(mesh_.points_of_convex(maybe_element[i_el])[ipt][idim]);
//       std::cout<< "Points dimesnion "<< p.size()<<std::endl;
      if (is_in(pt, {p[0 + shift_p0],p[1 + shift_p0],p[2 + shift_p0]} , 
		  {p[0 + shift_p1],p[1 + shift_p1],p[2 + shift_p1]} , 
		  {p[0 + shift_p2],p[1 + shift_p2],p[2 + shift_p2]})){
                  element=maybe_element[i_el];
// 		  std::cout<<"Found element "<<   maybe_element[i_el]<<std::endl;
		  }
  }
    return element;
}

//////////////////////////////////////////////////////////////

bool horizon::is_in (std::vector<double>  pt, std::vector<double>  v1, std::vector<double>  v2, std::vector<double>  v3)
{
    bool b1, b2, b3;
    
    b1 = sign(pt, v1, v2) < 0.0f;
    b2 = sign(pt, v2, v3) < 0.0f;
    b3 = sign(pt, v3, v1) < 0.0f;
    

    return ((b1 == b2) && (b2 == b3));
}



  int horizon::up_down(std::vector<double>& pt, int el){
    std::vector<double> p; // p1x p1y p1z...p2x p2y p2z
    for (int ipt=0; ipt< mesh_.points_of_convex(el).size(); ipt++)
      for (int idim=0; idim < mesh_.points_of_convex(el)[ipt].size(); idim++) 
	p.push_back(mesh_.points_of_convex(el)[ipt][idim]);
      //normal
    std::vector<double> n={0,0,0};
    double shiftp1=0;double shiftp2=3;double shiftp3=6;
    n[0] =     (p[1 + shiftp2] - p[1 + shiftp1])*(p[2 + shiftp3] - p[2 + shiftp1]) -
               (p[2 + shiftp2] - p[2 + shiftp1])*(p[1 + shiftp3] - p[1 + shiftp1]) ;
    n[1] =     (p[0 + shiftp2] - p[0 + shiftp1])*(p[2 + shiftp3] - p[2 + shiftp1]) -
               (p[2 + shiftp2] - p[2 + shiftp1])*(p[0 + shiftp2] - p[0 + shiftp1]) ;
    n[2] =     (p[0 + shiftp2] - p[0 + shiftp1])*(p[1 + shiftp3] - p[1 + shiftp1]) -
               (p[1 + shiftp2] - p[1 + shiftp1])*(p[0 + shiftp3] - p[0 + shiftp1]) ;
    double modn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);       
    n[0]/=modn;
    n[1]/=modn;
    n[2]/=modn;
    double d=-(n[0]*p[0] + n[1]*p[1]+n[2]*p[2]);
   std::cout<< "normale "<< n[0]<< " "<< n[1]<< " "<< n[2]<< std::endl;  
//    if(n[0]*pt[0] + n[1]*pt[1]+n[2]*pt[2] + d < 0)
//     return 1;
//    else
//     return 2;
   
   if(pt[2]< p[2])
    return 1;
   else
    return 2;
  }