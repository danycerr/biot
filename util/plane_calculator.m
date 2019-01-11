% module load gcc-glibc/7
% module load octave/4.2.1

clear all 
p1=[15499 , 3000 , 4858 ]
p2=[15499 , 7000 , 4858 ]
p3=[19501 , 3000 , 4144 ]

p1p2 = p2 - p1;
p1p3 = p3 - p1;

coeff = cross(p1p2,p1p3);

a = coeff(1)
b = coeff(2)
c = coeff(3)
d = a*p1(1) + b*p1(2) + c*p1(3)
