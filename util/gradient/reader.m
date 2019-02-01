clear all
%  for (istep=0:189)
for (istep=0:10)
filename = strcat("cell_grad/cel_gradi.", num2str(istep),".csv")
%  filename = "cell_grad/cel_gradi.2.csv";
M = csvread (filename);

for(j=1:3)
for(i=1:3)
a(j,i)= M(2,1+i+(j-1)*3) +eye(3)(j,i);
endfor
endfor
C=a'*a;
[V, LAMBDA] = eig (C);

n1=V(:,1);
n2=V(:,2);
n3=V(:,3);

D=LAMBDA(1,1) * n1*n1' +LAMBDA(2,2) * n2*n2' + LAMBDA(3,3) * n3*n3';

U=sqrt(LAMBDA(1,1)) * n1*n1' +sqrt(LAMBDA(2,2)) * n2*n2' + sqrt(LAMBDA(3,3)) * n3*n3';

R=a*inv(U);

pt1=[0.157279,-0.0895207, 0.983459];
pt2=[0.156757,-0.0892232, -0.3475];

pt1d=a*pt1';
pt2d=a*pt2';

vec1=pt2-pt1;

vec2=pt2d-pt1d;

x = atan2(R(3,2), R(3,3));
y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
z = atan2(R(2,1), R(1,1));

rtot(istep+1,:)=[R(1,:),R(2,:),R(3,:)];
angle_tot(istep+1,:)=[x,y,z];
endfor

filename_rot = "rotation_matrix.txt";
fid = fopen (filename_rot, "w");
fputs (fid, "# du1dx1 du1dx2 du1dx3 du2dx1 du2dx2 du2dx3 du3dx1 du3dx2 du3dx3\n");
fclose (fid);
tmp=size(rtot)(1);
save("-ascii","-append","rotation_matrix.txt","tmp")
save("-ascii","-append","rotation_matrix.txt","rtot")
save("-ascii","angle_xyz.txt","angle_tot")


