clear all
matrix_base_name = "lk_ls_ovp.kp";
iter =  int2str(0);
matrix_name_b=strcat(matrix_base_name,iter)
% load /u/archive/agip/cerroni/software/mygetfem/biot/resu/lk_ls_ovp.kp.0.mm
% Kf=full(spconvert(K));
% load Ku_s
% load Ku_x
% ks=full(spconvert(Ku_s));
% kx=full(spconvert(Ku_x));
% P=[ks, zeros(length(ks),length(kx)); zeros(length(kx),length(ks)),diag(diag(kx))];
% cond(Kf)
% cond(inv(P)*Kf)
