clear all
load K.mm
Kf=full(spconvert(K));
load Ku_s
load Ku_x
ks=full(spconvert(Ku_s));
kx=full(spconvert(Ku_x));
P=[ks, zeros(length(ks),length(kx)); zeros(length(kx),length(ks)),diag(diag(kx))];
cond(Kf)
cond(inv(P)*Kf)
