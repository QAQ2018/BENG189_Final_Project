function a=alpham17(v)
%filename: alphah.m
k1=0; k2=15.5; k3=-5; k4=-12.08;
a=k1+k2/(1+exp((v+k3)/k4));