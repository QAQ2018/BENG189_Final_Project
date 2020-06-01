function b=betas17(v)
%filename: alphah.m
k1=132.05; k2=-132.05; k3=-384.9; k4=28.5;
b=k1+k2/(1+exp((v+k3)/k4));