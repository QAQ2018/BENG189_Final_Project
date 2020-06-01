function a=alphas17(v)
%filename: alphah.m
k1=3e-5; k2=9.2e-4; k3=93.9; k4=16.6;
kv=150;
a=k1+k2/(1+exp((kv*v+k3)/k4));