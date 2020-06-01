function a=alphah17(v)
%filename: alphah.m
k1=0; k2=0.38685; k3=122.35; k4=15.29;
kv=150;
a=k1+k2/(1+exp((kv*v+k3)/k4));