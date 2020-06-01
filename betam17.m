function b=betam17(v)
%filename: alphah.m
k1=0; k2=35.2; k3=72.7; k4=16.7;
kv=150;
b=k1+k2/(1+exp((kv*v+k3)/k4));