function b=betam18(v)
%filename: alphah.m
k1=0; k2=7.6205; k3=46.463; k4=8.8289;
kv=150;
b=k1+k2/(1+exp((kv*v+k3)/k4));