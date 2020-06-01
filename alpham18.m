function a=alpham18(v)
%filename: alphah.m
k1=2.85; k2=-2.839; k3=-1.159; k4=13.95;
kv=150;
a=k1+k2/(1+exp((kv*v+k3)/k4));