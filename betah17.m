function b=betah17(v)
%filename: alphah.m
k1=-0.00283; k2=2.00283; k3=5.5266; k4=-12.70195;
kv=150;
b=k1+k2/(1+exp((kv*v+k3)/k4));