function b=betan1(v)
%filename: betan.m
theta = (v+55)/2.5;
kv=150;
b=k1+k2/(1+exp((kv*v+k3)/k4));