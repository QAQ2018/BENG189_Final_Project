function a=alphan1(v)
%filename: alphan.m
kv=150;
theta=(kv*v+14.273)/10;
if(theta==0)   %check for case that gives 0/0
  a=0.01265;  %in that case use L'Hospital's rule
else
  a=0.01265*theta/(1-exp(-theta));
end
