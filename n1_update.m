function n1 = n1_update(n1_old, v, dt)
theta=(v+14.273)/10;
if(theta==0)   %check for case that gives 0/0
  a=0.01265;  %in that case use L'Hospital's rule
else
  a=0.01265*theta/(1-exp(-theta));
end
theta = (v+55)/2.5;
b=0.125*exp(-theta);

n1_inf = 1/(1+exp(-(v+14.62)/18.38));
tau = 1/(a+b)+1;
n1 = (n1_old*tau+n1_inf*dt)/(tau+dt);