function P = NMDA_probability(V, t)
tau1=1.5;
tau2=152;
Pmax=0.99;
T=5;
if t<T
    P = 0;
else
    G = 1/(1+10/3.57*exp(-V/16.13));
    P = Pmax*exp(-t/tau2)*G;
end
