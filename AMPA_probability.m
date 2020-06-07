function P = AMPA_probability(t)
tau1=0.2;
tau2=1.3;
Pmax = 1;
T=1; % -log(1-Pmax)*tau1
if t<=T
    P=Pmax*(1-exp(-t/tau1));
else
    P = Pmax*exp(-t/tau2);
end