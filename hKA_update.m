function hKA = hKA_update(hKA_old, v, dt)
hKA_inf = 1/(1+exp((v+49.9)/4.6));
tau = 20+50*exp(-(v+40)^2/(2*40^2));
if tau < 5
    tau=5;
end
hKA = (hKA_old*tau+hKA_inf*dt)/(tau+dt);