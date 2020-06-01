function hKA = hKA_update(hKA_old, v, dt)
kv=150;
hKA_inf = 1/(1+exp((kv*v+49.9)/4.6));
tau = 20+50*exp(-(kv*v+40)^2/(2*40^2));
if tau < 5
    tau=5
end
hKA = hKA_old+(hKA_inf-hKA_old)*dt/tau;