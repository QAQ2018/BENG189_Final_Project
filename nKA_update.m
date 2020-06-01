function nKA = nKA_update(nKA_old, v, dt)
kv=150;
nKA_inf = 1/(1+exp(-(kv*v+5.4)/16.4))^4;
tau = 0.25+10.04*exp(-(kv*v+24.69)^2/(2*34.8^2));
nKA = nKA_old+(nKA_inf-nKA_old)*dt/tau;