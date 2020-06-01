function nKA = nKA_update(nKA_old, v, dt)
nKA_inf = (1+exp(-(v+5.4)/16.4))^-4;
tau = 0.25+10.04*exp(-(v+24.67)^2/(2*34.8^2));
nKA = (nKA_old*tau+nKA_inf*dt)/(tau+dt);