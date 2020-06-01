function h18 = h18_update(h18_old, v, dt)
kv=150;
h18_inf = 1/(1+exp(kv*v+32.2/4));
tau = 1.218+42.043*exp(-(kv*v+38.1)^2/(2*15.19^2));
h18 = h18_old+(h18_inf-h18_old)*dt/tau;