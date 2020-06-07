%filename HH.m
%numerical solution of the space-clamped Hodgkin-Huxley equations
clear all; clf; close all;
global check;
in_HH_TTX
%% Defining the neurons
nociceptors = {};interneuron = {};motorNeuron = {};
noci_ids = 1:20; % The indexes of all nociceptors
v1=vhold;
nociceptors.v(1,noci_ids)=v1; %The sensory neuron
nociceptors.m17(1,noci_ids)=alpham17(v1)/(alpham17(v1)+betam17(v1));
nociceptors.h17(1,noci_ids)=alphah17(v1)/(alphah17(v1)+betah17(v1));
nociceptors.s17(1,noci_ids)=alphas17(v1)/(alphas17(v1)+betas17(v1));
nociceptors.m18(1,noci_ids)=alpham18(v1)/(alpham18(v1)+betam18(v1));
nociceptors.h18(1,noci_ids)=1/(1+exp(v1+32.2/4));
nociceptors.n1(1,noci_ids)=1/(1+exp((-v1+14.62)/18.38));
nociceptors.nKA(1,noci_ids)=1/(1+exp(-(v1+5.4)/16.4))^4;
nociceptors.hKA(1,noci_ids)=1/(1+exp((v1+49.9)/4.6));

% initialize the neurotransmitter release times
neurotransmitter_release 

v2=vhold;
interneuron.v(1)=v2; %The interneuron in the spine
interneuron.m(1)=alpham(v2)/(alpham(v2)+betam(v2));
interneuron.h(1)=alphah(v2)/(alphah(v2)+betah(v2));
interneuron.n(1)=alphan(v2)/(alphan(v2)+betan(v2));
v3=vhold;
motorNeuron.v(1)=v3; %The efferent neuron
motorNeuron.m(1)=alpham(v3)/(alpham(v3)+betam(v3));
motorNeuron.h(1)=alphah(v3)/(alphah(v3)+betah(v3));
motorNeuron.n(1)=alphan(v3)/(alphan(v3)+betan(v3));

%% Dunning the simulation
%initialize parameters of subsequent current pulse:
t1p(noci_ids)=20;       %starting time (ms)
t2p(noci_ids)=180;       %stopping time (ms)
ip(noci_ids)=(0:1:19)+8;        %current applied (muA)
for klok=1:klokmax
  t=klok*dt;                      %note time
  %Neuron 1: The sensory neuron
  for id = noci_ids
  v1=nociceptors.v(klok,id);
  nociceptors.m17(klok+1,id)=snew(nociceptors.m17(klok,id),alpham17(v1),betam17(v1),dt); %update m1.7
  nociceptors.h17(klok+1,id)=snew(nociceptors.h17(klok,id),alphah17(v1),betah17(v1),dt); %update h1.7
  nociceptors.s17(klok+1,id)=snew(nociceptors.s17(klok,id),alphas17(v1),betas17(v1),dt); %update s1.7
  nociceptors.m18(klok+1,id)=snew(nociceptors.m18(klok,id),alpham18(v1),betam18(v1),dt); %update m1.8
  nociceptors.h18(klok+1,id)=h18_update(nociceptors.h18(klok,id), v1, dt); %update m1.8
  nociceptors.n1(klok+1,id)=n1_update(nociceptors.n1(klok,id), v1, dt);    %update n
  nociceptors.nKA(klok+1,id)=nKA_update(nociceptors.nKA(klok,id), v1, dt); %update nKA
  nociceptors.hKA(klok+1,id)=hKA_update(nociceptors.hKA(klok,id), v1, dt); %
  gNabar17 = 0.01;
  gNabar18= 1;
  gNa17=gNabar17*(nociceptors.m17(klok+1,id)^3)*nociceptors.h17(klok+1,id)*nociceptors.s17(klok+1,id);
  gNa18=gNabar18*nociceptors.m18(klok+1,id)*nociceptors.h18(klok+1,id);
  gNa=gNa17+gNa18;
  gK1=gKbar1*nociceptors.n1(klok+1,id);    %potassium conductance
  gKA=gKbar1A*nociceptors.nKA(klok+1,id)*nociceptors.hKA(klok+1,id);
  gK=gK1+gKA;
  g=gNa+gK+gLbar1;         %total conductance
  gE=gNa*ENa1+gK*EK1+gLbar1*EL1;         %gE=g*E
  %update v:
  v1=(v1+(dt/C1)*(gE+izero(t1p(id), t2p(id), ip(id), t)))/(1+(dt/C1)*g);
  nociceptors.v(klok+1,id)=v1;
  
  %When the action potential arrives at the terminal and Glutamate is not
  %released in the last 8 ms, there will be a release of Glutamate in the 
  %next 5ms, lasting for 1ms to activate AMPA.
    if nociceptors.v(klok,id) >= +7 && t > Glu_release_start(id) + 8
      Glu_release_start(id)=t+5;
      Glu_release_end(id)=Glu_release_start(id)+1;
    end
  end
  

  %Neuron 2: The interneuron
  gAMPA=0;  gNMDA=0;
  for id = noci_ids
    %if t>= Glu_release_start(id)+dt && t <= Glu_release_end(id)+20
    if t>= Glu_release_start(id)+dt && Glu_release_end(id) > 0
      gAMPA = gAMPA+gAMPA_max*AMPA_probability(t-Glu_release_start(id));
      gNMDA = gNMDA+gNMDA_max*NMDA_probability(v2, t-Glu_release_start(id));
      gAMPA = gAMPA;
      gNMDA = gNMDA;
    else
      gAMPA = gAMPA+0;
      gNMDA = gNMDA+0;
    end
  end
  
  % Weaken the synnaptic strength
  gAMPA = gAMPA/20;
  gNMDA = gNMDA/20;
  interneuron.m(klok+1)=snew(interneuron.m(klok),alpham(v2),betam(v2),dt); %update m
  interneuron.h(klok+1)=snew(interneuron.h(klok),alphah(v2),betah(v2),dt); %update h
  interneuron.n(klok+1)=snew(interneuron.n(klok),alphan(v2),betan(v2),dt); %update n
  Isyn(klok+1) = gNMDA*ECa2+gAMPA*ENa2;
  %Isyn(klok+1) = gAMPA*ENa2;
  gNa=gNabar2*(interneuron.m(klok+1)^3)*interneuron.h(klok+1);%sodium conductance
  gNa=gNa;
  gK =gKbar2*(interneuron.n(klok+1)^4);%potassium conductance
  g=gNa+gK+gLbar2+gNMDA+gAMPA;         %total conductance
  gE=gNa*ENa2+gK*EK2+gLbar2*EL2;         %gE=g*E
  %update v:
  v2=(v2+(dt/C2)*(gE+Isyn(klok+1)))/(1+(dt/C2)*g);
  interneuron.v(klok+1)=v2;
  
  %When the action potential arrives at the terminal and ACh is not
  %released in the last 5 ms, there will be a release of ACh in the 
  %next 5ms, lasting for 1ms, to activate nAChR.
  if v2 >= +0 && t > ACh_release_start + 5
      ACh_release_start=t+5;
      ACh_release_end=ACh_release_start+1;
  end
  
  %Neuron 3: The efferent neuron
  motorNeuron.m(klok+1)=snew(motorNeuron.m(klok),alpham(v3),betam(v3),dt); %update m
  motorNeuron.h(klok+1)=snew(motorNeuron.h(klok),alphah(v3),betah(v3),dt); %update h
  motorNeuron.n(klok+1)=snew(motorNeuron.n(klok),alphan(v3),betan(v3),dt); %update n
  %Check whether nAChR opens.
  if t>= ACh_release_start && t <= ACh_release_end
      gSyn = gnAChR_max*normpdf(t,ACh_release_start+0.5,0.25)^10;
  else
      gSyn = 0;
  end
  gNa=gNabar3*(motorNeuron.m(klok+1)^3)*motorNeuron.h(klok+1)+gSyn;
  gNa=gNa;
  gK =gKbar3*(motorNeuron.n(klok+1)^4);    %potassium conductance
  g=gNa+gK+gLbar3+gSyn;         %total conductance
  gE=gNa*ENa3+gK*EK3+gLbar3*EL3;         %gE=g*E
  %update v:
  v3=(v3+(dt/C3)*gE)/(1+(dt/C3)*g);
  motorNeuron.v(klok+1)=v3; 
  

end

%%
%subplot(2,1,1),plot(t_plot,v_plot);hold on
%subplot(2,1,2),plot(t_plot,mhn_plot);legend('m', 'h', 'n')
t_plot = 0:dt:tmax;
subplot(3,1,1)
plot(t_plot, nociceptors.v);
legend('nociceptor');ylabel('voltage (mV)');
subplot(3,1,2)
plot(t_plot, interneuron.v); 
legend('interneuron');ylabel('voltage (mV)');
subplot(3,1,3)
plot(t_plot, motorNeuron.v);
legend('motor neuron');ylabel('voltage (mV)');xlabel('time (ms)');
%legend('sensory', 'interneuron', 'efferent');


figure()
subplot(3,1,1)
plot(t_plot, sum(nociceptors.v>0, 2));
legend('nociceptor');ylabel('No. firing');
subplot(3,1,2)
plot(t_plot, interneuron.v>0); 
legend('interneuron');ylabel('No. firing');
subplot(3,1,3)
plot(t_plot, motorNeuron.v>0);
legend('motor neuron');ylabel('No. firing');xlabel('time (ms)');
%legend('sensory', 'interneuron', 'efferent');
%figure
%plot(t_plot,[nociceptors.m17; nociceptors.h17; nociceptors.s17])