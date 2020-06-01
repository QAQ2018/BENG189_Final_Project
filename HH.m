%filename HH.m
%numerical solution of the space-clamped Hodgkin-Huxley equations
clear all; clf; close all;
global check;
global t1p t2p ip; %parameters for the function izero(t)
in_HH
in_mhnv
receptor_activate_time
%%
neuron1 = {};neuron2 = {};neuron3 = {};
v1=vhold;
neuron1.v(1)=v1; %The sensory neuron
neuron1.m17(1)=alpham17(v1)/(alpham17(v1)+betam17(v1));
neuron1.h17(1)=alphah17(v1)/(alphah17(v1)+betah17(v1));
neuron1.s17(1)=alphas17(v1)/(alphas17(v1)+betas17(v1));
neuron1.m18(1)=alpham18(v1)/(alpham18(v1)+betam18(v1));
neuron1.h18(1)=1/(1+exp(v1+32.2/4));
neuron1.n1(1)=1/(1+exp((-v+14.62)/18.38));
neuron1.nKA(1)=1/(1+exp(-(v+5.4)/16.4))^4;
neuron1.hKA(1)=1/(1+exp((v1+49.9)/4.6));
v2=vhold;
neuron2.v(1)=v2; %The interneuron in the spine
neuron2.m(1)=alpham(v2)/(alpham(v2)+betam(v2));
neuron2.h(1)=alphah(v2)/(alphah(v2)+betah(v2));
neuron2.n(1)=alphan(v2)/(alphan(v2)+betan(v2));
v3=vhold;
neuron3.v(1)=v3; %The efferent neuron
neuron3.m(1)=alpham(v3)/(alpham(v3)+betam(v3));
neuron3.h(1)=alphah(v3)/(alphah(v3)+betah(v3));
neuron3.n(1)=alphan(v3)/(alphan(v3)+betan(v3));
for klok=1:klokmax
  t=klok*dt;                      %note time
  %Neuron 1: The sensory neuron
  neuron1.m17(klok+1)=snew(neuron1.m17(klok),alpham17(v1),betam17(v1),dt); %update m
  neuron1.h17(klok+1)=snew(neuron1.h17(klok),alphah17(v1),betah17(v1),dt); %update h
  neuron1.s17(klok+1)=snew(neuron1.s17(klok),alphas17(v1),betas17(v1),dt); %update s
  neuron1.m18(klok+1)=snew(neuron1.m18(klok),alpham18(v1),betam18(v1),dt); %update m
  neuron1.h18(klok+1)=h18_update(neuron1.h18(klok), v1, dt); %update m
  neuron1.n1(klok+1)=n1_update(neuron1.n1(klok), v1, dt); %update n
  neuron1.nKA(klok+1)=nKA_update(neuron1.nKA(klok), v1, dt);
  neuron1.hKA(klok+1)=hKA_update(neuron1.hKA(klok), v1, dt);
  gNa17=gNabar17*(neuron1.m17(klok+1)^3)*neuron1.h17(klok+1)*neuron1.s17(klok+1);    %sodium conductance
  gNa18=gNabar18*neuron1.m18(klok+1)*neuron1.h18(klok+1);
  gNa=gNa17+gNa18;
  gK1=gKbar1*neuron1.n1(klok+1);    %potassium conductance
  gKA=gKbar1A*neuron1.nKA(klok+1)*neuron1.hKA(klok+1);
  gK=gK1+gKA;
  g=gNa+gK+gLbar1;         %total conductance
  gE=gNa*(v1-ENa1)+gK*(v2-EK1)+gLbar1*(v3-EL1);         %gE=g*E
  %gE=gNa*ENa2+gK*EK2+gLbar2*EL2;         %gE=g*E

  %save old value of v for checking purposes:
  %update v:
  v1=(v1+(dt/C1)*(gE+izero(t)))/(1+(dt/C1)*g);
  neuron1.v(klok+1)=v1;
  
  %Neuron 2: The interneuron
  neuron2.m(klok+1)=snew(neuron2.m(klok),alpham(v2),betam(v2),dt); %update m
  neuron2.h(klok+1)=snew(neuron2.h(klok),alphah(v2),betah(v2),dt); %update h
  neuron2.n(klok+1)=snew(neuron2.n(klok),alphan(v2),betan(v2),dt); %update n
  if t>= AMPA_activate_start && t <= AMPA_activate_end
      gNa_addition = gAMPA*normpdf(t,AMPA_activate_start+0.5,0.25)^10;
  else
      gNa_addition = 0;
  end
  gNa=gNabar2*(neuron2.m(klok+1)^3)*neuron2.h(klok+1)+gNa_addition;%sodium conductance
  gK =gKbar2*(neuron2.n(klok+1)^4);%potassium conductance
  g=gNa+gK+gLbar2;         %total conductance
  gE=gNa*ENa2+gK*EK2+gLbar2*EL2;         %gE=g*E
  %update v:
  v2=(v2+(dt/C2)*gE)/(1+(dt/C2)*g);
  neuron2.v(klok+1)=v2;
  
  
  %Neuron 3: The efferent neuron
  neuron3.m(klok+1)=snew(neuron3.m(klok),alpham(v3),betam(v3),dt); %update m
  neuron3.h(klok+1)=snew(neuron3.h(klok),alphah(v3),betah(v3),dt); %update h
  neuron3.n(klok+1)=snew(neuron3.n(klok),alphan(v3),betan(v3),dt); %update n
  %Check whether nAChR opens.
  if t>= nAChR_activate_start && t <= nAChR_activate_end
      gNa_addition = gnAChR*normpdf(t,nAChR_activate_start+0.5,0.25)^10;
  else
      gNa_addition = 0;
  end
  gNa=gNabar3*(neuron3.m(klok+1)^3)*neuron3.h(klok+1)+gNa_addition;    %sodium conductance
  gK =gKbar3*(neuron3.n(klok+1)^4);    %potassium conductance
  g=gNa+gK+gLbar3;         %total conductance
  gE=gNa*ENa3+gK*EK3+gLbar3*EL3;         %gE=g*E
  %update v:
  v3=(v3+(dt/C3)*gE)/(1+(dt/C3)*g);
  neuron3.v(klok+1)=v3; 
  
  %When the action potential arrives at the terminal and AMPA is not
  %activated in the last 5 ms, there will be a release of Glutamate in the 
  %next 1ms, lasting for 1ms to activate AMPA.
  if v1 >= +7 && t > AMPA_activate_start + 5
      AMPA_activate_start=t+1;
      AMPA_activate_end=AMPA_activate_start+1;
  end
  %When the action potential arrives at the terminal and ACh is not
  %activated in the last 5 ms, there will be a activate of ACh in the 
  %next 1ms, lasting for 1ms, to activate nAChR.
  if v2 >= +7 && t > nAChR_activate_start + 5
      nAChR_activate_start=t+1;
      nAChR_activate_end=nAChR_activate_start+1;
  end
end

%%
%subplot(2,1,1),plot(t_plot,v_plot);hold on
%subplot(2,1,2),plot(t_plot,mhn_plot);legend('m', 'h', 'n')
t_plot = 0:dt:tmax;
plot(t_plot, neuron1.v); hold on;
plot(t_plot, neuron2.v); plot(t_plot, neuron3.v);
legend('sensory', 'interneuron', 'efferent');

figure
plot(t_plot,[neuron1.m17; neuron1.h17; neuron1.s17])