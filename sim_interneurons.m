%filename HH.m
%numerical solution of the space-clamped Hodgkin-Huxley equations
clear all; clf; close all;
global check;
in_HH
in_mhnv
%% Defining the neurons

v2=vhold;
interneuron.v(1,1:2)=v2; %The interneuron in the spine
interneuron.m(1,1:2)=alpham(v2)/(alpham(v2)+betam(v2));
interneuron.h(1,1:2)=alphah(v2)/(alphah(v2)+betah(v2));
interneuron.n(1,1:2)=alphan(v2)/(alphan(v2)+betan(v2));

% initialize the neurotransmitter release times
Glu_release_start(1:2) = -100;
Glu_release_end(1:2) = -100;
Glu_release(1:2) = false;
v2(1:2) = vhold;
gAMPA(1:1001,1:2)=0;
%% Running the simulation
%initialize parameters of subsequent current pulse:
t1p(1:2)=10;       %starting time (ms)
t2p(1:2)=20;       %stopping time (ms)
ip(1)=25;        %current applied (muA)
ip(2)=0;
for klok=1:klokmax
  t=klok*dt;                      %note time
  
  %When the action potential arrives at the terminal and Glutamate is not
  %released in the last 8 ms, there will be a release of Glutamate in the 
  %next 10ms, lasting for 1ms to activate AMPA.
  for id=1:2
    if interneuron.v(klok,id) >= 0 && t > Glu_release_start(id) + 5
      Glu_release_start(id)=t+5;
      Glu_release_end(id)=Glu_release_start(id)+1;
      Glu_release(id)=true;
    else
      Glu_release(id)=false;
    end
  end
  for id = 1:2
    %Neuron 2: The interneuron
    gNMDA=0;
    if t>= Glu_release_start(3-id)+dt && Glu_release_start(3-id)>=0
      gAMPA(klok+1,id) = gAMPA_max*AMPA_probability(t-Glu_release_start(3-id));
      gNMDA = gNMDA+gNMDA_max*NMDA_probability(v2(id), t-Glu_release_start(3-id));
      %gNMDA=0;
    else
      gAMPA(klok+1,id) = 0;
      gNMDA = 0;
    end
    interneuron.m(klok+1,id)=snew(interneuron.m(klok,id),alpham(v2(id)),betam(v2(id)),dt); %update m
    interneuron.h(klok+1,id)=snew(interneuron.h(klok,id),alphah(v2(id)),betah(v2(id)),dt); %update h
    interneuron.n(klok+1,id)=snew(interneuron.n(klok,id),alphan(v2(id)),betan(v2(id)),dt); %update n
    Isyn(klok+1,id) = gNMDA*ECa2+gAMPA(klok+1,id)*ENa2;
    %Isyn(klok+1) = gAMPA*ENa2;
    gNa=gNabar2*(interneuron.m(klok+1,id)^3)*interneuron.h(klok+1,id);%sodium conductance
    gK =gKbar2*(interneuron.n(klok+1,id)^4);%potassium conductance
    g=gNa+gK+gLbar2+gNMDA+gAMPA(klok+1,id);         %total conductance
    gE=gNa*ENa2+gK*EK2+gLbar2*EL2;         %gE=g*E
    %update v:
    v2(id)=(v2(id)+(dt/C2)*(gE+Isyn(klok+1,id)+izero(t1p(id), t2p(id), ip(id), t)))/(1+(dt/C2)*g);
    disp(izero(t1p(id), t2p(id), ip(id), t))
    interneuron.v(klok+1,id)=v2(id);
  end

end

%%
%subplot(2,1,1),plot(t_plot,v_plot);hold on
%subplot(2,1,2),plot(t_plot,mhn_plot);legend('m', 'h', 'n')
t_plot = 0:dt:tmax;
%plot(t_plot, nociceptors.v); hold on;
subplot(2,1,1);
plot(t_plot, interneuron.v); %plot(t_plot, motorNeuron.v);
legend('interneuron1', 'interneuron2');
ylabel('voltage (mV)')
%figure
%plot(t_plot,[nociceptors.m17; nociceptors.h17; nociceptors.s17])
%% Defining the neurons

v2=vhold;
interneuron.v(1,1:2)=v2; %The interneuron in the spine
interneuron.m(1,1:2)=alpham(v2)/(alpham(v2)+betam(v2));
interneuron.h(1,1:2)=alphah(v2)/(alphah(v2)+betah(v2));
interneuron.n(1,1:2)=alphan(v2)/(alphan(v2)+betan(v2));

% initialize the neurotransmitter release times
Glu_release_start(1:2) = -100;
Glu_release_end(1:2) = -100;
Glu_release(1:2) = false;
v2(1:2) = vhold;
gAMPA(1:1001,1:2)=0;
%% Running the simulation
%initialize parameters of subsequent current pulse:
t1p(1:2)=10;       %starting time (ms)
t2p(1:2)=20;       %stopping time (ms)
ip(1)=25;        %current applied (muA)
ip(2)=0;
for klok=1:klokmax
  t=klok*dt;                      %note time
  
  %When the action potential arrives at the terminal and Glutamate is not
  %released in the last 8 ms, there will be a release of Glutamate in the 
  %next 10ms, lasting for 1ms to activate AMPA.
  for id=1:2
    if interneuron.v(klok,id) >= 0 && t > Glu_release_start(id) + 5
      Glu_release_start(id)=t+5;
      Glu_release_end(id)=Glu_release_start(id)+1;
      Glu_release(id)=true;
    else
      Glu_release(id)=false;
    end
  end
  for id = 1:2
    %Neuron 2: The interneuron
    gNMDA=0;
    if t>= Glu_release_start(3-id)+dt && Glu_release_start(3-id)>=0
      gAMPA(klok+1,id) = gAMPA_max*AMPA_probability(t-Glu_release_start(3-id));
      gNMDA = gNMDA+gNMDA_max*NMDA_probability(v2(id), t-Glu_release_start(3-id));
      gNMDA=0;
    else
      gAMPA(klok+1,id) = 0;
      gNMDA = 0;
    end
    interneuron.m(klok+1,id)=snew(interneuron.m(klok,id),alpham(v2(id)),betam(v2(id)),dt); %update m
    interneuron.h(klok+1,id)=snew(interneuron.h(klok,id),alphah(v2(id)),betah(v2(id)),dt); %update h
    interneuron.n(klok+1,id)=snew(interneuron.n(klok,id),alphan(v2(id)),betan(v2(id)),dt); %update n
    Isyn(klok+1,id) = gNMDA*ECa2+gAMPA(klok+1,id)*ENa2;
    %Isyn(klok+1) = gAMPA*ENa2;
    gNa=gNabar2*0.16*(interneuron.m(klok+1,id)^3)*interneuron.h(klok+1,id);%sodium conductance
    gK =gKbar2*(interneuron.n(klok+1,id)^4);%potassium conductance
    g=gNa+gK+gLbar2+gNMDA+gAMPA(klok+1,id);         %total conductance
    gE=gNa*ENa2+gK*EK2+gLbar2*EL2;         %gE=g*E
    %update v:
    v2(id)=(v2(id)+(dt/C2)*(gE+Isyn(klok+1,id)+izero(t1p(id), t2p(id), ip(id), t)))/(1+(dt/C2)*g);
    disp(izero(t1p(id), t2p(id), ip(id), t))
    interneuron.v(klok+1,id)=v2(id);
  end

end

%%
%subplot(2,1,1),plot(t_plot,v_plot);hold on
%subplot(2,1,2),plot(t_plot,mhn_plot);legend('m', 'h', 'n')
t_plot = 0:dt:tmax;
%plot(t_plot, nociceptors.v); hold on;
subplot(2,1,2);
plot(t_plot, interneuron.v); %plot(t_plot, motorNeuron.v);
legend('interneuron1', 'interneuron2');
xlabel('time (ms)');
ylabel('voltage (mV)')
%figure
%plot(t_plot,[nociceptors.m17; nociceptors.h17; nociceptors.s17])