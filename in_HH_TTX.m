%filename: in_HH.m
%Units:
%voltage is in millivolts (mV)
%current is in microamperes (muA)
%length is in centimeters (cm)
%time is in milliseconds (ms)
%note that microfarads =
%microamperes*milliseconds/millivolt
%= muA*ms/mV
%
%initialize membrane parameters:
%Constants for neuron 1
%membrane capacitance per unit area:
C1=0.9361      %(muF/cm^2)
%max possible Na+ conductance per unit area:
gNabar17=18/100 %((muA/mV)/cm^2)
%max possible Na+ conductance per unit area:
gNabar18=7 %((muA/mV)/cm^2)
%max possible K+ conductance per unit area:
gKbar1=4.78   %((muA/mV)/cm^2)
%max possible K+ conductance per unit area:
gKbar1A=8.33   %((muA/mV)/cm^2)
%leakage conductance per unit area:
gLbar1=0.0575  %((muA/mV)/cm^2)
%Na+ equilibrium potential:
ENa1 = 67.1   %(mV)
%K+ equilibrium potential:
EK1 = -84.7   %(mV)
%leakage channel reversal potential:
EL1 = -58.91   %(mV)
%parameters for neuron 2
%membrane capacitance per unit area:
C2=2      %(muF/cm^2)
%max possible Na+ conductance per unit area:
gNabar2=10/10 %((muA/mV)/cm^2)
%max possible K+ conductance per unit area:
gKbar2=5   %((muA/mV)/cm^2)
%leakage conductance per unit area:
gLbar2=0.1  %((muA/mV)/cm^2)
%Na+ equilibrium potential:
ENa2 = 55   %(mV)
%K+ equilibrium potential:
EK2 = -80   %(mV)
%leakage channel reversal potential:
EL2 = -68   %(mV)
%Ca++ equilibrium potential:
ECa2 = 137   %(mV)

%parameters for neuron 3
%membrane capacitance per unit area:
C3=1      %(muF/cm^2)
%max possible Na+ conductance per unit area:
gNabar3=120/10 %((muA/mV)/cm^2)
%max possible K+ conductance per unit area:
gKbar3=80   %((muA/mV)/cm^2)
%leakage conductance per unit area:
gLbar3=0.1  %((muA/mV)/cm^2)
%Na+ equilibrium potential:
ENa3 = 55   %(mV)
%K+ equilibrium potential:
EK3 = -75   %(mV)
%leakage channel reversal potential:
EL3 = -68   %(mV)

%nAChR conductance per unit area:
gnAChR_max = 2 %((muA/mV)/cm^2)
%AMPA conductance per unit area:
gAMPA_max = 7 %((muA/mV)/cm^2)
%AMPA conductance per unit area:
gNMDA_max = 3.5 %((muA/mV)/cm^2)

%initialize time step and experiment duration:
dt=0.1     %time step duration (ms)
tmax=200    %duration of experiment (ms)
%total number of time steps in the experiment:
klokmax=ceil(tmax/dt)
%
%initialize arrays that hold data for plotting:
mhn_plot=zeros(3,klokmax);
v_plot=zeros(1,klokmax);
t_plot=zeros(1,klokmax);
%
%initialize parameters that define the experiment:
%the neuron is at rest (v= -70 mV) prior to t=0;
%at t=0 a current shock is applied after which v= -55 mV;
%then a subsequent 15 muA current pulse of 1 ms duration
%is applied beginning at t=10 ms.
%voltage prior to t=0:
vhold=  -70  %(mV)

%
%initialize checking parameter
check=0      %set check=1 to enable self-checking
             %set check=0 to disable self-checking
