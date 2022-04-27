%filename: in_HH.m
%Units:
%voltage is in millivolts (mV)
%current is in microamperes (muA)
%length is in centimeters (cm)
%time is in milliseconds (ms)
%note that microfarads =
%microamperes*milliseconds/millivolt
%= muA*ms/mV
%sodium more likely, potassium less likely
%initialize membrane parameters:
%membrane capacitance per unit area:
C=1.0      %(muF/cm^2)
%max possible Na+ conductance per unit area:
gNabar=120 %((muA/mV)/cm^2)
%max possible K+ conductance per unit area:
gKbar=36   %((muA/mV)/cm^2)
%leakage conductance per unit area:
gLbar=0.3  %((muA/mV)/cm^2)

gSynbar = 0.3; %0.33 min

%Na+ equilibrium potential:
ENa = 45   %(mV)
%K+ equilibrium potential:
EK = -82   %(mV)
%leakage channel reversal potential:
EL = -59   %(mV)

%excitatory >Vrest, shunting inhibition =Vrest, inhibitory <Vrest
ESyn = 0; %0, -55, -100 (?)

%initialize time step and experiment duration:
dt=0.1     %time step duration (ms)
tmax=25    %duration of experiment (ms)
%total number of time steps in the experiment:
klokmax=ceil(tmax/dt)
%
%initialize arrays that hold data for plotting:

%
%initialize parameters that define the experiment:
%the neuron is at rest (v= -70 mV) prior to t=0;
%at t=0 a current shock is applied after which v= -55 mV;
%then a subsequent 15 muA current pulse of 1 ms duration
%is applied beginning at t=10 ms.
%voltage prior to t=0:
vhold=  -70  %(mV)
%voltage just after t=0:
vstart= -70  %(mV)
%(change in v is result of current shock applied at t=0)
%
%initialize parameters of subsequent current pulse:
t1p=4      %starting time (ms)
t2p=6       %stopping time (ms)
ip=18      %current applied (muA), 28 min for action potential
%
%initialize checking parameter
check=0      %set check=1 to enable self-checking
             %set check=0 to disable self-checking
