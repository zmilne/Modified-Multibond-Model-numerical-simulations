%Call MultibondCriticalStretchRungKuttaPeriodicPotential.m Written by Zachary Banks Milne,
%University of Pennsylvania
%Copyright 2018, Zachary Banks Milne
close all
clear all

ns=6e-6;
nc=6e-6;
Z=0.5;%The noise multiplier 'zeta'
v=1e-3*ones(1,1);%Or, one speed can be tested

ksub=.5;kcant=1;%The substrate and cantilever spring constants.
n=1;%THe number of interaction sites
timeStep=5e-9;%The time step. 5e-10 is the smallest time step I have used
a=.1e-9;%Critical stretch length
% velocity=linspace(1e-6,1e-3,10);%Several speeds can be tested on a linear
% scale velocity=logspace(-8,-5,10);%Several speeds can be tested on a log
% scale

gammaSub=[ns*ones(1,length(v))];%Substrate damping constant(s)
gammaCant=[nc*ones(1,length(v))];%Cantilever damping constant(s)
aTimes=round(2500*(v).^.5/(1e-3)^.5);%Used to either lengthen or shorten the time to run the simulation
TotalTimeIndices=round(aTimes*a./(v)/timeStep);%Total number of time indices to use
Ender=round(TotalTimeIndices*3/4);%How many indices, counting back from
%the last index, to use for the average Ff calculation. This should not go
%so far back as the initial stick, but should extend well into the kinetic
%(quasi-equilibrium) friction.
Temp=[273*ones(1,length(v))];%Temperature(s)

PeriodicAmplitudeOffset=1e-9; %Add to force to make the periodic force never go negative 





for i=1:length(v)%Runs for each velocity
        [FF MeanFf(i) MaxFf(i) StdFf(i) NoiseParamSub(i) NoiseParamCant(i) t]=mMB_RK_noMass_Periodic(v(i),Temp(i),gammaSub(i),gammaCant(i),Ender(i),Z,n,ksub,kcant,timeStep,aTimes,PeriodicAmplitudeOffset);
end
noiseMult=Z;velocity=v;
CalcResultsPeriodicPotential%Calculate results Periodic Potential

SelectRegion=1;
if SelectRegion
   t(end)=[];
    ValueStart=ginput(1); TimeValue=ValueStart(1); TimeToAvg=t(t>=TimeValue);
    AvgFf=mean(FF(t>=TimeValue));
    StdFf=std(FF(t>=TimeValue));
    StatsFf=[AvgFf StdFf];
end
