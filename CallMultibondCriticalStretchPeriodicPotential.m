%Call MultibondCriticalStretchRungKuttaPeriodicPotential.m Written by Zachary Banks Milne,
%University of Pennsylvania
%Copyright 2018, Zachary Banks Milne
close all
clear all

n=1;%THe number of interaction sites
timeStep=5e-10;%The time step. 5e-10 is the smallest time step I have used
a=.2e-9;%Critical stretch length
% velocity=linspace(1e-6,1e-3,10);%Several speeds can be tested on a linear
% scale velocity=logspace(-8,-5,10);%Several speeds can be tested on a log
% scale
velocity=1e-5*ones(1,1);%Or, one speed can be tested
gammaSub=[6e-6*ones(1,length(velocity))];%Substrate damping constant(s)
gammaCant=[6e-6*ones(1,length(velocity))];%Cantilever damping constant(s)
aTimes=6;%Used to either lengthen or shorten the time to run the simulation
TotalTimeIndices=round(aTimes*a./(velocity)/timeStep);%Total number of time indices to use
Ender=round(TotalTimeIndices*3/4);%How many indices, counting back from
%the last index, to use for the average Ff calculation. This should not go
%so far back as the initial stick, but should extend well into the kinetic
%(quasi-equilibrium) friction.
Temp=[273*ones(1,length(velocity))];%Temperature(s)
noiseMult=1e4;%The noise multiplier 'zeta'
ksub=1.3;kcant=10;%The substrate and cantilever spring constants.






for i=1:length(velocity)%Runs for each velocity
        [MeanFf(i) MaxFf(i) StdFf(i) NoiseParamSub(i) NoiseParamCant(i)]=MultibondCriticalStretchRungeKuttaPeriodicPotential(velocity(i),Temp(i),gammaSub(i),gammaCant(i),Ender(i),noiseMult,n,ksub,kcant,timeStep,aTimes);
end
CalcResultsPeriodicPotential%Calculate results Periodic Potential