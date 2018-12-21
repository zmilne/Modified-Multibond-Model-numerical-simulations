%Call MultibondCriticalStretchRungKuttaHarmonicPotential.m Written by Zachary Banks Milne,
%University of Pennsylvania
%Copyright 2018, Zachary Banks Milne

close all
clear all

n=2;%THe number of interaction sites
timeStep=5e-9;%The time step. 5e-10 is the smallest time step I have used
a=.2e-9;%Critical stretch length
% velocity=linspace(1e-6,1e-3,10);%Several speeds can be tested on a linear
% scale velocity=logspace(-8,-5,10);%Several speeds can be tested on a log
% scale
velocity=1e-7*ones(1,1);%Or, one speed can be tested
V2=velocity;%The velocity which the code uses halfway through the total
%time it is set to run for. This can be changed for velocity-stepping
%experiments such as performed in rate-and-state research
gammaSub=[6e-6*ones(1,length(velocity))];%Substrate damping constant(s)
gammaCant=[0*ones(1,length(velocity))];%Cantilever damping constant(s)
aTimes=2;%Used to either lengthen or shorten the time to run the simulation
TotalTimeIndices=round(aTimes*a./(velocity)/timeStep);%Total number of time indices to use
Ender=round(TotalTimeIndices*3/4);%How many indices, counting back from
%the last index, to use for the average Ff calculation. This should not go
%so far back as the initial stick, but should extend well into the kinetic
%(quasi-equilibrium) friction.
Temp=[300*ones(1,length(velocity))];%Temperature(s)
noiseMult=5e4;%The noise multiplier 'zeta'
ksub=1.3;kcant=1;%The substrate and cantilever spring constants.






for i=1:length(velocity)%Runs for each velocity
    [MeanFf(i) MaxFf(i) StdFf(i) tNotBonded(i) NoiseParamSub(i) NoiseParamCant(i) FirstSlipForce(i)]=MultibondCriticalStretchRungeKuttaHarmonicPotential(velocity(i),V2,Temp(i),gammaSub(i),gammaCant(i),Ender(i),noiseMult,n,ksub,kcant,timeStep,aTimes);
end
CalcResultsHarmonicPotential%Calculate results Harmonic Potential