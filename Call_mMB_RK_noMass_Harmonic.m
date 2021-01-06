 %Call MultibondCriticalStretchRungKuttaHarmonicPotential.m Written by Zachary Banks Milne,
%University of Pennsylvania
%Copyright 2018, Zachary Banks Milne
clc
close all
clearvars -except FfTemp04112019 FfTempNM5e5n1nc6n6ns6n6 FfTempNM5e4n1n1nc6n6ns6n6 FfTempNM5e3n1nc6n6ns6n6 FfTempNM5e3n1nc0ns6n6 FfTempNM5e4n1nc0ns6n6 FfTempNM5e4n1nc6n6ns6n6 FfTemp04102019
%%
doCorrelatedStickSlip=1; %Bonds will only form whenever a site from the upper surface is close to a site from the lower surface
ns=1e-6;
nc=1e-6;
Z=1e4;%The noise multiplier 'zeta'
n=1;%THe number of interaction sites
v=1e-4*ones(1,1);%Or, one speed can be tested
timeStep=5e-10;%The time step. 5e-10 is the smallest time step I have used
% ksub=40;kcant=1;%The substrate and cantilever spring constants.timeStep=5e-9;%The time step. 5e-10 is the smallest time step I have used
a=.1e-9;%Critical stretch length
% velocity=linspace(1e-6,1e-3,10);%Several speeds can be tested on a linear
% scale velocity=logspace(-8,-5,10);%Several speeds can be tested on a log
% scale
V2=v;
gammaSub=[ns*ones(1,length(v))];%Substrate damping constant(s)
gammaCant=[nc*ones(1,length(v))];%Cantilever damping constant(s)
aTimes=round(2500*(v).^.5/(1e-3)^.5);%Used to either lengthen or shorten the time to run the simulation
TotalTimeIndices=round(aTimes*a./(v)/timeStep);%Total number of time indices to use
Ender=round(TotalTimeIndices*3/4);%How many indices, counting back from
%the last index, to use for the average Ff calculation. This should not go
%so far back as the initial stick, but should extend well into the kinetic
%(quasi-equilibrium) friction.
Temp=[300*ones(1,length(v))];%Temperature(s)
ksub=20;kcant=10;%The substrate and cantilever spring constants.
% FfTemp04112019
% figure
% hold on
for i=1:length(v)%Runs for each velocity
    [FF t MeanFf(i) MaxFf(i) StdFf(i) tNotBonded(i) NoiseParamSub(i) NoiseParamCant(i) FirstSlipForce(i)]=mMB_RK_noMass_Harmonic(v(i),V2(i),Temp(i),gammaSub(i),gammaCant(i),Ender(i),Z,n,ksub,kcant,timeStep,aTimes,doCorrelatedStickSlip);
end
% legend('large tStep','small tStep')
CalcResultsHarmonicPotential%Calculate results Harmonic Potential
%%
gatherdata=0;
if gatherdata
FfTemp04112019=[FfTemp04112019;C];
% close all
figure
scatter(FfTemp04102019(:,1),FfTemp04102019(:,15),'filled')
hold on
scatter(FfTemp04112019(:,1),FfTemp04112019(:,15),'filled')
end

SelectRegion=1;
if SelectRegion
   t(end)=[];
    ValueStart=ginput(1); TimeValue=ValueStart(1); TimeToAvg=t(t>=TimeValue);
    AvgFf=mean(FF(t>=TimeValue))
    StdFf=std(FF(t>=TimeValue))
    StatsFf=[AvgFf StdFf];
end