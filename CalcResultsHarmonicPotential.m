%Calculates results from MultibondCriticalStretchRungKutta.m for a *harmonic*
%potential for the substrate interaction
%Written by Zachary Banks Milne, University of Pennsylvania
%Copyright 2018, Zachary Banks Milne
format short
velocity=v;
noiseMult=Z;
if length(noiseMult)==1

C=[Temp' noiseMult*ones(length(Temp),1) n*ones(length(Temp),1) gammaSub' gammaCant' gammaSub'+gammaCant' gammaCant'./(gammaSub'+gammaCant')  ksub*ones(length(Temp),1) kcant*ones(length(Temp),1)  zeros(length(Temp),1) noiseMult*ones(length(Temp),1).*sqrt(2*1.38e-23*Temp'.*(gammaSub'+gammaCant')) 50*0.0000000003./velocity timeStep*ones(length(Temp),1) velocity MeanFf(1:length(Temp))' StdFf(1:length(Temp))' tNotBonded(1:length(Temp))' MaxFf' FirstSlipForce']
meanFirstSlipForce=mean(C(:,end))
StdFirstSlipForce=std(C(:,end))
% meanOfMeanFf=mean(MeanFf)
% stdOfMeanFf=std(MeanFf)
% close all
else
    C=[Temp' noiseMult n*ones(length(Temp),1) gammaSub' gammaCant' gammaSub'+gammaCant' gammaCant'./(gammaSub'+gammaCant')  ksub*ones(length(Temp),1) kcant*ones(length(Temp),1)  zeros(length(Temp),1) noiseMult.*sqrt(2*1.38e-23*Temp'.*(gammaSub'+gammaCant')) 50*0.0000000003./velocity timeStep*ones(length(Temp),1) velocity MeanFf(1:length(Temp))' StdFf(1:length(Temp))' tNotBonded(1:length(Temp))' MaxFf' FirstSlipForce']
meanFirstSlipForce=mean(C(:,end))
StdFirstSlipForce=std(C(:,end))
MaxFf=MaxFf';
FirstSlipForce=FirstSlipForce';
% close all
end 
