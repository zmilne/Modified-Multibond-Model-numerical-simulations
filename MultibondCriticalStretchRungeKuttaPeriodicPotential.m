%MultibondCriticalStretchRungKutta.m for a *periodic* (cosine) substrate
%interaction potential
%Written by Zachary Banks Milne, University of Pennsylvania 
%Copyright 2018, Zachary Banks Milne

function [meanFf maxFf stdFf NoiseParamSub NoiseParamCant]=MultibondCriticalStretchPeriodicPotential(v,T,GamSub,GamCant,ender,NoiseMult,n,ksub,kcant,timeStep,aTimes)
clc
clear X x MeanFf BondState tFree Ff StretchAtDebonding

m=.0000000000001;%The mass should be very very small
k=ksub;K=kcant;%Cantilever and substrate spring constants
a=.2E-9;%The critical stretch length
kB=1.38E-23;w0=1E10;

gam=GamSub+GamCant;%A convenient damping parameter
ThetaT=GamCant/gam;%Another convenient damping parameter

NoiseParamSub=sqrt(2*GamSub*T*kB);%Noise amplitude from fluctuation-dissipation theorem for substrate damping
NoiseParamCant=sqrt(2*GamCant*T*kB)%Noise amplitude from fluctuation-dissipation theorem for cantilever damping
NoiseParamIndep=0;%If external noise is independent of damping

xc=2e-10; %critical stretch length
Eon=1.5E-20;%Activation energy (J) to form a bond

tStep=timeStep; %This is the working tStep
t=0:tStep:aTimes*a/(v); %Good for 100 slips

N=1;%number of particles (sites) should be one for periodic potential
x=zeros(length(t),length(v),N);%Particle positions

X(1)=0;%Initial mass position set to 0
dXdt(1)=0;%Initial mass velocity set to 0

for i=1:length(t)-1
 
    Ff(i)=K*X(i);%Friction force
    
    NoiseTot=NoiseMult*normrnd(0,NoiseParamCant)+NoiseMult*normrnd(0,NoiseParamSub)+NoiseMult*normrnd(0,NoiseParamIndep);%The noise amplitude comes from a random distribution with standard deviations determined by the noise parameters
     %4th-order Runge-Kutta method
    k1=tStep*1/m*(NoiseTot+1e-9*sin(2*pi*(x(i,:)-X(i))/xc)-gam*dXdt(i)-K*X(i)+GamSub*v);
    l1=tStep*dXdt(i);
    k2=tStep*1/m*(NoiseTot+1e-9*sin(2*pi*(-(X(i)+l1/2)+x(i,:))/xc)-gam*(dXdt(i)+k1/2)-K*(X(i)+l1/2)+GamSub*v);
    l2=tStep*(dXdt(i)+k1/2);
    k3=tStep*1/m*(NoiseTot+1e-9*sin(2*pi*(-(X(i)+l2/2)+x(i,:))/xc)-gam*(dXdt(i)+k2/2)-K*(X(i)+l2/2)+GamSub*v);
    l3=tStep*(dXdt(i)+k2/2);
    k4=tStep*1/m*(NoiseTot+1e-9*sin(2*pi*(-(X(i)+l3/2)+x(i,:))/xc)-gam*(dXdt(i)+k3)-K*(X(i)+l3)+GamSub*v);
    l4=tStep*(dXdt(i)+k3);
    dXdt(i+1)=dXdt(i)+1/6*(k1+2*k2+2*k3+k4);
    X(i+1)=X(i)+1/6*(l1+2*l2+2*l3+l4);
    
    
    for l=1:N
        x(i+1,l)=x(i,l)+v*tStep;
                
    end
    
    
    
end

figure
plot(t,K*(X))
title(['Ff vs puller position N=' num2str(N) 'T=' num2str(T) ' v=' num2str(v) ' mps'] )
savename = ['Ff_' num2str(T) 'K_N' num2str(N) 'Speed' num2str(v) 'nms' 'gamma' num2str(gam) 'RandMechNoiseNoBondNoise' 'tStep=' num2str(tStep) '.fig'];
saveas(gca,savename)

maxFf=max(Ff);
meanFf=mean(Ff(end-ender:end));
stdFf=std(Ff(end-ender:end));
end
