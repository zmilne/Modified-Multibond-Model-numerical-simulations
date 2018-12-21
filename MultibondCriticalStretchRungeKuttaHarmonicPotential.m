%MultibondCriticalStretchRungKutta.m for a *harmonic* (spring) substrate
%interaction potential
%Written by Zachary Banks Milne, University of Pennsylvania 
%Copyright 2018, Zachary Banks Milne

function [meanFf maxFf stdFf TNotBonded NoiseParamSub NoiseParamCant firstSlipforce]=MultibondCriticalStretchRungeKuttaHarmonicPotential(v,v2,T,GamSub,GamCant,ender,NoiseMult,n,ksub,kcant,timeStep,aTimes)
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

N=n;%number of particles (sites)
x=zeros(length(t),length(v),N);%Site positions
BondState=round(rand(1,N));%Start sites with random bonding state
tFree=tStep*round(rand(1,N));%Start all sites with a random time unbonded
tFreeCounter=zeros(1,N);%Counts how long a site has been NOT interacting
tBound=tStep*round(rand(1,N));%Start all sites with a random time bonded
tBoundCounter=zeros(1,N);%Counts how long a site has been interacting

kon=w0*exp(-Eon/(kB*T));%The rate of bonding
IsBonded=zeros(length(t),1);%Tracks whether a site is bonded
DoFirstFf=1;%If you want to only track MaxFf, the Ff at the first slip


X(1)=0;%Initial mass position set to 0
dXdt(1)=0;%Initial mass velocity set to 0

for i=1:length(t)-1
    if i>length(t)/2;
        v=v2; %Can change the speed halfway through if you want to do velocity-stepping simulation for rate and state theory
    end
    Ff(i)=K*X(i); %Friction force
    
    NoiseTot=NoiseMult*normrnd(0,NoiseParamCant)+NoiseMult*normrnd(0,NoiseParamSub)+NoiseMult*normrnd(0,NoiseParamIndep);%The noise amplitude comes from a random distribution with standard deviations determined by the noise parameters
    %4th-order Runge-Kutta method
    k1=tStep*1/m*(NoiseTot+sum(k*(-X(i)+x(i,:)))-gam*dXdt(i)-K*X(i)+GamSub*v);
    l1=tStep*dXdt(i);
    k2=tStep*1/m*(NoiseTot+sum(k*(-(X(i)+l1/2)+x(i,:)))-gam*(dXdt(i)+k1/2)-K*(X(i)+l1/2)+GamSub*v);
    l2=tStep*(dXdt(i)+k1/2);
    k3=tStep*1/m*(NoiseTot+sum(k*(-(X(i)+l2/2)+x(i,:)))-gam*(dXdt(i)+k2/2)-K*(X(i)+l2/2)+GamSub*v);
    l3=tStep*(dXdt(i)+k2/2);
    k4=tStep*1/m*(NoiseTot+sum(k*(-(X(i)+l3)+x(i,:)))-gam*(dXdt(i)+k3)-K*(X(i)+l3)+GamSub*v);
    l4=tStep*(dXdt(i)+k3);
    dXdt(i+1)=dXdt(i)+1/6*(k1+2*k2+2*k3+k4);
    X(i+1)=X(i)+1/6*(l1+2*l2+2*l3+l4);
    
    
    for l=1:N
        
        RandBond=exprnd(1/kon,1,1);%Monte Carlo number to determine if a bond should form (per Barel and Urbakh)
        if BondState(l)==1 %If bonded, update the site position to match the puller position
            x(i+1,l)=x(i,l)+v*tStep;
        end
        
        if (BondState(l)==0||abs(X(i+1)-x(i+1,l))>=xc)&tFree(l)<=RandBond %If not bonded or if the stretch length exceeds Xc, set the bond position equal to the tip position
            if abs(X(i+1)-x(i+1,l))>=xc&&DoFirstFf
                firstSlipforce=K*X(i);
                DoFirstFf=0;
                
            end
            
            tBound(l)=0;
            BondState(l)=0;
            x(i+1,l)=X(i+1);
            tFree(l)=tFree(l)+tStep;
            tFreeCounter(1,l)=tFreeCounter(1,l)+1;
        end
        
        if tFree(l)>RandBond %If the rest time is longer than the calculated Arrhenius time, create a bond
            tFree(l)=0;
            BondState(l)=1;
            x(i+1,l)=X(i+1);
        end
        
    end
    
    
    
end




doplots=1;
if doplots %If doplots=1 Ff vs time will be plotted and saved with a name that contains the parameters used
    figure
    plot(t,K*(X))
%     title(['Ff vs puller position N=' num2str(N) 'T=' num2str(T) ' v=' num2str(v) ' mps'] )
%     savename = ['Ff_' num2str(T) 'K_N' num2str(N) 'Speed' num2str(v) 'nms' 'gamma' num2str(gam) 'tStep=' num2str(tStep) '.fig'];
%     saveas(gca,savename)
end

% figure plot(v*t,IsBonded)


% figure
tFreeAvg=tFreeCounter/length(t);%Calculate the average time not bonded for each site
TNotBonded=mean(tFreeAvg);%Calculate the average time not bonded for all sites
% scatter(1:1:N,tFreeAvg) savename = ['tFree_' num2str(T) 'K_N' num2str(N)
% 'Speed' num2str(v) ' mps' 'gamma' num2str(gam) 'RandMechNoiseNoBondNoise'
% 'tStep=' num2str(tStep) '.fig']; saveas(gca,savename)
maxFf=max(Ff);%The maximum Ff
meanFf=mean(Ff(end-ender:end));% The average Ff during quasi-equilibrium (kinetic Ff)
stdFf=std(Ff(end-ender:end));% The standard deviation of Ff during quasi-equilibrium (kinetic Ff)

end
