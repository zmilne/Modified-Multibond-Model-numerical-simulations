%mMB_RK_noMass_Harmonic.m for a *harmonic* (spring) substrate
%interaction potential, removed mass so the ODE is first order.
%Written by Zachary Banks Milne, University of Pennsylvania 
%Copyright 2019, Zachary Banks Milne

function [Ff t meanFf maxFf stdFf TNotBonded NoiseParamSub NoiseParamCant firstSlipforce stdNoiseTot]=mMB_RK_noMass_Harmonic(v,v2,T,GamSub,GamCant,ender,NoiseMult,n,ksub,kcant,timeStep,aTimes,DoCorrelatedStickSlip)
clc
clear X x MeanFf BondState tFree Ff StretchAtDebonding

k=ksub;K=kcant;%Cantilever and substrate spring constants
a=.2E-9;%The critical stretch length
b=.7e-9;%lattice constant
kB=1.38E-23;w0=1E10;

gam=GamSub+GamCant;%A convenient damping parameter
ThetaT=GamCant/gam;%Another convenient damping parameter
tStep=timeStep; %This is the working tStep
% NoiseParamSub=sqrt(2*GamSub*T*kB);%(Original) Noise amplitude from fluctuation-dissipation theorem for substrate damping
% NoiseParamCant=sqrt(2*GamCant*T*kB)%(Original) Noise amplitude from fluctuation-dissipation theorem for cantilever damping
NoiseParamSub=sqrt(2*GamSub*T*kB*5e-10/tStep);%Noise amplitude NOT from fluctuation-dissipation theorem for substrate damping
NoiseParamCant=sqrt(2*GamCant*T*kB*5e-10/tStep);%Noise amplitude NOT from fluctuation-dissipation theorem for cantilever damping

NoiseParamIndep=0;%If external noise is independent of damping

xc=2e-10; %critical stretch length
% Eon=1.5E-20;%Activation energy (J) to form a bond (1.5e-20 J commonly used)
Eon=15E-20;%Activation energy (J) to form a bond, edit this

% tStep=timeStep; %This is the working tStep
t=0:tStep:aTimes*a/(v); %Good for 100 slips

N=n;%number of particles (sites)
x=zeros(length(t),length(v),N);%Site positions
x2=x;x3=x;
BondState=round(rand(1,N));%Start sites with random bonding state
tFree=tStep*round(rand(1,N));%Start all sites with a random time unbonded
tFreeCounter=zeros(1,N);%Counts how long a site has been NOT interacting
tBound=tStep*round(rand(1,N));%Start all sites with a random time bonded
tBoundCounter=zeros(1,N);%Counts how long a site has been interacting

kon=w0*exp(-Eon/(kB*T));%The rate of bonding
IsBonded=zeros(length(t),1);%Tracks whether a site is bonded
DoFirstFf=0;%If you want to only track MaxFf, the Ff at the first slip
if DoFirstFf==0
   firstSlipforce=111; 
end
X=zeros(length(t),1);
X(1)=0;%Initial mass position set to 0
dXdt(1)=0;%Initial mass velocity set to 0
syms Xhold1
Xhold2=0;
RandBond=exprnd(1/kon,length(t),N);
%Original noise multiplier (a.k.a. Zeta or NM)VVV
NoiseTot=(NoiseMult*normrnd(0,NoiseParamCant,length(t),1)+NoiseMult*normrnd(0,NoiseParamSub,length(t),1)+NoiseMult*normrnd(0,NoiseParamIndep,length(t),1));%The noise amplitude comes from a random distribution with standard deviations determined by the noise parameters
% stdNoiseTot=std(NoiseTot);
% NoiseTot=tStep*sqrt(5e-10)/(sqrt(tStep)*gam)*(NoiseMult*normrnd(0,NoiseParamCant,length(t),1)+NoiseMult*normrnd(0,NoiseParamSub,length(t),1)+NoiseMult*normrnd(0,NoiseParamIndep,length(t),1));%The noise amplitude comes from a random distribution with standard deviations determined by the noise parameters


%Experimental noise multiplier (a.k.a. Zeta or NM) VVV
% NoiseTot=(timeStep/5e-10)^-1*(NoiseMult*normrnd(0,NoiseParamCant,length(t),1)+NoiseMult*normrnd(0,NoiseParamSub,length(t),1)+NoiseMult*normrnd(0,NoiseParamIndep,length(t),1));%The noise amplitude comes from a random distribution with standard deviations determined by the noise parameters
%     NoiseTot=(timeStep/5e-10)*(NoiseMult*normrnd(0,NoiseParamCant,length(t),1)+NoiseMult*normrnd(0,NoiseParamSub,length(t),1)+NoiseMult*normrnd(0,NoiseParamIndep,length(t),1));%The noise amplitude comes from a random distribution with standard deviations determined by the noise parameters
% NoiseTot=(sqrt(timeStep/5e-10))^-1*(NoiseMult*normrnd(0,NoiseParamCant,length(t),1)+NoiseMult*normrnd(0,NoiseParamSub,length(t),1)+NoiseMult*normrnd(0,NoiseParamIndep,length(t),1));%The noise amplitude comes from a random distribution with standard deviations determined by the noise parameters
% i=0;
% while i<=length(t)-1&DoFirstFf
%     i=i+1;
for i=1:length(t)-1
    
%     if i>length(t)/2;
%         v=v2; %Can change the speed halfway through if you want to do velocity-stepping simulation for rate and state theory
%     end
    Ff(i)=K*X(i); %Friction force
    
    %NoiseTot=NoiseMult*normrnd(0,NoiseParamCant)+NoiseMult*normrnd(0,NoiseParamSub)+NoiseMult*normrnd(0,NoiseParamIndep);%The noise amplitude comes from a random distribution with standard deviations determined by the noise parameters
    %The original 4th-order Runge-Kutta method for the first-order ODE
%     k1=tStep/gam*(NoiseTot(i)+sum(k*(-X(i)+x(i,:)))-K*X(i)+GamSub*v);
%     k2=tStep/gam*(NoiseTot(i)+sum(k*(-(X(i)+k1/2)+x(i,:)))-K*(X(i)+k1/2)+GamSub*v);
%     k3=tStep/gam*(NoiseTot(i)+sum(k*(-(X(i)+k2/2)+x(i,:)))-K*(X(i)+k2/2)+GamSub*v);
%     k4=tStep/gam*(NoiseTot(i)+sum(k*(-(X(i)+k3)+x(i,:)))-K*(X(i)+k3)+GamSub*v);
%     X(i+1)=X(i)+1/6*(k1+2*k2+2*k3+k4);

% The corrected 4th-order Runge-Kutta method for the first-order ODE

if gam==0
    Xhold2=solve(tStep*(NoiseTot(i)+sum(k*(-Xhold1+x(i,:)))-K*Xhold1+GamSub*v),Xhold1);
    X(i+1)=Xhold2;
else
    k1=tStep/gam*(NoiseTot(i)+sum(k*(-X(i)+x(i,:)))-K*X(i)+GamSub*v);
    k2=tStep/gam*(NoiseTot(i)+sum(k*(-(X(i)+k1/2)+x2(i,:)))-K*(X(i)+k1/2)+GamSub*v);
    k3=tStep/gam*(NoiseTot(i)+sum(k*(-(X(i)+k2/2)+x2(i,:)))-K*(X(i)+k2/2)+GamSub*v);
    k4=tStep/gam*(NoiseTot(i)+sum(k*(-(X(i)+k3)+x3(i,:)))-K*(X(i)+k3)+GamSub*v);
    X(i+1)=X(i)+1/6*(k1+2*k2+2*k3+k4);
end
    
% %Take noise out of force and put in position
%  k1=tStep/gam*(sum(k*(-X(i)+x(i,:)))-K*X(i)+GamSub*v);
%     k2=tStep/gam*(sum(k*(-(X(i)+k1/2)+x2(i,:)))-K*(X(i)+k1/2)+GamSub*v);
%     k3=tStep/gam*(sum(k*(-(X(i)+k2/2)+x2(i,:)))-K*(X(i)+k2/2)+GamSub*v);
%     k4=tStep/gam*(sum(k*(-(X(i)+k3)+x3(i,:)))-K*(X(i)+k3)+GamSub*v);
%     X(i+1)=X(i)+1/6*(k1+2*k2+2*k3+k4)+4e-6*NoiseTot(i);

% %Remove noise completely
%  k1=tStep/gam*(sum(k*(-X(i)+x(i,:)))-K*X(i)+GamSub*v);
%     k2=tStep/gam*(sum(k*(-(X(i)+k1/2)+x2(i,:)))-K*(X(i)+k1/2)+GamSub*v);
%     k3=tStep/gam*(sum(k*(-(X(i)+k2/2)+x2(i,:)))-K*(X(i)+k2/2)+GamSub*v);
%     k4=tStep/gam*(sum(k*(-(X(i)+k3)+x3(i,:)))-K*(X(i)+k3)+GamSub*v);
%     X(i+1)=X(i)+1/6*(k1+2*k2+2*k3+k4);
    
    for l=1:N
        
        %RandBond=exprnd(1/kon,1,1);%Monte Carlo number to determine if a bond should form (per Barel and Urbakh)
        if BondState(l)==1 %If bonded, update the site position to match the puller position
%             x(i+1,l)=x(i,l)+v*tStep;
%             x(i+1,l)=v*t(i);
%             x2(i+1,l)=v*(t(i)+tStep/2);
%             x3(i+1,l)=v*(t(i)+tStep);
            x(i+1,l)=x(i,l)+v*tStep;
            x2(i+1,l)=x2(i,l)+v*tStep/2;
            x3(i+1,l)=x3(i,l)+v*tStep;
        end
        
        if (BondState(l)==0||abs(X(i+1)-x(i+1,l))>=xc)&tFree(l)<=RandBond(i,l) %If not bonded or if the stretch length exceeds Xc, set the bond position equal to the tip position
           
                        
            if abs(X(i+1)-x(i+1,l))>=xc&&DoFirstFf
                firstSlipforce=K*X(i);
                DoFirstFf=0;
               
                
            end
           
            tBound(l)=0;
            BondState(l)=0;
            x(i+1,l)=X(i+1);
            x2(i+1,l)=X(i+1);
            x3(i+1,l)=X(i+1);
            tFree(l)=tFree(l)+tStep;
            tFreeCounter(1,l)=tFreeCounter(1,l)+1;
        end
        
        if DoCorrelatedStickSlip==1
            if (tFree(l)>RandBond(i,l)&&mod(X(i+1),b)<=a) %If the rest time is longer than the calculated Arrhenius time, create a bond
            tFree(l)=0;
            BondState(l)=1;
            x(i+1,l)=X(i+1);
            x2(i+1,l)=X(i+1);
            x3(i+1,l)=X(i+1);
            end
        else
        if tFree(l)>RandBond(i,l) %If the rest time is longer than the calculated Arrhenius time, create a bond
            tFree(l)=0;
            BondState(l)=1;
            x(i+1,l)=X(i+1);
            x2(i+1,l)=X(i+1);
            x3(i+1,l)=X(i+1);
        end
        end
        
    end
    
    
    
end




doplots=1;
if doplots %If doplots=1 Ff vs time will be plotted and saved with a name that contains the parameters used
%     scatter(x(:,1),X)
%     xlabel('x');ylabel('X');
    plot(t,K*(X))
    title(['Ff vs puller position N=' num2str(N) 'T=' num2str(T) ' v=' num2str(v) ' mps'] )
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
% meanFf=mean(Ff(end-ender:end));% (original) The average Ff during quasi-equilibrium (kinetic Ff)
% stdFf=std(Ff(end-ender:end));% (original) The standard deviation of Ff during quasi-equilibrium (kinetic Ff)
meanFf=mean(Ff);% (New) The average Ff during quasi-equilibrium (kinetic Ff)
stdFf=std(Ff);% (New) The standard deviation of Ff during quasi-equilibrium (kinetic Ff)

end
