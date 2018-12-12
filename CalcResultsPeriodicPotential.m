%Calculates results from MultibondCriticalStretchRungKutta.m for a
%*periodic* potential for the substrate interaction Written by Zachary
%Banks Milne, University of Pennsylvania Copyright 2018, Zachary Banks
%Milne
C=[Temp' noiseMult*ones(length(Temp),1) n*ones(length(Temp),1) gammaSub' gammaCant' gammaSub'+gammaCant' gammaCant'./(gammaSub'+gammaCant')  ksub*ones(length(Temp),1) kcant*ones(length(Temp),1)  zeros(length(Temp),1) noiseMult*ones(length(Temp),1).*sqrt(2*1.38e-23*Temp'.*(gammaSub'+gammaCant')) 50*0.0000000003./velocity' timeStep*ones(length(Temp),1) velocity' MeanFf(1:length(Temp))' StdFf(1:length(Temp))' MaxFf']