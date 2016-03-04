function [Fva,F,Fr,F0,Vcat,L,A,vtot0] = SteadyState1Tube(I)
% I = [ ethylene, acetic, water, ch4, P, T, Tube, Volume, ID]

% Species order
% 1-ethylene, 2-oxygen, 3-acetic acid, 4-water, 5-CH4, 6-VAM, 7-CO2, 8-Eth,
%   9-Argon, 10 - N2 11 -Pressure

%Input order
% 1-ethylene 2-acetic acid 3-water 4-ch4 5-P, 6-T, 7-Tube #, 8-Volume cat max, 9-ID
% 
safety= 0.8;
purge=0.3;
recoveryAA=0.5;
Rspec = [I(5) I(6)]; % Rspec = [ P T Tube Volume ID ]
Po= Rspec(1);%psia
T= Rspec(2); %Kelvin

%Densities at 349 F and 180 psig
density=[0.641861,0.717219,2.1532,0.436983,0.360817,2.50445,1.01649,0.692685,0.895527,0.625938];

global MM;
MM=[28.0532,31.9988,60.052,18.0153, 16.04, 86.0892,44.0095,30.069,39.948,28.0134];


%Calculation of inerts using volume %
FvolE = I(1).*MM(1)/453.59237./density(1);
FvolCH4=I(4).*MM(5)/453.59237./density(5);

%WTF
Feth = (FvolE*245/1000000+FvolCH4*800/1000000)*density(8)*453.59237/MM(8); 
FO2 = 0;
% FO2 = I(4);
FAr = 0;
FN2 = 0;
%     Oxyinerts = [1 -1 -1; -0.0011 1 0; -0.0005 0 1];
%     Totalinerts = [FO2;0;0];
%     ArandN2 = Oxyinerts\Totalinerts;
%     FAr = ArandN2(2)*density(9)*453.59237/MM(9);
%     FN2 = ArandN2(3)*density(10)*453.59237/MM(10);
% 
% %Iteratetion for getting methane
% Fmeth = 0;
% for i=1:100000
%         Drysum = I(1)+Fmeth+Feth+FO2+FAr+FN2;
%     % This is where bugs go to fuck and make more bugs
%     % Do you want bugs cause this is how you get bugs
%     PandT = 20.5-0.02*Po-0.02*(9/5*(T - 273) + 32);
%     Ethane = 0.08*Feth/Drysum*100;
%     Ethylene= 0.01*I(1)/Drysum*100;
%     Fmeth = (FO2*100/Drysum+3-PandT-Ethane-Ethylene)/0.04*100*Drysum;
%         Drynew = I(1)+I(4)+Feth+Fmeth+FAr+FN2;
%     Check = abs(Drysum-Drynew);
% 
%     if Check <= 10
%         break
%     end
%  end
% Iteration for getting Oxygen flow properly
for i = 1:100
    Drysum = I(1)+I(4)+Feth+FO2+FAr+FN2;

    % This is where bugs go to fuck and make more bugs
    % Do you want bugs cause this is how you get bugs
    PandT = 20.5-0.02*Po-0.02*(9/5*(T - 273) + 32);
    Methane = 0.04*I(4)/Drysum*100;
    Ethane = 0.08*Feth/Drysum*100;
    Ethylene= 0.01*I(1)/Drysum*100;
    FLO2=((PandT+Methane+Ethane+Ethylene)-3)*(Drysum)/100;

    FO2=FLO2*safety;
    FO2vol=FO2*MM(2)/453.59237/density(2);
    Oxyinerts = [1 -1 -1; -0.0011 1 0; -0.0005 0 1];
    Totalinerts = [FO2vol;0;0];
    ArandN2 = Oxyinerts\Totalinerts;
    FAr = ArandN2(2)*density(9)*453.59237/MM(9);
    FN2 = ArandN2(3)*density(10)*453.59237/MM(10);
    Drynew = I(1)+I(4)+Feth+FO2+FAr+FN2;
    Check = abs(Drysum-Drynew);

    if Check <= 0.001
        break
    end
end

F0=[I(1) FO2 I(2) I(3) I(4) 0 0 Feth FAr FN2];
Fr=zeros(1,10);


for n=1:100
    F=F0+Fr;
    New=Solver1Tube(F,Rspec);
    Fr=New(end,1:end-1);
    Fr([4 6])=0; %taking out heavy ends
    Fr(3)=Fr(3)*recoveryAA; %recycled AA from azeo column
    Fr(7)=Fr(7)*0.01; %CO2 treatment
    Fr=Fr*(1-purge); %purge stream
    Fnew=Fr+F0;
    sum(F-Fnew);
    if abs(sum(F-Fnew))<10
        break
    end
end

%F0;
%Fr;
[Fend,F,Fva, Vcat,L,A,vtot0]=Solver1Tube(Fnew,Rspec);

end

