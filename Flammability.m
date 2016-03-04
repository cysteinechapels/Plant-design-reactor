function [selectivity, O2perc,i]=Flammability(I)
Po=180;
T=446;
safety=0.8;

MM=[28.0532,31.9988,60.052,18.0153, 16.04, 86.0892,44.0095,30.069,39.948,28.0134];
density=[0.650502,0.725865,1.73532,0.435758,0.3637,2.48941,1.01799,0.702875,0.905922,0.633822];

%Calculation of inerts using volume
FvolE = I(1).*MM(1)/453.59237./density(1);
FvolCH4=I(4).*MM(5)/453.59237./density(5);

%WTF
Feth = (FvolE*245/1000000+FvolCH4*800/1000000)*density(8)*453.59237/MM(8)
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
    Check = abs(Drysum-Drynew)

    if Check <= 0.001
        break
    end
end
i;
F0=[I(1) FO2 I(2) I(3) I(4) 0 0 Feth FAr FN2];

p=F0./sum(F0)*Po;
rx= 1409*exp(-3674/T)*p(2)*p(1)*p(3)*(1+1.17*p(4))/((1+0.583*p(2)*(1+1.7*p(4)))*(1+6.8*p(3)))/3600;
ry=9625*exp(-4250/T)*p(2)*(1+0.68*p(4))/(1+0.76*p(2)*(1+0.68*p(4)))/3600;


O2perc = FO2/sum(F0);

selectivity = rx/ry;