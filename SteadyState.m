function [Fva,F,Fr,F0,Vcat,L,A,vtot0] = SteadyState(I)
%I don't know what I'm doing
%   1-ethylene, 2-oxygen, 3-acetic acid, 4-water, 5-CH4, 6-VAM, 7-CO2, 8-Eth,
%   9-Argon, 10 - N2 11 -Pressure

%6 - P, 7- T, 8 -
% Tube #, 9-Volume cat max, 10 - ID
safety= 0.8;
purge=1;
recoveryAA=0.6;
Rspec=I(5:9);
Po= Rspec(1);%psia
T= Rspec(2); %Kelvin

density=[0.650502,0.725865,1.73532,0.435758,0.3637,2.48941,1.01799,0.702875,0.905922,0.633822];

MM=[28.0532,31.9988,60.052,18.0153, 16.04, 86.0892,44.0095,30.069,39.948,28.0134];


%Calculation of inerts using volume %
FvolE = I(1).*MM(1)/453.59237./density(1);
FvolCH4=I(4).*MM(5)/453.59237./density(5);

Feth = (FvolE*245/1000000+FvolCH4*800/1000000)*density(8)*453.59237/MM(8); 

Drysum = I(1)+I(4);

FLO2=((20.5-0.02*Po-0.02*(T-459.67)/(5/9)+0.04*I(4)/(Drysum+I(4))*100+0.08*Feth/(Drysum+I(4))*100+0.01*100*I(1)/(Drysum+I(4)))-3)*(Drysum+I(4))/100;

FO2=FLO2*safety;

Oxyinerts = [1 -1 -1; -0.11 1 0; -0.05 0 1];
Totalinerts = [FO2;0;0];
ArandN2 = Oxyinerts\Totalinerts;
FAr = ArandN2(2)*density(9)*453.59237/MM(9);
FN2 = ArandN2(3)*density(10)*453.59237/MM(10);

F0=[I(1) FO2 I(2) I(3) I(4) 0 0 Feth FAr FN2];
Fr=zeros(1,10);

% 
% for n=1:100
%     F=F0+Fr;
%     New=Solver(F,Rspec);
%     Fr=New(end,1:end-1);
%     Fr([4 6])=0; %taking out heavy ends
%     Fr(3)=Fr(3)*recoveryAA; %recycled AA from azeo column
%     Fr(7)=Fr(7)*0.01; %CO2 treatment
%     Fr=Fr*(1-purge); %purge stream
%     Fnew=Fr+F0;
%     sum(F-Fnew);
%     if abs(sum(F-Fnew))<10
%         break
%     end
% end
% 
% F0;
% Fr;
[Fend,F,Fva, Vcat,L,A,vtot0]=Solver(F0,Rspec);

end

