function [Fend,F,Fva,Vcat,L,A,vtot0] = Solver (F0,Rspec)
%   1-ethylene, 2-oxygen, 3-acetic acid, 4-water, 5-CH4, 6-VAM, 7-CO2, 8-Eth,
%   9-Argon, 10 - N2 11 -Pressure

Volcat=[0 Rspec(4)]; %volume of catalyst, liters
T= Rspec(2); %Kelvin
D= Rspec(5)/12; %feet
A=D^2/4*pi()*Rspec(3); %total cross-sectional area of reactor
Po= Rspec(1);%psia


%Density at 180psig nd 340.1F]
density=[0.650502,0.725865,1.73532,0.435758,0.3637,2.48941,1.01799,0.702875,0.905922,0.633822];
global MM   
%MM=[28.0532,31.9988,60.052,18.0153, 16.04, 86.0892,44.0095,30.069,39.948,28.0134];

%Calculation of inerts using volume %
% FvolE = I(1).*MM(1)/453.59237./density(1);
% FvolAA = I(2).*MM(3)/453.59237./density(3);
% FvolW = I(3).*MM(4)/453.59237./density(4);
% FvolCH4=I(4).*MM(5)/453.59237./density(5);
% 
% Feth = (FvolE(1)*245/1000000+FvolCH4*800/1000000)*density(8)*453.59237/MM(8); 
% 
% Drysum = FvolE+FvolCH4;
% 
% FLO2=((20.5-0.02*Po-0.02*(T-459.67)/(5/9)+0.04*I(4)/(Drysum+I(4))*100+0.08*Feth/(Drysum+I(4))*100+0.01*100*I(1)/(Drysum+I(4)))-3)*(Drysum+I(4))/100;
% 
% FO2=FLO2*safety;

% Oxyinerts = [1 -1 -1; -0.11 1 0; -0.05 0 1];
% Totalinerts = [FO2;0;0];
% ArandN2 = Oxyinerts\Totalinerts;
% FAr = ArandN2(2)*density(9)*453.59237/MM(9);
% FN2 = ArandN2(3)*density(10)*453.59237/MM(10);

% %Broken Methane calcs
% % Drysum = sum(I(1:3))+Feth+FAr+FN2;
% % 
% % flame = @(FCH4) Flammability(FCH4,FO2,I(1),Feth,Drysum,T,Po);
% % FCH40 = 100;
% % FCH4 = fminsearch(flame, FCH40);

% F0=[I(1) FO2 I(2) I(3) I(4) 0 0 Feth FAr FN2];

Ft0=sum(F0);

%Convert gmol to lbmol flow for velocity calculations
F0lbmol=F0./453.59237;
Massflow=F0lbmol.*MM;
appdensity= sum(Massflow.*density)/sum(Massflow);
vtot0=sum(Massflow)/appdensity./A; %average velocity

F0=[F0 Po];

func = @(V,F) myReactor(V,F,A,T,Po,vtot0,Ft0);
opts=odeset('Events',@events);

[Vcat,F] = ode45(func,Volcat,F0,opts);

%Events used to stop integration
    function[value,isterminal,direction]= events(V,F)
        value = [F(1),F(2), F(3), F(4), F(5), F(6), F(7), F(8), F(9), F(10),F(11)-(Po-40)]; %detect when flow is negative or pressure drop greater than 40 
        isterminal =[1,1,1,1,1,1,1,1,1,1,1]; % stop integration
        direction = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]; % only find when pressure is dropping
    end
Vcat(end);
L=Vcat(end)/A;
Fva = F(end,6)/453.59237.*MM(6); %mass flow of Va in lb
Fend=F(end,:);
end
