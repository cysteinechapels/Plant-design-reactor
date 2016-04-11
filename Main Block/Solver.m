function [Fend,F,Fva,Vcat,L,A,vtot0] = Solver (F0,Rspec)
% This program manipulates all the necessary input variables into the
% format/unit the reactor differential equation solver needs. Then it calls
% the reactor code in a ODE solver to determine flows through the reactor
%   1-ethylene, 2-oxygen, 3-acetic acid, 4-water, 5-CH4, 6-VAM, 7-CO2, 8-Eth,
%   9-Argon, 10 - N2 11 -Pressure

T= Rspec(2); %Kelvin
D= 1.4/12; %feet
A=D^2/4*pi()*Rspec(3); %total cross-sectional area of reactor
Po= Rspec(1);%psia
Volmax = A*Rspec(4)*28.316847;
Volcat=[0 Volmax]; %volume of catalyst, liters

%Density at 180psig nd 340.1F]
density=[0.650502,0.725865,1.73532,0.435758,0.3637,2.48941,1.01799,0.702875,0.905922,0.633822];
global MM   
MM=[28.0532,31.9988,60.052,18.0153, 16.04, 86.0892,44.0095,30.069,39.948,28.0134];

Ft0=sum(F0);

%Convert gmol to lbmol flow for velocity calculations
F0lbmol=F0./453.59237;
Massflow=F0lbmol.*MM;
Massfrac = Massflow/sum(Massflow);
appdensity= sum(Massfrac.*density);
vtot0=sum(Massflow)/appdensity./A; %average velocity

F0=[F0 Po];

%ODE solver used to determine flows throughout reactor
func = @(V,F) myReactor(V,F,A,T,Po,vtot0,Ft0);
opts=odeset('Events',@events);

[Vcat,F] = ode45(func,Volcat,F0,opts);

%Events used to stop integration
%Integration is stopped if dP exceeds 40 or flows become negative
    function[value,isterminal,direction]= events(V,F)
        value = [F(1),F(2), F(3), F(4), F(5), F(6), F(7), F(8), F(9), F(10),F(11)-(Po-40)]; %detect when flow is negative or pressure drop greater than 40 
        isterminal =[1,1,1,1,1,1,1,1,1,1,1]; % stop integration
        direction = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]; % only find when pressure is dropping
    end
Vcat(end);
L=Vcat(end)/28.316847/A;
%Fva=F(end,6);
Fva = F(end,6)/453.59237.*MM(6); %mass flow of Va in lb
Fend=F(end,:);
end
