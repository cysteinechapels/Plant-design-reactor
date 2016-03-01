function [Fva] = Solvertest (I)
%   1-ethylene, 2-oxygen, 3-acetic acid, 4-water, 5-VAM, 6-CO2, 7-Eth,
%   8-Argon, 9 - N2 10 -Pressure
Volcat=[0 I(8)]; %volume of catalyst, liters
T= I(6); %Kelvin
D= 1.2/12; %feet
A=D^2/4*pi()*I(7); %total cross-sectional area of reactor
Po= I(5); %psia

%Density at 180psig and 340.1F]
density=[0.650502,0.725865,1.73532,0.435758,2.48941,1.01799,0.702875,0.905922,0.633822];
   
MM=[28.0532,31.9988,60.052,18.0153, 86.0892,44.0095,30.069,39.948,28.0134];

%NEED to calc ethane and argon and N2 from feed
Fvol = I(1:4).*MM(1:4)/453.59237/density(1:4);
Feth = Fvol(1)*245/1000000;
Farg = F(;
Fn2;

%All flows in reaction must be gmol
F0=[I(1:4) 0 0 0 0 0];

Ft0=sum(F0);

%Convert gmol to lbmol flow for velocity calculations
F0lbmol=F0./453.59237;
Massflow=F0lbmol.*MM;
v0=Massflow/density./A; %velocity of each component
vtot0=sum(v0.*Massflow)/sum(Massflow); %total velocity?? weighted average?????


F0 = [F0 Po];
func = @(V,F) myReactor(V,F,A,T,Po,vtot0,Ft0);
opts=odeset('Events',@events);

[Vcat,F] = ode45(func,Volcat,F0,opts);

MMM = MM;
for n=1:size(F,1)-1
    MMM=vertcat(MMM,MM);
end

Flb = F(:,1:9)/453.59237.*MMM;

 plot(Vcat,Flb(:,1:6))

    function[value,isterminal,direction]=events(V,F)
        value = [F(1),F(2), F(3), F(4), F(5), F(6), F(7), F(8), F(9), F(10)-(Po-40)]; %detect when flow is negative or pressure drop greater than 40
        isterminal =[1,1,1,1,1,1,1,1,1,1]; % stop integration
        direction = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]; % only find when pressure is dropping
    end
L= Vcat(end)/A;

Fva = F(end,5)/453.59237.*MM(5); %mass flow of Va in lb
end
