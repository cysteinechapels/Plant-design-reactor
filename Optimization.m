function Optimization
%   1-ethylene, 2-oxygen, 3-acetic acid, 4-water, 5-CH4, 6-VAM, 7-CO2, 8-Eth,
%   9-Argon, 10 - N2 11 -Pressure
Pmin = 150+14.69;
Pmax=180+14.69;
Tmin=(335+459.67)*(5/9);
Tmax=(350+459.67)*(5/9);
Tubemin = 1000; %minimum number of tubes
Tubemax = 8000; %maximum number of tubes
Lengthmin = 0;
Lengthmax = 40;

C2H4min=1;
C2H4max = 700;
AAmin =1;
AAmax = 400;
H2Omin = 0;
H2Omax = 0;
CH4min = 1;
CH4max = 1500;

% 1-ethylene, 2-acetic acid, 3-water, 4-CH4, 5 - P, 6- T, 7 -
% Tube #, 8-Volume cat max, 9 - ID
LB = [C2H4min AAmin H2Omin CH4min  Pmin Tmin Tubemin Lengthmin];
UB = [C2H4max AAmax H2Omax CH4max Pmax Tmax Tubemax Lengthmax];

Recovery = 0.8;

% desired vam per hour = yearlytarget * tons/gram / days/year / hours/day /
% seconds/hour / grams/pound / fudge factor
product = 300000*1000000/350/24/3600/453.59/Recovery

% density=[0.650502,0.725865,1.73532,0.435758,0.3637,2.48941,1.01799,0.702875,0.905922,0.633822];
global MM;   
MM=[28.0532,31.9988,60.052,18.0153, 16.04, 86.0892,44.0095,30.069,39.948,28.0134];

%New optimization function that currently takes into account desired VAM
%and velocity

Fprice1 = zeros(1,10);
Fprice2 = zeros(1,10);
    function error = goal(x)
        [Fva, F, Fr, F0, Vcat, L, A,vo]=SteadyState(x);
        error = (product - Fva)^4;
    end


S = fmincon(@(x) goal(x),[1 1 0 1 Pmin Tmin 4000 20],[],[],[],[],LB,UB);
% S=fmincon(@(x) goal(x),0.5,[],[],[],[],0,1)
% % 
%  S = [3900 240 0 140 192.26 449 39000 120];

%product
% 1-ethylene, 2-acetic acid, 3-water, 4-CH4, 5 - P, 6- T, 7-Tube #, 8-Volume cat max, 9 - ID

[Fva, F, Fr, F0, Vcat, L, A,vo, n]=SteadyState(S);
S

%%%==============================================================
% Remainder takes outputs of SteadyState using inputs determined by fmincon
% and displays them

% MMM = [ MM
%         MM
%         MM
%         MM
%       ]
MMM = ones(size(F,1),10);
for n=1:size(F,1)
    MMM(n,:)=MM;
end
Recyclemol=Fr
Freshmol=F0
Recycle=Fr/453.59237.*MM*3600
Fresh = F0/453.59237.*MM*3600
% Convert gmol/s flow to lb/hr flow
Flb = F(:,1:10)/453.59237.*MMM*3600;
totallbout=Flb(end,1:10)
totalmolout=F(end,1:10)
% Find %error in VAM production
error= (Fva-product)/product*100;

% flow of VAM in lb/s
Fva; 

Flb = F(:,1:10)/453.59237.*MMM*3600;
Fton1 = Fresh(:,1:10)*24*350*0.0005;
Fton2 = Flb(end,6)*24*350*0.0005;
Fprice1 = Fton1(1,1)*1300+Fton1(1,2)*200;
Fprice2 = Fton2*1400;


% percent CH4 of inflow to reactor
percentCH4=Flb(1,5)/sum(Flb(1,:))*100;

% percent of inerts in reactor output
percentinerts=sum(Flb(end,8:10))/sum(Flb(end,:))*100;


VolumeCatalyst=Vcat(end); %total volume of catalyst
L; %end length of reactor
O2fin = Flb(end,2); %final O2 amount
FeedtoRecycle=sum(F0)/sum(Fr)*100;
vo;
dP=F(1,11)-F(end,11);
dPcalc = 0.25*vo^2*L;
Ntubes=S(7);

Flbplot= [Flb(1,1:4) Flb(1,8:10)];

Larray=Vcat./28.31/A;



conversion= (sum(F(1,1:10))-sum(F(end,1:10)))/sum(F(1,1:10))*100;
yield=F(end,6)/F(1,1)*100;

f = figure('Position',[440 500 800 120]);

% create the data
d = [Fva error dP dPcalc percentCH4 percentinerts FeedtoRecycle VolumeCatalyst L vo Ntubes conversion yield];

% Create the column and row names in cell arrays 
cnames = {'VAM outflow','Error','dP','dPcacl','%CH4','%inerts','%Feed/Recycle', 'VolCat', 'Length','velocity','#tubes', 'conversion', 'yield'};

% Create the uitable
t = uitable(f,'Data',d,...
            'ColumnName',cnames);

% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);

figure
subplot(2,2,1)
    plot(Vcat, Flb(:,1:7))
    title('Process Flows')
    legend('ethylene','O2','AA','H20','CH4','VAM','CO2')
    xlabel('Volume of Cat (L)')
    ylabel('Flow Rate (lb/hr)')
subplot(2,2,2)
    bar(Recycle)
    title('Components of Recycle')
    xlabel('Recycle component')
    ylabel('lb/hr')
subplot(2,2,3)
    plot(Larray,F(:,11))
    title('Pressure profile')
    xlabel('Length of Reactor (ft)')
    ylabel('Pressure psia')
subplot(2,2,4)
    bar(Fresh)
    title('Components of Feed')
    xlabel('Feed component')
    ylabel('lb/hr')
Costcheck = 0;
Costcheck = Fprice2-Fprice1
end