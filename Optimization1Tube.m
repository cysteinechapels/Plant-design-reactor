function Optimization1Tube
%   1-ethylene, 2-oxygen, 3-acetic acid, 4-water, 5-CH4, 6-VAM, 7-CO2, 8-Eth,
%   9-Argon, 10 - N2 11 -Pressure
Pmin = 150+14.69;
Pmax=180+14.69;
Tmin=(335+459.67)*(5/9);
Tmax=(350+459.67)*(5/9);
Ethylenemin = 0;
Ethylenemax = 1;
AAmin=0;
AAmax = 1;
H20min=0;
H20max = 1;
CH4min=0;
CH4max=1;
% 1-ethylene, 2-acetic acid, 3-water, 4-CH4, 5 - P, 6- T
LB = [Ethylenemin AAmin H20min CH4min  Pmin Tmin];
UB = [Ethylenemax AAmax H20max CH4max Pmax Tmax];

Recovery = 0.8;

% desired vam per hour = yearlytarget * tons/gram / days/year / hours/day /
% seconds/hour / grams/pound / fudge factor
product = 300000*1000000/350/24/3600/453.59/Recovery;

% density=[0.650502,0.725865,1.73532,0.435758,0.3637,2.48941,1.01799,0.702875,0.905922,0.633822];
global MM;   
MM=[28.0532,31.9988,60.052,18.0153, 16.04, 86.0892,44.0095,30.069,39.948,28.0134];


%New optimization function that currently takes into account desired VAM
%and velocity
    function error = goal(x)
        [Fva, F, Fr, F0, Vcat, L, A,vo]=SteadyState1Tube(x);
        error = -Fva*100+(20-L)^2;
    end


S = fmincon(@(x) goal(x),[0.2 0.13 0 0.08 Pmax Tmax],[],[],[],[],LB,UB);
S
 
%product
% 1-ethylene, 2-acetic acid, 3-water, 4-CH4, 5 - P, 6- T, 7-Tube #, 8-Volume cat max, 9 - ID

[Fva, F, Fr, F0, Vcat, L, A,vo]=SteadyState1Tube(S);

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

% Convert gmol/s flow to lb/s flow
Flb = F(:,1:10)/453.59237.*MMM*3600;

% flow of VAM in lb/s
Fva 


% percent CH4 of inflow to reactor
percentCH4=Flb(1,5)/sum(Flb(1,:))*100;

% percent of inerts in reactor output
percentinerts=sum(Flb(end,8:10))/sum(Flb(end,:))*100;

VolumeCatalyst=Vcat(end); %total volume of catalyst
L; %end length of reactor
O2fin = Flb(end,2); %final O2 amount
FeedtoRecycle=sum(F0)/sum(Fr)*100;
vo;
dP=F(1,11)-F(end,11)
dPcalc = 0.25*vo^2*L;
Ntubes=product/Fva;



Larray=Vcat./A;


conversion= (sum(F(1,1:10))-sum(F(end,1:10)))/sum(F(1,1:10))*100;
yield=F(end,6)/F(1,1)*100;

f = figure('Position',[440 500 800 120]);

% create the data
d = [Fva dP dPcalc percentCH4 percentinerts FeedtoRecycle VolumeCatalyst vo Ntubes conversion yield];

% Create the column and row names in cell arrays 
cnames = {'VAM outflow','dP','dPcacl','%CH4','%inerts','%Feed/Recycle', 'VolCat', 'velocity','#tubes', 'conversion', 'yield'};

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
    plot(Vcat, Flb(:,8:10))
    title('Inert Flows')
    legend('ethane','Ar','N2')
    xlabel('Volume of Cat (L)')
    ylabel('Flow Rate (lb/hr)')
subplot(2,2,3)
    plot(Larray,F(:,11))
    title('Pressure profile')
    xlabel('Length of Reactor (ft)')
    ylabel('Pressure psia')
subplot(2,2,4)
    bar(Flb(1,1:6))
    title('Components of Feed')
    xlabel('Feed component')
    ylabel('lb/hr')
end