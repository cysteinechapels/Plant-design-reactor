%   1-ethylene, 2-oxygen, 3-acetic acid, 4-water, 5-CH4, 6-VAM, 7-CO2, 8-Eth,
%   9-Argon, 10 - N2 11 -Pressure
Pmin = 150+14.69;
Pmax=180+14.69;
Tmin=(335+459.67)*(5/9);
Tmax=(350+459.67)*(5/9);
Tubemin = 1000; %minimum number of tubes
Tubemax = 10000; %maximum number of tubes
Volmin = 1000; %minimum volume
Volmax = 10000; %maximum volume
IDmin = 1;
IDmax = 1.5;
% 1-ethylene, 2-acetic acid, 3-water, 4-CH4, 5 - P, 6- T, 7 -
% Tube #, 8-Volume cat max, 9 - ID
LB = [100 100 100 100  Pmin Tmin Tubemin Volmin IDmin];
UB = [10000 1000 1000 50000 Pmax Tmax Tubemax Volmax IDmax];
Recovery = 0.8;

product = 300000*1000000/350/24/3600/453.59/Recovery;

density=[0.650502,0.725865,1.73532,0.435758,0.3637,2.48941,1.01799,0.702875,0.905922,0.633822];
   
MM=[28.0532,31.9988,60.052,18.0153, 16.04, 86.0892,44.0095,30.069,39.948,28.0134];


goal = @(x) (product - subsref(SteadyState(x),struct('type','()','subs',{{1}})))^4;
S = fmincon(@(x) goal(x),[1000 300 200 50000 Pmin Tmin 8001 1000 1.2],[],[],[],[],LB,UB);
S

product
% 1-ethylene, 2-acetic acid, 3-water, 4-CH4, 5 - P, 6- T, 7-Tube #, 8-Volume cat max, 9 - ID
%Guess = [1000 605.3 200 1812.9 190.69 445.4 10000 1000000 1.5];
[Fva, F, Fr, F0, Vcat, L, A,vo]=SteadyState(S);

error= (Fva-product)/product*100;
Fva; % flow in lb/s
percentCH4=F(1,5)/sum(F(1,:))*100;
percentinerts=sum(F(1,8:10))/sum(F(1,:))*100;
VolumeCatalyst=Vcat(end); %total volume of catalyst
L; %end length of reactor
O2fin = F(end,2); %final O2 amount
FeedtoRecycle=sum(F0)/sum(Fr)*100;
vo;
dP=F(1,11)-F(end,11);


MMM = MM;
for n=1:size(F,1)-1
    MMM=vertcat(MMM,MM);
end

Larray=Vcat./A;

Flb = F(:,1:10)/453.59237.*MMM*3600; 
conversion= (sum(F(1,1:10))-sum(F(end,1:10)))/sum(F(1,1:10))*100;
yield=F(end,6)/F(1,1)*100;

f = figure('Position',[440 500 800 120]);

% create the data
d = [Fva error dP percentinerts percentCH4 FeedtoRecycle VolumeCatalyst L vo conversion yield];

% Create the column and row names in cell arrays 
cnames = {'VAM outflow','Error','dP','%CH4','%inerts','%Feed/Recycle', 'VolCat', 'Length','velocity', 'conversion', 'yield'};

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