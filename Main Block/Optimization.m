function Optimization
% This function solves the SteadyState function using the range of input
% values specified and tries to minimize "error" the difference between how
% much VAM is produced and how much Vam is required. It also prints and
% presents this data

%   1-ethylene, 2-oxygen, 3-acetic acid, 4-water, 5-CH4, 6-VAM, 7-CO2, 8-Eth,
%   9-Argon, 10 - N2 11 -Pressure
Pmin = 150+14.69; %minimum pressure
Pmax=180+14.69; %maximum pressure
Tmin=(335+459.67)*(5/9); %minimum temperature
Tmax=(350+459.67)*(5/9); %maximum temperature
Tubemin = 1000; %minimum number of tubes
Tubemax = 6000; %maximum number of tubes
Lengthmin = 0; %minimum length
Lengthmax = 20; %maximum length
Purgemin = 0;
Purgemax = 0.05;

%start up fresh feed ranges, in mol/s
C2H4min=1;
C2H4max = 1500;
AAmin =1;
AAmax =1000;
H2Omin = 0;
H2Omax = 0;
CH4min = 1;
CH4max = 800;

% 1-ethylene, 2-acetic acid, 3-water, 4-CH4, 5 - P, 6- T, 7 -
% Tube #, 8-Volume cat max, 9 - ID
LB = [C2H4min AAmin H2Omin CH4min  Pmin Tmin Tubemin Lengthmin Purgemin];
UB = [C2H4max AAmax H2Omax CH4max Pmax Tmax Tubemax Lengthmax Purgemax];

Recovery = 0.95; %estimated recovery of vinyl acetate
O2conversion = 90;
AAconversion = 90;


% desired vam in lb per second = yearlytarget * tons/gram / days/year / hours/day /
% seconds/hour / grams/pound / fudge factor
product = 300000*1000000/350/24/3600/453.59/Recovery

% density=[0.650502,0.725865,1.73532,0.435758,0.3637,2.48941,1.01799,0.702875,0.905922,0.633822];
MM=[28.0532,31.9988,60.052,18.0153, 16.04, 86.0892,44.0095,30.069,39.948,28.0134];

%New optimization function that currently takes into account desired VAM
Fprice1 = zeros(1,10);
Fprice2 = zeros(1,10);
Costcheck = 0;
spec=0;
cost = 0;
    function error = goal(x)
        [Fva, F, Fr, F0, Vcat, L, A,vo]=SteadyState(x);
        % conversions for cost check
        MMM = ones(size(F,1),10);
            for n=1:size(F,1)
            MMM(n,:)=MM;
            end
        Flb = F(:,1:10)/453.59237.*MMM*3600;
        Fresh = F0/453.59237.*MM*3600;
        Fton1 = Fresh(:,1:10)*24*350*0.0005;
        Fton2 = Flb(end,6)*24*350*0.0005;
        Fprice1 = Fton1(1,1)*1300+Fton1(1,2)*200;
        Fprice2 = Fton2*1400*Recovery;
        Fprice3 = Fton1(1,3)*850;
        Costcheck = Fprice2-Fprice1-Fprice3;
        
        %Limit O2 conversion 
        convO2= (F(1,2)-F(end,2))/F(1,2)*100;
        O2check = (convO2-O2conversion)^2;
        
        %Limit AA conversion
        convAA= (F(1,3)-F(end,3))/F(1,3)*100;
        AAcheck = (convAA-AAconversion)^2;
        
       
        %set production to goal value
        spec = (product - Fva)^2;
        
        %maximize cost
        cost = ((4.2E8-Costcheck)/1E9)^2;
        error = spec+cost+O2check+AAcheck;
    end


S = fmincon(@(x) goal(x),[1 1 0 1 Pmin Tmin 4000 20 0.005],[],[],[],[],LB,UB);

S
S(9)
%product
% Index key for flows: 1-ethylene, 2-acetic acid, 3-water, 4-CH4, 5 - P, 6- T, 7-Tube #, 8-Volume cat max, 9 - ID
[Fva, F, Fr, F0, Vcat, L, A,vo, n]=SteadyState(S);

spec
cost
%%%==============================================================
% Remainder takes outputs of SteadyState using inputs determined by fmincon
% and displays them

%generates a matrix of molecular masses; required for matrix math
MMM = ones(size(F,1),10);
for n=1:size(F,1)
    MMM(n,:)=MM;
end

%calculate the recycle flows in lb/hr
Recycle=Fr/453.59237.*MM*3600;

%Calculate fresh flows in lb/hr
Fresh = F0/453.59237.*MM*3600;

%Total allowable amount of VAM in vaporizer
Fppm =(sum(Recycle)+sum(Fresh))*0.0001/0.9999;

% Convert gmol/s flow to lb/hr flow
Flb = F(:,1:10)/453.59237.*MMM*3600;
Outputlbs=Flb(end,1:10)
F(end,1:10);
% Find %error in VAM production
error= (Fva-product)/product*100;

% flow of VAM in lb/s
Fva; 
MMM = ones(size(F,1),10);
for n=1:size(F,1)
       MMM(n,:)=MM;
end

%Raw material costs
Flb = F(:,1:10)/453.59237.*MMM*3600;
Fton1 = Fresh(:,1:10)*24*350*0.0005;
Fton2 = Flb(end,6)*24*350*0.0005;
Fprice1 = Fton1(1,1)*1300+Fton1(1,2)*200;
Fprice2 = Fton2*1400*Recovery;
Fprice3 = Fton1(1,3)*850;


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


%%%==============================================================
% Remainder takes outputs and presents it as a table and a series of graphs

convC2H4= (F(1,1)-F(end,1))/F(1,1)*100;
convO2= (F(1,2)-F(end,2))/F(1,2)*100;
convAA= (F(1,3)-F(end,3))/F(1,3)*100;

f = figure('Position',[440 500 800 120]);

% create the data
d = [Fva error dP dPcalc percentCH4 percentinerts FeedtoRecycle VolumeCatalyst L vo Ntubes convC2H4 convO2 convAA];

% Create the column and row names in cell arrays 
cnames = {'VAM outflow','Error','dP','dPcacl','%CH4','%inerts','%Feed/Recycle', 'VolCat', 'Length','velocity','#tubes', 'convC2H4', 'convO2','convAA'};

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
    bar(Flbplot)
    title('Components of Feed')
    xlabel('Feed component')
    ylabel('lb/hr')
Costcheck = 0;
Costcheck = Fprice2-Fprice1-Fprice3

Vamproduced= Fva*Recovery*453.59*3600*24*350/1000000
% 
% Titles={'Stream', 'ethylene', 'oxygen', 'acetic acid', 'water', 'CH4', 'VAM','CO2', 'Ethane','Argon','N2'};
% sheet= 1;
% filename= 'Reactor.xlsx';
% TitleRange='B2';
% 
% xlswrite(filename,Titles,sheet,TitleRange)
% 

end