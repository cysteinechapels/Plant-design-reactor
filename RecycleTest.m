P = 180+14.69; %pressure psia
T = 446; %Kelvin
NTube= 5000;
Volmax = 100;
ID = 1.3;



Fethyl=0.08; %flow of ethylene in gmol/s
FAA = 0.05; %flow acetic acid in gmol/s
Fh20= 0.001; %flow of water in gmol/s
FCH4 = 0.04; %flow of methane in gmol/s CONTROLS O2 flow

S=[Fethyl FAA Fh20 FCH4 P T];
[Fva, F, Fr, F0, Vcat, L, A,vo]=SteadyState1Tube(S);

Recovery = 0.8; %estimation of downstream recovery

product = 300000*1000000/350/24/3600/453.59/Recovery; %amount of VAM that must be produced
%

% error= (Fva-product)/product*100

Fva; % flow in lb/s
percentCH4=F(1,5)/sum(F(1,:))*100;
percentinerts=sum(F(1,8:10))/sum(F(1,:))*100;
VolumeCatalyst=Vcat(end); %total volume of catalyst
L; %end length of reactor
O2fin = F(end,2); %final O2 amount
FeedtoRecycle=sum(F0)/sum(Fr)*100;
vo;
dP=F(1,11)-F(end,11);

global MM
MM=[28.0532,31.9988,60.052,18.0153, 16.04, 86.0892,44.0095,30.069,39.948,28.0134];


MMM = MM;
for n=1:size(F,1)-1
    MMM=vertcat(MMM,MM);
end


Ntube=product/Fva
Larray=Vcat./A;

Flb = F(:,1:10)/453.59237.*MMM*3600; 
conversion= (sum(F(1,1:10))-sum(F(end,1:10)))/sum(F(1,1:10))*100;
yield=F(end,6)/F(1,1)*100;

f = figure('Position',[440 500 800 120]);

% create the data
d = [Fva dP percentinerts percentCH4 FeedtoRecycle VolumeCatalyst L vo Ntube conversion yield];

% Create the column and row names in cell arrays 
cnames = {'VAM outflow','dP','%CH4','%inerts','%Feed/Recycle', 'VolCat', 'Length','velocity','Number of Tubes', 'conversion', 'yield'};

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