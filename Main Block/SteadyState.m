function [Fva,F,Fr,F0,Vcat,L,A,vtot0, n] = SteadyState(I)
%This function determines the steady state flows through the reactor using
%a for loop. The purification section is approximated and the purge is a
%specified value. 

% I = [ ethylene, acetic, water, ch4, P, T, Tube, Volume, ID]

% Species order
% 1-ethylene, 2-oxygen, 3-acetic acid, 4-water, 5-CH4, 6-VAM, 7-CO2, 8-Eth,
%   9-Argon, 10 - N2 11 -Pressure

%Input order
% 1-ethylene 2-acetic acid 3-water 4-ch4 5-P, 6-T, 7-Tube #, 8-Volume cat max, 9-ID


%safety factor for flammability calculations
safety= 0.8;
%Approximated water recovery from the purification section
waterfactor = 0.14;
%Selected purge percent
purge=0.01;
%Approximated acetic acid recovery from the purification section
recoveryAA=0.99;

Rspec = I(5:8); % Rspec = [ P T Tube Volume ID ]
Po= Rspec(1);%psia
T= Rspec(2); %Kelvin

MM=[28.0532,31.9988,60.052,18.0153, 16.04, 86.0892,44.0095,30.069,39.948,28.0134];
%Densities at respective feed T and P in lb/ft^3 (T=80F and P-400psig for
%ethylene, ethane, oxygen, argon, N2; methane at 80F and 135psig; AAc at
%80F 100psig; 0 signifies not fed to process
feeddensity=[2.21607,2.34002,2.4017,0,0.422338,0,0,2.18847,2.29067,2.02514];
%separate density for the ethane that comes in with methane
ethinmethdensity=0.842457;

%Calculation of inerts using volume
FvolE = I(1).*MM(1)/453.59237./feeddensity(1);
FvolCH4=I(4).*MM(5)/453.59237./feeddensity(5);

%WTF
Feth = (FvolE*245/1000000*feeddensity(8)+FvolCH4*800/1000000*ethinmethdensity)*453.59237/MM(8);
FO2 = 0;
% FO2 = I(4);
FAr = 0;
FN2 = 0;


% % Iterative calculation used determine safe (FLO2) and operating (FO2)
% oxygen concentration based on input values for ethylene and methane
for i = 1:100
    Drysum = I(1)+I(4)+Feth+FO2+FAr+FN2;

    PandT = 20.5-0.02*Po-0.02*(9/5*(T - 273) + 32);
    Methane = 0.04*I(4)/Drysum*100;
    Ethane = 0.08*Feth/Drysum*100;
    Ethylene= 0.01*I(1)/Drysum*100;
    FLO2=((PandT+Methane+Ethane+Ethylene)-3)*(Drysum)/100;

    FO2=FLO2*safety;
    FO2vol=FO2*MM(2)/453.59237/feeddensity(2);
    Oxyinerts = [1 -1 -1; -0.0011 1 0; -0.0005 0 1];
    Totalinerts = [FO2vol;0;0];
    ArandN2 = Oxyinerts\Totalinerts;
    FAr = ArandN2(2)*feeddensity(9)*453.59237/MM(9);
    FN2 = ArandN2(3)*feeddensity(10)*453.59237/MM(10);
    Drynew = I(1)+I(4)+Feth+FO2+FAr+FN2;
    Check = abs(Drysum-Drynew);

    if Check <= 0.001
        break
    end
end
i;
F0=[I(1) FO2 I(2) I(3) I(4) 0 0 Feth FAr FN2];
F0star=[I(1) FO2 I(2) I(3) I(4) 0 0 Feth FAr FN2];
Fr=zeros(1,10);

%Outflows from a single pass through the reactor are modified to simulate
%the purification section, then sent back to reactor along with fresh feed.
%This is repeated until there is no change between the previous iteration
%outfow and the current outflow, signifying steady state
for n=1:500
    Foriginal = F0star;
    F=F0+Fr;
    New=Solver(F,Rspec);
    Fr=New(end,1:end-1);
    SumF = sum(abs(Fr));
    Fr(1)=Fr(1)*(1-purge);%ethylene purge
    Fr(2)=Fr(2)*(1-purge);%Oxy purge
    Fr(5)=Fr(5)*(1-purge);%Methane purge
    Fr(8)=Fr(8)*(1-purge);%ethane purge
    Fr(9)=Fr(9)*(1-purge);%argon purge
    Fr(10)=Fr(10)*(1-purge);%nitrogen purge
    Fr(4)=Fr(4)*waterfactor; %water recovery
    Fr(6)=SumF*0.0001/0.9999*(1-purge);%100ppm VAM in recycle
    Fr(3)=Fr(3)*recoveryAA; %recycled AA from azeo column
    Fr(7)=(Fr(7)/2+Fr(7)*0.01/2)*(1-purge); %CO2 treatment
    Fnew=Fr+F0;
%     TF = abs(F);
%     TFnew = abs(Fnew);
%     Test = sum(abs(TF-TFnew));
%      if Test<1
%          break
%      end
    Tinert = F(8)+F(9)+F(10);
    Tinertnew = Fnew(8)+Fnew(9)+Fnew(10);
    Test = (sum(abs(Tinert-Tinertnew)))^2;
        if ((n>1) && (Test<0.01))
            break
        end
    F0 = Foriginal;
    F0 = F0-Fr;
    F0(8)=Feth;  
    F0(9)=FAr;  
    F0(10)=FN2;  
    F0(4)=0;
    F0(6)=0;
    F0(7)=0;
    F0;
    Fr(8);
end
n;
[Fend,F,Fva, Vcat,L,A,vtot0]=Solver(Fnew,Rspec);

end

