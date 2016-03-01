function [ Error ] = Flammability(FCH4,FO2,Fethyl,Fethan,Drysum, T,P)
%Flammability calculation from hell
%   Run while you can
FLO2=((20.5-0.02*P-0.02*(T-459.67)/(5/9)+0.04*FCH4/(Drysum+FCH4)+0.08*Fethan/(Drysum+FCH4)+0.01*Fethyl/(Drysum+FCH4))-3)*(Drysum+FCH4);
if (FLO2-FO2)>0.05*FO2 && (FLO2-FO2)<0.1*FO2
    Error=(FLO2-FO2)^2;
else
    Error=10000;
end

end

