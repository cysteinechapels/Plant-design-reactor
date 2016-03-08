function [ dFdV ] = myReactor( V,F,A,T,Po,vo,Fto)
%Coupled differential equations that describe PBR
%   1-ethylene, 2-oxygen, 3-acetic acid, 4-water, 5-CH4, 6-VAM, 7-CO2, 8-Eth,
%   9-Argon, 10 - N2 11 -Pressure
%   rx = r1 from problem statement to avoid confusion, ry=r2
%   all flows in gmol

%L=V/1000*(3.28^3)/A
Ft=sum(F)-F(11);
vtube=vo*(Ft/Fto)*(Po/F(11));
%vtube = vo;

% dp=-0.25*vtube^2*L;

p=zeros(1,10); %partial pressures, indexed as above
for n=1:10
    p(n)=F(n)/Ft*F(11);
end
p;
r=zeros(1,11); %reaction rates, indexed as above

% This is how satan gets ahold of men's hearts
rx= 1409*exp(-3674/T)*(p(2)*p(1)*p(3)*(1+1.7*p(4)))/((1+0.583*p(2)*(1+1.7*p(4)))*(1+6.8*p(3)))/3600;
ry=9625*exp(-4250/T)*p(2)*(1+0.68*p(4))/(1+0.76*p(2)*(1+0.68*p(4)))/3600;


r(1)=-rx-1/2*ry;
r(2)=-1/2*rx-3/2*ry;
r(3)=-rx;
r(4)=rx+ry;
r(5)=0;
r(6)=rx;
r(7)=ry;
r(8)=0;
r(9)=0;
r(10)=0;
%r(11)=0;
r(11)=-0.25*vtube^2/(A)/28.316847;

dFdV=ones(1,11);
for n=1:11
    dFdV(n)=r(n);
end

dFdV = dFdV';

end

