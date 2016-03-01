function Reactor
%Optimization
A = [0,1,1,1];
b=1;
LB = [441.4833,0.001,0.001,0.001];
UB = [449.8167,1,1,1];


r1 = @(x) 1409*exp(-3674/x(1)*((x(3)*x(2)*x(4)*(1+1.7*0))/((1+0.583*x(3)*(1+1.7*0))*(1+6.8*x(4)))));
r2 = @(x) 9625*exp(-4250/x(1))*((x(3)*(1+0.68*0))/(1+0.76*x(3)*(1+0.68*0)));

S = fmincon(@(x) -r1(x),[445,1,1,1],A,b,[],[],LB,UB);
S

%Opt end

%Solve Reactor
T = S(1); % Kelvin | Conversion for ref. (F+459.67)*(5/9)
P0 = 180;
x0 = [S(2) S(3) S(4)];
w = linspace(0,5,500);

[w,x] = ode23s(@f,[w],[0 0]);



    function dxdw = f(w,x)
    D = 1.6/12; %Tube ID
    rho = 50;   %Catalyst bulk density
    v=50;       %velocity
    epsilon = x0(1)*-.5;    
    Fao = 10;
   
    %Concentrations
    P = P0+x(2);
    pA = x0(1)*(1-x(1))/(1+epsilon*x(1))*(P/P0); %Ethylene
    pB = x0(1)*(x0(2)/x0(1)-(7/4)*x(1))/(1+epsilon*x(1))*(P/P0);%Oxygen
    pC = x0(1)*(x0(3)/x0(1)-0.5*x(1))/(1+epsilon*x(1))*(P/P0);%Acetic Acid
    pD = x0(1)*(0.5*x(1))/(1+epsilon*x(1))*(P/P0);%Vinyl Acetate
    pE = x0(1)*(1.5*x(1))/(1+epsilon*x(1))*(P/P0);%Water
    pF = x0(1)*(x(1))/(1+epsilon*x(1))*(P/P0);%Carbon Dioxide
    
    %Rates
    r1 = 1409*exp(-3674/T)*((pB*pA*pC*(1+1.7*pE))/((1+0.583*pB*(1+1.7*pE))*(1+6.8*pC)));
    r2 = 9625*exp(-4250/T)*((pB*(1+0.68*pE))/(1+0.76*pB*(1+0.68*pE)));

    %Conversion
    dxdw(1) = -(-r1-0.5*r2)/Fao;

        
    %Pressure Drop
    dxdw(2) = 0.25*v^2*((4*w)/(rho*pi*D^2));

    dxdw = dxdw';
    end
end 