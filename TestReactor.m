function TestReactor
Po=180;
F0=[1000 700 1000 10 20 0 0 10 2 4 180];
T=446;
D = 1.4/12;
L=20;
A=D^2/4*pi()*2000;
Volmax=A*L
vo=2;
Ft0=sum(F0);

Volcat = [0 Volmax];
func = @(V,F) myReactor(V,F,A,T,Po,vo,Ft0);
opts=odeset('Events',@events);

[Vcat,F] = ode45(func,Volcat,F0,opts);

%Events used to stop integration
    function[value,isterminal,direction]= events(V,F)
        value = [F(1),F(2), F(3), F(4), F(5), F(6), F(7), F(8), F(9), F(10),F(11)-(Po-40)]; %detect when flow is negative or pressure drop greater than 40 
        isterminal =[1,1,1,1,1,1,1,1,1,1,1]; % stop integration
        direction = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]; % only find when pressure is dropping
    end
F
Vcat(end)
end
