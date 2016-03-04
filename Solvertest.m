F0= [500 200 500 10 100 0 0 20 10 20];
Rspec= [180 446 5000 500 1.4];

[Fend,F,Fva,Vcat,L,A,vtot0] = Solver(F0,Rspec);
vtot0
A
L
Fva
Larray = Vcat/A;

figure
subplot(2,2,1)
    plot(Vcat, F(:,1:7))
    title('Process Flows')
    legend('ethylene','O2','AA','H20','CH4','VAM','CO2')
    xlabel('Volume of Cat (L)')
    ylabel('Flow Rate (lb/hr)')
subplot(2,2,2)
    plot(Vcat, F(:,8:10))
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
    bar(F(1,1:6))
    title('Components of Feed')
    xlabel('Feed component')
    ylabel('lb/hr')