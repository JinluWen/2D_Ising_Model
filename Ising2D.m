function [Calculation]=Ising2D(T,L)

Lattice=ones(L,L);%Initial all equal 1 
Tmcs=1*10^5;%Mcs for thermalization
Mmcs=3*10^5/10*numel(Lattice);%Mcs for measurement

Lattice=Thermalize(Lattice,L,T,Tmcs);%Thermalize
[Ms,Es,Xs,Cs]=Measurement(Lattice,L,T,Mmcs);%Measurement

%Data output
Calculation=[T,L,Ms,Es,Xs,Cs];

end