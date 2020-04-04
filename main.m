%2D Ising Model using Metropolis Algorithm

clear;

tic

Data=zeros(900,6);

for i=1:3
    L=10*i;
    for j=1:300
        T=0.015*j;
        Data(j+300*(i-1),:)=Ising2D(T,L);
    end
end

%save data
save Data.txt -ascii Data

toc
