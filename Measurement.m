function [Ms,Es,Xs,Cs]=Measurement(Lattice,L,T,mcs)
    
    Mmean=zeros(1,mcs);
    Emean=zeros(1,mcs);
    
for i=1:mcs %Metropolis algorithm

    %generate two random numbers m & n here to swap (m,n)'th random spin
    ri=randi(L,1,2); %2 random integers between 1 & L
    m=ri(1); n=ri(2); %save column (m) and row (n) indices
    
    %find its nearest neighbors (periodic boundary conditions)
    above = mod(n - 1 - 1, size(Lattice,1)) + 1;
    below = mod(n + 1 - 1, size(Lattice,1)) + 1;
    left  = mod(m - 1 - 1, size(Lattice,2)) + 1;
    right = mod(m + 1 - 1, size(Lattice,2)) + 1; 
    neighbors=[Lattice(right,n);Lattice(left,n);Lattice(m,above);Lattice(m,below)]; 
    
    %calculate the energy component if that spin is flipped
    dE=2*Lattice(m,n)*sum(neighbors);
    if dE<=0 %flip the spin if that reduces the energy
        Lattice(m,n)=-Lattice(m,n);
    else %if that increases the energy, then flip the spin with probability,supressed by the corresponding boltzmann factor p=exp(-dE/T)
        p=exp(-dE/T); %update Boltzmann factor (probability)
        r=rand; %uniform random variable
        if r<p
            Lattice(m,n)=-Lattice(m,n);
        else
            Lattice(m,n)=Lattice(m,n);
        end
    end
    
    %Calculate magnetisation per spin
    Mmean(1,i)=mean(Lattice(:))/numel(Lattice);
    
    %Calculate energy per spin
    Sum_Neighbors=...
          circshift(Lattice,[ 0  1]) ...
        + circshift(Lattice,[ 0 -1]) ...
        + circshift(Lattice,[ 1  0]) ...
        + circshift(Lattice,[-1  0]);
    Em=-Lattice.*Sum_Neighbors;
    E=0.5*sum(Em(:));
    Emean(1,i)=E/numel(Lattice);
    
end

    Ms=abs(mean(Mmean));%all above zero
    Es=mean(Emean);
    Xs=(mean(Mmean.^2)-mean(Mmean)^2);
    Cs=(mean(Emean.^2)-mean(Emean)^2)/T^2;
    
end