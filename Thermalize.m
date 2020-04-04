function [Lattice]=Thermalize(Lattice,L,T,mcs)

for step=1:mcs %Metropolis algorithm
    
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
end
end