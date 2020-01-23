%Define a function M2SystemSolver that takes as its inputs a list of
%values of b and N in bValues and NValues respectively, where b is the 
%spacing between successive points and N is the total number of points
%The output of the function is timeMat, which is a matrix containing the
%times taken for the system to be solved for each value of b and N
%Note that this function will solve the system (PMP')(PM) = (PV) using the
%ideas discussed in Question 2
function [normF,timeMat] = M2SystemSolver(bValues,NValues)

%Initialise a matrix normF which will store the norms of each vector F
%found when solving MF = V for each value of b and N
normF = zeros(length(NValues),length(bValues));

%Initialise a matrix timeMat which will store the times taken for the LU
%decomposition to take place for each value of b and N
timeMat = zeros(length(NValues),length(bValues));

%Loop through each value on NValues, calling the current value N
for i=1:length(NValues)
    N = NValues(i);
    
    %Initialise the vector V2 of the y-components of the velocities of each
    %particle in the linear system MF = V so that V2 is a column of ones
    %Note this is all we need to consider when solving the system for M2
    V2 = ones(N,1);
    
    figure() 
    
    %Loop through each value of bValues, calling the current value b
    for j=1:length(bValues)
        b = bValues(j);

        %Obtain the matrix M using the provided function Msetup with our
        %current values of b and N as its inputs
        M = Msetup(b,N);

        %Start timing how long the system takes to solve using LU
        %decomposition
        tic;
        
        %Obtain the matrix M2 by taking the elements of M which occur in
        %the rows and columns of M with even indexes
        M2 = M(2:2:end,:);
        M2 = M2(:,2:2:end);

        %Obtain the LU decomposition of M2 using the provided function
        %parpivgelim, and also the projection matrix such that 
        %(P2)(M2) = (L2)(U2) - these will all be NxN matrices
        [L2, U2, P2] = parpivgelim(M2);

        %Since V2 is a vector of ones, permuting it will not change it so 
        %just use V2 in the calculations below since (P2)(V2) = V2

        %Since V2 = (P2)(M2)(F2) = (L2)(U2)(F2), solve the lower triangular
        %system (L2)y = V2 (where y = (U2)(F2)) using forward substitution
        y = zeros(N,1);
        y(1) = V2(1)/L2(1,1);
        for k=2:N
            y(k) = ( V2(k) - L2(k,1:k-1)*y(1:k-1) )/L2(k,k);
        end

        %Find F2 by solving the upper triangular system (U2)(F2) = y using
        %backward substitution
        F2 = zeros(N,1);
        F2(N) = y(N)/U2(N,N);
        for k=N-1:-1:1
            F2(k) = ( y(k) - U2(k,k+1:N)*F2(k+1:N) )/U2(k,k);
        end

        %Form the solution F to the entire system MF = V, since we know F1
        %is a vector of N zeros since V1 = [0 ... 0]' which are the odd 
        %entries of F, and the even entries of F are the entries of F2
        F = zeros(2*N,1);
        for k=1:N
            F(2*k) = F2(k);
        end

        %Store the time taken for the LU decomposition and the solution to
        %be computed, storing it in timeMat
        timeMat(i,j) = toc;  
        
        %Compute the norm of the vector F, and store it in the matrix normF
        normF(i,j) = norm(F2);

        %Plot the y-components of F for each particle i against xi/xN
        xPlot = linspace(0,(N-1)*b,N)/((N-1)*b);
        hold on
        plot(xPlot,F2)
        title(sprintf('Graph Of F_{y,i} Against x_{i}/x_{%s} For N={%s}',num2str(N),num2str(N)))
        xlabel(sprintf('Relative position of particle, x_{i}/x_{%s}',num2str(N)))
        ylabel('y-component of force, F_{y,i}')
        
    end 
    
    hold off
end

end