%Define a function systemSolver that takes as its inputs a list of values
%of b and N in bValues and NValues respectively, where b is the spacing
%between successive points and N is the total number of points
%The outputs normF, timeMat residualLU, normLU and normM are all matrices
%normF contains the norms of each solution F to MF = V
%timeMat contains the times taken for the LU decomposition calculations and
%to solve the system
%residualLU contians the values of the backwards errors 
%normLU contains the values of ||L||.||U|| for each LU factorisation
%normM contains the values of ||M|| for each M
function [normF,timeMat,residualLU,normLU,normM] = systemSolver(bValues,NValues)

%Initialise a matrix normF which will store the norms of each vector F
%found when solving MF = V for each value of b and N
normF = zeros(length(NValues),length(bValues));

%Initialise a matrix timeMat which will store the times taken for the LU
%decomposition to take place for each value of b and N
timeMat = zeros(length(NValues),length(bValues));

%Initialise a matrix residualLU which will store backwards error of each LU 
%decomposition, which is ||P'LU - M||/(||L||.||U||)
residualLU = zeros(length(NValues),length(bValues));

%Initialise a matrix normLU which will store each value of ||L||.||U||
normLU = zeros(length(NValues),length(bValues));

%Initialise a matrix normM which will store each value of ||M||
normM = zeros(length(NValues),length(bValues));

%Loop through each value on NValues, calling the current value N
for i=1:length(NValues)
    N = NValues(i);
    
    %Initialise the vector V for our linear system MF = V so that V is a
    %vector of length 2N consisting of [0,1]' repeated N times
    V = zeros(2*N,1);
    for k = 1:N
        V(2*k) = 1;
    end
    
    figure() 
    
    %Loop through each value of bValues, calling the current value b
    for j=1:length(bValues)
        b = bValues(j);

        %Obtain the matrix M using the provided function Msetup with our
        %current values of b and N as its inputs
        M = Msetup(b,N);
        
        %Start timing how long the LU decomposition calculations take
        tic;
        
        %Obtain the LU decomposition of M using the provided function
        %parpivgelim, and also the projection matrix P such that PM = LU
        [L, U, P] = parpivgelim(M);

        %In order to solve MF = V, multiply both sides from the left by P
        %to obtain PMF = PV, and store the right-hand side as PV
        PV = P*V;
        
        %Since PV = PMF = LUF, solve the lower triangular system Ly = PV
        %(where y = UF) using forward substitution
        y = zeros(2*N,1);
        y(1) = PV(1)/L(1,1);
        for k=2:2*N
            y(k) = ( PV(k) - L(k,1:k-1)*y(1:k-1) )/L(k,k);
        end

        %Find F by solving the upper triangular system UF = y using
        %backward substitution
        F = zeros(2*N,1);
        F(2*N) = y(2*N)/U(2*N,2*N);
        for k=2*N-1:-1:1
            F(k) = ( y(k) - U(k,k+1:2*N)*F(k+1:2*N) )/U(k,k);
        end
        
        %Store the time taken for the LU decomposition and the solution to
        %be computed, storing it in timeMat
        timeMat(i,j) = toc;
        
        %Compute the norm of the vector F, and store it in the matrix normF
        normF(i,j) = norm(F);    
        
        %Compute the backwards error ||P'LU - M||/(||L||.||U||), and store 
        %it in residualLU
        residualLU(i,j) = norm(P'*L*U - M)/(norm(L)*norm(U));
        normLU(i,j) = norm(L)*norm(U);
        normM(i,j) = norm(M);
        
        %Plot the y-components of F for each particle i against xi/xN
        xPlot = linspace(0,(N-1)*b,N)/((N-1)*b);
        fPlot = F(2 : 2 : end);
        hold on
        plot(xPlot,fPlot)
        title(sprintf('Graph Of F_{y,i} Against x_{i}/x_{%s} For N={%s}',num2str(N),num2str(N)))
        xlabel(sprintf('Relative position of particle, x_{i}/x_{%s}',num2str(N)))
        ylabel('y-component of force, F_{y,i}')
        
    end 
    
    legend('b=2','b=4','b=10')
    hold off
end

end