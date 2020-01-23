%Define a function levinsonSystemSolver that takes as its inputs a list of
%values of b and N in bValues and NValues respectively, where b is the 
%spacing between successive points and N is the total number of points
%The output of the function is timeMat, which is a matrix containing the
%times taken for the system to be solved for each value of b and N when
%using the Levinson algorithm from Chapter 4 of Golub and Van Loan
function [normF,timeMat] = levinsonSystemSolver(bValues,NValues)

%Initialise a matrix normF which will store the norms of each vector F
%found when solving MF = V for each value of b and N
normF = zeros(length(NValues),length(bValues));

%Initialise a matrix timeMat which will store the times taken for the
%system to be solved using Levinson's algorithm for each value of b and N
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
        
        %Start timing how long the system takes to solve using the Levinson
        %algorithm
        tic;

        %Obtain the matrix M2 by taking the elements of M which occur in
        %the rows and columns of M with even indexes
        M2 = M(2:2:end,:);
        M2 = M2(:,2:2:end);

        %Since the version of Levinson's algorithm provided requires a
        %leading diagonal of ones in the matrix, normalise both M2 and V2
        %by dividing all elements by this value
        V2 = V2/M2(1,1);
        M2 = M2/M2(1,1);
        
        %Compute the solution using the Levinson algorithm as provided in
        %Chapter 4 of Matrix Computations by Golub and Van Loan
        %Note that this is not my own algorithm
        r = M2(1,2:end)';
        y = zeros(N-1,1);
        y(1) = -r(1);
        F2 = zeros(N,1);
        F2(1) = V2(1);
        beta = 1;
        alpha = -r(1);

        for k=1:N-1
            beta = (1 - alpha^2)*beta;
            mu =(V2(k+1) - r(1:k)'*F2(k:-1:1))/beta;
            F2(1:k) = F2(1:k) + mu*y(k:-1:1);
            F2(k+1) = mu;

            if k < N-1
                alpha = (-r(k+1) - r(1:k)'*y(k:-1:1))/beta;
                y(1:k) = y(1:k) + alpha*y(k:-1:1);
                y(k+1) = alpha;
            end
        end

        %Form the solution F to the entire system MF = V, since we know F1
        %is a vector of N zeros since V1 = [0 ... 0]' which are the odd 
        %entries of F, and the even entries of F are the entries of F2
        F = zeros(2*N,1);
        for k=1:N
            F(2*k) = F2(k);
        end
        
        %Store the time taken for the Levinson algorithm to compute the 
        %solution to MF = V, storing it in timeMat
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