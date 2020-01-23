function [M] = Msetup(b,N)
% This function sets up the matrix M for N equally spaced particles in a
% line alon the x-direction.  The spacing between the particles is b.
%
% Inputs: N: number of particles.
%         b: spacing between the particles.
%
% Output: M: 2Nx2N mobility matrix M.
%

x = (0:1:N-1)'*b;
M = zeros(2*N);

for n = 1:N
    Xn = [x(n) 0]';
    indn = 2*(n - 1) + 1; 
    for m = n:N
        indm = 2*(m - 1) + 1;
        Xm = [x(m) 0]';
        Mnm = RPY2D(Xn, Xm, 1);
        M(indn:indn+1, indm:indm+1) = Mnm;
        M(indm:indm+1, indn:indn+1) = Mnm';
    end
end

