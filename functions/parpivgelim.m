function [L, U, P] = parpivgelim(A)

% This function performs standard Gaussian elimination with pivoting.
% It assumes the input matrix A is square. The result is the decomposition
% PA = LU.
%
% Input: A, m x m matrix.
%
% Outputs: L, m x m unit lower triangular matrix
%          U, m x m upper triangular matrix
%          P, m x m permutation matrix


M = size(A,1); % get dimensions of A

% Initialize U, L, and P.
U = A; 
L = eye(M);
P = eye(M);

for k = 1:M-1
   
    [~, i] = max(abs(U(k:M,k))); % search column for pivot
    i = i + k - 1; % correct index
    
    U([i k],k:M) = U([k i],k:M); % exchange rows in U
    L([i k],1:k-1) = L([k i],1:k-1); % exchange rows in L
    P([i k],:) = P([k i],:); % exchange row in P
    
    ukk = U(k,k); %set the pivot.
    
    for j = k+1:M
        Ljk = U(j,k)/ukk; %obtain entry L_ik 
        U(j,k:M) = U(j,k:M) - Ljk*U(k,k:M); % multiplication of jth row by matrix L_k
        L(j,k) = Ljk;
    end
end