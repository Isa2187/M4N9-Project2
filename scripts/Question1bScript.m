bValues1b = [4];
NValues1b = [100 200 400 800 1600];
[normsF1b,times1b,residuals1b,normsLU1b,normsM1b] = systemSolver(bValues1b,NValues1b);
times1b = times1b';

%Find the coefficients of the least squares regression line such that
%log(times1b) = coeff(1)*log(NValues1b) + coeff(2)
coeff1b = polyfit(log(NValues1b),log(times1b),1);

coeff1b

%Obtain a log-log plot of t against N, and also plot the corresponding
%curve found by polyfit
figure()
loglog(NValues1b,times1b)
hold on
loglog(NValues1b,exp(coeff1b(2)).*NValues1b.^coeff1b(1))
title('log-log Graph Of Time Against N Using LU Factorisation On M')
xlabel('Number of particles, N')
ylabel('Time, t')
legend('Recorded Data','t = e^c N^k')
hold off
xlim([100 1600])