bValues3a= [4];
NValues3a = [100 200 400 800 1600];
[normsF3a,times3a] = M2SystemSolver(bValues3a,NValues3a);
times3a = times3a';

%Assuming that the script Question1bScript has alredy been run and the
%vector normsF1b stored, check that the difference between the respective
%norms is smaller than tolerance = 10^(-13)
tolerance3a = 10^(-13);
checker3a = zeros(length(normsF3a),1);
for i=1:length(normsF3a)
    if abs(normsF1b(i) - normsF3a(i)) < tolerance3a
        checker3a(i) = 1;
    end
end

checker3a

%Find the coefficients of the least squares regression line such that
%log(times3a) = coeff(1)*log(NValues3a) + coeff(2)
coeff3a = polyfit(log(NValues3a),log(times3a),1);

coeff3a

%Obtain a log-log plot of time against N, and also plot the corresponding
%curve found using polyfit
figure()
loglog(NValues3a,times3a)
hold on
loglog(NValues3a,exp(coeff3a(2)).*NValues3a.^coeff3a(1))
title('log-log Graph Of Time Against N Using LU Factorisation On M_{2}')
xlabel('Number of particles, N')
ylabel('Time, t')
legend('Recorded Data','t = e^c N^k')
hold off
xlim([100 1600])