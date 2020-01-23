bValues3b= [4];
NValues3b = [100 200 400 800 1600];
[normsF3b,times3b] = levinsonSystemSolver(bValues3b,NValues3b);
times3b = times3b';

%Assuming that the script Question1bScript has alredy been run and the
%vector normsF1b stored, check that the difference between the respective
%norms is smaller than tolerance = 10^(-15)
tolerance3b = 10^(-13);
checker3b = zeros(length(normsF3b),1);
for i=1:length(normsF3b)
    if abs(normsF1b(i) - normsF3b(i)) < tolerance3b
        checker3b(i) = 1;
    end
end

checker3b

%Find the coefficients of the least squares regression line such that
%log(times3a) = coeff(1)*log(NValues3a) + coeff(2), but only looking at the
%values N=400, 800 and 1600
coeff3b = polyfit(log(NValues3b(3:end)),log(times3b(3:end)),1);

coeff3b

%Obtain a log-log plot of time against N, and also plot the corresponding
%curve found using polyfit
figure()
loglog(NValues3b,times3b)
hold on
loglog(NValues3b,exp(coeff3b(2)).*NValues3b.^coeff3b(1))
title('log-log Graph Of Time Against N Using The Levinson Algorithm On M_{2}')
xlabel('Number of particles, N')
ylabel('Time, t')
legend('Recorded Data','t = e^c N^k')
hold off
xlim([100 1600])