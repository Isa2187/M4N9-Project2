bValues1c = [4];
NValues1c = [100 200 400];
[normsF1c,times1c,residuals1c,normsLU1c,normsM1c] = systemSolver(bValues1c,NValues1c);
residuals1c
normsLU1c
normsM1c