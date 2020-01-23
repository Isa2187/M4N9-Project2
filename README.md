# M4N9-Project2
Second project of Computational Linear Algebra (M4N9) module taken in 4th year. (Grade = 82.5%)

All code was done in MATLAB, and further details of the task can be found in the folder Project_Files.

By modelling the flow of equally-spaced particles through a viscous fliud, a system of linear equations can be used. By first using LU decomposition, the system was solved for different parameters e.g. number of particles and spacing between them, to reveal how this affected the tme taken to solve the system.

Since this matrix takes on a particular structure, it could be manipulated using a permutation matrix to give two N-by-N systems, rather than one 2N-by-2N system, halfing the compuational cost. The impact of changing parameters was then investigated using the permuted form of the system.
