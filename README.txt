The initalization of data is performer using the files "init.txt" and "energies.txt".

In the file "init.txt" three int numbers should be given. The first one corresponds to the total number of levels given by N=\sum d_j, where d_j corresponds to the degeneracy of the j-th level. The second value gives the total number of the Cooper pairs. The third number corresponds to the type of problem (inital energy values) according to the following scheme:
0 - two levels model with initial energies given by -1 and +1.
1 - non-degenerate picket-fence model with the initial energies given by 2i+1, where i = 0,1,2,...,N
2 - non-degenerate model with energies proportional to (i+1)^2 (i = 0,1,2,...,N). 
3 - non-degenerate model with energies proportional to -(i+1)^(-2) (i = 0,1,2,...,N).
4 - degenerate 2D infinite square-potential (2D ISP) model for five levels (e = {2,5,8,10,13})
5 - degenerate picket-fence model for five levels (e = {2,4,6,8,10})

Example:
16 16 0 - corresponds to the 16 Cooper pairs distributed among two levels +1 and -1
12 6 1 - corresponds to the 6 Cooper pairs distributed among 12 levels in the picket-fence model
12 8 4 - corresponds to the 8 Cooper pairs distributed among 5 levels in the 2D ISP model such that the total degeneracy of the unoccupied (or not fully occupied) levels is equal to 12-8 = 4.

In the file "energies.txt" the first line corresponds to the occupation number for each level.

Example (in the case of "init.txt": 16 16 0):
10 6 - corresponds to 10 pairs occupying the level -1 and 6 pairs occupying the level +1
Example (in the case of "init.txt": 12 6 5):
1 2 3 0 0 - corresponds to 1 pair occupying the level 2, 2 pairs occupying the level 4, etc.

For non-degenerate cases (1,2,3) this line should be ommited.

The second line is responsible for the occupation (initial values of Lambda). The total number of entries should be equal to the number of levels while each number should be either 1 or 0. The total sum of the line should be equal to the number of pairs. In the degenerate cases, occupation should correspond to the values given in the first line.

Example (for "init.txt" 12 6 5 and the first line: 1 2 3 0 0)
1 1 1 1 1 1 0 0 0 0 0 0

In the non-degenerate cases, this line is responsible for the occupation.

Example (for "init.txt" 12 6 1):
1 1 1 1 1 1 0 0 0 0 0 0 - pairs occupy the following levels: {1, 3, 5, 7, 9, 11}
1 0 1 0 1 0 1 0 1 0 1 0 - pairs occupy the following levels: {1, 5, 9, 13, 17, 21}
1 1 1 1 1 0 1 0 0 0 0 0 - pairs occupy the following levels: {1, 3, 5, 7, 9, 13}
