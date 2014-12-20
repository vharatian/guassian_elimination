guassian_elimination
====================

multi thread implementation of guassian elemination

code can be compiled with the following command on linux distributions

$ gcc guassianElimination.c -fopenmp -pthread -o guassianElimination

after compiling executable file can be runned with the following command

$ ./guassianElimination -ds 10 1

which run the programm on debug mode '-d' with a serial method '-s'
the 10 number specify size of matrix and the 1 specify number of threads which is used for OpenMp '-o' and Pthread 'p' methods.


