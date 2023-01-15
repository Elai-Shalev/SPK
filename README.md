    # SPK
    SPK - The Normalized Spectral Clustering Algorithm 
    Implementation in C and Pythin, utilizing C/Python-API 

    Algorithm and method based on:
    1. Andrew Ng, Michael Jordan, and Yair Weiss. On spectral clustering: Analysis and an algorithm.
    Advances in neural information processing systems, 14:849–856, 2001.
    2. Ulrike Von Luxburg. A tutorial on spectral clustering. Statistics and computing, 17(4):395–416, 2007.

    === The Algorithm ===
    1: Form the weighted adjacency matrix W from X
    2: Compute the normalized graph Laplacian Lnorm
    3: Determine k and obtain the largest k eigenvectors u1, . . . , uk of Lnorm
    4: Let U ∈ Rn×k be the matrix containing the vectors u1, . . . , uk as columns
    5: Form the matrix T ∈ Rn×k from U by renormalizing each of U’s rows to have unit length,
    that is set tij = uij / sum(uij ** 2)**0.5)
    6: Treating each row of T as a point in Rk, cluster them into k clusters via the K-means algorithm

    == FILES ===

    --- spkmeans.c ---
    This c program provides the main functionality of the algorith. 
    1. This program reads user CMD input, including a "goal" (specified below), and an input file name.
    2. The program will calculate and output according to the selected goal
    3. This file provides the C implementation for all algorithm goals,
    providing service to the python interface to interact with the c-Python API 
    and reduce computation time and memory. 

    The goals:
    wam: Calculate and output the Weighted Adjacency Matrix.
    ddg: Calculate and output the Diagonal Degree Matrix.
    lnorm: Calculate and output the Normalized Graph Laplacian.
    jacobi: Calculate and output the eigenvalues and eigenvectors.

    ---spkmeans.py ---
    1. This program reads user CMD input, including a parameter "k", a "goal" (specified below),
    and an input file name. 
    2. With goal=spk, The program will implement the full spk-means algorithm. 
    3. Interfacing with the C extension spkmeansmodule, all implementations of the different goals
    will be performed by calling the C extension and calculated in C. 

    The goals as seen above in spkmeans.c.

    ---spkmeansmodule.py---
    In this file we define our C extension which will wrap the algorithm implemented in spkmeans.c.
    This uses only functions implemented  in spkmeans.c. 


    Code by Elai Shalev and Ido Meshulam, 2022
