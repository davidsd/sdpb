# SDPLIB 1.2, June 28, 1999.

> This GitHub repository is a "modern" presentation of the original SDPLIB 1.2
> library by  Brian Borchers hosted at
> http://euler.nmt.edu/~brian/sdplib/sdplib.html.

If you use SDPLIB and wish to cite it, please refer to:

> Borchers, B., SDPLIB 1.2, A Library of Semidefinite Programming Test Problems.
> Optimization Methods and Software. 11(1):683-690, 1999.

This is version 1.2 of SDPLIB.  The only change in this version of
SDPLIB is that problem `hinf37` has been removed- the problem was
removed because it was identical to problem `hinf9`.

The `data` directory contains the SDPLIB collection of semidefinite
programming test problems.  There are a total of 92 problems in the
collection.  All problems are stored in the SDPA sparse format [[5]](#ref5).

The following table describes the problems, and gives optimal objective
values that have been computed by SDPA and cross checked against results
from the problem originators.  Note that in some cases, very slight changes
in the optimal objective function value have occurred as a result of the
conversion into SDPA format.  Furthermore, different authors have adopted
different conventions for the primal and dual SDP problems.  The results
reported here are based on the SDPA conventions- thus some objective
function values have their signs changed, and primal or dual infeasibility
means infeasibility with respect to SDPA's primal and dual.

| Problem   |    m |    n | Optimal Objective Value | Notes |
| --------- | ---: | ---: | ----------------------: | :---: |
| arch0     |  174 |  335 |  5.66517e-01            |     1 |
| arch2     |  174 |  335 |  6.71515e-01            |     1 |
| arch4     |  174 |  335 |  9.726274e-01           |     1 |
| arch8     |  174 |  335 |  7.05698e+00            |     1 |
| control1  |   21 |   15 |  1.778463e+01           |     2 |
| control2  |   66 |   30 |  8.300000e+00           |     2 |
| control3  |  136 |   45 |  1.363327e+01           |     2 |
| control4  |  231 |   60 |  1.979423e+01           |     2 |
| control5  |  351 |   75 |  1.68836e+01            |     2 |
| control6  |  496 |   90 |  3.73044e+01            |     2 |
| control7  |  666 |  105 |  2.06251e+01            |     2 |
| control8  |  861 |  120 |  2.0286e+01             |     2 |
| control9  | 1081 |  135 |  1.46754e+01            |     2 |
| control10 | 1326 |  150 |  3.8533e+01             |     2 |
| control11 | 1596 |  165 |  3.1959e+01             |     2 |
| eqaulG11  |  801 |  801 |  6.291553e+02           |     3 |
| equalG51  | 1001 | 1001 |  4.005601e+03           |     3 |
| gpp100    |  101 |  100 | -4.49435e+01            |     4 |
| gpp124-1  |  125 |  124 | -7.3431e+00             |     4 |
| gpp124-2  |  125 |  124 | -4.68623e+01            |     4 |
| gpp124-3  |  125 |  124 | -1.53014e+02            |     4 |
| gpp124-4  |  125 |  124 | -4.1899e+02             |     4 |
| gpp250-1  |  250 |  250 | -1.5445e+01             |     4 |
| gpp250-2  |  250 |  250 | -8.1869e+01             |     4 |
| gpp250-3  |  250 |  250 | -3.035e+02              |     4 |
| gpp250-4  |  250 |  250 | -7.473e+02              |     4 |
| gpp500-1  |  501 |  500 | -2.53e+01               |     4 |
| gpp500-2  |  501 |  500 | -1.5606e+02             |     4 |
| gpp500-3  |  501 |  500 | -5.1302e+02             |     4 |
| gpp500-4  |  501 |  500 | -1.56702e+03            |     5 |
| hinf1     |   13 |   14 |  2.0326e+00             |     5 |
| hinf2     |   13 |   16 |  1.0967e+01             |     5 |
| hinf3     |   13 |   16 |  5.69e+01               |     5 |
| hinf4     |   13 |   16 |  2.74764e+02            |     5 |
| hinf5     |   13 |   16 |  3.63e+02               |     5 |
| hinf6     |   13 |   16 |  4.490e+02              |     5 |
| hinf7     |   13 |   16 |  3.91e+02               |     5 |
| hinf8     |   13 |   16 |  1.16e+02               |     5 |
| hinf9     |   13 |   16 |  2.3625e+02             |     5 |
| hinf10    |   21 |   18 |  1.09e+02               |     5 |
| hinf11    |   31 |   22 |  6.59e+01               |     5 |
| hinf12    |   43 |   24 |  2e-1                   |     5 |
| hinf13    |   57 |   30 |  4.6e+01                |     5 |
| hinf14    |   73 |   34 |  1.30e+01               |     5 |
| hinf15    |   91 |   37 |  2.5e+01                |     5 |
| infd1     |   10 |   30 |  dual infeasible        |     6 |
| infd2     |   10 |   30 |  dual infeasible        |     6 |
| infp1     |   10 |   30 |  primal infeasible      |     6 |
| infp2     |   10 |   30 |  primal infeasible      |     6 |
| maxG11    |  800 |  800 |  6.291648e+02           |     7 |
| maxG32    | 2000 | 2000 |  1.567640e+03           |     7 |
| maxG51    | 1000 | 1000 |  4.003809e+03           |     7 |
| maxG55    | 5000 | 5000 |  9.999210e+03           |     7 |
| maxG60    | 7000 | 7000 |  1.522227e+04           |     7 |
| mcp100    |  100 |  100 |  2.261574e+02           |     8 |
| mcp124-1  |  124 |  124 |  1.419905e+02           |     8 |
| mcp124-2  |  124 |  124 |  2.698802e+02           |     8 |
| mcp124-3  |  124 |  124 |  4.677501e+02           |     8 |
| mcp124-4  |  124 |  124 |  8.644119e+02           |     8 |
| mcp250-1  |  250 |  250 |  3.172643e+02           |     8 |
| mcp250-2  |  250 |  250 |  5.319301e+02           |     8 |
| mcp250-3  |  250 |  250 |  9.811726e+02           |     8 |
| mcp250-4  |  250 |  250 |  1.681960e+03           |     8 |
| mcp500-1  |  500 |  500 |  5.981485e+02           |     8 |
| mcp500-2  |  500 |  500 |  1.070057e+03           |     8 |
| mcp500-3  |  500 |  500 |  1.847970e+03           |     8 |
| mcp500-4  |  500 |  500 |  3.566738e+03           |     8 |
| qap5      |  136 |   26 | -4.360e+02              |     9 |
| qap6      |  229 |   37 | -3.8144e+02             |     9 |
| qap7      |  358 |   50 | -4.25e+02               |     9 |
| qap8      |  529 |   65 | -7.57e+02               |     9 |
| qap9      |  748 |   82 | -1.410e+03              |     9 |
| qap10     | 1021 |  101 | -1.093e+01              |     9 |
| qpG11     |  800 | 1600 |  2.448659e+03           |    10 |
| qpG51     | 1000 | 2000 |  1.181000e+03           |    10 |
| ss30      |  132 |  426 |  2.02395e+01            |     1 |
| theta1    |  104 |   50 |  2.300000e+01           |    11 |
| theta2    |  498 |  100 |  3.287917e+01           |    11 |
| theta3    | 1106 |  150 |  4.216698e+01           |    11 |
| theta4    | 1949 |  200 |  5.032122e+01           |    11 |
| theta5    | 3028 |  250 |  5.723231e+01           |    11 |
| theta6    | 4375 |  300 |  6.347709e+01           |    11 |
| thetaG11  | 2401 |  801 |  4.000000e+02           |    12 |
| thetaG51  | 6910 | 1001 |  3.49000e+02            |    12 |
| truss1    |    6 |   13 | -8.999996e+00           |    13 |
| truss2    |   58 |  133 | -1.233804e+02           |    13 |
| truss3    |   27 |   31 | -9.109996e+00           |    13 |
| truss4    |   12 |   19 | -9.009996e+00           |    13 |
| truss5    |  208 |  331 | -1.326357e+02           |    13 |
| truss6    |  172 |  451 | -9.01001e+02            |    13 |
| truss7    |   86 |  301 | -9.00001e+02            |    13 |
| truss8    |  496 |  628 | -1.331146e+02           |    13 |


### Table notes:

1. These truss topology design problems were contributed by Katsuki Fujisawa.
   They are originally from [[8]](#ref8).

2. These problems from control and system theory were contributed by Katsuki
   Fujisawa [[4]](#ref4).

3. These graph equipartition problems were supplied by Steve Benson
   [[2]](#ref2).  The random graphs were originally generated by Christoph Helmberg and Franz Rendl [[7]](#ref7).

4. These graph partitioning problems were contributed by Katsuki
   Fujisawa [[4]](#ref4)[[6]](#ref6).

5. These linear matrix inequalities from control systems engineering are taken
   from the SDPPACK web site [[1]](#ref1).  The problems were originally
   developed by P. Gahinet.

6. These infeasible problems were generated by a MATLAB procedure provided by
   Mike Todd.

7. These max cut problems were supplied by Steve Benson [[2]](#ref2).  The
   random graphs were originally generated by Christoph Helmberg and Franz
   Rendl [[7]](#ref7).

8. These max cut problems were contributed by Katsuki Fujisawa [[4]](#ref4).

9. These quadratic assignment problems were contributed by Katsuki Fujisawa
   [[6]](#ref6).

10. These SDP relaxations of box constrained quadratic programming problems
    were supplied by Steve Benson [[2]](#ref2).  The random graphs were
    originally generated by Christoph Helmberg and Franz Rendl [[7]](#ref7).

11. These Lovasz theta problems are taken from [[3]](#ref3).

12. These Lovasz theta problems were contributed by Steve Benson [[2]](#ref2).
    The random graphs were originally generated by Christoph Helmberg and Franz
    Rendl [[7]](#ref7).

13. These truss topology design problems are taken from the SDPPACK web site
    [[1]](#ref1).  The problems were originally developed by A. Nemirovski.



## SDPA sparse format

The SDP problems in this directory are all encoded in the SDPA sparse
format [[5]](#ref5).  This file provides a description of the file format.

### The SDP Problem:

We work with a semidefinite programming problem that has been written in the
following standard form:

    (P)    min   c1*x1+c2*x2+...+cm*xm
           s.t.  F1*x1+F2*x2+...+Fm*xm - F0 = X
                                         X >= 0

The dual of the problem is:

    (D)    max   tr(F0*Y)
           s.t.  tr(Fi*Y) = ci    i = 1,2,...,m
                       Y >= 0

Here all of the matrices `F0`, `F1`, ..., `Fm`, `X`, and `Y` are assumed to be
symmetric of size `n` by `n`.  The constraints `X >= 0` and `Y >= 0` mean that
`X` and `Y` must be positive semidefinite.

Note that several other standard forms for SDP have been used by a number of
authors- these can generally be translated into the SDPA standard form with
little effort.

### File Format:

The file consists of six sections:

1. Comments. The file can begin with arbitrarily many lines of comments.
   Each line of comments must begin with `"` or `*`.

2. The first line after the comments contains `m`, the number of constraint
   matrices.  Additional text on this line after `m` is ignored.

3. The second line after the comments contains `nblocks`, the number of
   blocks in the block diagonal structure of the matrices.  Additional text
   on this line after `nblocks` is ignored.

4. The third line after the comments contains a vector of numbers that
   give the sizes of the individual blocks.  The special characters
   `,`, `(`, `)`, `{`, and `}` can be used as punctuation and are ignored.
   Negative numbers may be used to indicate that a block is actually a diagonal
   submatrix.  Thus a block size of `-5` indicates a 5 by 5 block in which
   only the diagonal elements are nonzero.

5. The fourth line after the comments contains the objective function
   vector `c`.

6. The remaining lines of the file contain entries in the constraint
   matrices, with one entry per line.  The format for each line is

       <matno> <blkno> <i> <j> <entry>

Here `<matno>` is the number of the matrix to which this entry belongs,
`<blkno>` specifies the block within this matrix, `<i>` and `<j>` specify a
location within the block, and `<entry>` gives the value of the entry in the
matrix.  Note that since all matrices are assumed to be symmetric, only
entries in the upper triangle of a matrix are given.

### Example:

    min   10 * x1 + 20 * x2

    s.t.  X = F1 * x1 + F2 * x2 - F0
          X >= 0

where

    F0 = [1 0 0 0
          0 2 0 0
          0 0 3 0
          0 0 0 4]

    F1 = [1 0 0 0
          0 1 0 0
          0 0 0 0
          0 0 0 0]

    F2 = [0 0 0 0
          0 1 0 0
          0 0 5 2
          0 0 2 6]

In SDPA sparse format, this problem can be written as:

    "A sample problem.
    2 =mdim
    2 =nblocks
    {2, 2}
    10.0 20.0
    0 1 1 1 1.0
    0 1 2 2 2.0
    0 2 1 1 3.0
    0 2 2 2 4.0
    1 1 1 1 1.0
    1 1 2 2 1.0
    2 1 2 2 1.0
    2 2 1 1 5.0
    2 2 1 2 2.0
    2 2 2 2 6.0



## References:

1. <a id="ref1"></a>
   F. Alizadeh, J.P. Haberly, M. V. Nayakkankuppam, M. L. Overton, and
   S. Schmieta.  **SDPpack user's guide- version 0.9 beta**, Technical Report
   1997-737, Courant Institute of Mathematical Sciences, NYU, New York NY,
   June 1997. https://cs.nyu.edu/media/publications/TR1997-737.pdf

2. <a id="ref2"></a>
   S. J. Benson, Y. Ye, and X. Zhang. **Solving Large-Scale Sparse Semidefinite
   Programs for Combinatorial Optimization**, SIAM J. Optim., 10-2 (2000),
   443–461.
   [DOI: 10.1137/S1052623497328008](https://doi.org/10.1137/S1052623497328008)

3. <a id="ref3"></a>
   B. Borchers. **CSDP, a C library for semidefinite programming**,
   Optimization Methods and Software, 11:1-4 (1999), 613-623,
   [DOI: 10.1080/1055678990880576](https://doi.org/10.1080/1055678990880576)

4. <a id="ref4"></a>
   K. Fujisawa, M. Fukuda, M. Kojima, and K. Nakata. **Numerical Evaluation of
   SDPA (Semidefinite Programming Algorithm)**. In: Frenk H., Roos K.,
   Terlaky T., Zhang S. (eds) **High Performance Optimization.** Applied
   Optimization, vol 33., Springer (2000), Boston, MA.
   [DOI: 10.1007/978-1-4757-3216-0_11](https://doi.org/10.1007/978-1-4757-3216-0_11)

5. <a id="ref5"></a>
   K. Fujisawa, M. Kojima, and K. Nakata. **SDPA (Semidefinite Programming
   Algorithm) User's Manual**, Technical Report B-308, Department of
   Mathematical and Computing Sciences, Tokyo Institute of Technology.
   Revised, May, 1998. http://www.is.titech.ac.jp/~kojima/articles/b-308.pdf

6. <a id="ref6"></a>
   K. Fujisawa, M. Kojima and K. Nakata. **Exploiting Sparsity in Primal-Dual
   Interior-Point Methods for Semidefinite Programming**, Mathematical
   Programming, 79(1997):235-253.
   [DOI: 10.1007/BF02614319](https://doi.org/10.1007/BF02614319)

7. <a id="ref7"></a>
   C. Helmberg and F. Rendl. **A spectral bundle method for semidefinite
   programming**, SIAM J. Optim., 10-3 (2000), 673–696.
   [DOI: 10.1137/S1052623497328987](https://doi.org/10.1137/S1052623497328987)

8. <a id="ref8"></a>
   T. Nakamura and M. Ohsaki. **A Natural Generator of Optimum Topology of
   Plane Trusses for Specified Fundamental Frequency**, Computer Methods in
   Applied Mechanics and Engineering 94(1992):113-129.
   [DOI: 10.1016/0045-7825(92)90159-H](https://doi.org/10.1016/0045-7825%2892%2990159-H)



## Old CHANGELOG

- June 28, 1999: Removed problem `hinf37`, which was identical to problem
  `hinf9`.

- July 24, 1998: Reduced the size of many problem files by removing extraneous
  and unneeded zeros. For example, `1.00000e+00` is shortened to `1.0`.

- July 27, 1998: Slightly more accurate optimal objective function values for
  some problems.

- September 18, 1998: Added 11 new problems provided by Steve Benson.
