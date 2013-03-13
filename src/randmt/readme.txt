% Mersenne Twister MT19937 pseudorandom number generator

WARNING:
    Do NOT use for cryptographic purposes. Read Internet RFC4086,   
        http://tools.ietf.org/html/rfc4086


# OVERVIEW

This software contains a high-quality pseudorandom number generator,
the Mersenne Twister MT19937 by Takuji Nishimura and Makoto Matsumoto, and
routines for sampling from several distributions.


# LICENSE

This software is distributed under the BSD license.  See the included file
license.txt for details.


# USAGE

To use randmt in a program, only the files randmt.h and randmt.c are needed;
all other included files are for documentation and testing purposes.  Include
randmt.h and call init_randmt_auto() at the beginning of the program to
initialize the pseudorandom number generator.  The initialization function
seeds the generator with the current time so that different numbers are 
produced on each run of the program.

Uniform random numbers are generated using the rand_unif() function.  Samples
from several other distributions can be generated as well, see below.

Example usage:

    #include <stdio.h>
    #include "randmt.h"
    
    int main(void)
    {
        int i;
        init_randmt_auto();
        for(i = 0; i < 20; i++)
            printf("%10.8f\n", rand_unif());
        return 0;
    }

This program can be compiled with GCC using the included makefile as

   make -f makefile.gcc example

or with MSVC,

   nmake -f makefile.vc example.exe

Pseudorandom numbers can be generated for several different distributions:

    rand_unif       the continuous uniform distribution on (0,1)
    rand_uint32     32-bit unsigned integers uniformly on [0,0xFFFFFFFF]
    rand_normal     the standard normal distribution
    rand_exp        exponential distribution
    rand_gamma      Gamma distribution
    rand_geometric  geometric distribution
    rand_poisson    Poisson distribution

Detailed documentation of these functions can be found online at
 
   http://www.getreuer.info/home/randmt

and in the comments in randmt.h.  Alternatively, HTML documentation can be
generated using Doxygen <http://www.stack.nl/~dimitri/doxygen/index.html> by
running the command

   doxygen doxygen.conf


# REENTRANT VERSIONS

For use in multithreaded applications, reentrant versions of the functions
are also included which have the same name suffixed with "_r." For these 
functions, the generator state is passed using a randmt_t object.

    rand_unif_r       the continuous uniform distribution on (0,1)
    rand_uint32_r     32-bit unsigned integers uniformly on [0,0xFFFFFFFF]
    rand_normal_r     the standard normal distribution
    rand_exp_r        exponential distribution
    rand_gamma_r      Gamma distribution
    rand_geometric_r  geometric distribution
    rand_poisson_r    Poisson distribution

A randmt_t object represents the state of an MT19937 pseudo-random number 
generator (PRNG).  Use the following functions to create, initialize, and
destroy randmt_t objects:

    new_randmt          create a new randmt_t
    init_randmt_auto_r  initialize randmt_t with time and address
    init_randmt_r       initialize randmt_t with a specified seed value
    free_randmt         free memory associated with a randmt_t

See the online documentation for details.


# TEST PROGRAM

To verify the distributions of the samplers, a test program is included. 
Compile the test program using

    make -f makefile.gcc test

or with MSVC,

    nmake -f makefile.vc test.exe

The test program applies the Kolmogorov-Smirnov and chi-squared tests to 
verify that the pseudorandom samplers produce the intended distributions. 
Note that the output of this program is different on each run.

Typical output is shown below:

For each random number generator, we sample N=1000000
values and compare the sample distribution to the theoretical
density function with the Kolmogorov-Smirnov test.

    Sampler              D
    rand_unif()          0.001095
    rand_normal()        0.001077
    rand_exp(1)          0.000811
    rand_gamma(0.2,1)    0.000516
    rand_gamma(  1,1)    0.000906
    rand_gamma(  2,1)    0.001051
    rand_gamma( 20,1)    0.001079

Supposing the distributions are correct, the D values should be
small with high probability:
    D < 0.001358 with probability 0.95
    D < 0.001627 with probability 0.99
    D < 0.001949 with probability 0.999

We apply the chi-squared test to verify the distributions of the
geometric and Poisson generators (the Kolmogorov-Smirnov test
applies only to continuous distributions).

    Sampler              p-value
    rand_geometric(0.1)  0.883482
    rand_geometric(0.5)  0.527853
    rand_geometric(0.9)  0.651309
    rand_poisson(0.2)    0.700759
    rand_poisson(  1)    0.551421
    rand_poisson(  2)    0.619268
    rand_poisson( 20)    0.656092
    rand_poisson(200)    0.422257

Supposing the distributions are correct, the p-values should be
above zero with high probability:
    p-value > 0.05 with probability 0.95
    p-value > 0.01 with probability 0.99
    p-value > 0.001 with probability 0.999


This material is based upon work supported by the National Science
Foundation under Award No. DMS-1004694. Any opinions, findings, and
conclusions or recommendations expressed in this material are those of
the author(s) and do not necessarily reflect the views of the National
Science Foundation.
