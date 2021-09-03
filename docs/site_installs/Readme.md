The *.org files in this directory contain notes on building and
running SDPB on various HPC machines.  Most of the systems use
modules, so there is typically a module invocation like

    module load gcc/8.3.0 openmpi/3.1.4 openblas/0.3.6

to load all necessary modules.  Generally, I tried to use as many of
the builtin modules and system libraries as possible.  This was not
always possible.  For example, `GMP` might be built without C++
support.  More commonly, `Boost` was built with the system compiler,
so it is incompatible with the C++ 14 compiler required for the rest
of the code.

At the end of the file, there is an example of running TTTT_small.  It
should take about 510 iterations and several hours, ending with primal
and dual optimal solution.  The values for the primal and dual objectives should be

    0.0062678518839620999534329032850645348339149554

to within 10^-30.
