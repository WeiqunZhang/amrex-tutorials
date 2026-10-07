amrex_tutorial_add(ForkJoin
    SOURCES main.cpp
    INPUTS inputs
    REQUIRES MPI LSOLVERS
    NPROCS 4)
