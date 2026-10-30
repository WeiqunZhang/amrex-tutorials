# ForkJoin ★★

Split the MPI ranks into independent tasks with `amrex::ForkJoin`, give each
task its own share of the data, run work on each task concurrently, and
gather the results back.

The program solves `ncomp` independent variable-coefficient problems
`a u - b div(beta grad u) = f` (one per MultiFab component) with MLMG.
`ForkJoin` splits the ranks into `ntasks` tasks; each task gets a subset of
the components of `rhs` and `phi` (`Strategy::split`) and a full copy of the
shared coefficients (`Strategy::duplicate`), and solves its problems on its
own MPI communicator. After the join, `phi` holds all solutions, and each is
compared with the exact solution.

**Key APIs:** `amrex::ForkJoin`, `ForkJoin::reg_mf`, `ForkJoin::reg_mf_vec`,
`ForkJoin::Strategy`, `ForkJoin::Intent`, `ForkJoin::fork_join`,
`ForkJoin::get_mf`, `ForkJoin::ComponentBounds`, `ParallelContext`,
`MLABecLaplacian`, `MLMG`

## Strategies and intents

| Strategy    | Each task gets                               |
|-------------|----------------------------------------------|
| `split`     | a contiguous subset of the components        |
| `duplicate` | a full copy                                  |
| `single`    | nothing, except the owner task (a full copy) |

`Intent::in` copies data into the tasks before the fork, `Intent::out`
copies it back after the join, and `Intent::inout` does both.
Inside `fork_join`, `ParallelContext` and `amrex::Print` refer to the task's
ranks, so AMReX functions (here MLMG) run on the task only. Set
`task_output_dir=<dir>` to send each task's output to its own file.

## Build and run

Requires MPI. CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Coupling
cmake --build build --target Coupling_ForkJoin -j 8
cd build/Coupling/ForkJoin && mpiexec -n 4 ./ForkJoin inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 4 ./main3d.gnu.MPI.ex inputs
```

Try `ntasks=4` or `ncomp=8`. The number of tasks is limited by the number
of ranks and components.

## Expected output

```
Solving 4 problems with 2 ForkJoin tasks on 4 ranks
  Task 0: 2 rank(s), components [0, 2)
  Task 1: 2 rank(s), components [2, 4)
Problem 0: max error of phi/(n+1) = 0.004394236913
...
```
