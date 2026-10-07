# HelloWorld ★

The smallest AMReX program: initialize and finalize AMReX, read a runtime
parameter, and print from one or all MPI ranks.

**Key APIs:** `amrex::Initialize`, `amrex::Finalize`, `amrex::Print`,
`amrex::AllPrint`, `amrex::ParmParse`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_HelloWorld -j
cd build/Base/HelloWorld && mpiexec -n 2 ./HelloWorld name=AMReX
```

GNU Make, with AMReX cloned next to this repository:

```
make -j
mpiexec -n 2 ./main3d.gnu.TPROF.MPI.ex name=AMReX
```

## Expected output

```
Hello AMReX from AMReX version 26.10
  MPI ranks:      2
  ...
```
