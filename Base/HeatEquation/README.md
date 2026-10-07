# HeatEquation ★

Solve the heat equation $\partial\phi/\partial t = D\nabla^2\phi$ with
explicit (forward Euler) time stepping on a periodic unit domain. The
initial condition $\phi = 1 + \prod_d \sin(2\pi x_d)$ decays as
$\exp(-4\pi^2 d D t)$ in $d$ dimensions, and the program compares the
numerical solution with this exact solution.

**Key APIs:** `MultiFab`, `Geometry`, `FillBoundary`, `ParallelFor(mf, ...)`,
`const_arrays`, `WriteSingleLevelPlotfile`

## Algorithm

Each step fills the ghost cells with `FillBoundary` and then computes
`phi_new = phi_old + dt * D * Laplacian(phi_old)` with the standard
second-order stencil. The time step is `dt = cfl * dx^2 / (2 * dim * D)`.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_HeatEquation -j
cd build/Base/HeatEquation && mpiexec -n 2 ./HeatEquation inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j
mpiexec -n 2 ./main3d.gnu.TPROF.MPI.ex inputs
```

Add `plot_int=10` to write a plotfile every 10 steps, and e.g.
`n_cell=128 stop_time=0.05` for a larger run.

## Expected output

```
Step 10, time = 0.00146484375
...
Finished 69 steps at time 0.01, mode amplitude = 0.3059442057, relative max error = 0.006383029421
PASSED
```
