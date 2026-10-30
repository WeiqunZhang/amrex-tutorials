.. _applications:

Applications
============

Complete solvers that combine several topics, meant to be copied as the
starting point of a new application. They may use a ``Source/`` directory
for the solver and ``Exec/<Problem>/`` directories for problem setups.
``CNS`` needs AMReX built with ``AmrLevel`` support (``AMReX_AMRLEVEL``,
on by default) and runs in 2D and 3D only. Build only this
topic with ``-DAMReXTutorials_TOPICS=Applications``; see
:doc:`GettingStarted`.

.. _applications-cns:

CNS
---

:Level: Blueprint
:Source: `Applications/CNS <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Applications/CNS>`__
:Key APIs: ``Amr``, ``AmrLevel``, ``LevelBld``, ``StateDescriptor``
   (``desc_lst``), ``DeriveList`` (``derive_lst``), ``AmrLevel::RK``,
   ``FluxRegister``, ``average_down``, ``AMRErrorTag``, ``GpuBndryFuncFab``,
   ``ParmParse::queryAdd``, ``ParallelFor``, ``ParReduce``

A compressible Euler solver built on ``amrex::Amr`` and ``amrex::AmrLevel``,
with adaptive mesh refinement, subcycling in time, refluxing and GPU
support. It uses a second-order Godunov method with characteristic slope
limiting and an approximate Riemann solver, advanced in time with
Runge-Kutta (``cns.rk_order`` = 2 or 3). It comes with two problems: a Sod
shock tube (``Exec/Sod``) and a Rayleigh-Taylor instability with gravity
(``Exec/RT``). At the end of a run, the code checks that the total mass
changed only by round-off, if no mass can enter or leave the domain.

Compared with :ref:`amr-advectionamrcore`, which manages level data itself
on top of ``AmrCore``, here ``AmrLevel`` provides the state data with old
and new time levels, ``FillPatch``, checkpoint/restart and plotfile I/O.

Code structure
^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - Contents
   * - ``Source/main.cpp``
     - Initialize, build ``Amr``, time loop, mass check
   * - ``Source/CNS.H``, ``Source/CNS.cpp``
     - The ``AmrLevel`` subclass, one object per level
   * - ``Source/CNS_setup.cpp``
     - State data, boundary conditions, derived variables, ``LevelBld``
   * - ``Source/CNS_advance.cpp``
     - ``advance()`` and the flux computation
   * - ``Source/CNS_hydro_K.H``
     - GPU kernels: primitive variables, slopes, Riemann solver
   * - ``Source/CNS_K.H``
     - GPU kernels: time step, temperature
   * - ``Source/CNS_derive.H``, ``Source/CNS_derive.cpp``
     - Derived plot variables (pressure, velocity)
   * - ``Source/CNS_index.H``, ``Source/CNS_parm.H``, ``Source/CNS_parm.cpp``
     - Component indices; physical parameters (EOS, gravity)
   * - ``Exec/<Problem>/cns_prob.H``
     - ``ProbParm`` and the initial condition ``cns_initdata``
   * - ``Exec/<Problem>/cns_prob.cpp``
     - Default parameters (``cns_prob_defaults``) and ``prob.*`` parameters

``Amr`` owns the grid hierarchy and calls the virtual functions of ``CNS``:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Function
     - Called
   * - ``variableSetUp`` (static)
     - Once, to define the state data and derived variables
   * - ``initData``
     - For each level at start-up
   * - ``computeInitialDt``, ``computeNewDt``
     - On level 0, to choose ``dt`` for all levels
   * - ``advance``
     - For each level and each (sub)step
   * - ``post_timestep``
     - After a level has caught up with the next finer level: reflux and
       average down
   * - ``errorEst``
     - When regridding, to tag cells for refinement
   * - ``init(old)``, ``init()``
     - When regridding, to fill new grids from old data or from the coarser
       level
   * - ``post_init``, ``postCoarseTimeStep``
     - After initialization and after each coarse step
   * - ``restart``
     - When restarting from a checkpoint

Driver and setup
^^^^^^^^^^^^^^^^

The driver builds ``Amr`` from a ``LevelBld`` and takes coarse time steps;
``Amr`` recursively advances the finer levels and writes plotfiles and
checkpoints according to the ``amr.*`` parameters:

.. literalinclude:: ../../Applications/CNS/Source/main.cpp
   :language: cpp
   :start-at: // Amr manages the hierarchy of CNS levels
   :end-before: // Write final checkpoint and plotfile
   :dedent:

The state is registered once in ``CNS::variableSetUp``: one cell-centered
``StateData`` with ``NUM_GROW`` ghost cells and conservative interpolation,
the component names with their ``BCRec`` boundary types and fill function,
and derived quantities computed on demand, e.g., for
``amr.derive_plot_vars``:

.. literalinclude:: ../../Applications/CNS/Source/CNS_setup.cpp
   :language: cpp
   :start-at: // One cell-centered StateData
   :end-at: derive_lst.addComponent("velocity"
   :dedent:

Time step and synchronization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``CNS::advance`` calls ``AmrLevel::RK``, which runs the Runge-Kutta stages
and fills the ghost cells of each stage, including coarse/fine boundaries
interpolated in space and time. ``compute_dSdt`` computes the fluxes with
``ParallelFor`` kernels and adds them to the flux registers:

.. literalinclude:: ../../Applications/CNS/Source/CNS_advance.cpp
   :language: cpp
   :start-at: // As the coarse side, this level adds
   :end-at: [&] (int /*stage*/, MultiFab& S) { computeTemp(S, 0); });
   :dedent:

After the finer level has caught up, ``post_timestep`` corrects the coarse
level with the fine fluxes and overwrites the covered cells with the
average of the fine data, which together keep the scheme conservative:

.. literalinclude:: ../../Applications/CNS/Source/CNS.cpp
   :language: cpp
   :start-at: if (level < parent->finestLevel()) {
   :end-before: // Resets the FillPatcher
   :dedent:

Making your own code
^^^^^^^^^^^^^^^^^^^^

- **New problem:** copy ``Exec/Sod`` to ``Exec/MyProblem``. Edit
  ``ProbParm`` and ``cns_initdata`` in ``cns_prob.H``, and the defaults and
  ``prob.*`` parameters in ``cns_prob.cpp``. Add the problem to the list in
  ``CMakeLists.txt``.
- **New physics:** add components in ``CNS_index.H`` and their names and
  boundary conditions in ``CNS_setup.cpp``, and change the kernels in
  ``CNS_hydro_K.H`` and ``compute_dSdt``. Inflow boundaries need values in
  ``CnsFillExtDir``.
- **New AmrLevel code:** keep the structure of ``main.cpp``,
  ``CNS_setup.cpp`` (``LevelBld``) and the ``AmrLevel`` overrides in
  ``CNS.cpp``, and replace ``advance``, the kernels and the tagging
  criteria. Kernels capture small parameter structs such as ``Parm`` by
  value; large tables belong in GPU memory, e.g., ``Gpu::DeviceVector``.

Building and running
^^^^^^^^^^^^^^^^^^^^

One executable is built per problem: targets ``Applications_CNS_Sod`` and
``Applications_CNS_RT``, binaries ``CNS_Sod`` and ``CNS_RT``. Both problems
have an ``inputs`` file, which is read from the source tree:

.. code-block:: console

   cmake --build build --target Applications_CNS_Sod -j 8
   cd build/Applications/CNS
   mpiexec -n 2 ./CNS_Sod ../../../Applications/CNS/Exec/Sod/inputs

After copying ``Applications/CNS`` out of the repository, its
``CMakeLists.txt`` builds it on its own and fetches AMReX (or uses a local
checkout given by ``-DFETCHCONTENT_SOURCE_DIR_AMREX``):

.. code-block:: console

   cmake -S . -B build -DAMReX_SPACEDIM=2 -DAMReX_GPU_BACKEND=CUDA
   cmake --build build -j 8

With GNU Make, build in a problem directory:

.. code-block:: console

   cd Applications/CNS/Exec/Sod
   make -j8
   mpiexec -n 2 ./CNS3d.gnu.MPI.ex inputs

Without an inputs file, each problem runs to its ``stop_time`` with the
defaults in ``cns_prob_defaults``. The ``inputs`` files stop after 20
steps; remove ``max_step`` and increase ``amr.n_cell`` and
``amr.max_level`` for a real run. Boundary types (``cns.lo_bc``,
``cns.hi_bc``) are 0 interior (periodic), 1 inflow (not implemented),
2 outflow, 3 symmetry, 4 slip wall and 5 no-slip wall. ``cns.v=2`` prints
the conserved totals after every coarse step.
