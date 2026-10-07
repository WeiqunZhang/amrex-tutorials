.. _coupling:

Coupling
========

These tutorials run several computations in one MPI job: concurrent tasks
on subsets of the ranks of one program (``ForkJoin``), and separate
programs that exchange ``MultiFab`` data (``MPMD``). Both require MPI;
ForkJoin also needs the linear solvers. Build only this
topic with ``-DAMReXTutorials_TOPICS=Coupling``; see :doc:`GettingStarted`.

.. _coupling-forkjoin:

ForkJoin
--------

:Level: ★★
:Source: `Coupling/ForkJoin <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Coupling/ForkJoin>`__
:Key APIs: ``amrex::ForkJoin``, ``ForkJoin::reg_mf``, ``ForkJoin::reg_mf_vec``,
           ``ForkJoin::Strategy``, ``ForkJoin::Intent``, ``ForkJoin::fork_join``,
           ``ForkJoin::get_mf``, ``ForkJoin::ComponentBounds``, ``ParallelContext``,
           ``MLABecLaplacian``, ``MLMG``

``amrex::ForkJoin`` splits the MPI ranks into independent tasks, gives each
task its share of the data, runs a function on all tasks concurrently, and
gathers the results back. Inside the task, ``ParallelContext`` and
``amrex::Print`` refer to the task's ranks only, so AMReX functions such as
an MLMG solve run on the task's communicator.

.. figure:: figs/fork_join_tasks.png
   :scale: 80 %
   :align: center
   :alt: Fork-join of MPI ranks into two tasks

   The ranks of the parent task are split into two tasks that run
   independently and then join.

The program solves ``ncomp`` independent variable-coefficient problems
:math:`a u - b \nabla \cdot (\beta \nabla u) = f`, one per component of a
``MultiFab``, with MLMG. ``rhs`` and ``phi`` are registered with
``Strategy::split``, so each task gets a contiguous subset of their
components; the coefficients are registered with ``Strategy::duplicate``,
so each task gets a full copy. ``Intent::in`` copies data into the tasks
before the fork, ``Intent::out`` copies it back after the join, and
``Intent::inout`` does both. ``Strategy::single`` gives a copy to one owner
task only.

.. figure:: figs/mf_remap_hires.png
   :scale: 7 %
   :align: center
   :alt: MultiFab registration strategies

   How a MultiFab is distributed to the tasks with the ``duplicate``,
   ``split`` and ``single`` strategies.

Registration only records how each MultiFab is to be copied:

.. literalinclude:: ../../Coupling/ForkJoin/main.cpp
   :language: cpp
   :start-at: // Split the ranks evenly into ntasks tasks.
   :end-before: amrex::Print() << "Solving "
   :dedent:

The function passed to ``fork_join`` runs on every task, with the task's
part of the registered data:

.. literalinclude:: ../../Coupling/ForkJoin/main.cpp
   :language: cpp
   :start-at: fj.fork_join([&] (ForkJoin& f)
   :end-at: // After the join, phi holds the solutions from all tasks.
   :dedent:

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target Coupling_ForkJoin -j
   cd build/Coupling/ForkJoin
   mpiexec -n 4 ./ForkJoin inputs

Try ``ntasks=4`` or ``ncomp=8``; the number of tasks is limited by the
numbers of ranks and components. ``task_output_dir=<dir>`` sends each
task's output to its own file.

.. code-block:: text

   Solving 4 problems with 2 ForkJoin tasks on 4 ranks
     Task 0: 2 rank(s), components [0, 2)
     Task 1: 2 rank(s), components [2, 4)
   Problem 0: max error of phi/(n+1) = 0.004394236913
   ...

.. _coupling-mpmd:

MPMD
----

:Level: ★★★
:Source: `Coupling/MPMD <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Coupling/MPMD>`__
:Key APIs: ``MPMD::Initialize``, ``MPMD::Finalize``, ``MPMD::Copier``,
           ``Copier::send``, ``Copier::recv``,
           ``amrex::Initialize(argc, argv, true, comm)``;
           Python: ``MPMD_Initialize_without_split``, ``initialize_when_MPMD``,
           ``MPMD_Copier``

Two separately built programs run in the same MPI job (multiple program,
multiple data) and exchange ``MultiFab`` data with ``amrex::MPMD``. The
main program (``main.cpp``) owns the domain :math:`[0, n)^D`, fills
component 0 with a Gaussian, sends it, and receives component 1. The agent,
in C++ (``agent.cpp``) or Python with pyAMReX (``agent.py``), owns the
upper half :math:`[n/2, n)^D` with its own boxes, receives component 0,
sets component 1 to 10 times component 0, and sends it back. The main
program checks the round trip.

``MPMD::Initialize`` splits ``MPI_COMM_WORLD`` into one communicator per
program, and AMReX is initialized on it:

.. literalinclude:: ../../Coupling/MPMD/main.cpp
   :language: cpp
   :start-at: // Split MPI_COMM_WORLD into one communicator per program
   :end-at: amrex::Initialize(argc, argv, true, comm);
   :dedent:

``MPMD::Copier`` exchanges the BoxArrays and DistributionMappings of the two
programs once; ``send`` and ``recv`` then copy components on the
intersection of the two BoxArrays. The programs may use different boxes and
numbers of ranks, but must use the same dimension and precision. With
``Copier(ba, dm, true)`` on one side and ``Copier(true)`` on the other, the
second program adopts the BoxArray of the first.

.. literalinclude:: ../../Coupling/MPMD/main.cpp
   :language: cpp
   :start-at: // The Copier exchanges BoxArray and DistributionMapping
   :end-at: copier.recv(mf, 1, 1);
   :dedent:

The agent matches each ``send`` with a ``recv`` and vice versa:

.. literalinclude:: ../../Coupling/MPMD/agent.cpp
   :language: cpp
   :start-at: MPMD::Copier copier(ba, dm);
   :end-at: copier.send(mf, 1, 1);
   :dedent:

The Python agent splits the communicator with mpi4py:

.. literalinclude:: ../../Coupling/MPMD/agent.py
   :language: python
   :start-at: # Split MPI_COMM_WORLD with mpi4py
   :end-at: amr.MPMD_Finalize()

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target Coupling_MPMD_Main Coupling_MPMD_Agent -j
   cd build/Coupling/MPMD
   mpiexec -n 2 ./MPMD_Main inputs : -n 2 ./MPMD_Agent inputs
   mpiexec -n 2 ./MPMD_Main inputs : -n 2 python3 agent.py inputs

The Python agent needs pyAMReX built with MPI, and mpi4py; add ``dim=2``
for a 2D build. With ``-DAMReXTutorials_PYTHON=ON``, CTest also runs the
Python agent. With GNU Make, build the two programs with ``make PROG=main``
and ``make PROG=agent``.

.. code-block:: text

   Main: 2 of 4 ranks in MPI_COMM_WORLD
   Agent: 2 of 4 ranks in MPI_COMM_WORLD
   Agent: received component 0 with min = 1, max = 1.92937572
   Main: max error of the round trip = 0
