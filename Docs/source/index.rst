AMReX Tutorials
===============

These tutorials show how to write codes with `AMReX
<https://github.com/AMReX-Codes/amrex>`__, a software framework for
block-structured adaptive mesh refinement. Each tutorial is a small,
self-contained program that runs on CPUs and GPUs, with MPI and OpenMP, and
builds with CMake or GNU Make. Most of them check their result against an
exact or reference answer, and all of them run as tests in CI.

The source is in the `amrex-tutorials repository
<https://github.com/AMReX-Codes/amrex-tutorials>`__. For the concepts and
the full API, see the `AMReX documentation
<https://amrex-codes.github.io/amrex/docs_html/>`__.

Organization
------------

Tutorials are organized by topic. Each one lives in ``<Topic>/<Name>/``
with a ``README.md``, an ``inputs`` file, a ``CMakeLists.txt`` and a
``GNUmakefile``. Each has a level:

- ★ introductory
- ★★ intermediate
- ★★★ advanced
- **Blueprint**: a complete code meant to be copied as the starting point
  of an application

New users should start with :doc:`GettingStarted`, then
:ref:`base-helloworld`, :ref:`base-multifab` and :ref:`base-heatequation`
in :doc:`Base`.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Topic
     - Contents
   * - :doc:`Base`
     - Core data structures, kernels, reductions, I/O, ``ParmParse``
   * - :doc:`FFT`
     - ``amrex::FFT`` transforms and FFT-based solvers
   * - :doc:`LinearSolvers`
     - MLMG, GMRES and related solvers
   * - :doc:`AMR`
     - Adaptive mesh refinement with ``AmrCore``
   * - :doc:`EB`
     - Embedded boundaries
   * - :doc:`Particles`
     - Particle containers and particle-mesh operations
   * - :doc:`TimeIntegration`
     - ``amrex::TimeIntegrator`` and SUNDIALS
   * - :doc:`Coupling`
     - Concurrent tasks on subsets of ranks (ForkJoin) and multiple programs (MPMD)
   * - :doc:`Python`
     - pyAMReX
   * - :doc:`ML`
     - Machine-learning inference
   * - :doc:`Applications`
     - Complete solvers combining several topics
   * - :doc:`BuildSystems`
     - Standalone projects that build against AMReX

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   GettingStarted
   VisualizationAndDebugging

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   Base
   FFT
   LinearSolvers
   AMR
   EB
   Particles
   TimeIntegration
   Coupling
   Python
   ML
   Applications
   BuildSystems

Questions and help
------------------

Questions are welcome on `AMReX's GitHub Discussions
<https://github.com/AMReX-Codes/amrex/discussions>`__. To add a tutorial,
see `CONTRIBUTING.md
<https://github.com/AMReX-Codes/amrex-tutorials/blob/main/CONTRIBUTING.md>`__.

AMReX and these tutorials are distributed under the 3-clause BSD license;
see `LICENSE <https://github.com/AMReX-Codes/amrex-tutorials/blob/main/LICENSE>`__.
