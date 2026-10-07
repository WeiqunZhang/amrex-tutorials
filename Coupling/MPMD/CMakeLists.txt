# Two programs that exchange MultiFab data through amrex::MPMD. They are
# tested together below, not separately.
amrex_tutorial_add(MPMD_Main  SOURCES main.cpp  INPUTS inputs REQUIRES MPI NO_TEST)
amrex_tutorial_add(MPMD_Agent SOURCES agent.cpp INPUTS inputs REQUIRES MPI NO_TEST)

if(NOT TARGET Coupling_MPMD_Main)
    return()
endif()

file(COPY agent.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_package(MPI REQUIRED COMPONENTS CXX)
set(_main ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:Coupling_MPMD_Main> inputs)
set(_agent_args inputs)
if(AMReX_GPU_BACKEND STREQUAL "NONE")
    list(APPEND _main amrex.fpe_trap_invalid=1 amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1)
else()
    # Four ranks share the GPU, so start with small memory arenas.
    list(APPEND _main amrex.the_arena_init_size=268435456)
    list(APPEND _agent_args amrex.the_arena_init_size=268435456)
endif()

# mpiexec -n 2 MPMD_Main inputs : -n 2 MPMD_Agent inputs
add_test(NAME Coupling_MPMD
    COMMAND ${_main} : ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:Coupling_MPMD_Agent> ${_agent_args}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# mpiexec -n 2 MPMD_Main inputs : -n 2 python3 agent.py inputs
# The agent exchanges raw data, so pyAMReX and AMReX must have the same precision.
if(AMReXTutorials_PYTHON)
    if(AMReX_PRECISION STREQUAL "DOUBLE")
        find_package(Python3 REQUIRED COMPONENTS Interpreter)
        add_test(NAME Coupling_MPMD_Python
            COMMAND ${_main} : ${MPIEXEC_NUMPROC_FLAG} 2 ${Python3_EXECUTABLE} agent.py inputs
                    dim=${AMReXTutorials_DIM}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    else()
        message(STATUS "Skipping Coupling_MPMD_Python: pyAMReX uses double precision")
    endif()
endif()
