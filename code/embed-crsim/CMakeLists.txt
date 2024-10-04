project(AOS)
cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)

option(PARALLEL_AOS "Enable parallel AOS" OFF)

add_executable(${PROJECT_NAME}
    crsim_luts_mod.F
    crsim_mod.F
    extract_apar_wrf.F
    kwm_date_utilities.F
    module_access_wrf.F
    module_aircraft.F
    module_cfradial_output.F
    module_configuration.F
    module_crsim_wrapper.F
    module_external_attitude.F
    module_geometry.F
    module_llxy.F
    module_quicksort.F
    module_scanning.F
    module_timing_utilities.F
    phys_param_mod.F
    ReadConfParameters.F
    wrf_var_mod.F
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    Fortran_FORMAT FREE
)

target_compile_options(${PROJECT_NAME} PRIVATE $<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-none>)

# MacOS specific flags to prevent an "Illegal instruction" error from occurring.
# TODO: Investigate why this is necessary.
set_property(SOURCE "extract_apar_wrf.F" PROPERTY COMPILE_FLAGS $<$<PLATFORM_ID:Darwin>:"-Og">)
set_property(SOURCE "extract_apar_wrf.F" PROPERTY COMPILE_FLAGS $<$<PLATFORM_ID:Darwin>:"-ftree-sra">)

find_package(NetCDF REQUIRED COMPONENTS Fortran)
target_link_libraries(${PROJECT_NAME} PRIVATE NetCDF::NetCDF_Fortran NetCDF::NetCDF_C)

message(STATUS "Parallel AOS: ${PARALLEL_AOS}")
if(PARALLEL_AOS)
    if(NOT MPI_HOME AND NOT ENV{MPI_HOME})
        find_path(MPI_HOME "bin/mpiexec" PATH_SUFFIXES mpich openmpi PATHS /usr/lib /usr/lib64)
    endif()
    find_package(MPI
        REQUIRED
        COMPONENTS
            Fortran
    )
    target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_Fortran)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _PARALLEL_)
endif()