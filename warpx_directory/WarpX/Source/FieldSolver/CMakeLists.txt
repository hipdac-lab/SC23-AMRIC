target_sources(WarpX
  PRIVATE
    ElectrostaticSolver.cpp
    WarpXPushFieldsEM.cpp
    WarpX_QED_Field_Pushers.cpp
)

add_subdirectory(FiniteDifferenceSolver)
add_subdirectory(MagnetostaticSolver)
if(WarpX_PSATD)
    add_subdirectory(SpectralSolver)
endif()
