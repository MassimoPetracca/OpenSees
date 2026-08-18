#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        ASDTimoshenkoBeam2d.cpp
        ASDTimoshenkoBeam3d.cpp
    PUBLIC
        ASDTimoshenkoBeam2d.h
        ASDTimoshenkoBeam3d.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
