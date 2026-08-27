#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        ASDHinge.cpp
    PUBLIC
        ASDHinge.h
        ASDHingeKernel.h
        ASDHingeTransformation.h
        ASDHingeCorotationalTransformation.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
