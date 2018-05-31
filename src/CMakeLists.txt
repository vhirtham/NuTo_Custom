addLib(CustomIntegrands
    integrands/MoistureTransport.cpp
    integrands/MoistureTransportBoundary.cpp
    )


#target_link_libraries(CustomIntegrands
#    NuToBase
#    NuToMath
#    NuToMechanics
#    )

add_library(Tools
    tools/GNUPlot.cpp
    tools/PlotData.cpp
    tools/PlotEnum.cpp
    )
