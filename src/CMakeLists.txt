#add_library(CustomIntegrands
#    integrands/MoistureTransport.cpp
#    )


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
