create_symlink("ARAMIS3d.msh")
create_symlink("ARAMIS3d_coarse.msh")

addNuToApp(ShrinkageAramis3d
    structures/MultiPhysicsStructureNew.cpp
    CustomIntegrands
    Tools
    )

create_symlink("ARAMIS2d.msh")
addNuToApp(GradientDamageAramis2d
    structures/MultiPhysicsStructureNew.cpp
    CustomIntegrands
    Tools
    )


create_symlink("Cylinder.msh")
addNuToApp(GradientDamageCylinder3d
    structures/MultiPhysicsStructureNew.cpp
    CustomIntegrands
    Tools
    )


addNuToApp(MoistureDependentCylinderCompression3d
    structures/MultiPhysicsStructureNew.cpp
    CustomIntegrands
    Tools
    )

addNuToApp(MoistureTransportCylinder3d
    structures/MultiPhysicsStructureNew.cpp
    CustomIntegrands
    Tools
    )

addNuToApp(RandomTests)
