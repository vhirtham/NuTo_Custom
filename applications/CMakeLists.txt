create_symlink("ARAMIS3d.msh")
create_symlink("ARAMIS3d_coarse.msh")

addNuToApp(ShrinkageAramis3d
    structures/MultiPhysicsStructureNew.cpp
    CustomIntegrands
    Tools
    )

