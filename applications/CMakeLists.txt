create_symlink("ARAMIS3d.msh")
create_symlink("ARAMIS3d_coarse.msh")

addNuToApp(MultiPhysics3d
    structures/MultiPhysicsStructureNew.cpp
    Tools
    )

