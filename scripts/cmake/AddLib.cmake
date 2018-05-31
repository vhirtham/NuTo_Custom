#[==[
This function does all the necessary steps to set up a unit test
#]==]

function(addLib LibName)
    ### Add prefix "Test_" to avoid target naming conflicts (for example with benchmarks or other GDL targets)
    #set(LibName "Test_${LibName}")

    ### Separate sources from libs
    if(ARGN)
        foreach(input ${ARGN})
            string(FIND "${input}" ".cpp" IsCpp)
            if(${IsCpp} EQUAL "-1")
                set(AdditionalLibs
                    "${AdditionalLibs};${input}")
            else()
                set(Sources
                    "${Sources};${NUTO_APPS_ROOT_DIR}/src/${input}")
            endif()
        endforeach()
    endif()

    ### Create executable
    add_library(${LibName}
        ${Sources}
        )

    ### Link necessary libs
    add_dependencies(${LibName} NuTo_ext)
    target_link_Libraries(${LibName}
        NuTo
        ${AdditionalLibs}
        )

    ### Add source directory
    target_include_directories(${LibName}
        PUBLIC
            ${NUTO_APPS_ROOT_DIR}/src)


endfunction()
