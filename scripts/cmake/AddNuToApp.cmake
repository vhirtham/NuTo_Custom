#[==[
This function does all the necessary steps to set up a unit test
#]==]

function(addNuToApp AppName)
    ### Add prefix "Test_" to avoid target naming conflicts (for example with benchmarks or other GDL targets)
    #set(AppName "Test_${AppName}")

    ### Separate sources from libs
    if(ARGN)
        foreach(input ${ARGN})
            string(FIND "${input}" ".cpp" IsCpp)
            if(${IsCpp} EQUAL "-1")
                set(AdditionalLibs
                    "${AdditionalLibs};${input}")
            else()
                set(AdditionalSources
                    "${AdditionalSources};${NUTO_APPS_ROOT_DIR}/src/${input}")
            endif()
        endforeach()
    endif()

    ### Create executable
    add_executable(${AppName}
        "${AppName}.cpp"
        ${AdditionalSources}
        )

    ### Link necessary libs
    add_dependencies(${AppName} NuTo_ext)
    target_link_Libraries(${AppName}
        NuTo
        ${AdditionalLibs}
        )

    ### Add source directory
    target_include_directories(${AppName}
        PUBLIC
            ${NUTO_APPS_ROOT_DIR}/src)


endfunction()
