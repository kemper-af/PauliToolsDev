add_executable( iterative_cartan iterative_cartan.cc )
target_link_libraries( iterative_cartan PUBLIC 
    paulitools
    dla 
    models
    )


target_include_directories( iterative_cartan PRIVATE
    /opt/homebrew/include
    )


find_library( DLIB dlib )
