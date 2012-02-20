cmake_minimum_required( VERSION 2.6 )

PROJECT( BIMAP )

find_package( MPI )

include_directories(
    "include"
    "common/include/"
    "math/include"
    "mcmc/include"
    "eesampler/include"
    "containers/include"
)

add_subdirectory( src )
# add_subdirectory( common )
add_subdirectory( math )
add_subdirectory( mcmc )
add_subdirectory( containers )
add_subdirectory( eesampler )

enable_testing()

add_subdirectory(test)
