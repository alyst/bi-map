find_package( Boost REQUIRED COMPONENTS mpi )
find_package( MPI )
find_package( GSL )
find_package( GTest REQUIRED )

include_directories( "../../src" ${Boost_INCLUDE_DIR} ${GTEST_INCLUDE_DIRS} )

if ( MPI_FOUND )
    include_directories( ${MPI_INCLUDE_PATH} )

set( MPI_EESamplerTest_Sources
     "MPI_EESamplerTest.cpp"
)

add_executable( MPI_EESamplerTest ${MPI_EESamplerTest_Sources} )
    add_dependencies( MPI_EESamplerTest BIMAP )
    target_link_libraries( MPI_EESamplerTest BIMAP ${Boost_LIBRARIES} ${GTEST_BOTH_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} )
    set_target_properties( MPI_EESamplerTest PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}" )

endif( MPI_FOUND )

