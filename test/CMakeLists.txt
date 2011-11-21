find_package( GTest REQUIRED )
find_package( Boost REQUIRED COMPONENTS filesystem program_options )

include_directories( "../src" ${Boost_INCLUDE_DIR} ${GTEST_INCLUDE_DIRS} )

add_definitions( -DGTEST_HAS_TR1_TUPLE=0 )
add_definitions( -DGTEST_USE_OWN_TR1_TUPLE=0 )
#add_definitions( -DBOOST_HAS_TR1_TUPLE )

set( BIMAPTest_Sources
     "main.cpp"
     "OPADataTest.cpp"
     "ClusteringTest.cpp"
     "DistributionsTest.cpp"
     "MCMCTest.cpp"
     "ChessboardBiclusteringTest.cpp"
     "PitmanYorTest.cpp"
     "EESamplerTest.cpp"
     "BIMAPTest.cpp"
     "BIMAPWalkTest.cpp"
     "SemindexedMap.cpp"
     "ChessboardBiclusteringEvalTest.cpp"
     "array2dTest.cpp"
     "symmetric_array2dTest.cpp"
     "OrderedDistanceMatrixTest.cpp"
)

add_executable( BIMAP-test ${BIMAPTest_Sources} )
    add_dependencies( BIMAP-test BIMAP )
    target_link_libraries( BIMAP-test BIMAP ${GTEST_BOTH_LIBRARIES} ${Boost_LIBRARIES} )

add_test( BIMAP-test "${EXECUTABLE_OUTPUT_PATH}/BIMAP-test" )

add_subdirectory( mpisampler )

