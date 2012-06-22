set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON)
find_package( Boost REQUIRED COMPONENTS serialization iostreams filesystem program_options system log )

find_package( ZLIB )
find_package( BZip2 )
find_package( Threads )

find_package( GSL REQUIRED )

set( BIMAP_Sources
     "OPAData.cpp"
     "OPADataImportCSV.cpp"
     "ChessboardBiclustering.cpp"
     "DataSignalNoiseCache.cpp"
     "CellSignalLLHMaximizer.cpp"
     "CellsLLHEval.cpp"
     "CoOccurrenceGraph.cpp"
     "ChessboardBiclusteringCrossover.cpp"
     "ChessboardBiclusteringPriorEval.cpp"
     "ChessboardBiclusteringLLHEval.cpp"
     "ChessboardBiclusteringStructureLLHEval.cpp"
     "ChessboardBiclusteringGibbsHelper.cpp"
     "ChessboardBiclusteringGibbsSampler.cpp"
     "ChessboardBiclusteringFit.cpp"
     "ChessboardBiclusteringsIndexing.cpp"
     "ChessboardBiclusteringsPDFEval.cpp"
     "ObjectsPartition.cpp"
     "ProbesPartition.cpp"
     "PrecomputedData.cpp"
     "BIMAPSampler.cpp"
     "BIMAPWalk.cpp"
     "BIMAPResultsSerialize.cpp"
     "ParametersReader.cpp"
     "misc_utils.cpp"
)

set( BIMAP_Headers
     "BasicTypedefs.h"
     "ObjectSet.h"
     "ProbeSet.h"
     "OPAData.h"
     "ChessboardBiclustering.h"
     "ChessboardBiclusteringEval.h"
     "ChessboardBiclusteringGibbsHelper.h"
     "ChessboardBiclusteringGibbsSampler.h"
     "ChessboardBiclusteringsIndexing.h"
     "BIMAPSampler.h"
     "BIMAPWalk.h"
     "math/StatsMath.h"
)

include_directories( ${Boost_INCLUDE_DIRS} )

add_library( BIMAP SHARED ${BIMAP_Sources} )
    add_dependencies( BIMAP math eesampler mcmc )
    target_link_libraries( BIMAP math eesampler
                           ${GSL_LIBRARIES}
                           ${Boost_LIBRARIES}
                           ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES} 
                           ${CMAKE_THREAD_LIBS_INIT} )

set( BIMAPSampler_Sources
     "BIMAP-sampler.cpp"
)

add_executable( BIMAP-sampler ${BIMAPSampler_Sources} )
    add_dependencies( BIMAP-sampler BIMAP )
    target_link_libraries( BIMAP-sampler BIMAP
                           ${CMAKE_THREAD_LIBS_INIT}
                           ${Boost_LIBRARIES} )

add_subdirectory( mpisampler )

add_subdirectory( rcpp )
