set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON)
find_package( Boost REQUIRED COMPONENTS thread serialization iostreams filesystem program_options mpi )

find_package( MPI )

include_directories( ${Boost_INCLUDE_DIRS} )
add_definitions( ${MPI_COMPILE_FLAGS} )

set( BIMAPSamplerMPI_Sources
     "BIMAP-sampler-mpi.cpp"
)

if ( MPI_FOUND )
   include_directories( ${MPI_INCLUDE_PATH} )
   add_definitions( "-DBOOST_MPI_HOMOGENEOUS" ) # warning: would fail on non-homogeneous cluster

   add_executable( BIMAP-sampler-mpi ${BIMAPSamplerMPI_Sources} )
        add_dependencies( BIMAP-sampler-mpi BIMAP )
        target_link_libraries( BIMAP-sampler-mpi BIMAP ${Boost_LIBRARIES}
${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY} )
        set_target_properties( BIMAP-sampler-mpi PROPERTIES LINK_FLAGS
"${MPI_LINK_FLAGS}" )

endif( MPI_FOUND )
