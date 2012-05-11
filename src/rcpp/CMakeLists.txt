find_path( R_INCLUDE_DIR "Rinternals.h" 
           PATHS "/usr/include" "/usr/share" "/usr/local/share"
                 "/fast/opt/include" "/fast/opt/lib64"
           PATH_SUFFIXES "R/include" "R" )
find_path( Rcpp_INCLUDE_DIR "Rcpp.h" 
           PATHS "/usr/lib" "/usr/local/lib" "/usr/lib64" "/fast/opt/lib64"
           PATH_SUFFIXES "R/site-library/Rcpp/include" "R/library/Rcpp/include" )
find_library( Rcpp_LIB "Rcpp" 
              PATHS
              "/home/astukalov/R/x86_64-unknown-linux-gnu-library/2.11"
              "/usr/lib" "/usr/local/lib" "/usr/lib64" "/fast/opt/lib64"
              PATH_SUFFIXES "R/site-library/Rcpp/lib" "R/library/Rcpp/lib" "Rcpp/lib" )
find_library( R_LIB "R"
              PATHS "/usr/lib" "/usr/local/lib" "/usr/lib64" "/fast/opt/lib64"
              PATH_SUFFIXES "R/lib" )

if ( R_INCLUDE_DIR )
    include_directories( ${R_INCLUDE_DIR} )
    add_definitions( "-DR_USE_C99_IN_CXX" ) # for Rvprintf() in R 2.13+
else ( R_INCLUDE_DIR )
    message( "R-project headers not found" )
endif( R_INCLUDE_DIR )

if ( Rcpp_INCLUDE_DIR )
    include_directories( ${Rcpp_INCLUDE_DIR} )
else ( Rcpp_INCLUDE_DIR )
    message( "Rcpp headers not found" )
endif( Rcpp_INCLUDE_DIR )

set( RBIMAP_Sources
     "RBIMAP.cpp"
)

add_library( RBIMAP SHARED ${RBIMAP_Sources} )
    add_dependencies( RBIMAP BIMAP math eesampler mcmc )
    target_link_libraries( RBIMAP BIMAP math ${R_LIB} ${Rcpp_LIB} )

set( RBIMAPMath_Sources
     "ExtraStats.cpp"
     "PartitionMatching.cpp"
)

add_library( RBIMAP-math SHARED ${RBIMAPMath_Sources} )
    add_dependencies( RBIMAP-math math )
    target_link_libraries( RBIMAP-math math ${GSL_LIBRARIES} ${R_LIB} ${Rcpp_LIB} )
