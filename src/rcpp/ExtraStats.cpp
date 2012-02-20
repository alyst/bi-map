#include <cemm/bimap/BasicTypedefs.h>

#include <stdarg.h>
#include <cmath>

#include <R_ext/Rdynload.h>
#include <cemm/RUtils.h>

#ifndef NDEBUG
//#define DYNLOAD_DEBUG
#endif

#include <cemm/math/LagrangianPoisson.h>
#include <cemm/math/PitmanYorProcess.h>
#include <cemm/math/GenericDiscreteDistribution.h>

namespace cemm { namespace bimap {

gsl_rng* rndNumGen = NULL;

gsl_rng* getRndNumGen()
{
    if ( rndNumGen == NULL ) {
        Rprintf( "Initializing GSL random number generator...\n" );
        gsl_rng_env_setup();
        rndNumGen = gsl_rng_alloc( gsl_rng_default ); /// @todo: free on detach
    }
    return ( rndNumGen );
}

RcppExport SEXP LagrangianPoissonLnPdf(
    SEXP    xExp,
    SEXP    lnRateExp,
    SEXP    shapeExp
){
    BEGIN_RCPP
    LagrangianPoissonDistribution lp( Rf_asReal( lnRateExp ), Rf_asReal( shapeExp ) );
    return ( Rf_ScalarReal( lp.lnPdf( Rf_asReal( xExp ) ) ) );
    END_RCPP
}

RcppExport SEXP LagrangianPoissonLnCdfP(
    SEXP    xExp,
    SEXP    lnRateExp,
    SEXP    shapeExp
){
    BEGIN_RCPP
    LagrangianPoissonDistribution lp( Rf_asReal( lnRateExp ), Rf_asReal( shapeExp ) );
    return ( Rf_ScalarReal( lp.lnCdf_P( Rf_asReal( xExp ) ) ) );
    END_RCPP
}

RcppExport SEXP LagrangianPoissonLnCdfQ(
    SEXP    xExp,
    SEXP    lnRateExp,
    SEXP    shapeExp
){
    BEGIN_RCPP
    LagrangianPoissonDistribution lp( Rf_asReal( lnRateExp ), Rf_asReal( shapeExp ) );
    return ( Rf_ScalarReal( lp.lnCdf_Q( Rf_asReal( xExp ) ) ) );
    END_RCPP
}

RcppExport SEXP LagrangianPoissonRandom(
    SEXP    sizeExp,
    SEXP    lnRateExp,
    SEXP    shapeExp
){
    BEGIN_RCPP
    LagrangianPoissonDistribution lp( Rf_asReal( lnRateExp ), Rf_asReal( shapeExp ) );
    GenericDiscreteDistribution gdd( lp, 1E-3 );
    size_t  size = Rf_asInteger( sizeExp );
    Rcpp::IntegerVector res( size );
    for ( size_t i = 0; i < size; i++ ) {
        res[i] = gdd.random( getRndNumGen() );
    }
    return ( res );
    END_RCPP
}

RcppExport SEXP PitmanYorRandomSample(
    SEXP    nElementsExp,
    SEXP    concentrationExp,
    SEXP    discountExp
){
    BEGIN_RCPP
    size_t nElements = Rf_asInteger( nElementsExp );
    double concentration = Rf_asReal( concentrationExp );
    double discount = Rf_asReal( discountExp );
    PitmanYorProcess pyp( concentration, discount );
    PitmanYorSample pys = pyp.random( getRndNumGen(), nElements );
    Rcpp::IntegerVector clusterSizes( pys.clustersCount() );
    for ( size_t i = 0; i < pys.clustersCount(); i++ ) {
        clusterSizes[ i ] = pys.cluster( i ).size();
    }
    return ( clusterSizes );
    END_RCPP
}

#if 0
RcppExport SEXP CreateMemLeak(
){
    BEGIN_RCPP
    int*  mem = new int[ 10000000 ];
    mem[0] = 1;
    //sleep( 10 );
    return ( Rf_ScalarInteger( mem[0] ) );
    END_RCPP
}

RcppExport SEXP DontCreateMemLeak(
){
    BEGIN_RCPP
    int*  mem = new int[ 10000000 ];
    delete [] mem;
    int r = 4;
    //sleep( 10 );
    return ( R_NilValue );
    END_RCPP
}

RcppExport SEXP DoNothing(
){
    BEGIN_RCPP
    return ( R_NilValue );
    END_RCPP
}
#endif

} }