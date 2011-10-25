#include "dynamic_bitset_utils.h"
#include "math/Distributions.h"

#include "ChessboardBiclusteringPriorEval.h"

ChessboardBiclusteringPriorEval::ChessboardBiclusteringPriorEval(
    const OPAData&                  data,
    const ChessboardBiclusteringPriors&    priors,
    const ChessboardBiclustering&          clustering
) : _data( data )
  , _priors( priors ), _clustering( clustering )
{
}

size_t ChessboardBiclusteringPriorEval::countMeasurements( bool signal ) const
{
    size_t mesCnt = 0;
    for ( ChessboardBiclustering::const_cross_cluster_iterator cluIt = clustering().begin(); cluIt != clustering().end(); ++cluIt ) {
        if ( cluIt->isEnabled() ) {
            mesCnt += data().assaysCount( clustering().probesCluster( cluIt->probesClusterIndex() ).items() )
                    * clustering().objectsCluster( cluIt->objectsClusterIndex() ).size();
        }
    }
    if ( !signal ) {
        mesCnt = data().assaysCount() * data().objectsCount() - mesCnt;
    }
    return ( mesCnt );
}

log_prob_t ChessboardBiclusteringPriorEval::noiseParamsLPP() const
{
    return ( noisePrior()( clustering().noiseParams().successRate ) );
}

struct ObjectsClustersSizes {
    const ChessboardBiclustering& cc;
    ObjectsClustersSizes( const ChessboardBiclustering& cc )
    : cc( cc )
    {}

    size_t size() const {
        return ( cc.objectsClusters().size() );
    }
    size_t operator[]( object_clundex_t cluIx ) const {
        return ( cc.objectsCluster( cluIx ).size() );
    }
};

log_prob_t ChessboardBiclusteringPriorEval::objectsClusteringLPP() const
{
    return ( priors().objectClustering.lnP( ObjectsClustersSizes( clustering() ) ) );
}

struct ProbesClustersSizes {
    const ChessboardBiclustering& cc;
    ProbesClustersSizes( const ChessboardBiclustering& cc )
    : cc( cc )
    {}

    size_t size() const {
        return ( cc.probesClusters().size() );
    }
    size_t operator[]( probe_clundex_t cluIx ) const {
        return ( cc.probesCluster( cluIx ).size() );
    }
};

ChessboardBiclusteringPriorEval::CrossClusterEnablementPrior::CrossClusterEnablementPrior(
    const ChessboardBiclusteringPriors&            priors1,
    const ChessboardBiclusteringDerivedPriors&     priors2,
    signal_t signal
) : lppEnabled( priors1.crossClusterEnablementPrior( true ) + priors2.signalPrior( signal ) )
  , lppDisabled( priors1.crossClusterEnablementPrior( false ) )
{
}

log_prob_t ChessboardBiclusteringPriorEval::probesClusteringLPP() const
{
    return ( priors().probeClustering.lnP( ProbesClustersSizes( clustering() ) ) );
}

log_prob_t ChessboardBiclusteringPriorEval::objectsMultiplesLPP() const
{
    double res = 0;
    for ( object_clundex_t objCluIx = 0; objCluIx < clustering().objectsClusters().size(); objCluIx++ ) {
        const ObjectsCluster& objClu = clustering().objectsCluster( objCluIx );
        GeometricDistribution  prior = objectMultiplePrior( objClu.size() );
        for ( object_set_t::const_iterator objIt = objClu.items().begin(); objIt != objClu.items().end(); ++objIt ) {
            res +=  prior( clustering().objectMultiple( *objIt ) );
        }
    }
    return ( res );
}

log_prob_t ChessboardBiclusteringPriorEval::crossClusterLPP(
    object_clundex_t    objCluIx,
    probe_clundex_t     probeCluIx
) const {
    return ( crossClusterLPP( *_clustering.findCrossCluster( objCluIx, probeCluIx, false ) ) );
}

log_prob_t ChessboardBiclusteringPriorEval::crossClusterLPP( const ChessboardBiclustering::cross_cluster_proxy& cc ) const
{
    if ( cc.isEnabled() ) {
        GaussianDistribution signalPrior = clustering().derivedPriors().signalPrior;
        signal_t signal = cc.signal();
        if ( is_unset( signal ) ) THROW_RUNTIME_ERROR( "lpp(): signal for ("
            << cc.objectsClusterIndex() << ", " << cc.probesClusterIndex()
            << ") not set" );
        return ( signalPrior( signal ) + crossClusterEnablementPrior()( true ) );
    }
    else {
        return ( crossClusterEnablementPrior()( false ) );
    }
}

log_prob_t ChessboardBiclusteringPriorEval::lpp() const
{
    double      res = objectsClusteringLPP() 
                    + probesClusteringLPP() 
                    + objectsMultiplesLPP()
                    + noiseParamsLPP();
    // chessboard biclusterings priors
    {
        size_t enabledCnt = 0;
        GaussianDistribution signalPrior = clustering().derivedPriors().signalPrior;
        for ( ChessboardBiclustering::const_cross_cluster_iterator ccIt = clustering().begin(); ccIt != clustering().end(); ++ccIt ) {
            if ( ccIt->isEnabled() ) {
                enabledCnt++;
                signal_t signal = ccIt->signal();
                if ( is_unset( signal ) ) THROW_RUNTIME_ERROR( "lpp(): signal for ("
                    << ccIt->objectsClusterIndex() << ", "
                    << ccIt->probesClusterIndex()
                    << ") is not set" );
                res += signalPrior( signal );
            }
        }
        // cross-cluster probe prior
        res += crossClusterEnablementPrior()( true ) * enabledCnt
            + crossClusterEnablementPrior()( false ) * ( clustering().objectsClusters().size() * clustering().probesClusters().size() - enabledCnt );
    }

    BOOST_ASSERT( res < 0 );
    return ( res );
}
