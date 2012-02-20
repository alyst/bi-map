#include <cemm/containers/dynamic_bitset_foreach.h>
#include <cemm/math/Distributions.h>

#include "cemm/bimap/ChessboardBiclusteringPriorEval.h"

namespace cemm { namespace bimap {

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
    for ( ChessboardBiclustering::const_block_iterator cluIt = clustering().begin(); cluIt != clustering().end(); ++cluIt ) {
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

ChessboardBiclusteringPriorEval::BlockEnablementPrior::BlockEnablementPrior(
    const ChessboardBiclusteringPriors&            priors1,
    const ChessboardBiclusteringDerivedPriors&     priors2,
    signal_t signal
) : lppEnabled( priors1.blockEnablementPrior( true ) + priors2.signalPrior( signal ) )
  , lppDisabled( priors1.blockEnablementPrior( false ) )
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

log_prob_t ChessboardBiclusteringPriorEval::blockLPP(
    object_clundex_t    objCluIx,
    probe_clundex_t     probeCluIx
) const {
    return ( blockLPP( *_clustering.findBlock( objCluIx, probeCluIx, false ) ) );
}

log_prob_t ChessboardBiclusteringPriorEval::blockLPP( const ChessboardBiclustering::block_proxy& cc ) const
{
    if ( cc.isEnabled() ) {
        GaussianDistribution signalPrior = clustering().derivedPriors().signalPrior;
        signal_t signal = cc.signal();
        if ( is_unset( signal ) ) THROW_RUNTIME_ERROR( "lpp(): signal for ("
            << cc.objectsClusterIndex() << ", " << cc.probesClusterIndex()
            << ") not set" );
        return ( signalPrior( signal ) + blockEnablementPrior()( true ) );
    }
    else {
        return ( blockEnablementPrior()( false ) );
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
        for ( ChessboardBiclustering::const_block_iterator ccIt = clustering().begin(); ccIt != clustering().end(); ++ccIt ) {
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
        // block probe prior
        res += blockEnablementPrior()( true ) * enabledCnt
            + blockEnablementPrior()( false ) * ( clustering().objectsClusters().size() * clustering().probesClusters().size() - enabledCnt );
    }

    BOOST_ASSERT( res < 0 );
    return ( res );
}

} }