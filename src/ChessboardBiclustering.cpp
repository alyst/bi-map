#include "dynamic_bitset_utils.h"

#include <boost/format.hpp>

#include "statically_tracked.h"

#include "ChessboardBiclustering.h"

ENABLE_STATIC_TRACKING( ChessboardBiclusteringPriors )

VectorCache<signal_t> ObjectsClusterParams::SignalMapCache;
VectorCache<size_t> ObjectsClusterParams::MultipleMapCache;
VectorCache<signal_t> ProbesClusterParams::SignalMapCache;

ObjectsClusterParams::ObjectsClusterParams(
    const ChessboardBiclustering& clustering
)   : blocksMask( clustering.probesClusters().size() )
    , probeSignal( SignalMapCache.borrow( clustering.probesClusters().size(), unset() ) )
    , objectMultiple( MultipleMapCache.borrow( clustering.objectsCount(), 0 ) )
{
}

ObjectsClusterParams::ObjectsClusterParams(
    const ChessboardBiclustering::ObjectsClusterParamsProxy& paramsProxy
)   : blocksMask( paramsProxy.clus.probesClusters().size() )
    , probeSignal( SignalMapCache.borrow( paramsProxy.clus.probesClusters().size(), unset() ) )
    , objectMultiple( MultipleMapCache.borrow( paramsProxy.clus.objectsCount(), 0 ) )
{
    operator=( paramsProxy );
}

ObjectsClusterParams::ObjectsClusterParams(
    const ObjectsClusterParams& params
)  : blocksMask( params.blocksMask )
   , probeSignal( SignalMapCache.borrow( params.probeSignal ) )
   , objectMultiple( MultipleMapCache.borrow( params.objectMultiple ) )
{
}

ProbesClusterParams::ProbesClusterParams(
    const ChessboardBiclustering& clustering
) : blocksMask( clustering.objectsClusters().size() )
  , objectsSignal( SignalMapCache.borrow( clustering.objectsClusters().size(), unset() ) )
{
}

ProbesClusterParams::ProbesClusterParams(
    const ChessboardBiclustering::ProbesClusterParamsProxy& paramsProxy
) : blocksMask( paramsProxy.clus.objectsClusters().size() )
  , objectsSignal( SignalMapCache.borrow( paramsProxy.clus.objectsClusters().size(), unset() ) )
{
    operator=( paramsProxy );
}

ProbesClusterParams::ProbesClusterParams(
    const ProbesClusterParams& params
)  : blocksMask( params.blocksMask )
   , objectsSignal( SignalMapCache.borrow( params.objectsSignal ) )
{
}

ChessboardBiclustering::ChessboardBiclustering(
    size_t  objectsCount,
    size_t  probesCount
) : ChessboardBiclusteringData( ChessboardBiclusteringDerivedPriors(), signal_params_type(), noise_params_type() )
  , _objectToCluster( objectsCount, CLUSTER_NA )
  , _probeToCluster( probesCount, CLUSTER_NA )
  , _objectsClustersCleanupRequired( false )
  , _probesClustersCleanupRequired( false )
  , _objectMultiples( objectsCount, 1 )
  , _signals( objectsCount, probesCount, unset() )
{
}

ChessboardBiclustering::ChessboardBiclustering(
    const ChessboardBiclusteringDerivedPriors&     derivedPriors,
    const signal_params_type&               baselineSignal,
    const noise_params_type&                noiseParams,
    const PitmanYorSample&                  objectsClusters,
    const PitmanYorSample&                  probesClusters
) : ChessboardBiclusteringData( derivedPriors, baselineSignal, noiseParams )
  , _objectToCluster( objectsClusters.samplesCount(), CLUSTER_NA )
  , _probeToCluster( probesClusters.samplesCount(), CLUSTER_NA )
  , _objectsClustersCleanupRequired( false )
  , _probesClustersCleanupRequired( false )
  , _objectMultiples( objectsClusters.samplesCount(), 1 )
  , _signals( objectsClusters.samplesCount(), probesClusters.samplesCount(), unset() )
{
    for ( cluster_index_t i = 0; i < objectsClusters.clustersCount(); ++i ) {
        addObjectCluster( objectsClusters[ i ].begin(), objectsClusters[ i ].end() );
    }
    for ( cluster_index_t i = 0; i < probesClusters.clustersCount(); ++i ) {
        addProbeCluster( probesClusters[ i ].begin(), probesClusters[ i ].end() );
    }
}

void ChessboardBiclustering::setObjectMultiple(
    object_index_t  objIx,
    size_t          multiple
){
    if ( (int)objIx >= 0 && objIx < objectsCount() ) {
        if ( multiple == 0 ) throw std::range_error( "Object multiple cannot be 0" );
        if ( _objectMultiples[ objIx ] != multiple ) {
            _objectMultiples[ objIx ] = multiple;
            afterObjectMultipleChanged( objIx );
        }
    }
    else {
        throw std::range_error( "Incorrect object index" );
    }
}

void ChessboardBiclustering::setObjectCluster(
    object_index_t          objIx,
    const object_clundex_t  objCluIx
){
    const object_clundex_t oldCluIx = _objectToCluster[ objIx ];
    if ( oldCluIx != objCluIx ) {
        LOG_DEBUG2( "Putting object " << objIx << " to cluster " << objCluIx );
        if ( oldCluIx != CLUSTER_NA ) {
            ObjectsCluster& oldClu = _objectsClusters[ oldCluIx ];
            oldClu._items.erase( objIx );
            if ( oldClu._items.empty() ) {
                _objectsClustersCleanupRequired = true;
            } else {
                afterObjectsClusterChanged( oldCluIx );
            }
        }
        if ( objectsCluster( objCluIx ).size() > 0 ) {
            object_index_t newObjCluObjIx = objectsCluster( objCluIx ).representative();
            // copy signals
            _signals.section1( objIx ) = _signals.section1( newObjCluObjIx );
        }
        _objectsClusters[ objCluIx ]._items.insert( objIx );
        _objectToCluster[ objIx ] = objCluIx;

        afterObjectsClusterChanged( objCluIx );
    }
}

void ChessboardBiclustering::setProbeCluster(
    probe_index_t       probeIx,
    probe_clundex_t     probeCluIx
){
    const probe_clundex_t oldCluIx = _probeToCluster[ probeIx ];
    if ( oldCluIx != probeCluIx ) {
        LOG_DEBUG2( "Putting probe " << probeIx << " to cluster " << probeCluIx );
        if ( oldCluIx != CLUSTER_NA ) {
            ProbesCluster& oldClu = _probesClusters[ oldCluIx ];
            oldClu._items.set( probeIx, false );
            if ( oldClu._items.none() ) {
                _probesClustersCleanupRequired = true;
            } else {
                afterProbesClusterChanged( oldCluIx );
            }
        }
        if ( probesCluster( probeCluIx ).size() > 0 ) {
            probe_index_t newProbeCluProbeIx = probesCluster( probeCluIx ).representative();
            // copy signals
            _signals.section2( probeIx ) = _signals.section2( newProbeCluProbeIx );
        }
        _probesClusters[ probeCluIx ]._items.set( probeIx, true );
        _probeToCluster[ probeIx ] = probeCluIx;
        afterProbesClusterChanged( probeCluIx );
    }
}

/**
 *  Exchange given set of objects between 2 clusters.
 *  If second cluster index is OBJECT_NA, new cluster is created.
 */
object_clundex_t ChessboardBiclustering::exchangeObjects(
    object_clundex_t    clu1Ix,         /** i first cluster */
    object_clundex_t    clu2Ix,         /** i second cluster */
    const object_set_t& clu2Elements    /** i elements of updated clu2 cluster */
){
    object_clundex_t res = clu2Ix;
    // remember cluster signals
    std::vector<signal_t>   clu1Signals = _signals.section1( objectsCluster( clu1Ix ).representative() );
    std::vector<signal_t>   clu2Signals = clu2Ix != OBJECT_NA
                                        ? (std::vector<signal_t>)_signals.section1( objectsCluster( clu2Ix ).representative() )
                                        : clu1Signals;

    // move elements from clu1 to clu2
    for ( object_set_t::const_iterator it = clu2Elements.begin(); it != clu2Elements.end(); ++it ) {
        if ( clusterOfObject( *it ) == clu1Ix ) {
            if ( res == OBJECT_NA ) {
                res = addObjectCluster( *it );
                // copy block probes
                for ( probe_clundex_t i = 0; i < probesClusters().size(); i++ ) {
                    _blocksMask.set( res * _probesClusters.size() + i,
                                            _blocksMask.test( clu1Ix * _probesClusters.size() + i ) );
                }
            }
            else {
                // copy signals
                _signals.section1( *it ) = clu2Signals;
                _objectsClusters[ res ]._items.insert( *it );
                _objectsClusters[ clu1Ix ]._items.erase( *it );
                _objectToCluster[ *it ] = res;
            }
        }
    }
    // move elements from clu2 to clu1
    if ( clu2Ix != OBJECT_NA ) {
        // elements that are still in clu2, but should be moved to clu1
        object_set_t newClu1Items = objectsCluster( clu2Ix ).items();
        for ( object_set_t::const_iterator it = clu2Elements.begin(); it != clu2Elements.end(); ++it ) {
            newClu1Items.erase( *it );
        }
        // add newClu1Items to clu1
        for ( object_set_t::const_iterator it = newClu1Items.begin(); it != newClu1Items.end(); ++it ) {
            _objectsClusters[ clu1Ix ]._items.insert( *it );
            _objectsClusters[ clu2Ix ]._items.erase( *it );
            _objectToCluster[ *it ] = clu1Ix;
            // copy signals
            _signals.section1( *it ) = clu1Signals;
        }
        if ( _objectsClusters[ clu2Ix ]._items.empty() ) {
            _objectsClustersCleanupRequired = true;
        } else {
            afterObjectsClusterChanged( clu2Ix );
        }
    }
    if ( _objectsClusters[ clu1Ix ]._items.empty() ) {
        _objectsClustersCleanupRequired = true;
    } else {
        afterObjectsClusterChanged( clu1Ix );
    }
    BOOST_ASSERT( checkBlocks() );
    return ( res );
}

/**
 *  Exchange given set of elements between 2 clusters.
 *  If second cluster index is PROBE_NA, new cluster is created.
 */
probe_clundex_t ChessboardBiclustering::exchangeProbes(
    probe_clundex_t         clu1Ix, /** i first cluster */
    probe_clundex_t         clu2Ix, /** i second cluster, might be missing */
    const probe_bitset_t&   clu2Elements /** i elements of updated clu2 cluster */
){
    probe_clundex_t res = clu2Ix;
    // remember cluster signals
    std::vector<signal_t>   clu1Signals = _signals.section2( probesCluster( clu1Ix ).representative() );
    std::vector<signal_t>   clu2Signals = clu2Ix != PROBE_NA
                                        ? (std::vector<signal_t>)_signals.section2( probesCluster( clu2Ix ).representative() )
                                        : clu1Signals;

    // move elements from clu1 to clu2
    probe_bitset_t newClu2Items = clu2Elements &  probesCluster( clu1Ix ).items();
    foreach_bit( probe_index_t, probeIx, newClu2Items ) {
        if ( res == PROBE_NA ) {
            // create new cluster
            res = addProbeCluster( probeIx );
            // copy block probes
            for ( object_clundex_t i = 0; i < objectsClusters().size(); i++ ) {
                _blocksMask.set( i * _probesClusters.size() + res,
                                        _blocksMask.test( i * _probesClusters.size() + clu1Ix ) );
            }
        }
        else {
            // copy signals
            _signals.section2( probeIx ) = clu2Signals;
            _probesClusters[ res ]._items.set( probeIx, true );
            _probesClusters[ clu1Ix ]._items.set( probeIx, false );
            _probeToCluster[ probeIx ] = res;
        }
    }
    // move elements from clu2 to clu1
    if ( clu2Ix != PROBE_NA ) {
        probe_bitset_t newClu1Items = probesCluster( clu2Ix ).items();
        newClu1Items -= clu2Elements;
        foreach_bit( probe_index_t, probeIx, newClu1Items ) {
            _probesClusters[ clu1Ix ]._items.set( probeIx, true );
            _probesClusters[ clu2Ix ]._items.set( probeIx, false );
            _probeToCluster[ probeIx ] = clu1Ix;
            // copy signals
            _signals.section2( probeIx ) = clu1Signals;
        }
        if ( _probesClusters[ clu2Ix ]._items.none() ) {
            _probesClustersCleanupRequired = true;
        } else {
            afterProbesClusterChanged( clu2Ix );
        }
    }
    if ( _probesClusters[ clu1Ix ]._items.none() ) {
        _probesClustersCleanupRequired = true;
    } else {
        afterProbesClusterChanged( clu1Ix );
    }
    BOOST_ASSERT( checkBlocks() );
    return ( res );
}

ChessboardBiclustering::objects_cluster_index_remap ChessboardBiclustering::cleanupObjectsClusters()
{
    objects_cluster_index_remap old2newIx( _objectsClusters.size() );

    object_clundex_t oldIx = 0;
    for ( object_clundex_t newIx = 0; newIx < _objectsClusters.size(); oldIx++ ) {
        if ( _objectsClusters[ newIx ].items().empty() ) {
            LOG_DEBUG2( "Deleting objects cluster " << oldIx << "(" << newIx << ")" );
            beforeObjectsClusterRemoved( newIx );
            _objectsClusters.erase( _objectsClusters.begin() + newIx );
            old2newIx[ oldIx ] = CLUSTER_NA;
        }
        else {
            if ( oldIx != newIx ) {
                for ( object_set_t::const_iterator objIt = _objectsClusters[ newIx ]._items.begin(); 
                      objIt != _objectsClusters[ newIx ]._items.end(); ++objIt
                ){
                    _objectToCluster[ *objIt ] = newIx;
                }
            }
            old2newIx[ oldIx ] = newIx++;
        }
    }
    _objectsClustersCleanupRequired = false;
    return ( old2newIx );
}

ChessboardBiclustering::probes_cluster_index_remap ChessboardBiclustering::cleanupProbesClusters()
{
    probes_cluster_index_remap old2newIx( _probesClusters.size() );

    probe_clundex_t oldIx = 0;
    for ( probe_clundex_t newIx = 0; newIx < _probesClusters.size(); oldIx++ ) {
        if ( _probesClusters[ newIx ].items().none() ) {
            LOG_DEBUG2( "Deleting probes cluster " << oldIx << "(" << newIx << ")" );
            beforeProbesClusterRemoved( newIx );
            _probesClusters.erase( _probesClusters.begin() + newIx );
            old2newIx[ oldIx ] = CLUSTER_NA;
        }
        else {
            if ( oldIx != newIx ) {
                foreach_bit( probe_index_t, probeIx, _probesClusters[ newIx ]._items ) {
                    _probeToCluster[ probeIx ] = newIx;
                }
            }
            old2newIx[ oldIx ] = newIx++;
        }
    }
    _probesClustersCleanupRequired = false;
    return ( old2newIx );
}

/**
 *  Resizes mask to new dimensions, preserving data.
 *  @note
 *  Any existing BlockProxy instances after resizing would be invalid.
 */
void ChessboardBiclustering::resizeBlockMask(
    size_t objectClustersCount, 
    size_t probeClustersCount
){
    blocks_mask_t maskReindexed( objectClustersCount * probeClustersCount );
    for ( const_block_iterator cluIt = begin(); cluIt != end(); ++cluIt ) {
        object_clundex_t  objCluIx = cluIt->objectsClusterIndex();
        probe_clundex_t   probeCluIx = cluIt->probesClusterIndex();
        if ( objCluIx < objectClustersCount && probeCluIx < probeClustersCount ) {
            maskReindexed.set( objCluIx * probeClustersCount + probeCluIx );
        }
    }
    _blocksMask.swap( maskReindexed );
}

void ChessboardBiclustering::cleanupClusters()
{
    if ( !_objectsClustersCleanupRequired && !_probesClustersCleanupRequired ) return;
    // remember old dimensions
    size_t oldProbeClusters = probesClusters().size();

    objects_cluster_index_remap objRemap;
    if ( _objectsClustersCleanupRequired ) {
        LOG_DEBUG2( "Cleaning objects clusters..." );
        objRemap = cleanupObjectsClusters();
    }
    probes_cluster_index_remap probeRemap;
    if ( _probesClustersCleanupRequired ) {
        LOG_DEBUG2( "Cleaning probes clusters..." );
        probeRemap = cleanupProbesClusters();
    }
    // now the counts of object and probe clusters are updated
    blocks_mask_t maskReindexed( objectsClusters().size() * probesClusters().size() );
    // don't use BlockIterator, since dimensions are updated, but not the map itself
    foreach_bit( size_t, ccIndex, _blocksMask ) {
        object_clundex_t  oldObjIx = ccIndex / oldProbeClusters;
        probe_clundex_t   oldProbeIx = ccIndex % oldProbeClusters;
        object_clundex_t  newObjIx = !objRemap.empty() ? objRemap[ oldObjIx ] : oldObjIx;
        probe_clundex_t   newProbeIx = !probeRemap.empty() ? probeRemap[ oldProbeIx ] : oldProbeIx;
        if ( newObjIx != CLUSTER_NA && newProbeIx != CLUSTER_NA ) {
            LOG_DEBUG2( "Reindexing block (" << oldObjIx << "," << oldProbeIx 
                        << ") to (" << newObjIx << "," << newProbeIx << ")" );
            maskReindexed.set( newObjIx * probesClusters().size() + newProbeIx );
        }
    }
    _blocksMask.swap( maskReindexed );
    BOOST_ASSERT( !_objectsClustersCleanupRequired && !_probesClustersCleanupRequired );
}

/**
 *  Get 2D bitmap (probes x objects) of block enablement.
 */
blocks_mask_t ChessboardBiclustering::blocksMask(
    bool objectsMajor    /** i if true (default), rows are objects, probes are columns */
) const {
    if ( objectsMajor ) {
        return ( _blocksMask );
    } else {
        // transpose mask
        blocks_mask_t res( objectsClusters().size() * probesClusters().size() );
        // fill signals and mask
        for ( const_block_iterator ccit = begin(); ccit != end(); ++ccit ) {
            const BlockProxy& cc = *ccit;
            size_t  blockIx = objectsClusters().size() * cc.probesClusterIndex() + cc.objectsClusterIndex();
            res.set( blockIx );
        }
        return ( res );
    }
}

/**
 *  Sets block on or off.
 *  @note
 *  Block signal not set.
 */
ChessboardBiclustering::block_iterator ChessboardBiclustering::setBlock(
    object_clundex_t    objCluIx,
    probe_clundex_t     probeCluIx,
    bool                enable
){
    block_iterator cluIt = findBlock( objCluIx, probeCluIx, false );
    if ( cluIt->isEnabled() != enable ) {
        cluIt->setEnabled( enable );
        afterBlockFlipped( cluIt->objectsClusterIndex(), cluIt->probesClusterIndex() );
    }
    return ( cluIt );
}

void ChessboardBiclustering::setBlockSignal(
    object_clundex_t    objCluIx,
    probe_clundex_t     probeCluIx,
    signal_t            signal
){
    const object_set_t& objs = objectsCluster( objCluIx ).items();
    const probe_bitset_t& probes = probesCluster( probeCluIx ).items();
    for ( object_set_t::const_iterator oit = objs.begin(); oit != objs.end(); ++oit ) {
        object_index_t objIx = *oit;
        foreach_bit( probe_index_t, probeIx, probes ) {
            _signals( objIx, probeIx ) = signal;
        }
    }
    afterSignalChanged( objCluIx, probeCluIx );
}

void ChessboardBiclustering::updateObjectToClusterMap()
{
    if ( _objectToCluster.empty() ) {
        size_t objSize = 0;
        for ( object_clundex_t cluIx = 0; cluIx < _objectsClusters.size(); cluIx++ ) {
            objSize += _objectsClusters[ cluIx ].size();
        }
        _objectToCluster.resize( objSize, CLUSTER_NA );
    }
    for ( object_clundex_t cluIx = 0; cluIx < _objectsClusters.size(); cluIx++ ) {
        for ( object_set_t::const_iterator oit = _objectsClusters[ cluIx ].items().begin(); 
              oit != _objectsClusters[ cluIx ].items().end(); ++oit ) {
            _objectToCluster[ *oit ] = cluIx;
        }
    }
}

void ChessboardBiclustering::updateObjectsClustersFromMap()
{
    if ( _objectToCluster.empty() ) {
        _objectsClusters.clear();
        return;
    }
    object_clundex_t maxClu = *std::max_element( _objectToCluster.begin(), _objectToCluster.end() );
    LOG_DEBUG1( "Restoring objects clusters, " << maxClu );
    if ( maxClu == CLUSTER_NA ) {
        THROW_RUNTIME_ERROR( "Object to cluster map contained unassigned elements" );
    }
    size_t minSize = std::min( maxClu + 1, _objectsClusters.size() );
    for ( object_clundex_t cluIx = 0; cluIx < minSize; cluIx++ ) {
        _objectsClusters[cluIx]._items.clear();
    }
    _objectsClusters.resize( maxClu + 1, ObjectsCluster() );
    for ( object_index_t objIx = 0; objIx < _objectToCluster.size(); objIx++ ) {
        object_clundex_t cluIx = _objectToCluster[ objIx ];
        if ( cluIx != CLUSTER_NA ) _objectsClusters[ cluIx ]._items.insert( objIx );
    }
}

void ChessboardBiclustering::updateProbeToClusterMap()
{
    if ( !_probesClusters.empty() ) {
        _probeToCluster.resize( _probesClusters.front().items().size(), CLUSTER_NA );
        for ( probe_clundex_t cluIx = 0; cluIx < _probesClusters.size(); cluIx++ ) {
            foreach_bit( probe_index_t, probeIx, _probesClusters[ cluIx ].items() ) {
                _probeToCluster[ probeIx ] = cluIx;
            }
        }
    }
}

void ChessboardBiclustering::updateProbesClustersFromMap()
{
    if ( _probeToCluster.empty() ) {
        _probesClusters.clear();
        return;
    }
    probe_clundex_t maxClu = *std::max_element( _probeToCluster.begin(), _probeToCluster.end() );
    LOG_DEBUG1( "Restoring probes clusters, " << maxClu );
    if ( maxClu == CLUSTER_NA ) {
        THROW_RUNTIME_ERROR( "Probe to cluster map contained unassigned elements" );
    }
    size_t minSize = std::min( maxClu + 1, _probesClusters.size() );
    for ( probe_clundex_t cluIx = 0; cluIx < minSize; cluIx++ ) {
        _probesClusters[cluIx]._items.reset();
    }
    _probesClusters.resize( maxClu + 1, ProbesCluster( _probeToCluster.size() ) );
    for ( probe_index_t probeIx = 0; probeIx < _probeToCluster.size(); probeIx++ ) {
        probe_clundex_t cluIx = _probeToCluster[ probeIx ];
        if ( cluIx != CLUSTER_NA ) _probesClusters[ cluIx ]._items.set( probeIx );
    }
}

ChessboardBiclustering::ObjectsClusterParamsProxy& ChessboardBiclustering::ObjectsClusterParamsProxy::operator=(
    const ObjectsClusterParams& params
){
    BOOST_ASSERT( cluIx < clus.objectsClusters().size() );
    BOOST_ASSERT( params.blocksMask.size() == clus.probesClusters().size() );
    for ( probe_clundex_t probeCluIx = 0; probeCluIx < clus.probesClusters().size(); ++probeCluIx ) {
        bool  enabled = params.blocksMask.test( probeCluIx );
        block_iterator cluIt = clus.findBlock( cluIx, probeCluIx, false );
        cluIt->setEnabled( enabled );
        if ( enabled ) {
            signal_t signal = params.probeSignal[ probeCluIx ];
            if ( !is_unset( signal ) ) {
                LOG_DEBUG2( "Signal[" << cluIx << "," << probeCluIx << "]:=" << signal );
                cluIt->setSignal( signal );
            } else {
                LOG_WARN( "Signal[" << cluIx << "," << probeCluIx << "] not set" );
            }
        }
        BOOST_ASSERT( cluIt->check() );
        const object_set_t& effMask = effectiveMask();
        for ( object_set_t::const_iterator maskIt = effMask.begin(); maskIt != effMask.end(); ++maskIt ) {
            size_t mult = params.objectMultiple[ *maskIt ];
            if ( mult > 0 ) clus.setObjectMultiple( *maskIt, mult );
        }
    }
    return ( *this );
}

ChessboardBiclustering::ProbesClusterParamsProxy& ChessboardBiclustering::ProbesClusterParamsProxy::operator=(
    const ProbesClusterParams& params
){
    LOG_DEBUG2( "Setting probe cluster " << cluIx << " params" );
    BOOST_ASSERT( cluIx < clus.probesClusters().size() );
    BOOST_ASSERT( params.blocksMask.size() == clus.objectsClusters().size() );
    for ( object_clundex_t objCluIx = 0; objCluIx < clus.objectsClusters().size(); ++objCluIx ) {
        bool  enabled = params.blocksMask.test( objCluIx );
        block_iterator cluIt = clus.findBlock( objCluIx, cluIx, false );
        cluIt->setEnabled( enabled );
        if ( enabled ) {
            signal_t signal = params.objectsSignal[ objCluIx ];
            if ( !is_unset( signal ) ) {
                LOG_DEBUG2( "Signal[" << objCluIx << "," << cluIx << "]:=" << signal );
                cluIt->setSignal( signal );
            } else {
                LOG_WARN( "Signal[" << objCluIx << "," << cluIx << "] not found" );
            }
        }
        // @kludge: object clusters in new enabled block could have no signals defined
        // so no checking them at this moment, but the calling code should fill the missing params
        //BOOST_ASSERT( cluIt->check() );
    }
    return ( *this );
}

ObjectsClusterParams& ObjectsClusterParams::operator=(
    const ObjectsClusterParams& params
){
    blocksMask = params.blocksMask;
    probeSignal = params.probeSignal;
    objectMultiple = params.objectMultiple;
    return ( *this );
}

ObjectsClusterParams& ObjectsClusterParams::operator=(
    const ChessboardBiclustering::ObjectsClusterParamsProxy& proxy
){
    BOOST_ASSERT( proxy.cluIx < proxy.clus.objectsClusters().size() );
    // clear
    blocksMask.resize( proxy.clus.probesClusters().size() );
    probeSignal.resize( proxy.clus.probesClusters().size(), unset() );
    objectMultiple = proxy.clus.objectMultiples(); // just copy all multiples, it should be reasonably fast

    // fill signals and mask
    for ( probe_clundex_t cluIx = 0; cluIx < proxy.clus.probesClusters().size(); ++cluIx ) {
        probeSignal[ cluIx ] = proxy.clus.blockSignal( proxy.cluIx, cluIx );
        blocksMask.set( cluIx, proxy.clus.isBlockEnabled( proxy.cluIx, cluIx ) );
    }

    return ( *this );
}

ProbesClusterParams& ProbesClusterParams::operator=(
    const ProbesClusterParams& params
){
    blocksMask = params.blocksMask;
    objectsSignal = params.objectsSignal;
    return ( *this );
}

ProbesClusterParams& ProbesClusterParams::operator=(
    const ChessboardBiclustering::ProbesClusterParamsProxy& proxy
){
    BOOST_ASSERT( proxy.cluIx < proxy.clus.probesClusters().size() );

    // clear
    blocksMask.resize( proxy.clus.objectsClusters().size() );
    objectsSignal.resize( proxy.clus.objectsClusters().size(), unset() );

    // fill signals and mask
    const probe_bitset_t& probes = proxy.mask.any() ? proxy.mask : proxy.clus.probesCluster( proxy.cluIx ).items();
    if ( probes.none() ) {
        // nothing to test
        LOG_WARN( "Empty probes for ProbesClusterParams" );
        return ( *this );
    }
    probe_index_t probeIx = probes.find_first();

    for ( object_clundex_t objCluIx = 0; objCluIx < proxy.clus.objectsClusters().size(); ++objCluIx ) {
        blocksMask.set( objCluIx, proxy.clus.isBlockEnabled( objCluIx, proxy.cluIx ) );
        const object_index_t objIx = proxy.clus.objectsCluster( objCluIx ).representative(); // representative
        if ( blocksMask.test( objCluIx ) ) {
            signal_t signal = proxy.clus.cellSignal( objIx, probeIx );
            objectsSignal[ objCluIx ] = signal;
        }
    }
    return ( *this );
}

bool ChessboardBiclustering::checkObjectsPartition() const
{
    object_set_t    objs;
    for ( object_clundex_t objCluIx = 0; objCluIx < objectsClusters().size(); objCluIx++ ) {
        const ObjectsCluster& objClu = objectsCluster( objCluIx );
        size_t objsSize = objs.size();
        if ( objClu.size() == 0 ) THROW_RUNTIME_ERROR( "Empty objects cluster #" << objCluIx );
        objs.insert( objClu.items().begin(), objClu.items().end() );
        if ( objsSize + objClu.size() > objs.size() ) THROW_RUNTIME_ERROR( "Objects cluster #" << objCluIx << " overlaps with some other clusters" );
        for ( object_set_t::const_iterator it = objClu.items().begin(); it != objClu.items().end(); ++it ) {
            if ( clusterOfObject( *it ) != objCluIx ) {
                THROW_RUNTIME_ERROR( "Object-to-cluster[" << *it << "]=" 
                                     << clusterOfObject( *it ) << " != " << objCluIx );
            }
        }
    }
    if ( objs.size() < objectsCount() )  throw std::runtime_error( "Objects partition incomplete, some objects unassigned" );
    else if ( objs.size() > objectsCount() )  throw std::runtime_error( "Objects partition size incorrect" );
    return ( true );
}

bool ChessboardBiclustering::checkProbesPartition() const
{
    probe_bitset_t    probes( probesCount() );
    for ( probe_clundex_t probeCluIx = 0; probeCluIx < probesClusters().size(); probeCluIx++ ) {
        const ProbesCluster& probeClu = probesCluster( probeCluIx );
        if ( probeClu.size() == 0 ) THROW_RUNTIME_ERROR( "Empty probes cluster #" << probeCluIx );
        size_t probesSize = probes.count();
        probes |= probeClu.items();
        if ( probesSize + probeClu.size() > probes.count() ) THROW_RUNTIME_ERROR( "Probes cluster #" << probeCluIx << " overlaps with some other clusters" );
        foreach_bit( probe_index_t, probeIx, probeClu.items() ) {
            if ( clusterOfProbe( probeIx ) != probeCluIx ) {
                THROW_RUNTIME_ERROR( "Probe-to-cluster[" << probeIx << "]=" << clusterOfProbe( probeIx ) 
                << " != " << probeCluIx );
            }
        }
    }
    if ( probes.count() < probesCount() )  throw std::runtime_error( "Probes partition incomplete, some probes unassigned" );
    else if ( probes.count() > probesCount() )  throw std::runtime_error( "Probes partition size incorrect" );
    return ( true );
}

bool ChessboardBiclustering::BlockProxy::check() const
{
    if ( objectsClusterIndex() >= _clustering.objectsClusters().size() ) {
        throw std::out_of_range( "Biclustering block references incorrect objects cluster" );
    }
    if ( probesClusterIndex() >= _clustering.probesClusters().size() ) {
        throw std::out_of_range( "Biclustering block references incorrect probes cluster" );
    }
    if ( isEnabled() ) {
        const probe_bitset_t& probes = _clustering.probesCluster( probesClusterIndex() ).items();
        const object_set_t& objs = _clustering.objectsCluster( objectsClusterIndex() ).items();
        signal_t cluSignal = unset();
        foreach_bit( probe_index_t, probeIx, probes ) {
            for ( object_set_t::const_iterator oit = objs.begin(); oit != objs.end(); ++oit ) {
                signal_t objSignal = _clustering.cellSignal( *oit, probeIx );
                if ( is_unset( objSignal ) ) {
                    THROW_RUNTIME_ERROR( "Signal(" 
                                         << *oit << "@" << objectsClusterIndex() << ", "
                                         << probeIx << "@" << probesClusterIndex() << ") not set" );
                }
                else if ( !is_unset( cluSignal ) && ( cluSignal != objSignal ) ) {
                    THROW_RUNTIME_ERROR( "Signal(" 
                                         << *oit << "@" << objectsClusterIndex() << ", "
                                         << probeIx << "@" << probesClusterIndex() << ")="
                                         << objSignal << "!=" << cluSignal );
                }
                cluSignal = objSignal;
            }
        }
    }
    return ( true );
}

bool ChessboardBiclustering::checkBlocks() const
{
    for ( const_block_iterator cluIt = begin(); cluIt != end(); ++cluIt ) {
        cluIt->check();
    }
    return ( true );
}

bool ChessboardBiclustering::check() const
{
    checkObjectsPartition();
    checkProbesPartition();
    checkBlocks();
    return ( true );
}

void ChessboardBiclustering::outputSignals(
    std::ostream& out
) const {
    for ( object_index_t objIx = 0; objIx < objectsCount(); ++objIx ) {
        for ( probe_index_t probeIx = 0; probeIx < probesCount(); ++probeIx ) {
            out << boost::format("%.3f ") % _signals( objIx, probeIx );
        }
        out << "\n";
    }
}

std::ostream& operator<<( std::ostream& out, const ChessboardBiclustering& cc )
{
    out << "N_obj_clu=" << cc.objectsClusters().size()
        << " N_probe_clu=" << cc.probesClusters().size()
        << " N_cc=" << cc.enabledBlocksCount()
        << " Noise_rate=" << cc.clusteringData()._noiseParams.successRate;
    return ( out );
}
