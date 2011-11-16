#include "dynamic_bitset_utils.h"
#include "ChessboardBiclusteringFitInternal.h"

#include "ChessboardBiclusteringFit.h"

BOOST_CLASS_EXPORT( ChessboardBiclusteringFit );

ChessboardBiclusteringFit::ChessboardBiclusteringFit(
    const PrecomputedData&                  precomputed,
    const ChessboardBiclusteringPriors&     priors
) : ChessboardBiclustering( precomputed.data().objectsCount(), 
                            precomputed.data().probesCount() )
  , _precomputed( precomputed )
  , _priors( priors )
  , _blockIsSignalLLHValid( false )
  , _blockIsNoiseLLHValid( false )
{
    _baselineSignalParams._scShape = _precomputed.signalParams().scShape;
}

ChessboardBiclusteringFit::ChessboardBiclusteringFit(
    const PrecomputedData&          precomputed,
    const ChessboardBiclusteringPriors&    priors,
    const ChessboardBiclustering&          clus
) : ChessboardBiclustering( clus )
  , _precomputed( precomputed )
  , _priors( priors )
  , _objClusterLLH( objectsClusters().size(), LLHMetrics() )
  , _probeClusterLLH( probesClusters().size(), LLHMetrics() )
  , _objClusterLPP( objectsClusters().size(), unset() )
  , _blockLLH( objectsClusters().size(), probesClusters().size(), unset() )
  , _blockLPP( objectsClusters().size(), probesClusters().size(), unset() )
  , _blockIsSignalLLHValid( false )
  , _blockIsSignalLLH( objectsClusters().size(), probesClusters().size(), unset() )
  , _blockIsNoiseLLHValid( false )
  , _blockIsNoiseLLH( objectsClusters().size(), probesClusters().size(), unset() )
  , _blockToSample( objectsClusters().size(), probesClusters().size(), 0 )
{
    _baselineSignalParams._scShape = _precomputed.signalParams().scShape;
}

void ChessboardBiclusteringFit::updateMetrics() const
{
    const_cast<ChessboardBiclusteringFit&>( *this ).cleanupClusters();

    ChessboardBiclusteringPriorEval eval = PriorEval( *this );
    if ( is_unset( _metrics.lppBlocks ) ) {
        double res = 0;
        for ( object_clundex_t objCluIx = 0; objCluIx < objectsClusters().size(); objCluIx++ ) {
            log_prob_t& objcLPP = _objClusterLPP[ objCluIx ];
            if ( is_unset( objcLPP ) ) {
                double ccLnPSum = 0;
                for ( probe_clundex_t probeCluIx = 0; probeCluIx < probesClusters().size(); probeCluIx++ ) {
                    log_prob_t& ccLnP = _blockLPP( objCluIx, probeCluIx );
                    if ( is_unset( ccLnP ) ) {
                        ccLnP = eval.blockLPP( objCluIx, probeCluIx );
                        if ( is_unset( ccLnP ) ) THROW_RUNTIME_ERROR( "Block (" << objCluIx << ", " << probeCluIx << ") lpp eval error" );
                    }
                    ccLnPSum += ccLnP;
                }
                objcLPP = ccLnPSum;
                if ( is_unset( objcLPP ) ) THROW_RUNTIME_ERROR( "Object's cluster #" << objCluIx << " lpp eval error" );
                const object_set_t& objects = objectsCluster( objCluIx ).items();
                BOOST_ASSERT( objectMultiples().size() == objectsCount() );
                GeometricDistribution omp = eval.objectMultiplePrior( objects.size() );
                for ( object_set_t::const_iterator oit = objects.begin(); oit != objects.end(); ++oit ) {
                    objcLPP += omp( objectMultiple( *oit ) );
                }
            }
            res += objcLPP;
        }
        _metrics.lppBlocks = res + eval.noiseParamsLPP();
    }

    if ( is_unset( _metrics.lppObjClu ) ) {
        _metrics.lppObjClu = eval.objectsClusteringLPP();
        if ( is_unset( _metrics.lppObjClu ) ) THROW_RUNTIME_ERROR( "Objects' clustering lpp eval error" );
    }
    if ( is_unset( _metrics.lppProbesClu ) ) {
        _metrics.lppProbesClu = eval.probesClusteringLPP();
        if ( is_unset( _metrics.lppProbesClu ) ) THROW_RUNTIME_ERROR( "Probes' clustering lpp eval error" );
    }

    updateClustersLLH();
    if ( _metrics.llhObjs.is_unset() ) {
        _metrics.llhObjs.topo = _metrics.llhObjs.quant = _metrics.llhObjs.conf = 0;
        for ( object_clundex_t objCluIx = 0; objCluIx < objectsClusters().size(); objCluIx++ ) {
            BOOST_ASSERT( !is_unset( _objClusterLLH[ objCluIx ] ) );
            _metrics.llhObjs += _objClusterLLH[ objCluIx ];
        }
    }
    if ( _metrics.llhProbes.is_unset() ) {
        _metrics.llhProbes.topo = _metrics.llhProbes.quant = _metrics.llhProbes.conf = 0;
        for ( probe_clundex_t probeCluIx = 0; probeCluIx < probesClusters().size(); probeCluIx++ ) {
            BOOST_ASSERT( !is_unset( _probeClusterLLH[ probeCluIx ] ) );
            _metrics.llhProbes += _probeClusterLLH[ probeCluIx ];
        }
    }
}

void ChessboardBiclusteringFit::updateClustersLLH() const
{
    const_cast<ChessboardBiclusteringFit&>( *this ).cleanupClusters();
    ChessboardBiclusteringLLHEval eval = LLHEval( *this );
    ChessboardBiclusteringStructureLLHEval structEval = StructureLLHEval( *this );
    DataSignalNoiseCache& snCache = signalNoiseCache(); // implicitly updates cache
    // update objects clusters LLH
    for ( object_clundex_t objCluIx = 0; objCluIx < objectsClusters().size(); objCluIx++ ) {
        LLHMetrics& cluLLH = _objClusterLLH[ objCluIx ];
        const object_set_t& objects = objectsCluster( objCluIx ).items();
        if ( is_unset( cluLLH.quant ) ) {
            cluLLH.quant = 0.0;
            for ( probe_clundex_t probeCluIx = 0; probeCluIx < probesClusters().size(); probeCluIx++ ) {
                log_prob_t& blkLLH = _blockLLH( objCluIx, probeCluIx );
                if ( is_unset( blkLLH ) ) {
                    blkLLH = isBlockEnabled( objCluIx, probeCluIx )
                          ? eval.blockDataLLH( objCluIx, probeCluIx )
                          : snCache.noiseLLH( objects, probesCluster( probeCluIx ).items() );
                    if ( is_unset( blkLLH ) ) THROW_RUNTIME_ERROR( "Block[" << isBlockEnabled( objCluIx, probeCluIx ) << "](" 
                                                                  << objCluIx << ", " << probeCluIx << ") LLH eval error" );
                }
                cluLLH.quant += blkLLH;
            }
        }
        if ( is_unset( cluLLH.topo ) ) {
            cluLLH.topo = structEval.objectsClusterMismatchLLH( objects );
        }
        if ( is_unset( cluLLH.conf ) ) {
            cluLLH.conf = structEval.probeClustersPerObjectClusterLLH( boundProbesClusters( objects ).size() );
        }
        if ( cluLLH.is_unset() ) THROW_RUNTIME_ERROR( "Object's cluster #" << objCluIx << " LLH eval error" );
    }
    // update probes clusters LLH
    for ( probe_clundex_t probeCluIx = 0; probeCluIx < probesClusters().size(); probeCluIx++ ) {
        LLHMetrics& cluLLH = _probeClusterLLH[ probeCluIx ];
        if ( is_unset( cluLLH.conf ) ) {
            cluLLH.conf = structEval.objectClustersPerProbeClusterLLH( boundObjectsClusters( probeCluIx ).size() );
        }
        if ( is_unset( cluLLH.topo ) ) {
            cluLLH.topo = structEval.probesClusterMismatchLLH( probesClusters()[ probeCluIx ].items() );
        }
        if ( is_unset( cluLLH.quant ) ) {
            double ccLLHSum = 0;
            for ( object_clundex_t objCluIx = 0; objCluIx < objectsClusters().size(); objCluIx++ ) {
                const double& ccLLH = _blockLLH( objCluIx, probeCluIx );
                BOOST_ASSERT( !is_unset( ccLLH ) );
                ccLLHSum += ccLLH;
            }
            cluLLH.quant = ccLLHSum;
        }
        if ( cluLLH.is_unset() ) THROW_RUNTIME_ERROR( "Probes's cluster #" << probeCluIx << " llh eval error" );
    }
}

/**
 *  Objects clusters that contain baits of given probes.
 */
std::set<object_clundex_t> ChessboardBiclusteringFit::boundObjectsClusters(
    const probe_bitset_t& probes
) const {
    std::set<object_clundex_t> res;
    foreach_bit( probe_index_t, probeIx, probes ) {
        object_index_t  baitIx = data().probe( probeIx ).baitIndex();
        if ( baitIx != OBJECT_NA ) res.insert( clusterOfObject( baitIx ) );
    }
    return ( res );
}

/**
 *  Probe clusters, which baits are contained in given objects cluster.
 */
std::set<probe_clundex_t> ChessboardBiclusteringFit::boundProbesClusters(
    const object_set_t& objects
) const {
    std::set<probe_clundex_t> res;
    for ( object_set_t::const_iterator oit = objects.begin(); oit != objects.end(); ++oit ) {
        std::pair<OPAData::const_bait_to_probe_iterator, OPAData::const_bait_to_probe_iterator> probesRange = data().probesOfBait( *oit );
        for ( OPAData::const_bait_to_probe_iterator stit = probesRange.first; stit != probesRange.second; ++stit ) {
            res.insert( clusterOfProbe( stit->second ) );
        }
    }
    return ( res );
}

void ChessboardBiclusteringFit::updateBlocksProbesLLH() const
{
    // object block probe LLH
    for ( probe_clundex_t probeCluIx = 0; probeCluIx < probesClusters().size(); probeCluIx++ ) {
        for ( object_clundex_t objCluIx = 0; objCluIx < objectsClusters().size(); objCluIx++ ) {
            log_prob_t& noiseLLH = _blockIsNoiseLLH( objCluIx, probeCluIx );
            if ( is_unset( noiseLLH ) ) {
                noiseLLH = BlockEnablementDataLLH( signalNoiseCache(),
                                                          objectsCluster( objCluIx ).items(),
                                                          probesCluster( probeCluIx ).items() )( false );
                if ( is_unset( noiseLLH ) ) THROW_RUNTIME_ERROR( "Block (" << objCluIx << ", " << probeCluIx << ") noise llh eval error" );
            }
            log_prob_t& signalLLH = _blockIsSignalLLH( objCluIx, probeCluIx );
            if ( is_unset( signalLLH ) ) {
                // P(N<=n|noise)P(N<=|min_signal) if enabled
                signalLLH = BlockEnablementDataLLH( signalNoiseCache(), 
                                                           objectsCluster( objCluIx ).items(), 
                                                           probesCluster( probeCluIx ).items() )( true );
                if ( is_unset( signalLLH ) ) THROW_RUNTIME_ERROR( "Block (" << objCluIx << ", " << probeCluIx << ") signal llh eval error" );
            }
        }
    }
    _blockIsNoiseLLHValid = _blockIsSignalLLHValid = true;
}

ChessboardBiclusteringFit& ChessboardBiclusteringFit::operator=( const ChessboardBiclustering& that )
{
    ChessboardBiclustering::operator=( that );
    if ( _signalNoiseCache ) {
        _signalNoiseCache->setChessboardBiclusteringData( that.clusteringData() );
    }
    _blockToSample.reset( objectsClusters().size(),
                                 probesClusters().size(), 1 );
    resetAllCaches();
    return ( *this );
}

ChessboardBiclusteringFit& ChessboardBiclusteringFit::operator=( const StaticChessboardBiclustering& that )
{
    BOOST_ASSERT( that.clustering );
    ChessboardBiclustering::operator=( *that.clustering );
    if ( _signalNoiseCache ) {
        _signalNoiseCache->setChessboardBiclusteringData( that.clustering->clusteringData() );
    }
    if ( that.blocksToSample ) {
        BOOST_ASSERT( that.blocksToSample->size1() == objectsClusters().size() );
        BOOST_ASSERT( that.blocksToSample->size2() == probesClusters().size() );
        _blockToSample = *that.blocksToSample;
    } else {
        _blockToSample.reset( objectsClusters().size(),
                                     probesClusters().size(), 1 );
    }
    resetAllCaches();
    return ( *this );
}

ChessboardBiclusteringFit& ChessboardBiclusteringFit::operator=(const ChessboardBiclusteringFit& that) {
    BOOST_ASSERT( &_precomputed == &that._precomputed );
    BOOST_ASSERT( &_priors == &that._priors );
    _metrics = that._metrics;
    _probeClusterLLH = that._probeClusterLLH;
    _signalNoiseCache = that._signalNoiseCache;
    _blockLLH = that._blockLLH;
    _blockLPP = that._blockLPP;
    _blockToSample = that._blockToSample;
    _objClusterLLH = that._objClusterLLH;
    _objClusterLPP = that._objClusterLPP;
    _blockIsNoiseLLHValid = that._blockIsNoiseLLHValid;
    _blockIsNoiseLLH = that._blockIsNoiseLLH;
    _blockIsSignalLLHValid = that._blockIsSignalLLHValid;
    _blockIsSignalLLH = that._blockIsSignalLLH;
    ChessboardBiclustering::operator=( that );
    return ( *this );
}

void ChessboardBiclusteringFit::resetAllCaches() const
{
    _blockLLH.reset( objectsClusters().size(),
                            probesClusters().size(), unset() );
    _blockLPP.reset( objectsClusters().size(),
                            probesClusters().size(), unset() );
    _objClusterLLH.assign( objectsClusters().size(), LLHMetrics() );
    _objClusterLPP.assign( objectsClusters().size(), unset() );
    _probeClusterLLH.assign( probesClusters().size(), LLHMetrics() );
    _blockIsNoiseLLHValid = false;
    _blockIsNoiseLLH.reset( objectsClusters().size(),
                                   probesClusters().size(), unset() );
    _blockIsSignalLLHValid = false;
    _blockIsSignalLLH.reset( objectsClusters().size(),
                                    probesClusters().size(), unset() );
    _metrics.unset( true, true );
}

void ChessboardBiclusteringFit::resetObjectsClusterCache(
    object_clundex_t    objCluIx,
    bool                contentsChanged
) const {
    _blockLLH.fill1( objCluIx, unset() );
    _blockLPP.fill1( objCluIx, unset() );
    if ( contentsChanged ) {
        _blockIsNoiseLLHValid = false;
        _blockIsSignalLLHValid = false;
        _blockIsSignalLLH.fill1( objCluIx, unset() );
        _blockIsNoiseLLH.fill1( objCluIx, unset() );
    }
    _objClusterLLH[ objCluIx ].unset();
    _probeClusterLLH.assign( _probeClusterLLH.size(), LLHMetrics() );
    _objClusterLPP[ objCluIx ] = unset();
    _metrics.unset( true, false );
}

void ChessboardBiclusteringFit::resetProbesClusterCache(
    probe_clundex_t probeCluIx
) const {
    _blockLLH.fill2( probeCluIx, unset() );
    _blockLPP.fill2( probeCluIx, unset() );
    _blockIsNoiseLLHValid = false;
    _blockIsSignalLLHValid = false;
    _blockIsSignalLLH.fill2( probeCluIx, unset() );
    _blockIsNoiseLLH.fill2( probeCluIx, unset() );
    _probeClusterLLH[ probeCluIx ].unset();
    _objClusterLLH.assign( _objClusterLLH.size(), LLHMetrics() );
    _objClusterLPP.assign( _objClusterLPP.size(), unset() );
    _metrics.unset( false, true );
}

void ChessboardBiclusteringFit::resetCachesAfterLoading()
{
    _objClusterLLH.resize( objectsClusters().size(), LLHMetrics() );
    _objClusterLLH.assign( _objClusterLLH.size(), LLHMetrics() );
    _objClusterLPP.resize( objectsClusters().size(), unset() );
    _objClusterLPP.assign( _objClusterLPP.size(), unset() );
    _probeClusterLLH.resize( probesClusters().size(), LLHMetrics() );
    _probeClusterLLH.assign( probesClusters().size(), LLHMetrics() );
    // _metrics is serialized and restored, so it should not be reset
}

void ChessboardBiclusteringFit::setObjectsClusterSamples( object_clundex_t objCluIx, size_t counts )
{
    for ( probe_clundex_t cluIx = 0; cluIx < probesClusters().size(); cluIx++ ) {
        setBlockSamples( objCluIx, cluIx, counts );
    }
}

void ChessboardBiclusteringFit::setProbesClusterSamples( probe_clundex_t probeCluIx, size_t counts )
{
    for ( object_clundex_t cluIx = 0; cluIx < objectsClusters().size(); cluIx++ ) {
        setBlockSamples( cluIx, probeCluIx, counts );
    }
}

void ChessboardBiclusteringFit::setAllBlockSamples(size_t counts)
{
    for ( object_clundex_t objCluIx = 0; objCluIx < objectsClusters().size(); objCluIx++ ) {
        for ( probe_clundex_t probeCluIx = 0; probeCluIx < probesClusters().size(); probeCluIx++ ) {
            setBlockSamples( objCluIx, probeCluIx, counts );
        }
    }
}

void ChessboardBiclusteringFit::beforeObjectsClusterRemoved(object_clundex_t cluIx) const
{
    ChessboardBiclustering::beforeObjectsClusterRemoved(cluIx);
    const_cast<block_counts_matrix_type&>( _blockToSample ).remove1( cluIx );
    _blockLLH.remove1( cluIx );
    _blockLPP.remove1( cluIx );
    _blockIsNoiseLLH.remove1( cluIx );
    _blockIsSignalLLH.remove1( cluIx );
    _objClusterLLH.erase( _objClusterLLH.begin() + cluIx );
    _objClusterLPP.erase( _objClusterLPP.begin() + cluIx );
    _probeClusterLLH.assign( _probeClusterLLH.size(), LLHMetrics() );
    _metrics.unset( true, false );
}

void ChessboardBiclusteringFit::beforeProbesClusterRemoved(probe_clundex_t cluIx) const
{
    ChessboardBiclustering::beforeProbesClusterRemoved(cluIx);
    const_cast<block_counts_matrix_type&>( _blockToSample ).remove2( cluIx );
    _blockLLH.remove2( cluIx );
    _blockLPP.remove2( cluIx );
    _blockIsNoiseLLH.remove2( cluIx );
    _blockIsSignalLLH.remove2( cluIx );
    _probeClusterLLH.erase( _probeClusterLLH.begin() + cluIx );
    _objClusterLLH.assign( _objClusterLLH.size(), LLHMetrics() );
    _objClusterLPP.assign( _objClusterLPP.size(), unset() );
    _metrics.unset( false, true );
}

void ChessboardBiclusteringFit::afterObjectsClusterInserted(object_clundex_t cluIx) const
{
    ChessboardBiclustering::afterObjectsClusterInserted(cluIx);
    const_cast<block_counts_matrix_type&>( _blockToSample ).insert1( cluIx );
    _blockLLH.insert1( cluIx, unset() );
    _blockLPP.insert1( cluIx, unset() );
    _blockIsNoiseLLHValid = false;
    _blockIsSignalLLHValid = false;
    _blockIsNoiseLLH.insert1( cluIx, unset() );
    _blockIsSignalLLH.insert1( cluIx, unset() );
    _objClusterLLH.insert( _objClusterLLH.begin() + cluIx, LLHMetrics() );
    _objClusterLPP.insert( _objClusterLPP.begin() + cluIx, unset() );
    _probeClusterLLH.assign( _probeClusterLLH.size(), LLHMetrics() );
    _metrics.unset( true, false );
}

void ChessboardBiclusteringFit::afterProbesClusterInserted(probe_clundex_t cluIx) const
{
    ChessboardBiclustering::afterProbesClusterInserted(cluIx);
    const_cast<block_counts_matrix_type&>( _blockToSample ).insert2( cluIx );
    _blockLLH.insert2( cluIx, unset() );
    _blockLPP.insert2( cluIx, unset() );
    _blockIsNoiseLLHValid = false;
    _blockIsSignalLLHValid = false;
    _blockIsNoiseLLH.insert2( cluIx, unset() );
    _blockIsSignalLLH.insert2( cluIx, unset() );
    _probeClusterLLH.insert( _probeClusterLLH.begin() + cluIx, LLHMetrics() );
    _objClusterLLH.assign( _objClusterLLH.size(), LLHMetrics() );
    _objClusterLPP.assign( _objClusterLPP.size(), unset() );
    _metrics.unset( false, true );
}

void ChessboardBiclusteringFit::afterObjectsClusterChanged(object_clundex_t cluIx) const
{
    ChessboardBiclustering::afterObjectsClusterChanged(cluIx);
    resetObjectsClusterCache( cluIx, true );
}

void ChessboardBiclusteringFit::afterProbesClusterChanged(probe_clundex_t cluIx) const
{
    ChessboardBiclustering::afterProbesClusterChanged(cluIx);
    resetProbesClusterCache( cluIx );
}

void ChessboardBiclusteringFit::afterObjectMultipleChanged(object_index_t objIx) const
{
    ChessboardBiclustering::afterObjectMultipleChanged(objIx);
    resetObjectsClusterCache( clusterOfObject( objIx ), false );
}

void ChessboardBiclusteringFit::afterSignalChanged(
    object_clundex_t    objCluIx, 
    probe_clundex_t     probeCluIx
) const {
    ChessboardBiclustering::afterSignalChanged( objCluIx, probeCluIx );
    _blockLLH( objCluIx, probeCluIx ) = unset();
    _blockLPP( objCluIx, probeCluIx ) = unset();
    _objClusterLLH[ objCluIx ] = LLHMetrics();
    _probeClusterLLH[ probeCluIx ] = LLHMetrics();
    _objClusterLPP[ objCluIx ] = unset();
    _metrics.unset( false, false );
}

void ChessboardBiclusteringFit::afterBlockFlipped(
    object_clundex_t    objCluIx, 
    probe_clundex_t     probeCluIx
) const {
    ChessboardBiclustering::afterBlockFlipped(objCluIx, probeCluIx);
    _blockLLH( objCluIx, probeCluIx ) = unset();
    _blockLPP( objCluIx, probeCluIx ) = unset();
    _objClusterLLH[ objCluIx ] = LLHMetrics();
    _probeClusterLLH[ probeCluIx ] = LLHMetrics();
    _objClusterLPP[ objCluIx ] = unset();
    _metrics.unset( false, false );
}

void ChessboardBiclusteringFit::afterSignalPriorChanged() const
{
    ChessboardBiclustering::afterSignalPriorChanged();
    if ( !_signalNoiseCache.unique() ) {
        _signalNoiseCache.reset( new DataSignalNoiseCache( *_signalNoiseCache ) );
    }
    _signalNoiseCache->setChessboardBiclusteringData( clusteringData() );
    resetAllCaches();
}

void ChessboardBiclusteringFit::afterNoiseParamsChanged() const
{
    ChessboardBiclustering::afterNoiseParamsChanged();
    if ( !_signalNoiseCache.unique() ) {
        _signalNoiseCache.reset( new DataSignalNoiseCache( *_signalNoiseCache ) );
    }
    _signalNoiseCache->setChessboardBiclusteringData( clusteringData() );
    resetAllCaches();
}
