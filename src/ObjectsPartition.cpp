#include "dynamic_bitset_utils.h"

#include "ChessboardBiclusteringFitInternal.h"

#include "ObjectsPartition.h"

std::string ObjectsPartition::ObjectsClusterProxy::label() const
{
    std::ostringstream out;
    out << "O";
    for ( object_set_t::const_iterator objIt = cluster().items().begin(); objIt != cluster().items().end(); ++objIt ) {
        out << (int)( *objIt );
    }
    return ( out.str() );
}

void ObjectsPartitionEx::putToCluster(
    element_index_type  elmIx,
    cluster_index_type  cluIx
){
    unshare();
    wrapped().setObjectCluster( elmIx, cluIx );
}

/**
 *  Exchange given set of elements between 2 clusters.
 *  If second cluster index is PROBE_NA, new cluster is created.
 */
ObjectsPartitionEx::cluster_index_type ObjectsPartitionEx::exchangeElements(
    cluster_index_type clu1Ix, /** i first cluster */
    cluster_index_type clu2Ix, /** i second cluster */
    const element_index_set_type& clu2Elements /** i elements of updated clu2 cluster */
){
    unshare();
    object_clundex_t res = wrapped().exchangeObjects( clu1Ix, clu2Ix, clu2Elements );
    wrapped().cleanupClusters();
    return ( res );
}

LLHMetrics FixedObjectsPartitionStats::objectLLH(
    object_index_t          objIx,
    const params_type&      params
) const {
    object_set_t objs( &objIx, &objIx + 1 );
    return ( objectsLLH( objs, params ) );
}

LLHMetrics FixedObjectsPartitionStats::objectsLLH(
    const object_set_t&     objs,
    const params_type&      params
) const {
    LLHMetrics llh;

    // quantitative component
    llh.quant = 0;
    ChessboardBiclusteringLLHEval eval = LLHEval( clusFit );
    // enabled cells
    foreach_bit( probe_clundex_t, probeCluIx, params.blocksMask ) {
        const ProbesCluster&   cluster = clusFit.probesCluster( probeCluIx );
        BOOST_ASSERT( !is_unset( params.probeSignal[ probeCluIx ] ) );
        llh.quant += eval.cellsDataLLH( objs, cluster.items(), params.probeSignal[ probeCluIx ],
                                        params.objectMultiple );
        BOOST_ASSERT( is_finite( llh.quant ) );
    }

    // disabled cells
    probe_bitset_t notMask = ~params.blocksMask;
    const DataSignalNoiseCache& snCache = clusFit.signalNoiseCache();
    foreach_bit( probe_clundex_t, probeCluIx, notMask ) {
        llh.quant += snCache.noiseLLH( objs, clusFit.probesCluster( probeCluIx ).items() );
        BOOST_ASSERT( is_finite( llh.quant ) );
    }

    // other components
    ChessboardBiclusteringStructureLLHEval structEval = StructureLLHEval( clusFit );
    llh.topo = structEval.objectsClusterMismatchLLH( objs );
    llh.conf = structEval.probeClustersPerObjectClusterLLH( clusFit.boundProbesClusters( objs ).size() );

    return ( llh );
}

LLHMetrics FixedObjectsPartitionStats::llhDelta(
    const std::vector<ObjectsPartition::elements_set_proxy_type>& newClusters,    /** new clusters */
    const std::vector<params_type>&                     newParams,      /** params of new clusters */
    const ObjectsPartition::cluster_index_set_type&     oldIndexes      /** indicies of clusters,
                                                                            whose objects are in newClusters */
) const {
    typedef std::set<probe_clundex_t> probeclu_set_t;

    BOOST_ASSERT( newClusters.size() == newParams.size() );

    // calculate llh
    LLHMetrics llh( 0 );
    for ( size_t i = 0; i < newClusters.size(); i++ ) {
        // new cluster's internal LLH
        llh += objectsLLH( newClusters[i], newParams[i] );
    }
    // subtract llh of old clusters
    probeclu_set_t boundProbeClus;
    for ( ObjectsPartition::cluster_index_set_type::const_iterator oldCluIt = oldIndexes.begin(); oldCluIt != oldIndexes.end(); ++oldCluIt ) {
        // subtract block fitting LLH of old clusters
        cluster_index_t oldCluIx = *oldCluIt;
        llh -= clusFit.objectsClusterLLH( oldCluIx );
        probeclu_set_t probeClus = clusFit.boundProbesClusters( oldCluIx );
        boundProbeClus.insert( probeClus.begin(), probeClus.end() );
    }
    // update delta due to changed objects clusters bound to probe cluster
    ChessboardBiclusteringStructureLLHEval structEval = StructureLLHEval( clusFit );
    for ( probeclu_set_t::const_iterator scit = boundProbeClus.begin(); scit != boundProbeClus.end(); ++scit ) {
        object_set_t boundObjs = clusFit.data().probesToBaits( clusFit.probesCluster( *scit ).items() );
        std::set<object_clundex_t> boundObjClusNew;
        for ( object_set_t::const_iterator oit = boundObjs.begin(); oit != boundObjs.end(); ++oit ) {
            // try to search probe in new clusters
            bool found = false;
            for ( object_clundex_t i = 0; i < newClusters.size(); i++ ) {
                if ( newClusters[i].contains( *oit ) ) {
                    LOG_DEBUG3( "StructLLHDelta: found object " << *oit << " in new cluster " << i );
                    boundObjClusNew.insert( clusFit.objectsClusters().size() + i );
                    found = true;
                    break;
                }
            }
            // use unmodified clusters
            if ( !found ) {
                LOG_DEBUG3( "StructLLHDelta: for object " << *oit 
                            << " using old cluster " << clusFit.clusterOfObject( *oit ) );
                boundObjClusNew.insert( clusFit.clusterOfObject( *oit ) );
            }
        }
        size_t curClusCount = clusFit.boundObjectsClusters( *scit ).size();
        if ( curClusCount != boundObjClusNew.size() ) {
            LOG_DEBUG3( "Detected change of bound objects clusters in probes cluster #" << *scit
                        << ": from " << curClusCount << " to " << boundObjClusNew.size() );
            llh.conf += structEval.objectClustersPerProbeClusterLLH( boundObjClusNew.size() )
                     - structEval.objectClustersPerProbeClusterLLH( clusFit.boundObjectsClusters( *scit ).size() );
        }
    }
    return ( llh );
}

log_prob_t FixedObjectsPartitionStats::paramsLPP(
    const params_type&  params,
    const object_set_t& clusterObjects,
    const object_set_t& sampledObjects
) const {
    log_prob_t  lpp = 0;
    ChessboardBiclusteringPriorEval eval( clusFit.data(), clusFit.priors(), clusFit );
    for ( probe_clundex_t i = 0; i < params.blocksMask.size(); i++ ) {
        lpp += params.blocksMask.test( i )
             ? eval.blockEnablementPrior( params.probeSignal[i] )( true )
             : eval.blockEnablementPrior()( false );
    }
    object_set_t unionObjs = clusterObjects;
    unionObjs.insert( sampledObjects.begin(), sampledObjects.end() );
    for ( object_set_t::const_iterator oit = sampledObjects.begin(); oit != sampledObjects.end(); ++oit ) {
        if ( params.objectMultiple[ *oit ] > 0 ) {
            lpp += eval.objectMultipleLPP( unionObjs.size(), params.objectMultiple[ *oit ] );
        }
    }
    return ( lpp );
}

bool FixedObjectsParamsSampler::operator()(
    params_type&            params,
    const object_set_t&     clusterObjs,
    const object_set_t&     sampledObjects,
    bool                    overwrite,  /** i  overwrite parameters, 
                                               if false, only non-specified params are set */
    bool                    posterior   /** i  use current values (from partition, not params) to obtain posterior sample
                                               otherwise, get prior sample */
) const {
    LOG_DEBUG2( "Sampling objects cluster parameters, "
                << clusterObjs.size() << " in cluster, "
                << sampledObjects.size() << " sampled, "
                << (overwrite ? "overwrite" : "no-overwrite")
                << " " << (posterior ? "posterior" : "no-posterior" ) );
    bool randomized = false;
    object_set_t allObjects = clusterObjs;
    allObjects.insert( sampledObjects.begin(), sampledObjects.end() );
    if ( allObjects.empty() ) {
        THROW_EXCEPTION( std::invalid_argument, "Objects set for sampling is empty" );
    }

    // sample block enablement if told so, or new cluster
    LOG_DEBUG2( "Sampling probes-clusters enablement" );
    for ( probe_clundex_t probeCluIx = 0; probeCluIx < clusHelper->clustering().probesClusters().size(); ++probeCluIx ) {
        const ProbesCluster&    probeClu = clusHelper->clustering().probesCluster( probeCluIx );

        //bool wasEnabled = params.blocksMask.test( probeCluIx );
        bool isEnabled = params.blocksMask.test( probeCluIx );
        if ( overwrite && sampleBlockMask ) {
            isEnabled = posterior
                      ? clusHelper->sampleBlockEnablement( allObjects, probeClu.items(), isEnabled ).value
                      : PriorEval( clusHelper->clusteringFit() )
                        .blockEnablementPrior()
                        .generate( clusHelper->rndNumGen() );
            randomized = true;
        }
        params.blocksMask.set( probeCluIx, isEnabled );

        if ( isEnabled && sampleSignals ) {
            LOG_DEBUG2( "Block is enabled" );
            // if block wasn't enabled before, make sure all its probes would have signals
            // randomize signals for all cluster's probes
            randomized = true;
            object_clundex_t objCluIx = clusterObjs.empty() 
                                      ? CLUSTER_NA 
                                      : clusHelper->clustering().clusterOfObject( *clusterObjs.begin() ); /// @fixme non-perfect way to get cluster index
            LOG_DEBUG2( "Generating random signal for block (" << objCluIx << ", " << probeCluIx << ")" );
            signal_t signal = params.probeSignal[ probeCluIx ];
            if ( !is_unset( signal ) && !overwrite ) continue;
            const probe_bitset_t& probes = clusHelper->clustering().probesCluster( probeCluIx ).items();
            /// @todo keep track of logRatio for the whole params object
            signal_t newSignal = !is_unset( signal ) && ( objCluIx != CLUSTER_NA ) && posterior
                            ? clusHelper->sampleSignal( allObjects, probes, signal ).value
                            : clusHelper->initialSignal( allObjects, probes, NULL ).value;
            params.probeSignal[ probeCluIx ] = newSignal;
        }
        else {
            LOG_DEBUG2( "Cluster is disabled" );
        }
    }
    if ( sampleMultiples && overwrite ) {
        LOG_DEBUG2( "Sampling multiples" );
        GeometricDistribution omp = PriorEval( clusHelper->clusteringFit() ).objectMultiplePrior( clusterObjs.size() );
        const object_set_t& actualSampledObjects = sampledObjects.empty() ? allObjects : sampledObjects;
        for ( object_set_t::const_iterator objIt = actualSampledObjects.begin(); objIt != actualSampledObjects.end(); ++objIt ) {
            //const object_multiple_map_t& multMap = clusHelper->clustering().objectMultiples();
            //object_multiple_map_t::const_iterator   objMultIt = multMap.find( objIx );
            params.objectMultiple[ *objIt ] = posterior
                                            ? clusHelper->sampleObjectMultiple( *objIt ).value
                                            : omp.generate( clusHelper->rndNumGen() );
        }
        randomized = true;
    }
    return ( randomized );
}

bool ObjectsParamsSampler::operator()(
    const ObjectsPartition&     ptn, 
    params_type&                params, 
    const object_set_t&         clusterObjects, 
    const object_set_t&         sampledObjects, 
    bool                        overwrite,
    bool                        posterior
) const {
    ChessboardBiclusteringGibbsHelper helper = clusSampler->createGibbsHelper( ptn );
    return ( FixedObjectsParamsSampler( helper, sampleBlockMask, sampleSignals, sampleMultiples )
                                      ( params, clusterObjects, sampledObjects,
                                        overwrite, posterior ) );
}

log_prob_t ObjectsParamsSampler::transitionLP(
    const ObjectsPartition&                     before,
    const ObjectsPartition&                     after,
    object_clundex_t                            cluIx
) const {
    BOOST_ASSERT( ((const ChessboardBiclustering&)before).probesClusters().size()
                   == ((const ChessboardBiclustering&)after).probesClusters().size() );
    ChessboardBiclusteringGibbsHelper helperBefore = clusSampler->createGibbsHelper( before );
    ChessboardBiclusteringGibbsHelper helperAfter = clusSampler->createGibbsHelper( after );
    ChessboardBiclusteringPriorEval priorEval = PriorEval( helperBefore.clusteringFit() );
    ChessboardBiclusteringLLHEval   llhEvalBefore = LLHEval( helperBefore.clusteringFit() );
    ChessboardBiclusteringLLHEval   llhEvalAfter = LLHEval( helperAfter.clusteringFit() );
    const object_set_t& objsBefore = before.cluster( cluIx ).items();
    log_prob_t  lnProb = 0;

    for ( probe_clundex_t probeCluIx = 0; probeCluIx < helperBefore.clustering().probesClusters().size(); ++probeCluIx ) {
        ChessboardBiclustering::block_proxy cluBefore = *((const ChessboardBiclustering&)before).findBlock( cluIx, probeCluIx, false );
        ChessboardBiclustering::block_proxy cluAfter = *((const ChessboardBiclustering&)after).findBlock( cluIx, probeCluIx, false );
        bool enabledBefore = cluBefore.isEnabled();
        bool enabledAfter = cluAfter.isEnabled();
        if ( enabledBefore != enabledAfter ) {
            // block probe transition probability
            ChessboardBiclusteringGibbsHelper::BlockEnablementDataLLHCached llhBefore( helperBefore.energyEval(), helperBefore.clusteringFit(), cluIx, probeCluIx );
            ChessboardBiclusteringGibbsHelper::BlockEnablementDataLLHCached llhAfter( helperAfter.energyEval(), helperAfter.clusteringFit(), cluIx, probeCluIx );
            BernoulliDistribution prior = priorEval.blockEnablementPrior();
            log_prob_t  lnOddsRatio = prior( enabledAfter ) - prior( enabledBefore )
                                    + llhAfter( enabledAfter ) - llhBefore( enabledBefore );
            lnProb += ln_odds_to_prob( lnOddsRatio );
        }
        else if ( enabledBefore ) {
            // signals transition probability (both clusters enabled)
            const GaussianDistribution& signalPrior = priorEval.signalPrior();
            signal_t signalBefore = cluBefore.signal();
            signal_t signalAfter = cluAfter.signal();
            BOOST_ASSERT( !is_unset( signalBefore ) );
            BOOST_ASSERT( !is_unset( signalAfter ) );
            if ( signalBefore != signalAfter ) {
                log_prob_t lnOddsRatio = signalPrior( signalAfter ) - signalPrior( signalBefore )
                                        + helperAfter.clusteringFit().blockLLH( cluIx, probeCluIx )
                                        - helperBefore.clusteringFit().blockLLH( cluIx, probeCluIx );
                lnProb += ln_odds_to_prob( lnOddsRatio );
            }
        }
        BOOST_ASSERT( !is_unset( lnProb ) );
    }
    // object multiple transition probability (for before objects obly)
    GeometricDistribution multPrior = PriorEval( helperBefore.clusteringFit() ).objectMultiplePrior( objsBefore.size() );
    for ( object_set_t::const_iterator objIt = objsBefore.begin(); objIt != objsBefore.end(); ++objIt ) {
        object_index_t objIx = *objIt;
        size_t multBefore = ((const ChessboardBiclustering&)before).objectMultiple( objIx );
        size_t multAfter = ((const ChessboardBiclustering&)after).objectMultiple( objIx );
        if ( multBefore != multAfter ) {
            for ( probe_clundex_t probeCluIx = 0; probeCluIx < helperBefore.clustering().probesClusters().size(); ++probeCluIx ) {
                ChessboardBiclustering::block_proxy cluBefore = *((const ChessboardBiclustering&)before).findBlock( cluIx, probeCluIx, false );
                ChessboardBiclustering::block_proxy cluAfter = *((const ChessboardBiclustering&)after).findBlock( after.clusterIndex( objIx ), probeCluIx, false );
                bool enabledBefore = cluBefore.isEnabled();
                bool enabledAfter = cluAfter.isEnabled();
                if ( enabledBefore && enabledAfter ) {
                    signal_t signalAfter = cluAfter.signal();
                    BOOST_ASSERT( !is_unset( signalAfter ) );
                    log_prob_t lnOddsRatio = multPrior( multAfter ) - multPrior( multBefore )
                                        + llhEvalBefore.cellsDataLLH( objIx, cluAfter.probesCluster().items(),
                                                                signalAfter, multAfter ) 
                                        - llhEvalAfter.cellsDataLLH( objIx, cluAfter.probesCluster().items(),
                                                                signalAfter, multBefore );
                    BOOST_ASSERT( !is_unset( lnOddsRatio ) );
                    lnProb += ln_odds_to_prob( lnOddsRatio );
                    BOOST_ASSERT( !is_unset( lnProb ) );
                }
            }
        }
    }
    BOOST_ASSERT( !is_unset( lnProb ) );
    return ( lnProb );
}
