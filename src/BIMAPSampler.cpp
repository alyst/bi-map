#include "dynamic_bitset_utils.h"
#include "math/Permutation.h"

#include <boost/format.hpp>

#include "ChessboardBiclusteringFitInternal.h"

#include "BIMAPSampler.h"

const gsl_rng_type* BIMAPSamplerHelper::gslRngType = NULL;

BIMAPSampleCollectorParams::BIMAPSampleCollectorParams()
    : walkSamples( 1000 )
    , priorsStoragePeriod( 5 )
    , reportingPeriod( 100 )
    , maxReportingDelay( 60 )
{
}

gsl_rng* BIMAPSamplerHelper::CreateRndNumGen(
    size_t  rndNumGenSeed
){
    if ( gslRngType == NULL ) {
        gslRngType = gsl_rng_env_setup();
    }
    gsl_rng* res = gsl_rng_alloc( gslRngType );
    gsl_rng_set( res, rndNumGenSeed );
    return ( res );
}

BIMAPSamplerHelper::BIMAPSamplerHelper(
    const PrecomputedData&     precomputed,
    const ChessboardBiclusteringHyperPriors&   hyperpriors,
    const ChessboardBiclusteringPriors&        priors,
    const GibbsSamplerParams&           samplerParams,
    size_t                              rndNumGenSeed,
    TurbineCascadeExecutionMonitor*     pMonitor
) : precomputed( precomputed )
  , hyperpriors( hyperpriors )
  , priors( priors )
  , params( samplerParams )
  , rndNumGen( CreateRndNumGen( rndNumGenSeed ) )
  , dynamicCrossCluFactory( rndNumGen, precomputed, priors, hyperpriors, params )
  , crossoverGenerator( rndNumGen, precomputed, priors, hyperpriors )
{
}

BIMAPSamplerHelper::~BIMAPSamplerHelper()
{
    if ( rndNumGen ) gsl_rng_free( rndNumGen );
}

ChessboardBiclusteringDerivedPriors BIMAPSamplerHelper::randomPriors() const
{
    //Rprintf( "Generating random priors...\n" );
    return ( ChessboardBiclusteringDerivedPriors( hyperpriors.signalHyperprior.generate( rndNumGen ) ) );
}

GeometricDistribution BIMAPSamplerHelper::randomNoiseParams( 
    const ChessboardBiclustering&  clus
) const {
    #if 0
    //Rprintf( "Generating random signal params...\n" );
    double  signalSum = 0.0;
    size_t  signalsCount = 0;
    for ( object_index_t objectIx = 0; objectIx < _data.objectsCount(); objectIx++ ) {
        for ( probe_index_t probeIx = 0; probeIx < _data.probesCount(); probeIx++ ) {
            const OPAProbe& probe = _data.probe( probeIx );
            const bool cellCovered = clus.isBlockEnabled( clus.clusterOfObject( objectIx ), clus.clusterOfProbe( probeIx ) );
            if ( ( cellCovered && covered ) || ( !cellCovered && !covered ) ) {
                for ( assay_set_t::const_iterator assayIt = probe.assays().begin(); assayIt != probe.assays().end(); ++assayIt ) {
                    signalsCount++;
                    signalSum += _data.measurement( objectIx, *assayIt );
                }
            }
        }
    }
    /// @todo check is
    return ( OPACellSignalParams( signalsCount > 0 ? signalSum / signalsCount : 0.0, shape ) );
    #endif
    return ( GeometricDistribution::BySuccessRate( ChessboardBiclusteringPriorEval( data(), priors, clus ).noisePrior()
                                                   .generate( rndNumGen ), 0 ) );
}

#if 0
ChessboardBiclusteringEval BIMAPSampler::createEval( const ChessboardBiclustering& clus ) const
{
    return ( ChessboardBiclusteringEval( rndNumGen(), _data, _priors, clus ) );
}

ChessboardBiclusteringGibbsHelper BIMAPSampler::createGibbsHelper( const ChessboardBiclustering& clus ) const
{
    return ( ChessboardBiclusteringGibbsHelper( rndNumGen(), _data, _priors, _hyperpriors, clus, _signalsCache ) );
}
#endif

/**
    Generate random clustering and signal parameters.
 */
void BIMAPSamplerHelper::randomizeMissingData( ChessboardBiclustering& clus ) const
{
    if ( clus.objectsClusters().size() == 0 ) {
        PitmanYorSample ptn = priors.objectClustering.random( rndNumGen, data().objectsCount() );
        for ( cluster_index_t i = 0; i < ptn.clustersCount(); ++i ) {
            clus.addObjectCluster( ptn[ i ].begin(), ptn[ i ].end() );
        }
    }
    BOOST_ASSERT( clus.checkObjectsPartition() );
    if ( clus.probesClusters().size() == 0 ) {
        PitmanYorSample ptn = priors.probeClustering.random( rndNumGen, data().probesCount() );
        for ( cluster_index_t i = 0; i < ptn.clustersCount(); ++i ) {
            clus.addProbeCluster( ptn[ i ].begin(), ptn[ i ].end() );
        }
    }
    BOOST_ASSERT( clus.checkProbesPartition() );
    if ( clus.enabledBlocksCount() == 0 ) {
        for ( object_clundex_t objCluIx = 0; objCluIx < clus.objectsClusters().size(); objCluIx++ ){
            for ( probe_clundex_t probeCluIx = 0; probeCluIx < clus.probesClusters().size(); probeCluIx++ ){
                clus.setBlock( objCluIx, probeCluIx, gsl_ran_bernoulli( rndNumGen, priors.cellEnablementProb ) );
            }
        }
    }
    ChessboardBiclusteringFit clusFit( precomputed, priors, clus );

    for ( ChessboardBiclustering::const_block_iterator ccIt = clus.begin(); ccIt != clus.end(); ++ccIt ) {
        if ( ccIt->isEnabled() ) {
            generateMissingProbeSignals( clus, ccIt->objectsClusterIndex(), ccIt->probesClusterIndex() );
        }
    }
    BOOST_ASSERT( clus.objectMultiples().size() == clus.objectsCount() );

    // update random basic signals
    clus.setSignalPrior( randomPriors().signalPrior );
    clus.setNoiseParams( randomNoiseParams( clus ) );
    BOOST_ASSERT( clus.checkBlocks() );
}

ChessboardBiclustering BIMAPSamplerHelper::randomClustering() const
{
    ChessboardBiclustering res( data().objectsCount(), data().probesCount() );
    randomizeMissingData( res );
    return ( res );
}

/**
 *  Put it object and probe to separate cluster,
 *  enable bicluster, if non-zero data in the cell.
 */
ChessboardBiclustering BIMAPSamplerHelper::trivialClustering() const
{
    ChessboardBiclustering res( data().objectsCount(), data().probesCount() );
    LOG_DEBUG2( "Putting each probe to separate cluster..." );
    for ( probe_index_t probeIx = 0; probeIx < res.probesCount(); probeIx++ ) {
        res.addProbeCluster( probeIx );
    }
    LOG_DEBUG2( "Putting each object to separate cluster..." );
    for ( object_index_t objIx = 0; objIx < res.objectsCount(); objIx++ ) {
        res.addObjectCluster( objIx );
    }
    res.cleanupClusters();
    LOG_DEBUG2( "Setting initial blocks..." );
    for ( probe_clundex_t probeCluIx = 0; probeCluIx < res.probesCount(); probeCluIx++ ) {
        const OPAProbe& probe = data().probe( probeCluIx );
        for ( object_clundex_t objCluIx = 0; objCluIx < res.objectsCount(); objCluIx++ ) {
            for ( assay_container_t::const_iterator assayIt = probe.assayIndexes().begin();
                    assayIt != probe.assayIndexes().end(); ++assayIt
            ){
                if ( data().measurement( objCluIx, *assayIt ).sc > 0 ) {
                    res.setBlock( objCluIx, probeCluIx, true );
                    break;
                }
            }
        }
    }
    randomizeMissingData( res );
    return ( res );
}

void BIMAPSamplerHelper::generateMissingProbeSignals(
    ChessboardBiclustering&     clus,
    object_clundex_t            objCluIx,
    probe_clundex_t             probeCluIx
) const {
    ChessboardBiclusteringFit fit( precomputed, priors, clus );
    ChessboardBiclusteringGibbsHelper helper( rndNumGen, fit, SamplingTransform(), 0 );
    clus.setBlockSignal( objCluIx, probeCluIx,
                           helper.initialSignal( objCluIx, probeCluIx, NULL ).value ); 
}

BIMAPSampleCollector::BIMAPSampleCollector(
    ChessboardBiclusteringsIndexing&   chessboardBiclusteringsIndexing,
    const BIMAPSampleCollectorParams&  params
) : _walk( chessboardBiclusteringsIndexing )
  , _params( params )
  , _lastSampleReportTime( 0 )
{
}

bool BIMAPSampleCollector::storeSample( double time, turbine_ix_t originIx, const particle_type& particle )
{
    _walk.step( time, originIx, *particle.clustering, particle.metrics );
    if ( ( time > _lastSampleReportTime + _params.maxReportingDelay )
      || ( ( _walk.stepsCount() % _params.reportingPeriod ) == 0 )
    ){
        LOG_INFO( "BIMAP walk(" << ( boost::format( "%.0f" ) % time ) << "s): " << _walk.stepsCount() << " samples collected" );
        _lastSampleReportTime = time;
    }
    if ( ( _walk.stepsCount() % _params.priorsStoragePeriod ) == 0 ) {
        _walk.step( time, originIx, particle.clustering->derivedPriors() );
    }
    return ( _walk.stepsCount() >= _params.walkSamples );
}

ChessboardBiclusteringEnergyEval DynamicChessboardBiclusteringFactory::updateEnergyEval(
    const ChessboardBiclusteringEnergyEval& energyEval,
    std::vector<StatsMetrics>& energyLandscape
) const {
    /// @TODO
    return ( ChessboardBiclusteringEnergyEval() );
}

BIMAPWalk BIMAPSampler_run(
    const BIMAPSamplerHelper&  helper,
    ChessboardBiclusteringsIndexing&   chessboardBiclusteringsIndexing,
    const GibbsSamplerParams&   gibbsSamplerParams,
    const TurbineCascadeParams& eeCascadeParams,
    const BIMAPSampleCollectorParams&  collectorParams,
    const ChessboardBiclustering&      iniClus,
    const TurbineCascadeExecutionMonitor*   pMon
){
    typedef TurbineCascadeUnit<BIMAPSampleCollector, DynamicChessboardBiclusteringFactory, ChessboardBiclusteringCrossoverGenerator> equi_energy_sampler_type;

    BIMAPSampleCollector   collector( chessboardBiclusteringsIndexing, collectorParams );

    equi_energy_sampler_type eeSampler( 0, eeCascadeParams, 
            createTurbineCascadeStructure( eeCascadeParams.levelsCount, 
                                        eeCascadeParams.turbinesCount, 1 ),
                                        helper.rndNumGen,
                                        helper.dynamicCrossCluFactory,
                                        &helper.crossoverGenerator, &collector );
    {
        ChessboardBiclustering iniClu( iniClus );
        helper.randomizeMissingData( iniClu );
        StaticChessboardBiclustering iniClus;
        iniClus.clustering.reset( new ChessboardBiclustering( iniClu ) );
        eeSampler.init( iniClus );
    }
    eeSampler.run();

    return ( collector.walk() );
}
