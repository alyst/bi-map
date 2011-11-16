#include "dynamic_bitset_utils.h"

#include "mcmc/MetropolisHastingsStep.h"
#include "mcmc/SingleElementMembershipStep.h"

#include "ChessboardBiclusteringFitInternal.h"
#include "ObjectsPartition.h"
#include "ProbesPartition.h"

#include "ChessboardBiclusteringGibbsHelper.h"

ChessboardBiclusteringGibbsHelper::ChessboardBiclusteringGibbsHelper(
    const gsl_rng*                      rndNumGen,
    const ChessboardBiclusteringFit&           clusteringFit,
    const SamplingTransform&            samplingTransform,
    log_prob_t                          totalLnP
) : _rndNumGen( rndNumGen )
  , _fit( clusteringFit )
  , _objectsSetDistanceThreshold( 1 + clusteringFit.data().objectsCount() / 10 )
  , _samplingTransform( samplingTransform )
  , _totalLnP( is_unset( totalLnP ) ? _fit.totalLnP() : totalLnP )
{
}

log_prob_t ChessboardBiclusteringGibbsHelper::BlockEnablementDataLLHCached::operator()(
    bool isEnabled
) const {
    log_prob_t res = isEnabled
               ? clusFit.blockIsSignalLLH( objCluIx, probeCluIx )
               : clusFit.blockIsNoiseLLH( objCluIx, probeCluIx );
    LOG_DEBUG2( "CC_" << isEnabled << "[" << objCluIx << ", " << probeCluIx << "]=" << res );
    return ( res );
}

GibbsSample<bool> ChessboardBiclusteringGibbsHelper::sampleBlockEnablement(
    object_clundex_t    objCluIx,
    probe_clundex_t     probeCluIx
){
    ChessboardBiclustering::const_block_iterator ccIt = _fit.findBlock( objCluIx, probeCluIx, false );
    return ( MetropolisHastringsPosteriorSample<bool>(
        rndNumGen(), 
        BinaryTransition(),
        BlockEnablementDataLLHCached( _fit, objCluIx, probeCluIx ),
        PriorEval( _fit ).blockEnablementPrior( 
                ccIt->isEnabled() ? ccIt->signal()
                : _fit.precomputed().blockSignal( ccIt->objectsCluster().items(), 
                                                    ccIt->probesCluster().items(),
                                                    _fit.objectMultiples()
                                                  ) ),
        ccIt->isEnabled(), _totalLnP, _samplingTransform ) );
}

GibbsSample<bool> ChessboardBiclusteringGibbsHelper::sampleBlockEnablement(
    const object_set_t&     objects,
    const probe_bitset_t&   probes,
    bool                    curEnabled
){
    return ( MetropolisHastringsPosteriorSample<bool>( 
        rndNumGen(), 
        BinaryTransition(),
        BlockEnablementDataLLH( _fit.signalNoiseCache(), objects, probes ),
        PriorEval( _fit ).blockEnablementPrior( _fit.precomputed().blockSignal( objects, probes, _fit.objectMultiples() ) ),
        curEnabled, _totalLnP, _samplingTransform ) );
}

GibbsSample<object_clundex_t> ChessboardBiclusteringGibbsHelper::sampleClusterOfObject( 
    object_index_t                      objIx,
    const ClusterOfElementStepParams&   stepParams,
    ObjectsClusterParams*               stripeParamsNew,
    ObjectsClusterParams*               stripeParamsOld,
    const object_cluster_set_type&      cluCandidates
){
    return ( SampleClusterOfElement( rndNumGen(), 
                                       ObjectsPartition( _fit ), 
                                       ObjectsPartitionStats( rndNumGen(), _fit.data(), _fit.priors() ),
                                       objIx, 
                                       FixedObjectsParamsSampler( *this, true, true, true ),
                                       stepParams, _samplingTransform, _totalLnP,
                                       stripeParamsNew, stripeParamsOld, cluCandidates ) );
}

GibbsSample<probe_clundex_t> ChessboardBiclusteringGibbsHelper::sampleClusterOfProbe( 
    probe_index_t                       probeIx,
    const ClusterOfElementStepParams&   stepParams,
    ProbesClusterParams*                stripeParamsNew,
    ProbesClusterParams*                stripeParamsOld,
    const probe_cluster_set_type&       cluCandidates
){
    return ( SampleClusterOfElement( rndNumGen(),
                                       ProbesPartition( _fit ),
                                       ProbesPartitionStats( rndNumGen(), _fit.data(), _fit.priors() ),
                                       probeIx,
                                       FixedProbesParamsSampler( *this, true, true ),
                                       stepParams, _samplingTransform, _totalLnP,
                                       stripeParamsNew, stripeParamsOld, cluCandidates ) );
}

log_prob_t ChessboardBiclusteringGibbsHelper::ObjectMultipleDataLLH::operator()(
    size_t multiple
) const {
    params.objectMultiple[ objIx ] = multiple;
    return ( FixedObjectsPartitionStats( clusFit ).objectLLH( objIx, params ).total() );
}

GibbsSample<size_t> ChessboardBiclusteringGibbsHelper::sampleObjectMultiple(
    object_index_t  objIx
){
    return ( MetropolisHastringsPosteriorSample<size_t>( 
        rndNumGen(), 
        objectMultipleTransition(),
        ObjectMultipleDataLLH( _fit, objIx ),
        PriorEval( _fit ).objectMultiplePrior( _fit.objectsCluster( _fit.clusterOfObject( objIx ) ).size() ),
        _fit.objectMultiple( objIx ), _totalLnP, _samplingTransform ) );
}

GibbsSample<signal_t> ChessboardBiclusteringGibbsHelper::sampleSignal(
    const object_set_t&     objects,
    const probe_bitset_t&   probes,
    signal_t                curSignal
) const {
    BOOST_ASSERT( !is_unset( curSignal ) );
    GibbsSample<signal_t> res = MetropolisHastringsPosteriorSample<signal_t>( 
            rndNumGen(), 
            signalTransition(),
            LLHEval( _fit ).signalDataLLHEval( objects, probes ),
            PriorEval( _fit ).signalPrior(),
            curSignal, _totalLnP, _samplingTransform );

    return ( res );
}

GibbsSample<signal_t> ChessboardBiclusteringGibbsHelper::sampleSignal(
    const object_clundex_t  objCluIx,
    const probe_clundex_t   probeCluIx
) const {
    signal_t curSignal = clustering().blockSignal( objCluIx, probeCluIx );
    return ( !is_unset( curSignal )
             ? sampleSignal( clustering().objectsCluster( objCluIx ).items(),
                             clustering().probesCluster( probeCluIx ).items(),
                             curSignal )
             : initialSignal( objCluIx, probeCluIx, NULL ) );
}

GibbsSample<signal_t> ChessboardBiclusteringGibbsHelper::randomSignal() const
{
    double randomSignal = _fit.derivedPriors().signalPrior.generate( rndNumGen() );
    return ( GibbsSample<signal_t>( randomSignal, 0, _fit.derivedPriors().signalPrior( randomSignal ) ) );
}

/**
 *  Number of evidence samples used to obtain posterior signal distribution.
 */
#define FAKE_SAMPLES ((size_t)3u)

GibbsSample<signal_t> ChessboardBiclusteringGibbsHelper::initialSignal(
    const object_set_t&                 objects,
    const probe_bitset_t&               probes,
    const signal_t*                     pCurSignal
) const {
#if 1
    size_t samplesCnt;
    signal_t curSignal = std::numeric_limits<double>::quiet_NaN();
    if ( pCurSignal ) {
        curSignal = *pCurSignal;
        samplesCnt = FAKE_SAMPLES;
    }
    else if ( !objects.empty() ) {
        signal_t signal = _fit.precomputed().blockSignal( objects, probes, _fit.objectMultiples() );
        _fit.precomputed().maximizeSignalLLH( objects, probes, _fit.objectMultiples(), signal );
        //LOG_INFO( "Generated signal for (" << *objects.begin() << ", " << probeIx << ")=" << signal );
        return ( GibbsSample<signal_t>( signal, 0, _fit.derivedPriors().signalPrior( signal ) ) );
    }
    if ( !is_unset( curSignal ) ) {
        // generate new signal with reduced variability based on signal in cache
        // (kind of normal distribution posterior)
        GaussianDistribution signalPosterior = _fit.derivedPriors().signalPrior
            .posterior( curSignal, samplesCnt );
        signal_t signal = signalPosterior.generate( rndNumGen() );
        //LOG_INFO( "Generated signal for (" << *objects.begin() << ", " << probeIx << ")=" << signal );
        return ( GibbsSample<signal_t>( signal, 0, signalPosterior( signal ) ) );
    }
    else {
        return ( randomSignal() );
    }
#else
    return ( randomSignal() );
#endif
}

std::vector<signal_t> ChessboardBiclusteringGibbsHelper::Signals(
    const ChessboardBiclustering& clustering
){
    std::vector<signal_t> signals;
    signals.reserve( clustering.probesCount() * clustering.objectsClusters().size() / 2 );

    LOG_DEBUG2( "Signals for transform (E=" << _samplingTransform.maxLnP 
                << " T=" << _samplingTransform.temperature << ")" );
    for ( probe_clundex_t probeCluIx = 0; probeCluIx < clustering.probesClusters().size(); probeCluIx++ ) {
        for ( object_clundex_t objCluIx = 0; objCluIx < clustering.objectsClusters().size(); objCluIx++ ) {
            const ChessboardBiclustering::const_block_iterator cluIt = clustering.findBlock( objCluIx, probeCluIx );
            if ( cluIt != clustering.blockNotFound() ) {
                signal_t signal = clustering.blockSignal( objCluIx, probeCluIx );
                if ( !is_unset( signal ) ) {
                    LOG_DEBUG2( "Signal " << signal );
                    signals.push_back( signal );
                }
            }
        }
    }

    return ( signals );
}

std::vector<signal_t> ChessboardBiclusteringGibbsHelper::Measurements(
    const OPAData&         data,
    const ChessboardBiclustering& clustering,
    bool                   enabled,
    bool                   disabled
){
    std::vector<signal_t> measurements;
    measurements.reserve( data.assaysCount() * clustering.objectsCount() / 2 );

    for ( probe_clundex_t probeCluIx = 0; probeCluIx < clustering.probesClusters().size(); probeCluIx++ ) {
        const ProbesCluster& probesClu = clustering.probesCluster( probeCluIx );
        for ( object_clundex_t objCluIx = 0; objCluIx < clustering.objectsClusters().size(); objCluIx++ ) {
            const ChessboardBiclustering::const_block_iterator cluIt = clustering.findBlock( objCluIx, probeCluIx );
            bool isCCOn = cluIt != clustering.blockNotFound();
            if ( ( isCCOn && enabled ) || ( !isCCOn && disabled ) ) {
                foreach_bit( probe_index_t, probeIx, probesClu.items() ) {
                   const OPAProbe& probe = data.probe( probeIx );
                   for ( assay_index_t assayIx = probe.assayIndexes().front(); assayIx <= probe.assayIndexes().back(); ++assayIx ) {
                       for ( object_set_t::const_iterator oit = clustering.objectsCluster( objCluIx ).items().begin();
                             oit != clustering.objectsCluster( objCluIx ).items().end(); ++oit ) {
                             measurements.push_back( data.measurement( *oit, assayIx ).sc );
                       }
                   }
                }
            }
        }
    }

    return ( measurements );
}

GaussianDistribution ChessboardBiclusteringGibbsHelper::sampleSignalPrior(
    const ChessboardBiclusteringHyperPriors& hyperpriors
) const {
    return ( hyperpriors.signalHyperprior.posterior( Signals( clustering() ) ).generate( rndNumGen() ) );
}

GeometricDistribution ChessboardBiclusteringGibbsHelper::sampleNoiseParams() const
{
    std::vector<signal_t> signals = Measurements( _fit.data(), clustering(), false, true );
    // transform multiple signals to single to avoid degeneration of distribution
    for ( size_t i = 0; i < signals.size(); i++ ) {
        signals[ i ] = signals[ i ] > 0 ? 1.0 : 0.0;
    }
    return ( GeometricDistribution::BySuccessRate( _fit.priors().noise.posterior( signals, 0 )
                                                   .generate( rndNumGen() ) ) );
}
