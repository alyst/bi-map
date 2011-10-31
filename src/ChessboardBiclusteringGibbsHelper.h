#pragma once

#include "BasicTypedefs.h"

#include "array2d.h"
#include "ChessboardBiclusteringFit.h"

#include "mcmc/MetropolisHastingsStep.h"
#include "mcmc/SamlpingTransform.h"
#include "mcmc/TransitionDistributions.h"
#include "mcmc/SingleElementMembershipStep.h"
#include "mcmc/SplitMergeStep.h"

/**
    Transition operator for object counts.
 */
class ObjectMultipleTransition
{
protected:
    const size_t maxObjects;
    prob_t p( size_t value ) const
    {
        return ( (double)( std::min( value - 0.5, maxObjects - 1.0 ) ) / ( maxObjects - 1 ) );
    }
public:
    ObjectMultipleTransition( size_t maxObjects )
    : maxObjects( maxObjects )
    {}

    size_t generate( const gsl_rng* rng, size_t value ) const
    {
        return ( 1 + gsl_ran_binomial( rng, p( value ), maxObjects - 1 ) );
    }

    log_prob_t lnPdf( size_t oldValue, size_t newValue ) const
    {
        return ( log( gsl_ran_binomial_pdf( newValue - 1, p( oldValue ), maxObjects - 1 ) ) );
    }
};

/**
    Helper for making Gibbs steps in MCMC.

    @see ChessboardBiclusteringGibbsSampler
 */
class ChessboardBiclusteringGibbsHelper {
private:
    const gsl_rng*                      _rndNumGen;
    const ChessboardBiclusteringFit&           _fit;
    const size_t                        _objectsSetDistanceThreshold;
    const SamplingTransform             _samplingTransform;
    log_prob_t                          _totalLnP;

public:
    typedef ChessboardBiclustering::const_block_iterator const_block_iterator;
    typedef boost::unordered_set<object_clundex_t> object_cluster_set_type;
    typedef boost::unordered_set<probe_clundex_t> probe_cluster_set_type;

    static const signal_t NOISE_SIGNAL_TRANSITION_SIGMA = 0.1;
    static const signal_t BASELINE_SIGNAL_TRANSITION_SIGMA = 0.1;
    static const signal_t PROBE_SIGNAL_TRANSITION_SIGMA = 0.1;
    static const size_t MAX_OBJECT_COUNTS_PER_CLUSTER = 5;

    class BaseLLH {
    protected:
        const ChessboardBiclusteringFit&   clusFit;

        BaseLLH( const ChessboardBiclusteringFit& clusFit )
        : clusFit( clusFit )
        {}
    };

    class ObjectMultipleDataLLH: public BaseLLH {
        object_index_t              objIx;
        ObjectsClusterParams        params;

    public:
        ObjectMultipleDataLLH( const ChessboardBiclusteringFit& clusFit, object_index_t objIx )
        : BaseLLH( clusFit ), objIx( objIx )
        , params( clusFit.objectsClusterParams( clusFit.clusterOfObject( objIx ) ) )
        {}
        log_prob_t operator()( size_t multiple ) const;
    };

    /**
     *  Calculation of LLH for biclustering block,
     *  when the block is a part of current biclustering.
     */
    class BlockEnablementDataLLHCached: public BaseLLH {
        object_clundex_t    objCluIx;
        probe_clundex_t     probeCluIx;

    public:
        BlockEnablementDataLLHCached( const ChessboardBiclusteringFit& clusFit, object_clundex_t objCluIx, probe_clundex_t probeCluIx )
        : BaseLLH( clusFit ), objCluIx( objCluIx ), probeCluIx( probeCluIx )
        {}

        log_prob_t operator()( bool isEnabled ) const;
    };

    ChessboardBiclusteringGibbsHelper( const gsl_rng* rndNumGen,
                                const ChessboardBiclusteringFit& fit,
                                const SamplingTransform& samplingTransform,
                                log_prob_t totalLnP = unset() );

    const gsl_rng* rndNumGen() const {
        return ( _rndNumGen );
    }

    const ChessboardBiclustering& clustering() const {
        return ( (const ChessboardBiclustering&)_fit );
    }

    const ChessboardBiclusteringFit& clusteringFit() const {
        return ( _fit );
    }

    void probeToClusterProbs( probe_index_t probeIx, probability_vector_t& out );

    GibbsSample<signal_t> randomSignal() const;
    GibbsSample<signal_t> initialSignal( const object_set_t& objects,
                                         const probe_bitset_t& probes,
                                         const signal_t* pCurSignal ) const;
    GibbsSample<signal_t> initialSignal( object_clundex_t objCluIx,
                                         probe_clundex_t probeCluIx,
                                         const signal_t* pCurSignal ) const {
        return ( initialSignal( clustering().objectsCluster( objCluIx ).items(),
                                clustering().probesCluster( probeCluIx ).items(),
                                pCurSignal ) );
    }

    GibbsSample<object_clundex_t> sampleClusterOfObject( object_index_t objIx, const ClusterOfElementStepParams& stepParams,
                                            ObjectsClusterParams* stripeParamsNew = NULL, ObjectsClusterParams* stripeParamsOld = NULL,
                                            const object_cluster_set_type& cluCandidates = object_cluster_set_type() );
    GibbsSample<probe_clundex_t> sampleClusterOfProbe( probe_index_t probeIx, const ClusterOfElementStepParams& stepParams,
                                          ProbesClusterParams* stripeParamsNew = NULL, ProbesClusterParams* stripeParamsOld = NULL,
                                          const probe_cluster_set_type& cluCandidates = probe_cluster_set_type() );

    ObjectMultipleTransition objectMultipleTransition() const {
        return ( ObjectMultipleTransition( MAX_OBJECT_COUNTS_PER_CLUSTER ) );
    }
    GibbsSample<size_t> sampleObjectMultiple( object_index_t objIx );

    GaussianTransitionDistribution signalTransition() const {
        return ( GaussianTransitionDistribution( PROBE_SIGNAL_TRANSITION_SIGMA * _samplingTransform.temperature ) );
    }
    GibbsSample<signal_t> sampleSignal( const object_set_t& objects, const probe_bitset_t& probes, signal_t curSignal ) const;
    GibbsSample<signal_t> sampleSignal( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const;
    GibbsSample<signal_t> sampleSignal( const ChessboardBiclustering::block_proxy& block ) const {
        return ( sampleSignal( block.objectsClusterIndex(),
                               block.probesClusterIndex() ) );
    }

    GibbsSample<bool> sampleBlockEnablement( const object_set_t& objects, const probe_bitset_t& probes, bool curEnabled );
    GibbsSample<bool> sampleBlockEnablement( object_clundex_t objCluIx, probe_clundex_t probeCluIx );

    static std::vector<signal_t> Signals( const ChessboardBiclustering& clustering );
    static std::vector<signal_t> Measurements( const OPAData& data, const ChessboardBiclustering& clustering, bool enabled, bool disabled );

    GaussianDistribution sampleSignalPrior( const ChessboardBiclusteringHyperPriors& hyperpriors ) const;
    GeometricDistribution sampleNoiseParams() const;
};
