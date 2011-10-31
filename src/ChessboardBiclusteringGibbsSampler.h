#pragma once

#include "BasicTypedefs.h"

#include "mcmc/SamlpingTransform.h"

#include "ChessboardBiclustering.h"
#include "ChessboardBiclusteringLLHEval.h"
#include "ChessboardBiclusteringGibbsHelper.h"

/**
    Parameters for ChessboardBiclusteringGibbsSampler.
 */
struct GibbsSamplerParams {
    size_t    priorsUpdatePeriod;           /** period (# of iterations) at which sampling step is made for prior parameters */

    size_t    blockResamples;        /** # of sampling steps to make after block elements where modified */

    SplitMergeStepParams    objectsSplitMergeParams;
    ClusterOfElementStepParams  objectClusterParams;
    prob_t    objectsSplitMergeRate;        /** rate of applying split-merge sampling step to objects partition */
    prob_t    objectMembershipRate;         /** rate of sampling single object membership */

    SplitMergeStepParams    probesSplitMergeParams;
    ClusterOfElementStepParams  probeClusterParams;
    prob_t    probesSplitMergeRate;         /** rate of applying split-merge sampling step to probes partition */
    prob_t    probeMembershipRate;          /** rate of sampling single probe membership */

    prob_t    blockFlipRate;         /** rate of sampling on/off flag of block */

    prob_t    objectMultipleRate;           /** rate of sampling object's multiple */
    prob_t    signalRate;                   /** rate of sampling signal of enabled block */

    size_t    meanObjectRank;
    size_t    meanProbeRank;

    GibbsSamplerParams();

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( priorsUpdatePeriod );
        ar & BOOST_SERIALIZATION_NVP( blockResamples );
        ar & BOOST_SERIALIZATION_NVP( objectsSplitMergeParams );
        ar & BOOST_SERIALIZATION_NVP( objectClusterParams );
        ar & BOOST_SERIALIZATION_NVP( objectsSplitMergeRate );
        ar & BOOST_SERIALIZATION_NVP( objectMembershipRate );
        ar & BOOST_SERIALIZATION_NVP( probesSplitMergeParams );
        ar & BOOST_SERIALIZATION_NVP( probeClusterParams );
        ar & BOOST_SERIALIZATION_NVP( probesSplitMergeRate );
        ar & BOOST_SERIALIZATION_NVP( probeMembershipRate );
        ar & BOOST_SERIALIZATION_NVP( blockFlipRate );
        ar & BOOST_SERIALIZATION_NVP( objectMultipleRate );
        ar & BOOST_SERIALIZATION_NVP( signalRate );
        ar & BOOST_SERIALIZATION_NVP( meanObjectRank );
        ar & BOOST_SERIALIZATION_NVP( meanProbeRank );
    }
};

/**
    MCMC chessboard biclustering sampler with Gibbs steps.
    Invoked from Equi-Energy Sampler.
 */
class ChessboardBiclusteringGibbsSampler {
private:
    typedef boost::unordered_set<size_t> cluster_set_type;

    const gsl_rng*                  _rndNumGen;

    const PrecomputedData&          _precomputed;
    const ChessboardBiclusteringHyperPriors&   _hyperpriors;
    const GibbsSamplerParams&       _params;

    SamplingTransform               _samplingTransform;

    size_t                          _iteration;

public:
    ChessboardBiclusteringGibbsSampler( 
                   const gsl_rng* rndNumGen,
                   const PrecomputedData&   precomputed,
                   const ChessboardBiclusteringHyperPriors& hyperpriors, 
                   const GibbsSamplerParams& params,
                   const SamplingTransform& transform = SamplingTransform() );

    ChessboardBiclusteringGibbsHelper createGibbsHelper( const ChessboardBiclusteringFit& clusFit ) const;

    void setTransform( const SamplingTransform& transform )
    {
        _samplingTransform = transform;
    }

    probability_vector_t elementsPickupRates( const ChessboardBiclusteringFit& clusFit, bool useObjects, double MaxRateRatio = 10.0 ) const;

    object_index_t randomObjectToModify( const ChessboardBiclusteringFit& clusFit, double MaxRateRatio = 10.0 ) const;
    probe_index_t randomProbeToModify( const ChessboardBiclusteringFit& clusFit, double MaxRateRatio = 10.0 ) const;

    log_prob_t doSamplingStep( ChessboardBiclusteringFit& clusFit );

    const gsl_rng* rndNumGen() const {
        return ( _rndNumGen );
    }

    bool isPriorsUpdateRequired() const {
        return ( ( ( _iteration + 7 ) % _params.priorsUpdatePeriod ) == 0 );
    }
};
