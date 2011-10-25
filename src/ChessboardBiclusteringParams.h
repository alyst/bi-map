#pragma once

#include "BasicTypedefs.h"


#include "math/PitmanYorProcess.h"
#include "math/Distributions.h"
#include "Signals.h"

/**
    Parameters of prior chessboard biclustering distributions.
 */
class ChessboardBiclusteringPriors {
public:
    /**
        Parameter of the geometric distribution of the number of object of particular type in bicluster
     */
    PitmanYorProcess    objectClustering;
    PitmanYorProcess    probeClustering;

    prob_t  cellEnablementProb;     /** prior probability that element/projection cell contains signal */
    prob_t  probesCluOffObjectsCluRate; /** failure rate for the distribution of the number of probes clusters bound to given objects cluster via bait */
    prob_t  objectsCluOffProbesCluRate; /** failure rate for the distribution of the number of objects clusters bound to given probes cluster via bait */

    BetaDistribution    noise; /** noise distribution rate prior */

    prob_t    objectMultipleRate;   /** failure rate of object's multiple distribution */

    /** cached distributions */
    GeometricDistribution         multipleSingletonPrior;    /** prior for single object in cluster */
    GeometricDistribution         multipleInClusterPrior;    /** prior for clusters with >1 objects */
    BernoulliDistribution         crossClusterEnablementPrior;    /** cross-cluster enabled prior */

    ChessboardBiclusteringPriors( 
                    const PitmanYorProcess& objectClustering = PitmanYorProcess( 0.1, 0.0 ),
                    const PitmanYorProcess& probeClustering = PitmanYorProcess( 0.1, 0.0 ), 
                    prob_t  cellEnablementProb = 0.1,
                    prob_t  probesCluOffObjectsCluRate = 0.1,
                    prob_t  objectsCluOffProbesCluRate = 0.1,
                    prob_t  objectMultipleRate = 0.02,
                    const BetaDistribution& noise = BetaDistribution( 99, 1 )
    ) : objectClustering( objectClustering )
      , probeClustering( probeClustering )
      , cellEnablementProb( cellEnablementProb )
      , probesCluOffObjectsCluRate( probesCluOffObjectsCluRate )
      , objectsCluOffProbesCluRate( objectsCluOffProbesCluRate )
      , noise( noise )
      , objectMultipleRate( objectMultipleRate )
    {
        updateCachedDistributions();
    }

    void updateCachedDistributions()
    {
        crossClusterEnablementPrior = BernoulliDistribution( cellEnablementProb );
        multipleSingletonPrior = GeometricDistribution::ByFailureRate( 0.1 * objectMultipleRate, 1 );
        multipleInClusterPrior = GeometricDistribution::ByFailureRate( objectMultipleRate, 1 );
    }

private:
    friend class boost::serialization::access;

    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> BOOST_SERIALIZATION_NVP( objectClustering );
        ar >> BOOST_SERIALIZATION_NVP( probeClustering );
        ar >> BOOST_SERIALIZATION_NVP( cellEnablementProb );
        ar >> BOOST_SERIALIZATION_NVP( noise );
        ar >> BOOST_SERIALIZATION_NVP( objectMultipleRate );
        updateCachedDistributions();
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << BOOST_SERIALIZATION_NVP( objectClustering );
        ar << BOOST_SERIALIZATION_NVP( probeClustering );
        ar << BOOST_SERIALIZATION_NVP( cellEnablementProb );
        ar << BOOST_SERIALIZATION_NVP( noise );
        ar << BOOST_SERIALIZATION_NVP( objectMultipleRate );
    }
};

BOOST_CLASS_IMPLEMENTATION( ChessboardBiclusteringPriors, object_serializable )
//BOOST_IS_BITWISE_SERIALIZABLE( ChessboardBiclusteringPriors )

/**
    Parameters of prior chessboard biclustering distributions
    sampled using hyperpriors
 */
struct ChessboardBiclusteringDerivedPriors {
    GaussianDistribution signalPrior;

    ChessboardBiclusteringDerivedPriors( 
                    const GaussianDistribution& signalPrior = GaussianDistribution( 0.0, 1.0 )
    ) : signalPrior( signalPrior )
    {
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( signalPrior );
    }
};

BOOST_IS_BITWISE_SERIALIZABLE( ChessboardBiclusteringDerivedPriors )
BOOST_CLASS_VERSION( ChessboardBiclusteringDerivedPriors, 1 );

/**
    Parameters for hyperpriors.
 */
class ChessboardBiclusteringHyperPriors {
public:
    NormalScaledInverseGammaDistribution    signalHyperprior;

    ChessboardBiclusteringHyperPriors( 
        const NormalScaledInverseGammaDistribution& signalHyperprior
        = NormalScaledInverseGammaDistribution( 0.0, 10.0, 1.0, 0.1 )
    ) : signalHyperprior( signalHyperprior )
    {
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( signalHyperprior );
    }
};

BOOST_IS_BITWISE_SERIALIZABLE( ChessboardBiclusteringHyperPriors )

/**
    Basic parameters of the chessboard biclustering model.
 */
struct ChessboardBiclusteringData {
    typedef ObjectsClusterSignal    signal_params_type;
    typedef GeometricDistribution   noise_params_type;

    ChessboardBiclusteringDerivedPriors    _derivedPriors;         /**< priors derived from hyperpriors */
    signal_params_type              _baselineSignalParams;  /**< mean signal from cells (calibration constant + "variance") */
    noise_params_type               _noiseParams;           /**< signal from cells, uncovered by clustering -- i.e. true/false negatives */

    ChessboardBiclusteringData()
    {}

    ChessboardBiclusteringData( const ChessboardBiclusteringDerivedPriors& derivedPriors,
                         const signal_params_type& baselineSignalParams,
                         const noise_params_type& noiseParams
    ) : _derivedPriors( derivedPriors )
      , _baselineSignalParams( baselineSignalParams )
      , _noiseParams( noiseParams )
    {
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp( "derivedPriors", _derivedPriors );
        ar & boost::serialization::make_nvp( "baselineSignalParams", _baselineSignalParams );
        ar & boost::serialization::make_nvp( "noiseParams", _noiseParams );
    }
};

BOOST_CLASS_VERSION( ChessboardBiclusteringData, 1 )
//BOOST_IS_BITWISE_SERIALIZABLE( ChessboardBiclusteringData )
