#pragma once

#include "BasicTypedefs.h"

#include "math/Distributions.h"
#include "ChessboardBiclustering.h"
#include "OPAData.h"

/**
    Evaluations for chessboard biclustering.
    Evaluates likelihoods, priors etc.
 */
class ChessboardBiclusteringPriorEval {
private:
    const OPAData&                      _data;
    const ChessboardBiclusteringPriors&        _priors;
    const ChessboardBiclustering&              _clustering;

public:
    typedef ChessboardBiclustering::const_cross_cluster_iterator const_cross_cluster_iterator;

    ChessboardBiclusteringPriorEval( const OPAData& data,
                              const ChessboardBiclusteringPriors& priors,
                              const ChessboardBiclustering& clustering );

    log_prob_t lpp() const;

    const BetaDistribution& noisePrior() const {
        return ( _priors.noise );
    }

    const OPAData& data() const {
        return ( _data );
    }

    const ChessboardBiclustering& clustering() const {
        return ( _clustering );
    }

    const ChessboardBiclusteringPriors& priors() const {
        return ( _priors );
    }

    log_prob_t objectsClusteringLPP() const;
    log_prob_t probesClusteringLPP() const;

    log_prob_t objectsMultiplesLPP() const;
    log_prob_t objectMultipleLPP( size_t clusterSize, size_t multiple ) const {
        return ( objectMultiplePrior( clusterSize )( multiple ) );
    }
    log_prob_t noiseParamsLPP() const;

    log_prob_t signalParamsLPP( signal_t signal ) const {
        const GaussianDistribution& signalPrior = clustering().derivedPriors().signalPrior;
        return ( signalPrior( signal ) );
    }
    log_prob_t crossClusterLPP( const ChessboardBiclustering::cross_cluster_proxy& crossCluster ) const;
    log_prob_t crossClusterLPP( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const;

    /**
     *  Prior that also includes current signal prior.
     */
    struct CrossClusterEnablementPrior {
        log_prob_t  lppEnabled;
        log_prob_t  lppDisabled;

        CrossClusterEnablementPrior( 
            const ChessboardBiclusteringPriors&        priors1,
            const ChessboardBiclusteringDerivedPriors& priors2,
            signal_t                            signal );

        log_prob_t operator()( bool isEnabled ) const
        {
            return ( isEnabled ? lppEnabled : lppDisabled );
        }
    };

    CrossClusterEnablementPrior crossClusterEnablementPrior( signal_t signal ) const {
        return ( CrossClusterEnablementPrior( _priors, _clustering.derivedPriors(), signal ) );
    }

    const BernoulliDistribution& crossClusterEnablementPrior() const {
        return ( _priors.crossClusterEnablementPrior );
    }

    const GeometricDistribution& objectMultiplePrior( size_t objects ) const {
        return ( objects > 1
                 ? _priors.multipleInClusterPrior
                 : _priors.multipleSingletonPrior );
    }
    const GaussianDistribution& signalPrior() const {
        return ( clustering().derivedPriors().signalPrior );
    }
    size_t countMeasurements( bool signal ) const;
};
