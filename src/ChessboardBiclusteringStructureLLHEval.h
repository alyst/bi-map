#pragma once

#include "BasicTypedefs.h"

#include "math/Distributions.h"
#include "ChessboardBiclustering.h"
#include "PrecomputedData.h"

/**
    Evaluations for chessboard biclustering.
    Evaluates likelihoods.
 */
class ChessboardBiclusteringStructureLLHEval {
protected:
    const prob_t                        _lnPriorCCProbeRatio;   /** ratio of disabled/enabled probabilities */
    const GeometricDistribution         _probeCluPerObjClu;
    const GeometricDistribution         _objCluPerProbeClu;

    const PrecomputedData&              _precomputed;
    const ChessboardBiclustering&              _clustering;

public:
    typedef OPAData::hit_counts_type hit_counts_type;
    typedef array2d<prob_t> lnprob_matrix_type;
    typedef lnprob_matrix_type::section_type lnprob_section_type;
    typedef boost::dynamic_bitset<> observations_mask_type;

public:
    typedef ChessboardBiclustering::const_block_iterator const_block_iterator;

    ChessboardBiclusteringStructureLLHEval( const PrecomputedData& precomputed,
                                     prob_t priorCCEnabledProb,
                                     prob_t probeCluPerObjCluRate,
                                     prob_t objCluPerProbeCluRate,
                                     const ChessboardBiclustering& clustering );

    static log_prob_t ClustersPatternMismatchLLH(
        const observations_mask_type&   clu1mask,
        const observations_mask_type&   clu2mask
    );

    log_prob_t objectsClusterMismatchLLH( const object_set_t& objects ) const;
    log_prob_t probesClusterMismatchLLH( const probe_bitset_t& probes ) const;

    const OPAData& data() const {
        return ( _precomputed.data() );
    }

    const ChessboardBiclustering& clustering() const {
        return ( _clustering );
    }

    observations_mask_type compactProbesMask( const observations_mask_type& probesMask ) const;
    observations_mask_type compactObjectsMask( const observations_mask_type& objectsMask ) const;

    log_prob_t probeClustersPerObjectClusterLLH( size_t probeClustersCount ) const {
        return ( data().isMapBaitsToObjects() && probeClustersCount > 0
                 ? _probeCluPerObjClu( probeClustersCount ) : 0 );
    }
    log_prob_t objectClustersPerProbeClusterLLH( size_t objectClustersCount ) const {
        return ( data().isMapBaitsToObjects() && objectClustersCount > 0
                 ? _objCluPerProbeClu( objectClustersCount ) : 0 );
    }
};
