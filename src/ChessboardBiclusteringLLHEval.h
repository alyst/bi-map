#pragma once

#include "CellsLLHEval.h"
#include "DataSignalNoiseCache.h"

/**
    Likelihood evaluations for chessboard biclustering.
    Evaluates likelihoods of given signals for
    cell of objects-x-probes matrix or cross-clusters.
 */
class ChessboardBiclusteringLLHEval: public CellsLLHEval {
protected:
    const ChessboardBiclustering&        _clustering;
    DataSignalNoiseCache&         _cache;

public:
    typedef OPAData::hit_counts_type hit_counts_type;

    class BaseLLH {
    protected:
        const ChessboardBiclusteringLLHEval&   eval;

        BaseLLH( const ChessboardBiclusteringLLHEval& eval )
        : eval( eval )
        {}
    };

    class ClusterSignalDataLLH: public BaseLLH {
    protected:
        friend class ChessboardBiclusteringLLHEval;

        const object_set_t&     objects;
        const probe_bitset_t&   probes;

        ClusterSignalDataLLH( const ChessboardBiclusteringLLHEval& eval,
                              object_clundex_t objCluIx,
                              probe_clundex_t probeCluIx )
        : BaseLLH( eval )
        , objects( eval.clustering().objectsCluster( objCluIx ).items() )
        , probes( eval.clustering().probesCluster( probeCluIx ).items() )
        {}

        ClusterSignalDataLLH( const ChessboardBiclusteringLLHEval& eval,
                              const object_set_t&   objects,
                              const probe_bitset_t& probes )
        : BaseLLH( eval ), objects( objects ), probes( probes )
        {}

    public:
        log_prob_t operator()( signal_t signal ) const
        {
            return ( eval.cellsDataLLH( objects, probes, signal ) );
        }
    };

    class AllCellsDataLLH: public BaseLLH {
    protected:
        friend class ChessboardBiclusteringLLHEval;

        bool                        enabledCells;

        AllCellsDataLLH( const ChessboardBiclusteringLLHEval& eval, bool enabledCells )
        : BaseLLH( eval ), enabledCells( enabledCells )
        {}

    public:
        log_prob_t operator()( signal_t signal ) const;
    };

public:
    typedef ChessboardBiclustering::const_cross_cluster_iterator const_cross_cluster_iterator;

    ChessboardBiclusteringLLHEval( DataSignalNoiseCache& cache,
                            const ChessboardBiclustering& clustering );

    ClusterSignalDataLLH signalDataLLHEval( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const
    {
        return ( ClusterSignalDataLLH( *this, objCluIx, probeCluIx ) );
    }
    ClusterSignalDataLLH signalDataLLHEval( const object_set_t& objects, const probe_bitset_t& probes ) const
    {
        return ( ClusterSignalDataLLH( *this, objects, probes ) );
    }

    AllCellsDataLLH allCellsDataLLHEval( bool enabledCells ) const
    {
        return ( AllCellsDataLLH( *this, enabledCells ) );
    }

    log_prob_t cellsDataLLH( const object_set_t& objects, const probe_bitset_t& probes, signal_t signal, const multiple_map_t& objMultiples = multiple_map_t() ) const;
    log_prob_t cellsDataLLH( const object_set_t& objects, probe_index_t probe, signal_t signal, size_t objMultiple ) const;
    log_prob_t cellsDataLLH( object_index_t objectIx, const probe_bitset_t& probes, signal_t signal, size_t objMultiple ) const;

    log_prob_t crossClusterDataLLH( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const;

    log_prob_t allCellsDataLLH( signal_t baselineSignal, prob_t zeroRate, bool sumEnabledCrosses = true, bool sumDisabledCrosses = true ) const;
    log_prob_t allCellsDataLLH( bool sumEnabledCrosses = true, bool sumDisabledCrosses = true ) const;

#if 0
    OPACellSignalParams cellSignalParams( object_index_t objIx, signal_t signal, signal_t lnMultiplier ) const {
        return ( signalParams( baselineObjectSignalParams( objIx, clustering().objectMultiples().find( objIx )->second ), 
                               signal + lnMultiplier ) );
    }
#endif
    log_prob_t cellLLH( const signal_params_type& objectParams, object_index_t objIx, const OPAProbe& probe, signal_t signal ) const;
    log_prob_t cellLLH( object_index_t objIx, const OPAProbe& probe, signal_t signal ) const;
    log_prob_t cellNoiseLLH( const noise_params_type& noiseParams, object_index_t objIx, const OPAProbe& probe ) const;
    log_prob_t cellsNoiseLLH( const noise_params_type& noiseParams, const object_set_t& objects, const probe_bitset_t& probes ) const;

    const ChessboardBiclustering& clustering() const {
        return ( _clustering );
    }
};

class CrossClusterEnablementDataLLH {
    const DataSignalNoiseCache&     cache;
    const object_set_t&             objects;
    const probe_bitset_t&           probes;

public:
    CrossClusterEnablementDataLLH( const DataSignalNoiseCache& cache, const object_set_t& objects, const probe_bitset_t& probes )
    : cache( cache ), objects( objects ), probes( probes )
    {}

    log_prob_t operator()( bool isEnabled ) const;
};

