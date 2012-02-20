#pragma once

#include "BasicTypedefs.h"

#include <cemm/containers/LazyClone.h>

#include "ChessboardBiclusteringGibbsSampler.h"

namespace cemm { namespace bimap {

/**
    Realization of Partition concept for projections.
    Read-only part.

    @see ProbesPartitionEx, ObjectsPartition
 */
class ProbesPartition: public LazyClone<ChessboardBiclusteringFit> {
public:
    typedef probe_index_t element_index_type;
    typedef probe_bitset_t element_index_set_type;

    typedef probe_clundex_t cluster_index_type;
    typedef boost::unordered_set<cluster_index_type> cluster_index_set_type;

    typedef ProbesClusterParams cluster_params_type;
    typedef ChessboardBiclusteringFit::ProbesClusterParamsProxy cluster_params_proxy_type;

protected:
    typedef LazyClone<ChessboardBiclusteringFit> super_type;

public:
    /** Proxy (wrapper) for probe_bitset_t */
    class ProbesSetProxy: public LazyClone<probe_bitset_t> {
    protected:
        typedef LazyClone<probe_bitset_t> super_type;
    public:
        ProbesSetProxy( const probe_bitset_t& probes )
        : super_type( probes )
        {
        }
        ProbesSetProxy( const ProbesPartition& ptn, element_index_type elmIx = PROBE_NA )
        : super_type( boost::make_shared<probe_bitset_t>( ptn.elementsCount() ), false )
        {
            if ( elmIx != PROBE_NA ) wrapped().set( elmIx );
        }
        bool contains( element_index_type elmIx ) const {
            return ( wrapped().test( elmIx ) );
        }
        element_index_type operator[]( size_t index ) const
        {
            element_index_type elmIx = wrapped().find_first();
            if ( elmIx == probe_bitset_t::npos ) return ( PROBE_NA );
            while ( index > 0 ) {
                elmIx = wrapped().find_next( elmIx );
                if ( elmIx == probe_bitset_t::npos ) return ( PROBE_NA );
                index--;
            }
            return ( elmIx );
        }
        size_t size() const {
            return ( wrapped().count() );
        }
        ProbesSetProxy& operator+=( const ProbesSetProxy& that )
        {
            unshare();
            const_object_ptr_type thatpProbes = that.lock();
            wrapped() |= *thatpProbes;
            return ( *this );
        }
        ProbesSetProxy& operator+=( element_index_type elmIx )
        {
            unshare();
            wrapped().set( elmIx );
            return ( *this );
        }
        ProbesSetProxy& operator-=( const ProbesSetProxy& that )
        {
            unshare();
            const locked_ptr_type thatpProbes = that.lock();
            wrapped() -= *thatpProbes;
            return ( *this );
        }
        ProbesSetProxy& operator-=( element_index_type elmIx )
        {
            unshare();
            wrapped().set( elmIx, false );
            return ( *this );
        }
        ProbesSetProxy& operator&=( const ProbesSetProxy& that )
        {
            unshare();
            const locked_ptr_type thatpProbes = that.lock();
            wrapped() &= *thatpProbes;
            return ( *this );
        }
        bool operator==(const ProbesSetProxy& that ) const {
            const locked_ptr_type thatpProbes = that.lock();
            return ( wrapped() == *thatpProbes );
        }
        bool operator!=(const ProbesSetProxy& that ) const {
            const locked_ptr_type thatpProbes = that.lock();
            return ( wrapped() != *thatpProbes );
        }
    };

    typedef ProbesSetProxy elements_set_proxy_type;

    class ProbesClusterProxy {
    private:
        friend class ProbesPartition;
        friend class ProbesPartitionEx;

        ChessboardBiclusteringFit& clus;
        const probe_clundex_t cluIx;

        const ProbesCluster& cluster() const {
            return ( clus.probesCluster( cluIx ) );
        }

        ProbesClusterProxy( ChessboardBiclusteringFit& clus, probe_clundex_t cluIx )
        : clus( clus ), cluIx( cluIx )
        {
            BOOST_ASSERT( (int)cluIx >= 0 );
            BOOST_ASSERT( cluIx < clus.probesClusters().size() );
        }

    public:
        cluster_index_type index() const {
            return ( cluIx );
        }
        std::string label() const;
        size_t size() const {
            return ( cluster().size() );
        }
        const elements_set_proxy_type items() const {
            return ( ProbesSetProxy( cluster().items() ) );
        }
        const cluster_params_proxy_type params(
            const probe_bitset_t& mask = probe_bitset_t()
        ) const {
            return ( clus.probesClusterParams( cluIx, mask ) );
        }
        cluster_params_proxy_type params(
            const probe_bitset_t& mask = probe_bitset_t()
        ){
            return ( clus.probesClusterParams( cluIx, mask ) );
        }
    };

    typedef ProbesClusterProxy cluster_proxy_type;

    const static cluster_index_type ClusterNA = CLUSTER_NA;

    ProbesPartition( const ChessboardBiclusteringFit& clus )
    : super_type( clus )
    {}

    ProbesPartition( const shared_object_ptr_type& pClus )
    : super_type( pClus, true )
    {}

    virtual ~ProbesPartition()
    {}

    size_t clustersCount() const {
        return ( wrapped().probesClusters().size() );
    }

    size_t elementsCount() const {
        return ( wrapped().probesCount() );
    }

    const cluster_proxy_type cluster( cluster_index_type cluIx ) const {
        return ( cluster_proxy_type( const_cast<ChessboardBiclusteringFit&>( wrapped() ), cluIx ) );
    }

    cluster_index_type clusterIndex( element_index_type elmIx ) const {
        return ( wrapped().clusterOfProbe( elmIx ) );
    }

    operator const ChessboardBiclustering&() const {
        return ( (const ChessboardBiclustering&)wrapped() );
    }
};

/**
    Realization of Partition model for projections.
 */
class ProbesPartitionEx: public ProbesPartition {
public:
    ProbesPartitionEx( const ChessboardBiclusteringFit& clus )
    : ProbesPartition( clus )
    {}

    ProbesPartitionEx( const shared_object_ptr_type& pClus )
    : ProbesPartition( pClus )
    {}

    using ProbesPartition::cluster;

    cluster_proxy_type cluster( cluster_index_type cluIx ) {
        unshare();
        return ( cluster_proxy_type( wrapped(), cluIx ) );
    }

    void putToCluster( 
        element_index_type              elmIx,
        cluster_index_type              cluIx
    );

    cluster_index_type exchangeElements( 
        cluster_index_type              clu1Ix, 
        cluster_index_type              clu2Ix,
        const element_index_set_type&   clu2Elements 
    );

    void setClusterSamples( cluster_index_type cluIx, size_t samples ) {
        unshare();
        wrapped().setProbesClusterSamples( cluIx, samples );
    }
};

/**
    Sampler of objects' parameters for the fixed chessboard biclustering.

    Realization of ParamsSampler model.
 */
class FixedProbesParamsSampler {
public:
    typedef ProbesClusterParams params_type;

    /**
        Default parameters.
     */
    params_type defaults() const {
        return ( params_type( clusHelper->clustering() ) );
    }
    bool operator()( params_type& params, const probe_bitset_t& clusterProbes, const probe_bitset_t& sampledProbes, 
                     bool overwrite, bool posterior ) const;

    FixedProbesParamsSampler( ChessboardBiclusteringGibbsHelper& clusHelper, bool sampleBlockMask, bool sampleSignals )
    : clusHelper( &clusHelper ), sampleBlockMask( sampleBlockMask ), sampleSignals( sampleSignals )
    {}

private:
    ChessboardBiclusteringGibbsHelper* const   clusHelper;
    bool                                sampleBlockMask;
    bool                                sampleSignals;
};

/**
    Sampler of projections parameters for given chessboard biclustering.

    Realization of ParamsSampler model.
    Internally uses FixedObjectsParamsSampler.
*/
class ProbesParamsSampler {
public:
    typedef ProbesClusterParams params_type;

    params_type defaults( const ProbesPartition& ptn ) const {
        return ( params_type( ptn ) );
    }
    bool operator()( const ProbesPartition& ptn, params_type& params, 
                     const probe_bitset_t& clusterProbes, const probe_bitset_t& sampledProbes, 
                     bool overwrite, bool posterior ) const;

    log_prob_t transitionLP( const ProbesPartition& before,
                             const ProbesPartition& after,
                             probe_clundex_t cluIx ) const;

    ProbesParamsSampler( ChessboardBiclusteringGibbsSampler& clusSampler, bool sampleBlockMask, bool sampleSignals )
    : clusSampler( &clusSampler ), sampleBlockMask( sampleBlockMask ), sampleSignals( sampleSignals )
    {}

private:
    ChessboardBiclusteringGibbsSampler const*  clusSampler;
    bool                                sampleBlockMask;
    bool                                sampleSignals;
};


struct FixedProbesPartitionStats {
    typedef ProbesClusterParams params_type;

    /**
     *  Likelihood of objects cluster in the partition,
     *  not concerning its relation to other clusters
     *  (co-occurence, block probes and abundances, object multipliers).
     */
    LLHMetrics probesLLH( const probe_bitset_t& probes, const params_type& params ) const;

    /**
     *  Likelihood of modified projections clusters,
     *  (co-occurence of projections, relation to other clusters,
     *   block probes and abundances).
     */
    LLHMetrics llhDelta( const std::vector<ProbesPartition::elements_set_proxy_type>& newClusters,
                const std::vector<params_type>& newParams,
                const ProbesPartition::cluster_index_set_type& oldIndexes ) const;

    LLHMetrics llh() const {
        return ( clusFit.metrics().llhProbes );
    }
    log_prob_t paramsLPP( const params_type& params, const probe_bitset_t& clusterProbes, const probe_bitset_t& sampledProbes ) const;

    FixedProbesPartitionStats( const ChessboardBiclusteringFit& clusFit )
    : clusFit( clusFit )
    {}

    log_prob_t clusterAssignmentLPP( size_t clusterSize, size_t clustersCount, size_t totalElements ) const
    {
        return ( log( clusFit.priors().probeClustering.clusterAssignmentPrior( clusterSize, clustersCount, clusFit.probesCount() - 1 ) ) );
    }

private:
    const ChessboardBiclusteringFit&  clusFit;
};

struct ProbesPartitionStats {
    typedef ProbesClusterParams params_type;

    log_prob_t llhDelta( const ProbesPartition& ptn,
                const std::vector<ProbesPartition::elements_set_proxy_type>& newClusters,
                const std::vector<params_type>& newParams,
                const ProbesPartition::cluster_index_set_type& oldIndexes
    ) const {
        return ( FixedProbesPartitionStats( ptn ).llhDelta( newClusters, newParams, oldIndexes )( weights ) );
    }
    log_prob_t llh( const ProbesPartition& ptn ) const {
        const ChessboardBiclusteringFit& clusFit = ptn;
        return ( clusFit.metrics().llhProbes( weights ) );
    }
    log_prob_t lpp( const ProbesPartition& ptn ) const {
        const ChessboardBiclusteringFit& clusFit = ptn;
        return ( clusFit.lpp() );
    }
    log_prob_t paramsLPP( const ProbesPartition& ptn, const params_type& params, const probe_bitset_t& clusterProbes, const probe_bitset_t& sampledProbes ) const {
        return ( FixedProbesPartitionStats( ptn ).paramsLPP( params, clusterProbes, sampledProbes ) );
    }

    log_prob_t clusterAssignmentLPP( size_t clusterSize, size_t clustersCount, size_t totalElements ) const
    {
        return ( log( priors.probeClustering.clusterAssignmentPrior( clusterSize, clustersCount, totalElements - 1 ) ) );
    }

    ProbesPartitionStats(
        const gsl_rng*                  rndNumGen,
        const OPAData&                  data,
        const ChessboardBiclusteringPriors&    priors,
        const LLHPartitionWeights&      weights
    ) : rndNumGen( rndNumGen ), data( data )
    , priors( priors ), weights( weights )
    {}

private:
    const gsl_rng*                      rndNumGen;
    const OPAData&                      data;
    const ChessboardBiclusteringPriors& priors;
    const LLHPartitionWeights&          weights; 
};

} }