#pragma once

#include "BasicTypedefs.h"

#include "ChessboardBiclusteringGibbsSampler.h"
#include "LazyClone.h"

/**
    Realization of Partition concept for objects.
    Read-only part.

    @see ObjectsPartitionEx, ProbesPartition
 */
class ObjectsPartition: public LazyClone<ChessboardBiclusteringFit> {
public:
    typedef object_index_t element_index_type;
    typedef object_set_t element_index_set_type;

    typedef object_clundex_t cluster_index_type;
    typedef boost::unordered_set<cluster_index_type> cluster_index_set_type;

    typedef ObjectsClusterParams cluster_params_type;
    typedef ChessboardBiclusteringFit::ObjectsClusterParamsProxy cluster_params_proxy_type;

protected:
    typedef LazyClone<ChessboardBiclusteringFit> super_type;

public:
    /** Proxy (wrapper) for object_set_t */
    class ObjectsSetProxy: public LazyClone<object_set_t> {
    private:
        typedef LazyClone<object_set_t> super_type;

    public:
        ObjectsSetProxy( const object_set_t& objects )
        : super_type( objects )
        {
        }
        ObjectsSetProxy( const ObjectsPartition& ptn, element_index_type elmIx = OBJECT_NA )
        : super_type( boost::make_shared<object_set_t>(), false )
        {
            if ( elmIx != OBJECT_NA ) wrapped().insert( elmIx );
        }
        bool contains( element_index_type elmIx ) const {
            return ( wrapped().find( elmIx ) != wrapped().end() );
        }
        element_index_type operator[]( size_t index ) const
        {
            element_index_set_type::const_iterator it = wrapped().begin();
            std::advance( it, index );
            return ( it != wrapped().end() ? *it : OBJECT_NA );
        }
        size_t size() const {
            return ( wrapped().size() );
        }
        ObjectsSetProxy& operator+=( const ObjectsSetProxy& that )
        {
            unshare();
            const locked_ptr_type thatpObjs = that.lock();
            wrapped().insert( thatpObjs->begin(), thatpObjs->end() );
            return ( *this );
        }
        ObjectsSetProxy& operator+=( element_index_type elmIx )
        {
            unshare();
            wrapped().insert( elmIx );
            return ( *this );
        }
        ObjectsSetProxy& operator-=( const ObjectsSetProxy& that )
        {
            if ( that.lock()->empty() ) return ( *this );
            unshare();
            const locked_ptr_type thatpObjs = that.lock();
            element_index_set_type::iterator thisIt = wrapped().begin();
            element_index_set_type::iterator thatIt = thatpObjs->begin();
            do {
                if ( *thisIt < *thatIt ) {
                    thisIt = wrapped().lower_bound( *thatIt );
                }
                else if ( *thisIt > *thatIt ) {
                    thatIt = thatpObjs->lower_bound( *thisIt );
                }
                else /*if ( *thisIt == *thatIt )*/ {
                    wrapped().erase( thisIt++ );
                    ++thatIt;
                }
            } while ( thatIt != thatpObjs->end() && thisIt != wrapped().end() );
            return ( *this );
        }
        ObjectsSetProxy& operator-=( element_index_type elmIx )
        {
            unshare();
            wrapped().erase( elmIx );
            return ( *this );
        }
        ObjectsSetProxy& operator&=( const ObjectsSetProxy& that )
        {
            unshare();
            const locked_ptr_type thatpObjs = that.lock();
            if ( thatpObjs->empty() ) {
                wrapped().clear();
                return ( *this );
            }
            element_index_set_type::iterator thisIt = wrapped().begin();
            element_index_set_type::iterator thatIt = thatpObjs->begin();
            do {
                if ( *thisIt < *thatIt ) {
                    wrapped().erase( thisIt++ );
                }
                else if ( *thisIt > *thatIt ) {
                    thatIt = thatpObjs->lower_bound( *thisIt );
                }
                else /*if ( *thisIt == *thatIt )*/ {
                    ++thisIt;
                    ++thatIt;
                }
            } while ( thatIt != thatpObjs->end() && thisIt != wrapped().end() );
            wrapped().erase( thisIt, wrapped().end() );
            return ( *this );
        }
        bool operator==(const ObjectsSetProxy& that ) const {
            const locked_ptr_type thatpObjs = that.lock();
            return ( wrapped() == *thatpObjs );
        }
        bool operator!=(const ObjectsSetProxy& that ) const {
            const locked_ptr_type thatpObjs = that.lock();
            return ( wrapped() != *thatpObjs );
        }
    };

    typedef ObjectsSetProxy elements_set_proxy_type;

    /** 
     *  Proxy (wrapper) for objects cluster,
     *  providing access to its elements and parameters.
     */
    class ObjectsClusterProxy {
    private:
        friend class ObjectsPartition;
        friend class ObjectsPartitionEx;

        ChessboardBiclusteringFit& clus;
        const object_clundex_t cluIx;

        const ObjectsCluster& cluster() const {
            return ( clus.objectsCluster( cluIx ) );
        }

        ObjectsClusterProxy( ChessboardBiclusteringFit& clus, object_clundex_t cluIx )
        : clus( clus ), cluIx( cluIx )
        {
            BOOST_ASSERT( (int)cluIx >= 0 );
            BOOST_ASSERT( cluIx < clus.objectsClusters().size() );
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
            return ( ObjectsSetProxy( cluster().items() ) );
        }
        const cluster_params_proxy_type params(
            const object_set_t& mask = object_set_t()
        ) const {
            return ( clus.objectsClusterParams( cluIx, mask ) );
        }
        cluster_params_proxy_type params(
            const object_set_t& mask = object_set_t()
        ){
            return ( clus.objectsClusterParams( cluIx, mask ) );
        }
    };

    typedef ObjectsClusterProxy cluster_proxy_type;

    const static cluster_index_type ClusterNA = CLUSTER_NA;

    ObjectsPartition( const ChessboardBiclusteringFit& clus )
    : super_type( clus )
    {}

    ObjectsPartition( const shared_object_ptr_type& pClus )
    : super_type( pClus, true )
    {}

    virtual ~ObjectsPartition()
    {}

    size_t clustersCount() const {
        return ( wrapped().objectsClusters().size() );
    }

    size_t elementsCount() const {
        return ( wrapped().objectsCount() );
    }

    const cluster_proxy_type cluster( cluster_index_type cluIx ) const {
        return ( cluster_proxy_type( const_cast<ChessboardBiclusteringFit&>( wrapped() ), cluIx ) );
    }

    cluster_index_type clusterIndex( element_index_type elmIx ) const {
        return ( wrapped().clusterOfObject( elmIx ) );
    }

    operator const ChessboardBiclustering&() const {
        return ( (const ChessboardBiclustering&)wrapped() );
    }
};

/**
    Realization of Partition model for objects.
 */
class ObjectsPartitionEx: public ObjectsPartition {
public:
    ObjectsPartitionEx( const ChessboardBiclusteringFit& clus )
    : ObjectsPartition( clus )
    {}

    ObjectsPartitionEx( const shared_object_ptr_type& pClus )
    : ObjectsPartition( pClus )
    {}

    using ObjectsPartition::cluster;

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
        wrapped().setObjectsClusterSamples( cluIx, samples );
    }
};

/**
    Sampler of objects' parameters for the fixed chessboard biclustering.

    Realization of ParamsSampler model.
 */
class FixedObjectsParamsSampler {
public:
    typedef ObjectsClusterParams params_type;

    /**
        Default parameters.
     */
    params_type defaults() const {
        return ( params_type( clusHelper->clustering() ) );
    }
    bool operator()( params_type& params, const object_set_t& clusterObjects, const object_set_t& sampledObjects, 
                     bool overwrite, bool posterior ) const;

    FixedObjectsParamsSampler( ChessboardBiclusteringGibbsHelper& clusHelper, 
                                bool sampleBlockMask, 
                                bool sampleSignals, 
                                bool sampleMultiples )
    : clusHelper( &clusHelper ), sampleBlockMask( sampleBlockMask )
    , sampleSignals( sampleSignals ), sampleMultiples( sampleMultiples )
    {}

private:
    ChessboardBiclusteringGibbsHelper* const clusHelper;
    bool                            sampleBlockMask;
    bool                            sampleSignals;
    bool                            sampleMultiples;
};

/**
    Sampler of objects' parameters for given chessboard biclustering.

    Realization of ParamsSampler model.
    Internally uses FixedObjectsParamsSampler.
*/
class ObjectsParamsSampler {
public:
    typedef ObjectsClusterParams params_type;

    params_type defaults( const ObjectsPartition& ptn ) const {
        return ( params_type( ptn ) );
    }
    bool operator()( const ObjectsPartition& ptn, params_type& params, const object_set_t& clusterObjects, const object_set_t& sampledObjects, 
                     bool overwrite, bool posterior ) const;

    log_prob_t transitionLP( const ObjectsPartition& before,
                             const ObjectsPartition& after,
                             object_clundex_t cluIx ) const;

    ObjectsParamsSampler( ChessboardBiclusteringGibbsSampler& clusSampler, 
                          bool sampleBlockMask, 
                          bool sampleSignals, 
                          bool sampleMultiples )
    : clusSampler( &clusSampler ), sampleBlockMask( sampleBlockMask )
    , sampleSignals( sampleSignals ), sampleMultiples( sampleMultiples )
    {}

private:
    ChessboardBiclusteringGibbsSampler const* clusSampler;
    bool                            sampleBlockMask;
    bool                            sampleSignals;
    bool                            sampleMultiples;
};

struct FixedObjectsPartitionStats {
    typedef ObjectsClusterParams params_type;

    /**
     *  Likelihood of object parameters in the partition
     *  (multiplier, cluster).
     */
    LLHMetrics objectLLH( object_index_t objIx, const params_type& params ) const;

    /**
     *  Likelihood of objects cluster in the partition,
     *  not concerning its relation to other clusters
     *  (co-occurrence, block probes and abundances, object multipliers).
     */
    LLHMetrics objectsLLH( const object_set_t& objs, const params_type& params ) const;

    /**
     *  Likelihood of modified objects clusters,
     *  (co-occurrence of objects, relation to other clusters,
     *   block probes and abundances, object multipliers).
     */
    LLHMetrics llhDelta( const std::vector<ObjectsPartition::elements_set_proxy_type>& newClusters,
                const std::vector<params_type>& newParams,
                const ObjectsPartition::cluster_index_set_type& oldIndexes ) const;

    LLHMetrics llh() const {
        return ( clusFit.metrics().llhObjs );
    }
    log_prob_t paramsLPP( const params_type& params, const object_set_t& clusterObjects, const object_set_t& sampledObjects ) const;

    FixedObjectsPartitionStats( const ChessboardBiclusteringFit& clusFit )
    : clusFit( clusFit )
    {}

    log_prob_t clusterAssignmentLPP( size_t clusterSize, size_t clustersCount ) const
    {
        return ( log( clusFit.priors().objectClustering.clusterAssignmentPrior( clusterSize, clustersCount, clusFit.objectsCount() - 1 ) ) );
    }

private:
    const ChessboardBiclusteringFit& clusFit;
};

struct ObjectsPartitionStats {
    typedef ObjectsClusterParams params_type;

    LLHMetrics objectLLH( const ObjectsPartition& clus, object_index_t objIx, const params_type& params ) const
    {
        return ( FixedObjectsPartitionStats( clus ).objectLLH( objIx, params ) );
    }
    LLHMetrics llhDelta( const ObjectsPartition& ptn,
                const std::vector<ObjectsPartition::elements_set_proxy_type>& newClusters,
                const std::vector<params_type>& newParams,
                const ObjectsPartition::cluster_index_set_type& oldIndexes
    ) const {
        return ( FixedObjectsPartitionStats( ptn ).llhDelta( newClusters, newParams, oldIndexes ) );
    }
    LLHMetrics llh( const ObjectsPartition& ptn ) const {
        return ( ((const ChessboardBiclusteringFit&)ptn).llh() );
    }
    log_prob_t lpp( const ObjectsPartition& ptn ) const {
        return ( ((const ChessboardBiclusteringFit&)ptn).lpp() );
    }

    log_prob_t clusterAssignmentLPP( size_t clusterSize, size_t clustersCount, size_t totalElements ) const
    {
        return ( log( priors.objectClustering.clusterAssignmentPrior( clusterSize, clustersCount, totalElements - 1 ) ) );
    }
    log_prob_t paramsLPP( const ObjectsPartition& ptn, const params_type& params,
                          const object_set_t& clusterObjects, const object_set_t& sampledObjects ) const {
        return ( FixedObjectsPartitionStats( ptn ).paramsLPP( params, clusterObjects, sampledObjects ) );
    }

    ObjectsPartitionStats(
        const gsl_rng*                  rndNumGen,
        const OPAData&                  data,
        const ChessboardBiclusteringPriors&    priors
    ) : rndNumGen( rndNumGen ), data( data ), priors( priors )
    {}

private:
    const gsl_rng*                      rndNumGen;
    const OPAData&                      data;
    const ChessboardBiclusteringPriors&        priors;
};
