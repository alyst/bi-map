#pragma once

#include "BasicTypedefs.h"

#include <stack>
#include <vector>

#include <boost/unordered_map.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>

#include <boost/shared_ptr.hpp>

#include "unordered_serialization.h"
#include "dynamic_bitset_utils.h"

#include "ObjectSet.h"
#include "ProbeSet.h"

#include "ChessboardBiclusteringParams.h"

#include "array2d.h"
#include "VectorCache.h"

typedef size_t cluster_index_t;
typedef cluster_index_t probe_clundex_t;
typedef cluster_index_t object_clundex_t;

/// Non-defined cluster index
#define CLUSTER_NA (cluster_index_t)(-1)

typedef std::vector<size_t>     multiple_map_t;
typedef std::vector<signal_t>   signal_map_t;

typedef boost::dynamic_bitset<> cross_clusters_mask_t;

class ChessboardBiclustering;

/**
    Cluster of objects.
    Element of the objects partition.
 */
class ObjectsCluster {
private:
    friend class ChessboardBiclustering;
    friend class boost::serialization::access;

    object_set_t    _items;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp( "objects", _items );
    }

    /// for serialization
    ObjectsCluster()
    {}

public:

    ObjectsCluster( const object_set_t& items )
    : _items( items )
    {}

    const object_set_t& items() const {
        return ( _items );
    }
    size_t size() const {
        return ( _items.size() );
    }
    object_index_t representative() const {
        return ( *_items.begin() );
    }
};

BOOST_CLASS_IMPLEMENTATION( ObjectsCluster, object_serializable )

/**
    Cluster of probes.
    Element of the probes partition.
 */
class ProbesCluster {
private:
    friend class boost::serialization::access;
    friend class ChessboardBiclustering;

    probe_bitset_t  _items;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp( "probes", _items );
    }

    /// for serialization
    ProbesCluster()
    {
    }

public:
    ProbesCluster( const probe_bitset_t& items )
    : _items( items )
    {}

    ProbesCluster( size_t size )
    : _items( size )
    {}

    const probe_bitset_t& items() const {
        return ( _items );
    }
    size_t size() const {
        return ( _items.count() );
    }
    probe_index_t representative() const {
        return ( _items.find_first() );
    }
};

BOOST_CLASS_IMPLEMENTATION( ProbesCluster, object_serializable )

/**
    Chessboard biclustering of objects x probes matrix.
 */
class ChessboardBiclustering: protected ChessboardBiclusteringData {
private:
    friend class boost::serialization::access;

    typedef std::vector<ObjectsCluster> object_cluster_container_type;
    typedef std::vector<ProbesCluster> probe_cluster_container_type;
    typedef std::vector<object_clundex_t> objects_cluster_index_remap;
    typedef std::vector<probe_clundex_t>  probes_cluster_index_remap;

    object_cluster_container_type       _objectsClusters;
    probe_cluster_container_type        _probesClusters;
    std::vector<object_clundex_t>       _objectToCluster;
    std::vector<probe_clundex_t>        _probeToCluster;
    bool                                _objectsClustersCleanupRequired;
    bool                                _probesClustersCleanupRequired;
    multiple_map_t                      _objectMultiples;
    array2d<signal_t>                   _signals;   /** signals for object-probe pair, in clustering, 
                                                        all objects of the cluster should have the same signal */

    cross_clusters_mask_t               _crossClustersMask; /** 2D array mask of enabled cross clusters */

    objects_cluster_index_remap cleanupObjectsClusters();
    probes_cluster_index_remap cleanupProbesClusters();

    void updateObjectsClustersFromMap();
    void updateProbesClustersFromMap();
    void updateProbeToClusterMap();
    void updateObjectToClusterMap();

    /**
        Wrapper for bi-cluster, a cross product of objects cluster and probes cluster.
    */
    class CrossClusterProxy {
    private:
        friend class ChessboardBiclustering;
        friend class ChessboardBiclusteringIterator;

        ChessboardBiclustering&    _clustering;
        size_t              _ccIndex;       /** index of cross-cluster in cross-clusters mask */
        probe_clundex_t     _probesClusterIndex;

        CrossClusterProxy( ChessboardBiclustering& clustering,
                    size_t index = cross_clusters_mask_t::npos
        ) : _clustering( clustering ), _ccIndex( index )
        {
        }

    public:
        const ChessboardBiclustering& clustering() const {
            return ( _clustering );
        };

        ChessboardBiclustering& clustering() {
            return ( _clustering );
        };

        bool isEnabled() const {
            return ( _clustering._crossClustersMask.test( _ccIndex ) );
        }

        void setEnabled( bool enabled = true ) {
            bool oldEnabled = isEnabled();
            if ( oldEnabled != enabled ) {
                _clustering._crossClustersMask.set( _ccIndex, enabled );
                _clustering.afterCrossClusterFlipped( objectsClusterIndex(), probesClusterIndex() );
            }
        }

        signal_t signal() const {
            return ( _clustering.clusterSignal( objectsClusterIndex(),
                                                probesClusterIndex() ) );
        }

        void setSignal( signal_t signal ) {
            _clustering.setClusterSignal( objectsClusterIndex(), 
                                          probesClusterIndex(), signal );
        }

        object_clundex_t objectsClusterIndex() const {
            return ( _ccIndex != cross_clusters_mask_t::npos 
                    ? _ccIndex / _clustering.probesClusters().size()
                    : CLUSTER_NA );
        }
        const ObjectsCluster& objectsCluster() const {
            return ( _clustering.objectsCluster( objectsClusterIndex() ) );
        }
        probe_clundex_t probesClusterIndex() const {
            return ( _ccIndex != cross_clusters_mask_t::npos 
                    ? _ccIndex % _clustering.probesClusters().size()
                    : CLUSTER_NA );
        }
        const ProbesCluster& probesCluster() const {
            return ( _clustering.probesCluster( probesClusterIndex() ) );
        }
        bool check() const;
    };

    template<class Reference>
    class CrossClusterIterator {
    private:
        friend class ChessboardBiclustering;
        typedef CrossClusterIterator<Reference> self_type;

        CrossClusterProxy    _ccProxy;

        CrossClusterIterator( ChessboardBiclustering& clustering, bool init )
        : _ccProxy( clustering,
                    init ? clustering._crossClustersMask.find_first() 
                         : cross_clusters_mask_t::npos )
        {
        }

        template<class ThatReference>
        void check_clustering( const CrossClusterIterator<ThatReference>& that ) const
        {
            if ( &that._ccProxy._clustering != &_ccProxy._clustering ) {
                throw std::invalid_argument( "Cross-cluster iterator references another clustering" );
            }
        }

        void check_position() const
        {
            if ( _ccProxy._ccIndex == cross_clusters_mask_t::npos ) {
                throw std::runtime_error( "Cannot advance cross-cluster iterator: position not specified" );
            }
        }

    public:
        template<class ThatReference>
        CrossClusterIterator( const CrossClusterIterator<ThatReference>& it )
        : _ccProxy( it._ccProxy )
        {
        }

        template<class ThatReference>
        self_type& operator=( const CrossClusterIterator<ThatReference>& that )
        {
            check_clustering( that );
            _ccProxy._ccIndex = that._ccProxy._ccIndex;
            return ( *this );
        }
        template<class ThatReference>
        bool operator==( const CrossClusterIterator<ThatReference>& that ) const
        {
            check_clustering( that );
            return ( _ccProxy._ccIndex == that._ccProxy._ccIndex );
        }
        template<class ThatReference>
        bool operator!=( const CrossClusterIterator<ThatReference>& that ) const
        {
            check_clustering( that );
            return ( _ccProxy._ccIndex != that._ccProxy._ccIndex );
        }
        self_type operator++(int) const
        {
            check_position();
            self_type res( _ccProxy._clustering, false );
            res._ccProxy._ccIndex = res._ccProxy._clustering._crossClustersMask.find_next( res._ccProxy._ccIndex );
            return ( res );
        }
        self_type& operator++()
        {
            check_position();
            _ccProxy._ccIndex = _ccProxy._clustering._crossClustersMask.find_next( _ccProxy._ccIndex );
            return ( *this );
        }
        const Reference& operator*() const {
            check_position();
            return ( _ccProxy );
        }
        const Reference* operator->() const {
            check_position();
            return ( &_ccProxy );
        }
        Reference& operator*() {
            check_position();
            return ( _ccProxy );
        }
        Reference* operator->() {
            check_position();
            return ( &_ccProxy );
        }
    };

    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        LOG_DEBUG1( "Loading chessboard biclustering..." );
        ar >> boost::serialization::make_nvp( "clusteringData", 
                                              boost::serialization::base_object<ChessboardBiclusteringData>( *this ) );
        ar >> boost::serialization::make_nvp( "objectToCluster", _objectToCluster );
        ar >> boost::serialization::make_nvp( "probeToCluster", _probeToCluster );
        updateProbesClustersFromMap();
        updateObjectsClustersFromMap();
        _objectsClustersCleanupRequired = _probesClustersCleanupRequired = false;
        ar >> boost::serialization::make_nvp( "objectMultiples", _objectMultiples );
        ar >> boost::serialization::make_nvp( "signals", _signals );
        ar >> boost::serialization::make_nvp( "crossClustersMask", _crossClustersMask );
        BOOST_ASSERT( check() );
//        if ( !check() ) THROW_RUNTIME_ERROR( "check() failed for clustering just loaded" );
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        // cleanup clusters (remove empties) before serializing them
        LOG_DEBUG1( "Saving chessboard biclustering..." );
        const_cast<ChessboardBiclustering&>( *this ).cleanupClusters();
        BOOST_ASSERT( !_objectsClustersCleanupRequired
                      && !_probesClustersCleanupRequired );
        BOOST_ASSERT( check() );
//        if ( !check() ) THROW_RUNTIME_ERROR( "check() failed for clustering to be saved" );
        ar << boost::serialization::make_nvp( "clusteringData", 
                                              boost::serialization::base_object<ChessboardBiclusteringData>( *this ) );
        ar << boost::serialization::make_nvp( "objectToCluster", _objectToCluster );
        ar << boost::serialization::make_nvp( "probeToCluster", _probeToCluster );
        ar << boost::serialization::make_nvp( "objectMultiples", _objectMultiples );
        ar << boost::serialization::make_nvp( "signals", _signals );
        ar << boost::serialization::make_nvp( "crossClustersMask", _crossClustersMask );
    }

    void resizeCrossClusterMask( size_t objectClustersCount, size_t probeClustersCount );

public:
    typedef CrossClusterProxy cross_cluster_proxy;
    typedef CrossClusterIterator<cross_cluster_proxy> cross_cluster_iterator;
    typedef CrossClusterIterator<const cross_cluster_proxy> const_cross_cluster_iterator;
    typedef object_cluster_container_type::const_iterator const_object_cluster_iterator;
    typedef probe_cluster_container_type::const_iterator const_probe_cluster_iterator;
    typedef std::pair<object_clundex_t, probe_clundex_t> cluster_cell_key_type;

    class ObjectsClusterParamsProxy;
    friend class ObjectsClusterParamsProxy;

    class ProbesClusterParamsProxy;
    friend class ProbesClusterParamsProxy;

    ChessboardBiclustering( size_t objectsCount = 0, size_t probesCount = 0,
                            bool mapBaitsToObjects = true );

    ChessboardBiclustering( const ChessboardBiclusteringDerivedPriors& derivedPriors,
                     const signal_params_type& baselineSignal, const noise_params_type& noiseParams,
                     const PitmanYorSample& objectsClusters, const PitmanYorSample& probesClusters,
                     bool mapBaitsToObjects = true );

    static cluster_cell_key_type clusterCellKey( object_clundex_t objCluIx, probe_clundex_t probeCluIx )
    {
        return ( std::make_pair( objCluIx, probeCluIx ) );
    }

    const ChessboardBiclusteringData& clusteringData() const {
        return ( *this );
    }

    cross_cluster_iterator findCrossCluster( object_clundex_t objCluIx, probe_clundex_t probeCluIx, bool onlyEnabled = true ) {
        cross_cluster_iterator res( const_cast<ChessboardBiclustering&>( *this ), false );
        size_t ix = objCluIx * probesClusters().size() + probeCluIx;
        if ( !onlyEnabled || _crossClustersMask.test( ix ) ) res._ccProxy._ccIndex = ix;
        return ( res );
    }
    const_cross_cluster_iterator findCrossCluster( object_clundex_t objCluIx, probe_clundex_t probeCluIx, bool onlyEnabled = true ) const {
        const_cross_cluster_iterator res( const_cast<ChessboardBiclustering&>( *this ), false );
        size_t ix = objCluIx * probesClusters().size() + probeCluIx;
        if ( !onlyEnabled || _crossClustersMask.test( ix ) ) res._ccProxy._ccIndex = ix;
        return ( res );
    }

    cross_cluster_iterator begin() {
        return ( cross_cluster_iterator( *this, true ) );
    }
    cross_cluster_iterator end() {
        return ( cross_cluster_iterator( *this, false ) );
    }

    const_cross_cluster_iterator begin() const {
        return ( const_cross_cluster_iterator( const_cast<ChessboardBiclustering&>( *this ), true ) );
    }
    const_cross_cluster_iterator end() const {
        return ( const_cross_cluster_iterator( const_cast<ChessboardBiclustering&>( *this ), false ) );
    }

    cross_clusters_mask_t crossClustersMask( bool objectsMajor = true ) const;
    size_t enabledCrossClustersCount() const {
        return ( _crossClustersMask.count() );
    }

    dynamic_bitset_view objectsClusterSectionMask( object_clundex_t cluIx ) const {
        return ( dynamic_bitset_view( _crossClustersMask, cluIx * _probesClusters.size(), 1, 
                                      probesClusters().size() ) );
    }

    dynamic_bitset_view probesClusterSectionMask( probe_clundex_t cluIx ) const {
        return ( dynamic_bitset_view( _crossClustersMask, cluIx, _probesClusters.size(), 
                                      objectsClusters().size() ) );
    }

    bool isCrossClusterEnabled( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const {
        return ( _crossClustersMask.test( objCluIx * _probesClusters.size() + probeCluIx ) );
    }
    cross_cluster_iterator setCrossCluster( object_clundex_t objCluIx, probe_clundex_t probeCluIx, bool enable = true );

    const const_cross_cluster_iterator crossClusterNotFound() const {
        return ( const_cross_cluster_iterator( const_cast<ChessboardBiclustering&>( *this ), false ) );
    }

    size_t objectsCount() const {
        return ( _objectToCluster.size() );
    }

    size_t probesCount() const {
        return ( _probeToCluster.size() );
    }

    object_clundex_t clusterOfObject( object_index_t objIx ) const {
        return ( _objectToCluster[ objIx ] );
    };

    probe_clundex_t clusterOfProbe( probe_index_t probeIx ) const {
        return ( _probeToCluster[ probeIx ] );
    };

    const ObjectsCluster& objectsCluster( object_clundex_t objCluIx ) const {
        BOOST_ASSERT( (int)objCluIx >= 0 && objCluIx < _objectsClusters.size() );
        return ( _objectsClusters[ objCluIx ] );
    }

    const ProbesCluster& probesCluster( probe_clundex_t probeCluIx ) const {
        BOOST_ASSERT( (int)probeCluIx >= 0 && probeCluIx < _probesClusters.size() );
        return ( _probesClusters[ probeCluIx ] );
    }

    const object_cluster_container_type& objectsClusters() const {
        return ( _objectsClusters );
    };

    const probe_cluster_container_type& probesClusters() const {
        return ( _probesClusters );
    };

    void cleanupClusters();

    signal_t cellSignal( object_index_t objIx, probe_index_t probeIx ) const {
        return ( _signals( objIx, probeIx ) );
    }
    double clusterSignal( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const {
        return ( _signals( objectsCluster( objCluIx ).representative(), 
                           probesCluster( probeCluIx ).representative() ) );
    }
    void setClusterSignal( object_clundex_t objCluIx, probe_clundex_t probeCluIx, signal_t signal ) ;

    size_t objectMultiple( object_index_t objIx ) const {
        return ( _objectMultiples[ objIx ] );
    }
    void setObjectMultiple( object_index_t objIx, size_t multiple );

    const multiple_map_t& objectMultiples() const {
        return ( _objectMultiples );
    }

    void setObjectCluster( object_index_t objIx, const object_clundex_t objCluIx );

    void setProbeCluster( probe_index_t probeIx, const probe_clundex_t probeCluIx );

    object_clundex_t exchangeObjects(
        object_clundex_t clu1Ix, object_clundex_t clu2Ix,
        const object_set_t& newClu2Objs );

    probe_clundex_t exchangeProbes(
        probe_clundex_t clu1Ix, probe_clundex_t clu2Ix,
        const probe_bitset_t& newClu2Probes );

    object_clundex_t addObjectCluster( object_index_t objIx )
    {
        return ( addObjectCluster( &objIx, &objIx + 1 ) );
    }
    template<class Iterator>
    object_clundex_t addObjectCluster( const Iterator& begin, const Iterator& end )
    {
        resizeCrossClusterMask( _objectsClusters.size() + 1, _probesClusters.size() );
        object_clundex_t newCluIx = _objectsClusters.size();
        _objectsClusters.push_back( ObjectsCluster() );
        afterObjectsClusterInserted( newCluIx );
        for ( Iterator it = begin; it != end; ++it ) {
            setObjectCluster( *it, newCluIx );
        }
        return ( newCluIx );
    }
    probe_clundex_t addProbeCluster( probe_index_t probeIx )
    {
        return ( addProbeCluster( &probeIx, &probeIx + 1 ) );
    }
    template<class Iterator>
    probe_clundex_t addProbeCluster( const Iterator& begin, const Iterator& end )
    {
        resizeCrossClusterMask( _objectsClusters.size(), _probesClusters.size() + 1 );
        probe_clundex_t newCluIx = _probesClusters.size();
        _probesClusters.push_back( ProbesCluster( probesCount() ) );
        afterProbesClusterInserted( newCluIx );
        for ( Iterator it = begin; it != end; ++it ) {
            setProbeCluster( *it, newCluIx );
        }
        return ( newCluIx );
    }

    const signal_params_type& baselineSignalParams() const {
        return ( _baselineSignalParams );
    }
    const noise_params_type& noiseParams() const {
        return ( _noiseParams );
    }
    bool isMapBaitsToObjects() const {
        return ( _mapBaitsToObjects );
    }

    void setNoiseParams( const noise_params_type& noiseParams ) {
        if ( _noiseParams != noiseParams ) {
            _noiseParams = noiseParams;
            afterNoiseParamsChanged();
        }
    }
    const ChessboardBiclusteringDerivedPriors& derivedPriors() const {
        return ( _derivedPriors );
    }
    void setSignalPrior( const GaussianDistribution& signalPrior ) {
        if ( _derivedPriors.signalPrior != signalPrior ) {
            _baselineSignalParams = signalPrior.mean;
            _derivedPriors.signalPrior = signalPrior;
            afterSignalPriorChanged();
        }
    }
    bool checkObjectsPartition() const;
    bool checkProbesPartition() const;
    bool checkCrossClusters() const;
    bool check() const;

    const ObjectsClusterParamsProxy objectsClusterParams( object_clundex_t cluIx, const object_set_t& mask = object_set_t() ) const;
    ObjectsClusterParamsProxy objectsClusterParams( object_clundex_t cluIx, const object_set_t& mask = object_set_t() );

    const ProbesClusterParamsProxy probesClusterParams( probe_clundex_t cluIx, const probe_bitset_t& mask = probe_bitset_t() ) const;
    ProbesClusterParamsProxy probesClusterParams( probe_clundex_t cluIx, const probe_bitset_t& mask = probe_bitset_t() );

    void outputSignals( std::ostream& out ) const;

    friend std::ostream& operator<<( std::ostream& out, const ChessboardBiclustering& cc );

protected:
    virtual void beforeObjectsClusterRemoved( object_clundex_t cluIx ) const
    {}

    virtual void beforeProbesClusterRemoved( probe_clundex_t cluIx ) const
    {}

    virtual void afterObjectsClusterInserted( object_clundex_t cluIx ) const
    {}

    virtual void afterProbesClusterInserted( probe_clundex_t cluIx ) const
    {}

    virtual void afterObjectsClusterChanged( object_clundex_t cluIx ) const
    {}

    virtual void afterProbesClusterChanged( probe_clundex_t cluIx ) const
    {}

    virtual void afterObjectMultipleChanged( object_index_t objIx ) const
    {}

    virtual void afterSignalChanged( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const
    {}
    virtual void afterCrossClusterFlipped( object_clundex_t objCluIx, probe_clundex_t probeCluIx ) const
    {}
    virtual void afterSignalPriorChanged() const
    {}
    virtual void afterNoiseParamsChanged() const
    {}
};

/**
    Parameters of the stripe of cells of object cluster.
 */
struct ObjectsClusterParams {
    cross_clusters_mask_t   crossClustersMask;  /** enablement of cross clusters in object stripe */
    signal_map_t&           probeSignal;        /** signals of object for each probe cluster (defined only for enabled cells) */
    multiple_map_t&         objectMultiple;     /** objects' multiplicity */

    ObjectsClusterParams( const ChessboardBiclustering& clustering );
    ObjectsClusterParams( const ChessboardBiclustering::ObjectsClusterParamsProxy& paramsProxy );
    ObjectsClusterParams( const ObjectsClusterParams& params );
    ~ObjectsClusterParams()
    {
        MultipleMapCache.redeem( objectMultiple );
        SignalMapCache.redeem( probeSignal );
    }

    ObjectsClusterParams& operator=( const ObjectsClusterParams& params );
    ObjectsClusterParams& operator=( const ChessboardBiclustering::ObjectsClusterParamsProxy& paramsProxy );

private:
    static VectorCache<size_t> MultipleMapCache;
    static VectorCache<signal_t> SignalMapCache;
};

/**
    Proxy to get/set parameters for objects cluster stripe.
 */
class ChessboardBiclustering::ObjectsClusterParamsProxy {
private:
    friend class ChessboardBiclustering;
    friend class ObjectsClusterParams;

    ChessboardBiclustering&    clus;
    object_clundex_t    cluIx;
    const object_set_t& mask;

    ObjectsClusterParamsProxy( ChessboardBiclustering& clus, object_clundex_t cluIx, const object_set_t& mask )
    : clus( clus )
    , cluIx( cluIx )
    , mask( mask )
    {
        BOOST_ASSERT( (int)cluIx >= 0 && cluIx < clus.objectsClusters().size() );
    };

public:
    const object_set_t& effectiveMask() const {
        return ( mask.empty() ? clus.objectsCluster( cluIx ).items()
                              : mask );
    }

    bool isUsed( object_index_t objIx ) const {
        return ( clus.clusterOfObject( objIx ) == cluIx 
                 && ( mask.empty() || mask.find( objIx ) != mask.end() ) );
    }
    ObjectsClusterParamsProxy& operator=( const ObjectsClusterParams& params );
};

/**
    Parameters of the stripe of cells of object cluster.
 */
struct ProbesClusterParams {
    cross_clusters_mask_t   crossClustersMask;  /** enablement of cross clusters in probe stripe */
    signal_map_t&           objectsSignal;      /** signals of probe for each objects cluster (defined only for enabled cells) */

    ProbesClusterParams( const ChessboardBiclustering& clustering );
    ProbesClusterParams( const ChessboardBiclustering::ProbesClusterParamsProxy& paramsProxy );
    ProbesClusterParams( const ProbesClusterParams& params );
    ~ProbesClusterParams()
    {
        SignalMapCache.redeem( objectsSignal );
    }

    ProbesClusterParams& operator=( const ProbesClusterParams& params );
    ProbesClusterParams& operator=( const ChessboardBiclustering::ProbesClusterParamsProxy& paramsProxy );

private:
    static VectorCache<signal_t> SignalMapCache;
};

/**
    Proxy to get/set parameters for probes cluster stripe.
 */
class ChessboardBiclustering::ProbesClusterParamsProxy {
private:
    friend class ChessboardBiclustering;
    friend class ProbesClusterParams;

    ChessboardBiclustering&        clus;
    probe_clundex_t         cluIx;
    /** mask of probes to set */
    const probe_bitset_t&   mask;

    ProbesClusterParamsProxy( ChessboardBiclustering& clus, probe_clundex_t cluIx, const probe_bitset_t& mask )
    : clus( clus )
    , cluIx( cluIx )
    , mask( mask )
    {
        BOOST_ASSERT( (int)cluIx >= 0 && cluIx < clus.probesClusters().size() );
    };

public:
    /** is given probe enabled by the mask  */
    bool isUsed( probe_index_t probeIx ) const {
        return ( clus.clusterOfProbe( probeIx ) == cluIx
                 && ( mask.none() || mask.test( probeIx ) ) );
    }
    ProbesClusterParamsProxy& operator=( const ProbesClusterParams& params );
};


inline const ChessboardBiclustering::ObjectsClusterParamsProxy ChessboardBiclustering::objectsClusterParams(
    object_clundex_t        cluIx,
    const object_set_t&     mask
) const {
    return ( ObjectsClusterParamsProxy( const_cast<ChessboardBiclustering&>( *this ), cluIx, mask ) );
}

inline ChessboardBiclustering::ObjectsClusterParamsProxy ChessboardBiclustering::objectsClusterParams(
    object_clundex_t        cluIx,
    const object_set_t&     mask
) {
    return ( ObjectsClusterParamsProxy( *this, cluIx, mask ) );
}

inline const ChessboardBiclustering::ProbesClusterParamsProxy ChessboardBiclustering::probesClusterParams(
    probe_clundex_t         cluIx,
    const probe_bitset_t&   mask
) const {
    return ( ProbesClusterParamsProxy( const_cast<ChessboardBiclustering&>( *this ), cluIx, mask ) );
}

inline ChessboardBiclustering::ProbesClusterParamsProxy ChessboardBiclustering::probesClusterParams(
    probe_clundex_t         cluIx,
    const probe_bitset_t&   mask
) {
    return ( ProbesClusterParamsProxy( *this, cluIx, mask ) );
}
