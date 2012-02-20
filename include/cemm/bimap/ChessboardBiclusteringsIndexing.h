#pragma once

#include "BasicTypedefs.h"

#include <set>

#include <cemm/containers/dynamic_bitset_hash.h>
#include <boost/functional/hash.hpp>

#include <boost/serialization/version.hpp>
#include <boost/serialization/nvp.hpp>

#include <cemm/containers/EntityIndexing.h>
#include <cemm/containers/EntityCollectionIndexing.h>

#include "ObjectSet.h"
#include "ProbeSet.h"
#include "ChessboardBiclustering.h"

namespace cemm { namespace bimap {

struct ChessboardBiclusteringScaffold {
    typedef cemm::containers::EntityCollectionIndexing<object_set_t>     object_partition_indexing;
    typedef object_partition_indexing::collection_pointer_type indexed_object_partition_pointer;
    typedef object_partition_indexing::indexed_collection_type indexed_object_partition_type;
    typedef object_partition_indexing::collection_type::const_iterator const_object_cluster_iterator;

    typedef cemm::containers::EntityCollectionIndexing<probe_bitset_t>     probe_partition_indexing;
    typedef probe_partition_indexing::collection_pointer_type indexed_probe_partition_pointer;
    typedef probe_partition_indexing::indexed_collection_type indexed_probe_partition_type;
    typedef probe_partition_indexing::collection_type::const_iterator const_probe_cluster_iterator;

    typedef cemm::containers::default_serial_t objects_cluster_serial_type;
    typedef cemm::containers::default_serial_t probes_cluster_serial_type;

    typedef boost::dynamic_bitset<>         block_mask_t;
    typedef boost::dynamic_bitset<>         cells_mask_t;

    indexed_object_partition_pointer    pObjectsPartition;
    indexed_probe_partition_pointer     pProbesPartition;

    block_mask_t                blockMask;

    bool operator==( const ChessboardBiclusteringScaffold& that ) const {
        return ( pObjectsPartition == that.pObjectsPartition 
              && pProbesPartition == that.pProbesPartition 
              && blockMask == that.blockMask );
    }

    object_clundex_t objectsClusterIndex( objects_cluster_serial_type objectsCluSerial ) const;
    const_object_cluster_iterator objectClusterBegin() const {
        return ( pObjectsPartition->value().begin() );
    }
    const_object_cluster_iterator objectClusterEnd() const {
        return ( pObjectsPartition->value().end() );
    }
    probe_clundex_t probesClusterIndex( probes_cluster_serial_type probesCluSerial ) const;
    const_probe_cluster_iterator probeClusterBegin() const {
        return ( pProbesPartition->value().begin() );
    }
    const_probe_cluster_iterator probeClusterEnd() const {
        return ( pProbesPartition->value().end() );
    }

    bool isBlockEnabled( objects_cluster_serial_type objCluSerial, probes_cluster_serial_type probeCluSerial ) const;

    bool check() const;

    void getCellsMask( cells_mask_t& mask ) const;
};

// forward declaration for ChessboardBiclusteringsIndexing
struct ChessboardBiclusteringScaffoldSerializer;

} }

namespace boost {
inline std::size_t hash_value( const cemm::bimap::ChessboardBiclusteringScaffold& value )
{
    std::size_t seed = 0;
    boost::hash_combine(seed, value.pObjectsPartition->serial() );
    boost::hash_combine(seed, value.pProbesPartition->serial() );
    boost::hash_combine(seed, boost::hash_value( value.blockMask ) );
    return ( seed );
}
}

namespace cemm { namespace containers {
// workaround for namespace-related problem
// boost::hash_value(ns::T) overloads not recognized by gcc?
template<>
struct IndexedEntityValueHasher<cemm::bimap::ChessboardBiclusteringScaffold>
{
    std::size_t operator()( const cemm::bimap::ChessboardBiclusteringScaffold& value ) const {
        //using boost::hash_value;
        return ( boost::hash_value( value ) );
    }
};

} }

namespace cemm { namespace containers {

template<>
struct EntityIndexingTraits<cemm::bimap::ChessboardBiclusteringScaffold>
{
    typedef cemm::bimap::ChessboardBiclusteringScaffoldSerializer entity_serializer_type;
};

} }

namespace cemm { namespace bimap {

#if 0
class PSIClusterIndexed: private PSIClusterData {
private:
    typedef IndexedEntity<Bielements> indexed_type;
    typedef indexed_type::pointer_type  cluster_pointer_type;
    
    cluster_pointer_type     _pElements;
    
protected:
    friend class PSIClusteringIndexed;
    
    PSIClusterIndexed()
    : _pElements( NULL )
    {}

    PSIClusterIndexed( const cluster_pointer_type& pElements, const PSIClusterData& data )
    : PSIClusterData( data )
    , _pElements( pElements )
    {}
    
    const Bielements& bielements() const {
        return ( _pElements->value() );
    }
public:
    indexed_type::serial_type serial() const {
        return ( _pElements->serial() );
    }

    size_t objectsCount() const {
        return ( bielements().objects().size() );
    }
    const object_set_t& objects() const {
        return ( bielements().objects() );
    }
    size_t probesCount() const {
        return ( bielements().probes().count() );
    }
    const probe_bitset_t& probes() const {
        return ( bielements().probes() );
    }
    size_t objectMultiple( object_index_t objIx ) const {
        object_multiple_map_t::const_iterator objIt = _objectMultiples.find( objIx );
        return ( ( objects().find( objIx ) != objects().end() ) && objIt != _objectMultiples.end() ? objIt->second : 0 );
    }
    signal_t probeSignal( probe_index_t probeIx ) const {
        probe_signal_map_t::const_iterator probeIt = _probeSignals.find( probeIx );
        return ( probes().test( probeIx ) && probeIt != _probeSignals.end() ? probeIt->second : gsl_nan() );
    }
};
#endif

class ChessboardBiclusteringIndexed: private ChessboardBiclusteringData {
private: 
    typedef cemm::containers::EntityCollectionIndexing<object_set_t>    object_partition_indexing;
    typedef cemm::containers::EntityCollectionIndexing<probe_bitset_t>  probe_partition_indexing;
    typedef cemm::containers::IndexedEntity<ChessboardBiclusteringScaffold> indexed_type;
    typedef indexed_type::pointer_type  clustering_pointer_type;

public:
    typedef ChessboardBiclusteringScaffold::objects_cluster_serial_type objects_cluster_serial_type;
    typedef ChessboardBiclusteringScaffold::probes_cluster_serial_type probes_cluster_serial_type;
    typedef object_partition_indexing::collection_type objects_cluster_collection_type;
    typedef probe_partition_indexing::collection_type probes_cluster_collection_type;

public:
    typedef std::pair<objects_cluster_serial_type, probes_cluster_serial_type> block_key_type;
    typedef boost::unordered_map<block_key_type, signal_t> block_data_map_type;

private:
    clustering_pointer_type     _pScaffold;
    block_data_map_type         _blocksData;
    multiple_map_t              _objectsData;

protected:
    friend class ChessboardBiclusteringsIndexing;
    friend class BIMAPWalk;
    friend class boost::serialization::access;

    ChessboardBiclusteringIndexed( const clustering_pointer_type&      pScaffold,
                            const block_data_map_type&  blocksData,
                            const multiple_map_t&               objectsData,
                            const ChessboardBiclusteringData&          clusteringData
    ) : ChessboardBiclusteringData( clusteringData ), _pScaffold( pScaffold )
      , _blocksData( blocksData ), _objectsData( objectsData )
    {}

public:
    indexed_type::serial_type serial() const {
        return ( _pScaffold->serial() );
    }

    const ChessboardBiclusteringScaffold& scaffold() const {
        return ( _pScaffold->value() );
    }

    size_t scaffoldRefCount() const {
        return ( _pScaffold->refCount() );
    }

    indexed_type::serial_type objectsPartitionSerial() const {
        return ( _pScaffold->value().pObjectsPartition->serial() );
    }

    size_t objectsPartitionRefCount() const {
        return ( _pScaffold->value().pObjectsPartition->refCount() );
    }

    indexed_type::serial_type probesPartitionSerial() const {
        return ( _pScaffold->value().pProbesPartition->serial() );
    }

    size_t probesPartitionRefCount() const {
        return ( _pScaffold->value().pObjectsPartition->refCount() );
    }

    const objects_cluster_collection_type& objectsClusters() const {
        return ( _pScaffold->value().pObjectsPartition->value() );
    }

    const probes_cluster_collection_type& probesClusters() const {
        return ( _pScaffold->value().pProbesPartition->value() );
    }

    const multiple_map_t& objectsData() const {
        return ( _objectsData );
    }

    const block_data_map_type& blocksData() const {
        return ( _blocksData );
    }

    #if 0
    class ConstClusterIterator {
    private:
        friend class PSIClusteringIndexed;

        typedef BielementsCollection::const_iterator const_internal_iterator;

        const cluster_data_map_type&    _clustersData;
        const_internal_iterator         _it;
        mutable PSIClusterIndexed       _val;
        mutable bool                    _updateVal;

        ConstClusterIterator( const const_internal_iterator& it, const cluster_data_map_type& clustersData )
        : _clustersData( clustersData ), _it( it ), _updateVal( true )
        {}

    protected:
        void updateValue() const {
            if ( _updateVal ) {
                _val = PSIClusterIndexed( *_it, _clustersData.find( (*_it)->serial() )->second );
                _updateVal = false;
            }
        }

    public:
        const PSIClusterIndexed& operator*() const {
            updateValue();
            return ( _val );
        }
        const PSIClusterIndexed* operator->() const {
            updateValue();
            return ( &_val );
        }
        void operator++() {
            ++_it;
            _updateVal = true;
        }
        bool operator==(const ConstClusterIterator& that) const {
            return ( _it == that._it );
        }
        bool operator!=(const ConstClusterIterator& that) const {
            return ( _it != that._it );
        }
    };

    typedef ConstClusterIterator const_iterator;

    const_iterator begin() const {
        return ( const_iterator( clusters().begin(), _clustersData ) );
    }

    const_iterator end() const {
        return ( const_iterator( clusters().end(), _clustersData ) );
    }

    size_t clustersCount() const {
        return ( clusters().size() );
    }
#endif
    const ChessboardBiclusteringScaffold::block_mask_t& blockMask() const {
        return ( _pScaffold->value().blockMask );
    }

    bool isBlockEnabled( objects_cluster_serial_type objCluSerial, probes_cluster_serial_type probeCluSerial ) const {
        return ( _pScaffold->value().isBlockEnabled( objCluSerial, probeCluSerial ) );
    }

    const signal_params_type& baselineSignalParams() const {
        return ( _baselineSignalParams );
    }
    const noise_params_type& noiseParams() const {
        return ( _noiseParams );
    }

    bool check() const;
};

class ChessboardBiclusteringScaffoldSerializer;

class ChessboardBiclusteringsIndexing: public cemm::containers::EntityIndexing<ChessboardBiclusteringScaffold> {
public:
    typedef cemm::containers::EntityCollectionIndexing<object_set_t>    object_partition_indexing;
    typedef cemm::containers::EntityCollectionIndexing<probe_bitset_t>  probe_partition_indexing;

private:
    object_partition_indexing     _objectPartitionIndexing;
    probe_partition_indexing      _probePartitionIndexing;

    typedef object_partition_indexing::collection_pointer_type object_partition_pointer;
    typedef probe_partition_indexing::collection_pointer_type probe_partition_pointer;

    typedef cemm::containers::EntityIndexing<ChessboardBiclusteringScaffold>    base_type;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        LOG_DEBUG3("Serializing objects index");
        ar & boost::serialization::make_nvp( "objectsPartitionIndex", _objectPartitionIndexing );
        LOG_DEBUG3("Serializing probes index");
        ar & boost::serialization::make_nvp( "probesPartitionIndex", _probePartitionIndexing );
        LOG_DEBUG3("Serializing chessboard biclusterings index");
        ar & boost::serialization::make_nvp( "chessboardBiclusteringsIndex", boost::serialization::base_object<base_type>(*this) );
    }

public:
    ChessboardBiclusteringsIndexing(size_type minSize = MinIndexSize, size_type maxSize = MaxIndexSize);

    ChessboardBiclusteringIndexed index( const ChessboardBiclustering& clustering );

    object_partition_indexing& objectPartitionIndexing() {
        return ( _objectPartitionIndexing );
    }

    const object_partition_indexing& objectPartitionIndexing() const {
        return ( _objectPartitionIndexing );
    }

    probe_partition_indexing& probePartitionIndexing() {
        return ( _probePartitionIndexing );
    }

    const probe_partition_indexing& probePartitionIndexing() const {
        return ( _probePartitionIndexing );
    }

    /**
     *  Removes collections lacking outside references.
     */
    void remove_unreferenced() {
        base_type::remove_unreferenced();
        _objectPartitionIndexing.remove_unreferenced();
        _probePartitionIndexing.remove_unreferenced();
    }
};

struct ChessboardBiclusteringScaffoldSerializer
{
    typedef cemm::containers::EntityIndexing<ChessboardBiclusteringScaffold> entity_indexing_type;

    unsigned int entity_version;

    ChessboardBiclusteringsIndexing& indexing;

    ChessboardBiclusteringScaffoldSerializer( const entity_indexing_type& indexing, unsigned int entity_version )
    : entity_version( entity_version )
    , indexing( const_cast<ChessboardBiclusteringsIndexing&>( static_cast<const ChessboardBiclusteringsIndexing&>( indexing ) ) )
    {
    }

    template<class Archive>
    void save_entity(Archive & ar, const ChessboardBiclusteringScaffold& v) const
    {
        cemm::containers::default_serial_t objPtnSerial = v.pObjectsPartition->serial();
        cemm::containers::default_serial_t probePtnSerial = v.pProbesPartition->serial();

        ar << boost::serialization::make_nvp( "objectsPartitionSerial", objPtnSerial );
        ar << boost::serialization::make_nvp( "probesPartitionSerial", probePtnSerial );
        ar << boost::serialization::make_nvp( "blockMask", const_cast<ChessboardBiclusteringScaffold::block_mask_t&>( v.blockMask ) );
    }

    template<class Archive>
    ChessboardBiclusteringScaffold load_entity(Archive & ar) const
    {
        boost::serialization::collection_size_type objPtnSerial;
        boost::serialization::collection_size_type probePtnSerial;

        LOG_DEBUG3("Loading objects partition serial");
        ar >> boost::serialization::make_nvp( "objectsPartitionSerial", objPtnSerial );
        LOG_DEBUG3("Loading probes partition serial");
        ar >> boost::serialization::make_nvp( "probesPartitionSerial", probePtnSerial );

        ChessboardBiclusteringScaffold entity;

        entity.pObjectsPartition = *indexing.objectPartitionIndexing().iterator_to( objPtnSerial );
        entity.pProbesPartition = *indexing.probePartitionIndexing().iterator_to( probePtnSerial );
        LOG_DEBUG3("Loading block mask");
        ar >> boost::serialization::make_nvp( "blockMask", entity.blockMask );
        return ( entity );
    }
};

} }

BOOST_CLASS_VERSION( cemm::bimap::ChessboardBiclusteringsIndexing, 1 );
