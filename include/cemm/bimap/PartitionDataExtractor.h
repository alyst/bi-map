#pragma once

#include <cemm/containers/dynamic_bitset_iterator.h>

#include "BIMAPWalk.h"

namespace cemm { namespace bimap {

/**
 *  Model of PartitionDataExtractor for extracting data
 *  from partitions, used by IndexedPartitionsCollection.
 */
template<class Part>
struct PartitionDataExtractor
{
    typedef size_t partition_serial;
    typedef size_t partition_indexing;
    typedef size_t part_indexing;
    typedef size_t part_serial;
    typedef void const_partition_iterator;
    typedef void const_part_in_index_iterator;
    typedef void const_part_iterator;

    typedef size_t elements_container;
    typedef size_t element_type;
    typedef void const_element_iterator;

    static partition_serial PartitionSerial( const ChessboardBiclusteringIndexed& cc );

    static const partition_indexing& PartitionsIndexing( const BIMAPWalk& walk );
    static const_partition_iterator PartitionsBegin( const BIMAPWalk& walk )
    {
        return ( PartitionsIndexing( walk ).serialMap().begin() );
    }
    static const_partition_iterator PartitionsEnd( const BIMAPWalk& walk )
    {
        return ( PartitionsIndexing( walk ).serialMap().end() );
    }
    static const part_indexing& PartsIndexing( const BIMAPWalk& walk );
    static part_indexing& PartsIndexing( BIMAPWalk& walk );
    static const_part_iterator PartsBegin( const ChessboardBiclusteringIndexed& cc );
    static const_part_iterator PartsEnd( const ChessboardBiclusteringIndexed& cc );
    static part_serial PartSerial( const ChessboardBiclusteringIndexed& cc );

    static size_t Size( const elements_container& set );
    static bool IsEmpty( const elements_container& set );
    static bool IsSubset( const elements_container& superset, const elements_container& subset );
    static bool Contains( const elements_container& set, element_type elem );
    static bool IsIntersecting( const elements_container& a, const elements_container& b );
    static elements_container Intersect( const elements_container& a, const elements_container& b );

    /**
     *  Splits divisible (A) by divisor (B) into A&B and A/B.
     */
    static std::pair<elements_container, elements_container> Split( const elements_container& divisible,
                                                                    const elements_container& divisor );
    static const_element_iterator ElementsBegin( const elements_container& part );
    static const_element_iterator ElementsEnd( const elements_container& part );
    static void ElementInsert( elements_container& part, size_t element );
    static size_t TotalElementsCount( const BIMAPWalk& walk );
    static elements_container EmptyElementsSet( const BIMAPWalk& walk );
};

/**
 *  Realization of PartitionDataExtractor model for ObjectsCluster.
 */
template<>
struct PartitionDataExtractor<ObjectsCluster>
{
    typedef ChessboardBiclusteringsIndexing::object_partition_indexing partition_indexing;
    typedef partition_indexing::const_serial_iterator const_partition_iterator;
    typedef ChessboardBiclusteringsIndexing::object_partition_indexing::element_indexing_type part_indexing;
    typedef part_indexing::const_serial_iterator const_part_in_index_iterator;
    typedef partition_indexing::collection_type::const_iterator const_part_iterator;
    typedef object_set_t elements_container;
    typedef object_set_t::const_iterator const_element_iterator;
    typedef object_index_t element_type;

    static size_t PartitionSerial( const ChessboardBiclusteringIndexed& cc )
    {
        return ( cc.objectsPartitionSerial() );
    }

    static const_partition_iterator PartitionsBegin( const BIMAPWalk& walk )
    {
        return ( PartitionsIndexing( walk ).serialMap().begin() );
    }
    static const_partition_iterator PartitionsEnd( const BIMAPWalk& walk )
    {
        return ( PartitionsIndexing( walk ).serialMap().end() );
    }

    static const_part_iterator PartsBegin( const ChessboardBiclusteringIndexed& cc )
    {
        return ( cc.objectsClusters().begin() );
    }
    static const_part_iterator PartsEnd( const ChessboardBiclusteringIndexed& cc )
    {
        return ( cc.objectsClusters().end() );
    }

    static const partition_indexing& PartitionsIndexing( const BIMAPWalk& walk )
    {
        return ( walk.indexing().objectPartitionIndexing() );
    }
    static const part_indexing& PartsIndexing( const BIMAPWalk& walk )
    {
        return ( walk.indexing().objectPartitionIndexing().elementIndexing() );
    }
    static part_indexing& PartsIndexing( BIMAPWalk& walk )
    {
        return ( walk.indexing().objectPartitionIndexing().elementIndexing() );
    }
    static void ElementInsert( elements_container& part, size_t element )
    {
        part.insert( element );
    }
    static const_element_iterator ElementsBegin( const elements_container& part )
    {
        return ( part.begin() );
    }
    static const_element_iterator ElementsEnd( const elements_container& part )
    {
        return ( part.end() );
    }
    static size_t TotalElementsCount( const BIMAPWalk& walk ) {
        return ( walk.stepsBegin()->clustering.objectsData().size() );
    }
    static elements_container EmptyElementsSet( const BIMAPWalk& walk )
    {
        return ( elements_container() );
    }
    static size_t Size( const elements_container& set )
    {
        return ( set.size() );
    }
    static bool IsEmpty( const elements_container& set )
    {
        return ( set.empty() );
    }
    static bool IsSubset( const elements_container& sub, const elements_container& super )
    {
        return ( std::includes( super.begin(), super.end(), sub.begin(), sub.end() ) );
    }
    static bool Contains( const elements_container& set, const element_type& elem )
    {
        return ( set.count( elem ) > 0 );
    }
    static bool IsIntersecting( const elements_container& a, const elements_container& b )
    {
        for ( elements_container::const_iterator it = a.begin(); it != a.end(); ++it ) {
            if ( b.count( *it ) > 0 ) return ( true );
        }
        return ( false );
    }
    static elements_container Intersect( const elements_container& a, const elements_container& b )
    {
        elements_container res;
        std::set_intersection( a.begin(), a.end(), b.begin(), b.end(),
                               std::inserter( res, res.begin() ) );
        return ( res );
    }

    static std::pair<elements_container, elements_container> Split( const elements_container& divisible,
                                                                    const elements_container& divisor )
    {
        elements_container isect;
        std::set_intersection(
            divisible.begin(), divisible.end(),
            divisor.begin(), divisor.end(),
            std::inserter( isect, isect.begin() ) );
        elements_container diff;
        std::set_difference( divisible.begin(), divisible.end(),
                             divisor.begin(), divisor.end(),
                             std::inserter( diff, diff.begin() ) );
        return ( std::make_pair( isect, diff ) );
    }
};

/**
 *  Realization of PartitionDataExtractor model for ProbesCluster.
 */
template<>
struct PartitionDataExtractor<ProbesCluster>
{
    typedef ChessboardBiclusteringsIndexing::probe_partition_indexing partition_indexing;
    typedef partition_indexing::const_serial_iterator const_partition_iterator;
    typedef ChessboardBiclusteringsIndexing::probe_partition_indexing::element_indexing_type part_indexing;
    typedef part_indexing::const_serial_iterator const_part_in_index_iterator;
    typedef partition_indexing::collection_type::const_iterator const_part_iterator;
    typedef probe_bitset_t elements_container;
    typedef probe_index_t element_type;
    typedef dynamic_bitset_const_iterator const_element_iterator;

    static size_t PartitionSerial( const ChessboardBiclusteringIndexed& cc )
    {
        return ( cc.probesPartitionSerial() );
    }

    static const_partition_iterator PartitionsBegin( const BIMAPWalk& walk )
    {
        return ( PartitionsIndexing( walk ).serialMap().begin() );
    }
    static const_partition_iterator PartitionsEnd( const BIMAPWalk& walk )
    {
        return ( PartitionsIndexing( walk ).serialMap().end() );
    }

    static const_part_iterator PartsBegin( const ChessboardBiclusteringIndexed& cc )
    {
        return ( cc.probesClusters().begin() );
    }
    static const_part_iterator PartsEnd( const ChessboardBiclusteringIndexed& cc )
    {
        return ( cc.probesClusters().end() );
    }

    static const partition_indexing& PartitionsIndexing( const BIMAPWalk& walk )
    {
        return ( walk.indexing().probePartitionIndexing() );
    }
    static const part_indexing& PartsIndexing( const BIMAPWalk& walk )
    {
        return ( walk.indexing().probePartitionIndexing().elementIndexing() );
    }
    static part_indexing& PartsIndexing( BIMAPWalk& walk )
    {
        return ( walk.indexing().probePartitionIndexing().elementIndexing() );
    }
    static size_t TotalElementsCount( const BIMAPWalk& walk ) {
        return ( (*walk.indexing().probePartitionIndexing().elementIndexing().valueMap().begin() )->value().size() );
    }
    static void ElementInsert( elements_container& part, size_t element )
    {
        part.set( element );
    }
    static const_element_iterator ElementsBegin( const elements_container& part )
    {
        return ( const_element_iterator::Begin( part ) );
    }
    static const_element_iterator ElementsEnd( const elements_container& part )
    {
        return ( const_element_iterator::End( part ) );
    }
    static elements_container EmptyElementsSet( const BIMAPWalk& walk )
    {
        return ( elements_container( TotalElementsCount( walk ) ) );
    }
    static size_t Size( const elements_container& set )
    {
        return ( set.count() );
    }
    static bool IsEmpty( const elements_container& set )
    {
        return ( set.none() );
    }
    static bool IsSubset( const elements_container& subset, const elements_container& superset )
    {
        return ( subset.is_subset_of( superset ) );
    }
    static bool Contains( const elements_container& set, const element_type& elem )
    {
        return ( set.test( elem ) );
    }
    static bool IsIntersecting( const elements_container& a, const elements_container& b )
    {
        return ( a.intersects( b ) );
    }
    static elements_container Intersect( const elements_container& a, const elements_container& b )
    {
        return ( a & b );
    }
    static std::pair<elements_container, elements_container> Split( const elements_container& divisible,
                                                                    const elements_container& divisor )
    {
        return ( std::make_pair( divisible & divisor, divisible - divisor ) );
    }
};

} }