#pragma once

#include "BasicTypedefs.h"

#include <algorithm>

#include <boost/unordered/unordered_map.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/functional/hash/extensions.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include "PartitionStatistics.h"
#include "PartitionDataExtractor.h"

#include "BIMAPWalk.h"

namespace cemm { namespace bimap {

typedef boost::numeric::ublas::symmetric_matrix<size_t> pair_counts_matrix_t;

/**
 *  Collection of indexed partitions.
 * 
 */
template<class Part>
class IndexedPartitionsCollection
{
public:
    typedef size_t part_serial;
    typedef size_t partition_serial;
    typedef size_t component_index;
    typedef size_t subpartition_serial;
    typedef PartitionDataExtractor<Part> extractor_type;

    typedef boost::unordered_map<part_serial, component_index> part_component_map;
    typedef std::set<part_serial> component_parts;
    typedef std::vector<subpartition_serial> partition_composition;
    typedef boost::unordered_map<partition_serial, partition_composition> ptn_composition_map;
    typedef boost::unordered_map<part_serial, prob_t> part_prob_map;

private:
    typedef typename extractor_type::partition_indexing partition_indexing;
    typedef typename extractor_type::const_partition_iterator const_partition_iterator;

    typedef typename extractor_type::part_indexing part_indexing;
    typedef typename extractor_type::const_part_iterator const_part_iterator;
    typedef typename extractor_type::const_part_in_index_iterator const_part_in_index_iterator;

    typedef typename extractor_type::elements_container elements_container;
    typedef typename extractor_type::const_element_iterator const_elem_iterator;

    typedef EntityIndexing<component_parts> subpartition_indexing;

    BIMAPWalk&             _walk;
    counts_map_type         _partitionCounts;
    part_stats_map_type     _partStats;
    pair_counts_matrix_t    _pairCounts;
    std::vector<elements_container> _components;
    part_component_map      _partComponent;

    ptn_composition_map     _ptnComposition;
    boost::ptr_vector<subpartition_indexing> _subptnIndexes;
    boost::ptr_vector<counts_map_type>   _subptnCounts;

    void countPartitions();
    void countParts();
    void countElementPairs();
    std::vector<elements_container> coOccurenceSingleLinkageClusters( size_t threshold ) const;
    void classifyParts();
    void indexPartitionsComposition();
    void countSubpartitions();
    void countElementPairsAvgCooccurrencePerPart( bool onlyUnused );
    void countPartsInclusion( bool onlyUnused );
    size_t generateNeighboursIntersections( const gsl_rng* rng, const elements_container& clot,
                                          size_t maxIsections, size_t isectionTrials );
    size_t generateIntersectionsWithinClots( const gsl_rng* rng,  size_t clotThreshold,
                                           size_t maxIsections, size_t isectionTrials );

public:
    IndexedPartitionsCollection( BIMAPWalk& walk )
    : _walk( walk )
    {}

    const BIMAPWalk& walk() const {
        return ( _walk );
    }

    const std::vector<elements_container>& components() const {
        return ( _components );
    }

    const PartStats& partStats( size_t partSerial ) const  {
        part_stats_map_type::const_iterator pit = _partStats.find( partSerial );
        return ( pit != _partStats.end() ? pit->second : PartStats() );
    }

    const part_stats_map_type& partStats() const {
        return ( _partStats );
    }

    const counts_map_type& partitionCountsMap() const {
        return ( _partitionCounts );
    }

    size_t partitionCounts( size_t partitionSerial ) const  {
        counts_map_type::const_iterator pit = _partitionCounts.find( partitionSerial );
        return ( pit != _partitionCounts.end() ? pit->second : 0 );
    }

    const elements_container& partElements( size_t partSerial ) const {
        return ( extractor_type::PartsIndexing( _walk ).entity( partSerial )->value() );
    }

    void init( const gsl_rng* rng, size_t independenceThreshold, size_t clotThreshold );

    const ptn_composition_map& partitionCompositionMap() const {
        return ( _ptnComposition );
    }

    const boost::ptr_vector<counts_map_type>& subpartitionCountsMaps() const {
        return ( _subptnCounts );
    }
};

/**
 *  Count how many times partition part (cluster)
 *  had occurred in the whole collection of partitions.
 * 
 *  Requires countPartitions() to be called before.
 */
template<class Part>
void IndexedPartitionsCollection<Part>::countParts()
{
    std::set<size_t> ptnObserved;
    for ( BIMAPWalk::const_step_iterator sit = _walk.stepsBegin(); sit != _walk.stepsEnd(); ++sit ) {
        const ChessboardBiclusteringIndexed& cc = sit->clustering;
        const size_t ptnSerial = extractor_type::PartitionSerial( cc );
        // process clustering, if it was not done before
        if ( ptnObserved.find( ptnSerial ) == ptnObserved.end() ) {
            size_t ptnCounts = _partitionCounts.find( ptnSerial )->second;
            for ( const_part_iterator pit = extractor_type::PartsBegin( cc ); pit != extractor_type::PartsEnd( cc ); ++pit ){
                part_stats_map_type::iterator cit = _partStats.find( (*pit)->serial() );
                if ( cit != _partStats.end() ) {
                    cit->second.nsteps += ptnCounts;
                } else {
                    _partStats.insert( cit, std::make_pair( (*pit)->serial(), 
                                                            PartStats( extractor_type::Size( (*pit)->value() ),
                                                                       ptnCounts ) ) );
                }
            }
            ptnObserved.insert( ptnSerial ); // mark as processed
        }
    }
}

/**
 *  Count how many times specific partition
 *  occurred in the whole collection.
 */
template<class Part>
void IndexedPartitionsCollection<Part>::countPartitions()
{
    for ( BIMAPWalk::const_step_iterator sit = _walk.stepsBegin(); sit != _walk.stepsEnd(); ++sit ) {
        const ChessboardBiclusteringIndexed& cc = sit->clustering;
        const size_t ptnSerial = extractor_type::PartitionSerial( cc );
        counts_map_type::iterator pit = _partitionCounts.find( ptnSerial );
        if ( pit != _partitionCounts.end() ) {
            pit->second++;
        } else {
            _partitionCounts[ ptnSerial ] = 1;
        }
    }
}

/**
 *  Calculates the counts of element pairs co-occurrences in the same cluster,
 *  for a collection of partitions.
 *  Requires countParts() to be called before.
 */
template<class Part>
void IndexedPartitionsCollection<Part>::countElementPairs()
{
    const part_indexing& partsIndexing = extractor_type::PartsIndexing( _walk );
    _pairCounts.resize( extractor_type::TotalElementsCount( _walk ), false );
    _pairCounts.clear();

    for ( const_part_in_index_iterator pit = partsIndexing.serialMap().begin(); 
          pit != partsIndexing.serialMap().end(); ++pit
    ){
        const part_serial partSerial = (*pit)->serial();
        const size_t partCounts = _partStats.find( partSerial )->second.nsteps;
        const const_elem_iterator elEnd = extractor_type::ElementsEnd( (*pit)->value() );
        for ( const_elem_iterator e1it = extractor_type::ElementsBegin( (*pit)->value() );
              e1it != elEnd; ++e1it
        ){
            for ( const_elem_iterator e2it = e1it; e2it != elEnd; ++e2it ) {
                _pairCounts( *e1it, *e2it ) += partCounts;
            }
        }
    }
}

template<class Part>
struct is_element_contaned {
    typedef PartitionDataExtractor<Part> extractor_type;
    typedef typename extractor_type::element_type element_type;
    typedef typename extractor_type::const_part_in_index_iterator::value_type set_type;

    element_type    _elem;

    is_element_contaned( const element_type& elem )
    : _elem( elem )
    {}

    bool operator()( const set_type& set ) const
    {
        return ( extractor_type::Contains( set->value(), _elem ) );
    }
};

template<class Part>
struct is_intersecting {
    typedef PartitionDataExtractor<Part> extractor_type;
    typedef typename extractor_type::elements_container elements_container;
    typedef typename extractor_type::element_type element_type;
    typedef typename extractor_type::const_part_in_index_iterator::value_type set_type;

    elements_container    _set;

    is_intersecting( const elements_container& set )
    : _set( set )
    {}

    bool operator()( const set_type& indexedSet ) const
    {
        return ( extractor_type::IsIntersecting( indexedSet->value(), _set ) );
    }
};

/**
 *  Count how many times partition part (cluster)
 *  had occurred in the whole collection of partitions.
 * 
 *  Requires countParts() to be called before.
 */
template<class Part>
void IndexedPartitionsCollection<Part>::countPartsInclusion(
    bool onlyUnused /** i count inclusion only for sets not being a part of any partition */
){
    typedef boost::filter_iterator<is_element_contaned<Part>,
                       typename part_indexing::const_serial_iterator> filtered_part_iterator;

    const part_indexing& partsColl = extractor_type::PartsIndexing( _walk );

    for ( part_stats_map_type::iterator psit = _partStats.begin(); psit != _partStats.end(); ++psit ) {
        if ( onlyUnused && psit->second.nsteps > 0 ) continue;
        const_part_in_index_iterator pit = partsColl.serialMap().find( psit->first );
        if ( pit == partsColl.serialMap().end() ) throw std::runtime_error( "Part not found" );
        size_t& nStepsIncluded = psit->second.nstepsIncluded;
        nStepsIncluded = 0;

        is_element_contaned<Part> containsElem( *extractor_type::ElementsBegin( (*pit)->value() ) ); // use first element to reduce search space
        filtered_part_iterator pit2end( containsElem, partsColl.serialMap().end(), partsColl.serialMap().end() );
        for ( filtered_part_iterator pit2( containsElem, partsColl.serialMap().begin(), partsColl.serialMap().end() );
              pit2 != pit2end; ++pit2
        ){
            if ( extractor_type::IsSubset( (*pit)->value(), (*pit2)->value() ) ) {
                part_stats_map_type::const_iterator psit2 = _partStats.find( (*pit2)->serial() );
                nStepsIncluded += psit2->second.nsteps;
            }
        }
    }
}

template<class Part>
void IndexedPartitionsCollection<Part>::countElementPairsAvgCooccurrencePerPart(
    bool onlyUnused /** i count inclusion only for sets not being a part of any partition */
){
    const part_indexing& partsColl = extractor_type::PartsIndexing( _walk );
    part_prob_map partAvgFreqMap( partsColl.size() );

    for ( part_stats_map_type::iterator psit = _partStats.begin(); psit != _partStats.end(); ++psit ) {
        if ( onlyUnused && psit->second.nsteps > 0 ) continue;
        const_part_in_index_iterator pit = partsColl.serialMap().find( psit->first );
        if ( pit == partsColl.serialMap().end() ) throw std::runtime_error( "Part not found" );
        size_t cooccurSum = 0;
        const const_elem_iterator elEnd = extractor_type::ElementsEnd( (*pit)->value() );
        for ( const_elem_iterator e1it = extractor_type::ElementsBegin( (*pit)->value() );
              e1it != elEnd; ++e1it
        ){
            for ( const_elem_iterator e2it = e1it; e2it != elEnd; ++e2it ) {
                cooccurSum += _pairCounts( *e1it, *e2it );
            }
        }
        psit->second.avgPairCoOccurrence = 2 * cooccurSum / (double)( psit->second.size * ( psit->second.size + 1 ) );
    }
}

/**
 *  Clusters elements together based on their co-occurrence
 *  matrix and threshold
 */
template<class Part>
std::vector<typename IndexedPartitionsCollection<Part>::elements_container>
IndexedPartitionsCollection<Part>::coOccurenceSingleLinkageClusters(
    size_t threshold    /** co-occurrence threshold to consider pair of elements independent */
) const {
    const size_t ClusterNotSet = (size_t)(-1);

    size_t  nextCluIx = 0;
    std::vector<size_t>    elmToCluster( _pairCounts.size1(), ClusterNotSet );
    for ( size_t elm1Ix = 0; elm1Ix < _pairCounts.size1(); elm1Ix++ ) {
        // set the index of the cluster for elm1 -- either new or one of linked elements
        size_t& elm1CluIx = elmToCluster[ elm1Ix ];
        // search for cluster assignment for the elements, elm1 is linked to
        // searching only succeeding elements, as preceding elements should have been assigned elm1Ix cluster
        for ( size_t elm2Ix = elm1Ix + 1; elm2Ix < _pairCounts.size2(); elm2Ix++ ) {
            const size_t elm2CluIx = elmToCluster[ elm2Ix ];
            if ( ( elm2CluIx != ClusterNotSet ) 
                    && ( _pairCounts( elm1Ix, elm2Ix ) >= threshold )
            ){
                if ( elm1CluIx == ClusterNotSet ) {
                    elm1CluIx = elm2CluIx;
                }
                else if ( elm1CluIx != elm2CluIx ) {
                    // merge clusters, replace index by the minimal, update other indicies
                    size_t mergedCluIx = std::min( elm1CluIx, elm2CluIx );
                    size_t removedCluIx = std::max( elm1CluIx, elm2CluIx );
                    nextCluIx--;
                    for ( size_t elm3Ix = 0; elm3Ix < elmToCluster.size(); elm3Ix++ ) {
                        size_t& elm3CluIx = elmToCluster[ elm3Ix ];
                        if ( elm3CluIx == removedCluIx ) {
                            elm3CluIx = mergedCluIx;
                        }
                        else if ( ( elm3CluIx != ClusterNotSet ) && ( elm3CluIx > removedCluIx ) ) {
                            elm3CluIx--;
                        }
                    }
                }
            }
        }
        // no cluster assignments yet, create new cluster
        if ( elm1CluIx == ClusterNotSet ) {
            elm1CluIx = nextCluIx++;
        }
        // set cluster index of elements linked to elm1
        for ( size_t elm2Ix = elm1Ix + 1; elm2Ix < _pairCounts.size2(); elm2Ix++ ) {
            if ( _pairCounts( elm1Ix, elm2Ix ) >= threshold ) {
                size_t& elm2CluIx = elmToCluster[ elm2Ix ];
                if ( elm2CluIx != ClusterNotSet ) {
                    if ( elm2CluIx != elm1CluIx ) {
                        throw std::runtime_error( "Inconsistent single-linkage clustering" );
                    }
                }
                else {
                    // create new cluster
                    elm2CluIx = elm1CluIx;
                }
            }
        }
    }
    // create clusters from element-to-cluster map
    std::vector<elements_container> res( nextCluIx, extractor_type::EmptyElementsSet( _walk ) );
    for ( size_t elm1Ix = 0; elm1Ix < elmToCluster.size(); elm1Ix++ ) {
        extractor_type::ElementInsert( res[ elmToCluster[ elm1Ix ] ], elm1Ix );
    }
    return ( res );
}

/**
 *  Classify indexed parts by independent components to which they belong.
 */
template<class Part>
void IndexedPartitionsCollection<Part>::classifyParts()
{
    const typename extractor_type::part_indexing& pidx = extractor_type::PartsIndexing( _walk );
    for ( const_part_in_index_iterator pit = pidx.serialMap().begin();
          pit != pidx.serialMap().end(); ++pit
    ){
        for ( size_t compIx = 0; compIx < _components.size(); ++compIx ) {
            if ( extractor_type::IsSubset( (*pit)->value(), _components[ compIx ] ) ) {
                _partComponent[ (*pit)->serial() ] = compIx;
                break;
            }
        }
    }
}

/**
 *  Re-index parts of partitions within independent components.
 *  AKA calculate the marginal frequencies.
 */
template<class Part>
void IndexedPartitionsCollection<Part>::indexPartitionsComposition()
{
    partition_indexing ptnIndex = extractor_type::PartitionsIndexing( _walk );
    for ( const_partition_iterator ptnIt = ptnIndex.serialMap().begin(); ptnIt != ptnIndex.serialMap().end(); ++ptnIt ) {
        std::vector<component_parts>    subptn( _components.size() );
        for ( const_part_iterator pit = (*ptnIt)->value().begin(); pit != (*ptnIt)->value().end(); ++pit ) {
            size_t partSerial = (*pit)->serial();
            part_component_map::const_iterator compIt = _partComponent.find( partSerial );
            if ( compIt != _partComponent.end() ) {
                // short path, part if fully contained in the component
                subptn[ compIt->second ].insert( partSerial );
            } else {
                // long path, part is partially in the component, split it by the components
                elements_container partRemainder = (*pit)->value();
                for ( size_t compIx = 0; compIx < _components.size(); ++compIx ) {
                    std::pair<elements_container, elements_container> split = 
                        extractor_type::Split( partRemainder, _components[ compIx ] );
                    if ( !extractor_type::IsEmpty( split.first ) ) {
                        // put part's part to the component partition
                        subptn[ compIx ].insert( extractor_type::PartsIndexing( _walk ).index( split.first )->serial() );
                    }
                    if ( extractor_type::IsEmpty( split.second ) ) {
                        // part fully processed
                        break;
                    } else {
                        // leave the rest for further indexing
                        partRemainder = split.second;
                    }
                }
            }
        }
        // index the component's composition
        partition_composition ptnComponentSerials( _components.size() );
        for ( size_t compIx = 0; compIx < _components.size(); ++compIx ) {
            BOOST_ASSERT( !subptn[ compIx ].empty() );
            ptnComponentSerials[ compIx ] = _subptnIndexes[ compIx ].index( subptn[ compIx ] )->serial();
        }
        _ptnComposition[ (*ptnIt)->serial() ] = ptnComponentSerials;
    }
}

/**
 *  Count subpartition occurrences for each component.
 */
template<class Part>
void IndexedPartitionsCollection<Part>::countSubpartitions()
{
    for ( counts_map_type::const_iterator ptnIt = _partitionCounts.begin(); ptnIt != _partitionCounts.end(); ++ptnIt ) {
        partition_serial ptnSerial = ptnIt->first;
        size_t ptnCounts = ptnIt->second;
        ptn_composition_map::const_iterator compIt = _ptnComposition.find( ptnSerial );
        if ( compIt != _ptnComposition.end() ) {
            const partition_composition& comp = compIt->second;
            for ( size_t i = 0; i < comp.size(); i++ ) {
                counts_map_type::iterator countsIt = _subptnCounts[ i ].find( comp[ i ] );
                if ( countsIt == _subptnCounts[ i ].end() ) {
                    _subptnCounts[ i ][ comp[ i ] ] = ptnCounts;
                }
                else {
                    countsIt->second += ptnCounts;
                }
            }
        }
    }
}

template<class Part>
size_t IndexedPartitionsCollection<Part>::generateNeighboursIntersections(
    const gsl_rng* rng, 
    const elements_container& clot,
    size_t maxIsections,
    size_t isectionTrials
){
    size_t n_new_sets = 0;

    typedef boost::filter_iterator<is_intersecting<Part>, typename part_indexing::const_serial_iterator> filtered_part_iterator;

    std::vector<part_serial> isectingSets;

    const part_indexing& partsColl = extractor_type::PartsIndexing( _walk );
    is_intersecting<Part> intersectsClot( clot );
    filtered_part_iterator fpitEnd( intersectsClot, partsColl.serialMap().end(), partsColl.serialMap().end() );

    for ( filtered_part_iterator fpit( intersectsClot, partsColl.serialMap().begin(), partsColl.serialMap().end() );
          fpit != fpitEnd; ++fpit
    ){
        isectingSets.push_back( (*fpit)->serial() );
    }
    if ( isectingSets.size() < 2 ) return ( 0 ); // trivial cases -- already in indexing

    maxIsections = std::min( maxIsections, isectingSets.size() );
    std::vector<part_serial>    trialSet( maxIsections, 0 );
    for ( size_t n = 0; n < isectionTrials; ++n ) {
        // random shuffled subset of isectingSets
        gsl_ran_choose( rng, trialSet.data(), trialSet.size(),
                        isectingSets.data(), isectingSets.size(),
                        sizeof( part_serial ) );
        gsl_ran_shuffle( rng, trialSet.data(), trialSet.size(),
                         sizeof( part_serial ) );
        // put to indexing all non-trivial intersections
        elements_container isection;
        for ( size_t i = 0; i < trialSet.size(); i++ ) {
            const elements_container& set = ( *partsColl.serialMap().find( isectingSets[i] ) )->value();
            if ( i == 0 ) {
                isection = set; // initialize intersection
            }
            else {
                isection = extractor_type::Intersect( isection, set );
                if ( isection.size() > 1 ) {
                    // put intersection into indexing
                    typename const_part_in_index_iterator::value_type part = extractor_type::PartsIndexing( _walk ).index( isection );
                    std::pair<part_stats_map_type::const_iterator, bool> res = 
                        _partStats.insert( std::make_pair( part->serial(), 
                                                           PartStats( extractor_type::Size( part->value() ),
                                                           0 ) ) );
                    if ( res.second ) n_new_sets++;
                } else {
                    break; // trivial cases
                }
            }
        }
    }
    return ( n_new_sets );
}

template<class Part>
size_t IndexedPartitionsCollection<Part>::generateIntersectionsWithinClots(
    const gsl_rng* rng, 
    size_t clotThreshold,
    size_t maxIsections,
    size_t isectionTrials
){
    size_t n_new_sets = 0;
    std::vector<elements_container> clots = coOccurenceSingleLinkageClusters( clotThreshold );
    for ( size_t i = 0; i < clots.size(); i++ ) {
        n_new_sets += generateNeighboursIntersections( rng, clots[i], maxIsections, isectionTrials );
    }
    // update inclusion statistics for newly generated clusters
    countPartsInclusion( true );
    countElementPairsAvgCooccurrencePerPart( true );
    return ( n_new_sets );
}

template<class Part>
void IndexedPartitionsCollection<Part>::init(
    const gsl_rng*  rng,
    size_t          independenceThreshold,
    size_t          clotThreshold
){
    countPartitions();
    countParts();
    countElementPairs();
    countElementPairsAvgCooccurrencePerPart( false );
    countPartsInclusion( false );
    _components = coOccurenceSingleLinkageClusters( independenceThreshold );
    if ( rng && std::isfinite( clotThreshold ) ) {
        size_t n_new_sets = generateIntersectionsWithinClots( rng, clotThreshold, 20, 20 );
        LOG_DEBUG1( n_new_sets << " new sets around clots generated" );
    }
    classifyParts();
    _subptnIndexes.resize( _components.size() );
    indexPartitionsComposition();
    _subptnCounts.resize( _components.size() );
    countSubpartitions();
}

} }