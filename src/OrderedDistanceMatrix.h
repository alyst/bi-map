#pragma once

#include "BasicTypedefs.h"

#include "array2d.h"
#include "symmetric_array2d.h"

template<typename ElementIndex, typename Distance>
class OrderedDistanceMatrix {
public:
    typedef Distance dist_type;
    typedef ElementIndex elm_ix_type;
    typedef std::pair<elm_ix_type, dist_type>   dist_to_elm_type;       /** element and distance to it from given one */
    typedef symmetric_array2d<dist_type>        distance_matrix_type;
    typedef array2d<dist_to_elm_type>           distance_container_type;
    typedef OrderedDistanceMatrix<ElementIndex, Distance> ordered_distance_matrix_type;

    OrderedDistanceMatrix()
    {
    }

    OrderedDistanceMatrix( const distance_matrix_type& distances )
    : _data( distances.size(), distances.size() - 1, dist_to_elm_type( -1, unset() ) )
    {
        SortDistances( distances, _data );
    }

    ordered_distance_matrix_type& swap( ordered_distance_matrix_type& that )
    {
        std::swap( _data, that._data );
        return ( that );
    }

    ordered_distance_matrix_type& operator=( const distance_matrix_type& distances )
    {
        _data = ordered_distance_matrix_type( distances )._data;
        return ( *this );
    }

    static void SortDistances(
        const distance_matrix_type& distances,
        distance_container_type&    sortedDistances
    );

    const dist_to_elm_type& operator()( elm_ix_type elmIx, size_t rank ) const {
        return ( _data( elmIx, rank ) );
    }

    template<class Context>
    dist_to_elm_type rankedDistance(
        elm_ix_type elmIx, const Context& context,
        bool withinContext, size_t rank, bool ascending ) const;

    size_t size() const {
        return ( _data.size1() );
    }


private:
    friend class boost::serialization::access;

    struct DistToElmLess {
    bool operator()(
        const dist_to_elm_type& a,
        const dist_to_elm_type& b
    ) const {
        return ( a.second < b.second );
        }
    };

    distance_container_type    _data;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( _data );
    }
};

template<typename ElementIndex, typename Distance>
template<class Context>
typename OrderedDistanceMatrix<ElementIndex, Distance>::dist_to_elm_type
OrderedDistanceMatrix<ElementIndex, Distance>::rankedDistance(
    elm_ix_type     elmIx,
    const Context&  context,
    bool            withinContext,
    size_t          rank,
    bool            ascending
) const {
    BOOST_ASSERT( elmIx >= 0 && elmIx < _data.size1() );
    dist_to_elm_type res( (elm_ix_type)-1, unset() );

    size_t  curRank = 0;
    for ( int i = ascending ? 0 : (int)_data.size2() - 1;
          ascending ? i < (int)_data.size2() : i >= 0;
          ascending ? i++ : i--
    ){
        const dist_to_elm_type& d2e = _data( elmIx, i );
        if ( withinContext == context( d2e.first ) ){
            res = d2e;
            if ( curRank == rank ) return ( res );
            curRank++;
        }
    }
    return ( res );
}

/**
 *  Sort distances for each element.
 */
template<typename ElementIndex, typename Distance>
void OrderedDistanceMatrix<ElementIndex, Distance>::SortDistances(
    const distance_matrix_type& distances,
    distance_container_type&    sortedDistances
){
    sortedDistances.reset( distances.size(), distances.size() - 1, dist_to_elm_type( -1, unset() ) );
    std::vector<dist_to_elm_type> sortVec( distances.size() - 1 );
    for ( size_t elm1Ix = 0; elm1Ix < distances.size(); ++elm1Ix ) {
        sortVec.clear();
        for ( size_t elm2Ix = 0; elm2Ix < distances.size(); ++elm2Ix ) {
            if ( elm2Ix != elm1Ix ) {
                sortVec.push_back( dist_to_elm_type( elm2Ix, distances( elm1Ix, elm2Ix ) ) );
            }
        }
        // sort by distance ascending
        std::sort( sortVec.begin(), sortVec.end(), DistToElmLess() );
        for ( object_index_t i = 0; i < sortVec.size(); ++i ) {
            sortedDistances( elm1Ix, i ) = sortVec[ i ];
        }
    }
}
