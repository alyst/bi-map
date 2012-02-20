#pragma once

#include "BasicTypedefs.h"

#include <cemm/math/Distributions.h>
#include <cemm/containers/OrderedDistanceMatrix.h>
#include <cemm/containers/dynamic_bitset_foreach.h>

namespace cemm { namespace bimap {

template<typename ElementIndex>
class CoOccurenceDistanceMatrix: public OrderedDistanceMatrix<ElementIndex, log_prob_t>
{
public:
    typedef boost::dynamic_bitset<> obs_mask_type;
    typedef size_t obs_ix_type;
    typedef OrderedDistanceMatrix<ElementIndex, log_prob_t> super_type;
    typedef typename super_type::distance_matrix_type distance_matrix_type;
    typedef typename super_type::elm_ix_type elm_ix_type;

    using super_type::operator=;

    static obs_mask_type SpecificHits(
        const std::vector<obs_mask_type>&   observations,
        double                              freqThreshold
    );

    static log_prob_t ObservationsBias(
        const std::vector<obs_mask_type>&   observations,
        const obs_mask_type&                mask
    );

    static distance_matrix_type CoOccurenceDistances(
        const std::vector<obs_mask_type>&   observations,
        const obs_mask_type&                mask,
        log_prob_t                          observationsBias );

    CoOccurenceDistanceMatrix()
    {}

    CoOccurenceDistanceMatrix(
        const std::vector<obs_mask_type>&   observations,
        double                              freqThreshold,
        bool                                fixObservationsBias )
    : _observations( observations )
    , _specificObservations( SpecificHits( observations, freqThreshold ) )
    , _observationsBias( fixObservationsBias ? ObservationsBias( _observations, _specificObservations ) : 0 )
    {
        operator=( CoOccurenceDistances( _observations, _specificObservations, _observationsBias ) );
    }

    double observationsBias() const {
        return ( _observationsBias );
    }

    const obs_mask_type& specificObservations() const {
        return ( _specificObservations );
    }

    const std::vector<obs_mask_type>& allObservations() const {
        return ( _observations );
    }

    const obs_mask_type& observations( elm_ix_type elmIx ) const {
        return ( _observations[ elmIx ] );
    }

private:
    std::vector<obs_mask_type>  _observations;          /** observations per element */
    obs_mask_type               _specificObservations;  /** masks for specific observations */
    double                      _observationsBias;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp( "orderedDistances", boost::serialization::base_object<super_type>( *this ) );
        ar & BOOST_SERIALIZATION_NVP( _observations );
        ar & BOOST_SERIALIZATION_NVP( _specificObservations );
        ar & BOOST_SERIALIZATION_NVP( _observationsBias );
    }
};

/**
 *  Mask of frequent observations.
 */
template<typename ElementIndex>
typename CoOccurenceDistanceMatrix<ElementIndex>::obs_mask_type
CoOccurenceDistanceMatrix<ElementIndex>::SpecificHits(
    const std::vector<obs_mask_type>&       observations,
    double                                  freqThreshold
){
    // count total objects and probes occurrences
    std::vector<size_t> hitCounts( observations.front().size(), 0 );
    for ( elm_ix_type elmIx = 0; elmIx < observations.size(); elmIx++ ) {
        const obs_mask_type& elmHits = observations[ elmIx ];
        foreach_bit( obs_ix_type, hitIx, elmHits ) {
            hitCounts[ hitIx ]++;
        }
    }

    // set masks by the number of occurrences
    obs_mask_type specificHits( hitCounts.size() );
    for ( obs_ix_type obsIx = 0; obsIx < hitCounts.size(); obsIx++ ) {
        if ( hitCounts[ obsIx ] <= freqThreshold * observations.size() ) {
            specificHits.set( obsIx );
        }
    }
    return ( specificHits );
}

/**
 *  Sampling bias for Fisher noncentral hypergeometric (Fisher) distribution.
 *  A bias to simultaneously observe hit of the same element in different experiments
 *  (cleared of ratio of hits in the experiments) is calculated.
 */
template<typename ElementIndex>
log_prob_t CoOccurenceDistanceMatrix<ElementIndex>::ObservationsBias(
    const std::vector<obs_mask_type>&   observations,
    const obs_mask_type&                mask
){
    if ( mask.count() == 0 ) return ( 0 );
    size_t  nTotal = observations.size() * ( observations.size() - 1 ) * mask.count(); // number of all observations (multiplied across pairs of experiments)
    size_t  nHits = 0;
    size_t  nHitPairs = 0;
    for ( size_t i = 0; i < observations.size(); i++ ) {
        const obs_mask_type& iHits = observations[ i ];
        nHits += ( iHits & mask ).count();
        for ( size_t j = 0; j < observations.size(); j++ ) {
            if ( j != i ) {
                const obs_mask_type& jHits = observations[ j ];
                nHitPairs += ( iHits & jHits & mask ).count();
            }
        }
    }
    nHits *= observations.size() - 1;
    // observed rate of paired hits : nHitPairs/nHits
    // versus 0 hypothesis of hit independence : nHits / nTotal
    // versus
    // observed rate of not picking hit pair : (nHits - nHitPairs)/nHits
    // versus 0 hypothesis of hit independence : (nTotal - nHits)/ nTotal
    // add 0.01 to avoid degeneration
    return (   log( nHitPairs * ( nTotal + 0.01 - nHits ) )
             - log( nHits * ( nHits + 0.01 - nHitPairs ) ) );
}

/**
 *   Distances between vectors of observations of element pairs,
 *   based on Fisher hypergeometric p-value
 *   of elements co-observation.
 */
template<typename ElementIndex>
typename CoOccurenceDistanceMatrix<ElementIndex>::distance_matrix_type
CoOccurenceDistanceMatrix<ElementIndex>::CoOccurenceDistances(
    const std::vector<obs_mask_type>&   observations,       /** hit/miss array of observations */
    const obs_mask_type&                mask,               /** mask of observations */
    log_prob_t                          observationsBias    /** array of distances to fill */
){
    // evaluate matrix
    LOG_INFO( "Calculating co-occurrence distances..." );
    symmetric_array2d<log_prob_t>   distances( observations.size() );
    size_t maskedObservations = mask.count();
    for ( size_t elm1Ix = 0; elm1Ix < distances.size(); ++elm1Ix ) {
        obs_mask_type elm1Obs = observations[ elm1Ix ] & mask;
        size_t elm1ObsCount = elm1Obs.count();
        for ( size_t elm2Ix = elm1Ix; elm2Ix < distances.size(); ++elm2Ix ) {
            obs_mask_type elm2Obs = observations[ elm2Ix ] & mask;
            distances( elm1Ix, elm2Ix ) = FisherHypergeometricDistribution( elm1ObsCount, maskedObservations - elm1ObsCount,
                                                                            elm2Obs.count(), observationsBias )
                .lnCdf_Q( ( elm1Obs & elm2Obs ).count(), true );
        }
    }
    return ( distances );
}

} }