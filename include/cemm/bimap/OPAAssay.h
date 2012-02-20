#pragma once

#include "BasicTypedefs.h"

#include <boost/serialization/serialization.hpp>

namespace cemm { namespace bimap {

typedef boost::dynamic_bitset<> assay_bitset_t;

/**
    Assay (measurement) of all objects of system in specific probe.
    @todo rename to EPAAssay
 */
struct OPAAssay {
private:
    friend class OPAData;
    friend class boost::serialization::access;

    probe_index_t   _probeIndex;
    assay_index_t   _index;
    assay_label_t   _label;
    prob_t          _multiplier;
    log_prob_t      _lnMultiplier;

    OPAAssay()
    : _probeIndex( PROBE_NA ), _index( (assay_index_t)(-1) )
    , _multiplier( 1.0 ), _lnMultiplier( 0 )
    {}

    OPAAssay( probe_index_t probeIndex, assay_index_t index, const assay_label_t& label, prob_t multiplier = 1.0 )
    : _probeIndex( probeIndex )
    , _index( index ), _label( label )
    , _multiplier( multiplier )
    , _lnMultiplier( log( multiplier ) )
    {
        if ( _multiplier <= 0 ) {
            throw std::invalid_argument( "Multiplier must be positive" );
        }
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::make_nvp( "probeIndex", _probeIndex );
        ar & boost::serialization::make_nvp( "index", _index );
        ar & boost::serialization::make_nvp( "label", _label );
        ar & boost::serialization::make_nvp( "multiplier", _multiplier );
        _lnMultiplier = log( _multiplier );
    }

public:
    /**
     *  @todo rename to projectionIndex
     */
    probe_index_t probeIndex() const {
        return ( _probeIndex );
    }

    assay_index_t index() const {
        return ( _index );
    }

    assay_label_t label() const {
        return ( _label );
    }

    prob_t multiplier() const {
        return ( _multiplier );
    }
    log_prob_t lnMultiplier() const {
        return ( _lnMultiplier );
    }
};

} }

BOOST_CLASS_IMPLEMENTATION( cemm::bimap::OPAAssay, object_serializable )
BOOST_CLASS_TRACKING( cemm::bimap::OPAAssay, track_never )
