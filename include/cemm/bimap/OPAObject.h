#pragma once

#include "BasicTypedefs.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>

namespace cemm { namespace bimap {

/**
    Object of the object-probe-assay data.
    @todo rename to EPAElement
 */
struct OPAObject {
private:
    /// object unique index (internally assigned)
    object_index_t      _index;
    /// object unique label (externally assigned)
    object_label_t      _label;
    /// length of protein AA sequence
    size_t              _sequenceLength;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::make_nvp( "index", _index );
        ar & boost::serialization::make_nvp( "label", _label );
        ar & boost::serialization::make_nvp( "sequenceLength", _sequenceLength );
    }

    OPAObject()
    : _index( OBJECT_NA ), _sequenceLength( 0 )
    {}

protected:
    friend class OPAData;

    OPAObject( object_index_t index, const object_label_t& label, size_t sequenceLength )
    : _index( index ), _label( label ), _sequenceLength( sequenceLength )
    {}

public:
    const object_label_t& label() const {
        return ( _label );
    }
    size_t sequenceLength() const {
        return ( _sequenceLength );
    }
    object_index_t index() const {
        return ( _index );
    }
};

} }

BOOST_CLASS_IMPLEMENTATION( cemm::bimap::OPAObject, object_serializable )
BOOST_CLASS_TRACKING( cemm::bimap::OPAObject, track_never )
