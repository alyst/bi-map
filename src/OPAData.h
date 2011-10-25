#pragma once

#include "BasicTypedefs.h"

#include <exception>
#include <boost/unordered_map.hpp>

#include <boost/shared_ptr.hpp>

#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/ptr_container/serialize_ptr_vector.hpp>

#include "array2d.h"
#include "ObjectSet.h"
#include "ProbeSet.h"

#include "OPAObject.h"
#include "OPAAssay.h"

typedef std::vector<assay_index_t>      assay_container_t;

class OPAData;

class entity_error: public std::runtime_error {
protected:
    std::string compose_error_msg( const std::string& entityClass, const std::string& entityName, const std::string& message );

public:
    entity_error( const std::string& entityClass, const std::string& entityName, const std::string& message )
    : std::runtime_error( compose_error_msg( entityClass, entityName, message ) )
    {
    }
};

class entity_not_found: public entity_error {
public:
    entity_not_found( const std::string& entityClass, const std::string& entityName, const std::string& message = "not found" )
    : entity_error( entityClass, entityName, message )
    {
    }
};

class entity_already_exists: public entity_error {
public:
    entity_already_exists( const std::string& entityClass, const std::string& entityName, const std::string& message = "already exists" )
    : entity_error( entityClass, entityName, message )
    {
    }
};

/**
    Probe of object-probe-assay data.

    @remark
    The more generic approach is to decouple probe of the system
    and using specific object as bait.
    That is, there could be several samples of the same probe,
    but with different bait (each sample could have several assays).
 */
struct OPAProbe {
private:
    /// probe unique index (internally assigned)
    probe_index_t                   _index;
    /// probe unique label (externally assigned)
    probe_label_t                   _label;
    /// assays of the probe
    assay_container_t               _assayIndices;
    /// bait object of the probe
    object_index_t                  _baitIndex;

    friend class boost::serialization::access;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> boost::serialization::make_nvp( "index", _index );
        ar >> boost::serialization::make_nvp( "label", _label );
        object_index_t baitIndex;
        ar >> boost::serialization::make_nvp( "baitIndex", baitIndex );
        _baitIndex = baitIndex > 0 ? ( baitIndex - 1 ) : OBJECT_NA;
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << boost::serialization::make_nvp( "index", _index );
        ar << boost::serialization::make_nvp( "label", _label );
        object_index_t baitIndex = _baitIndex != OBJECT_NA ? _baitIndex + 1 : 0;
        ar << boost::serialization::make_nvp( "baitIndex", baitIndex );
    }

    OPAProbe()
    : _index( PROBE_NA ), _baitIndex( OBJECT_NA )
    {}

protected:
    friend class OPAData;

    void addAssay( assay_index_t assayIx ) {
        if ( !_assayIndices.empty() && _assayIndices.back() != (assayIx-1) ) {
            THROW_EXCEPTION( std::invalid_argument, "Assay #" << assayIx << " do not immediately follow " << _assayIndices.back() );
        }
        _assayIndices.push_back( assayIx );
    }

    OPAProbe( probe_index_t index, const probe_label_t& label, object_index_t baitIndex )
    : _index( index ), _label( label ), _baitIndex( baitIndex )
    {
    }

public:
    const probe_label_t& label() const {
        return ( _label );
    }
    const assay_container_t& assayIndexes() const {
        return ( _assayIndices );
    }

    /**
     *  Bait is the object
     *  used to construct the probe, it should
     *  be present in all sets, observed in this projection
     */
    object_index_t baitIndex() const {
        return ( _baitIndex );
    }
    probe_index_t index() const {
        return ( _index );
    }
};

BOOST_CLASS_IMPLEMENTATION( OPAProbe, object_serializable )
BOOST_CLASS_TRACKING( OPAProbe, track_never )

/**
    Object-Probe-Assay Input data.
 */
class OPAData {
public:
    struct MassSpectraData {
        size_t  sc; // spectra counts
        size_t  pc; // peptide counts

        MassSpectraData( size_t sc = 0, size_t pc = 0 )
        : sc( sc ), pc( pc )
        {
            if ( pc > sc ) throw std::invalid_argument( "Peptide count less than spectra count" );
        }

        MassSpectraData& operator=( const MassSpectraData& msData )
        {
            sc = msData.sc;
            pc = msData.pc;
            return ( *this );
        }
        MassSpectraData& operator=( size_t counts )
        {
            sc = counts;
            pc = 0;
            return ( *this );
        }
        operator size_t() const {
            return ( sc );
        }
        bool operator==( size_t that ) const {
            return ( pc == 0 && that == sc );
        }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version) {
            ar & BOOST_SERIALIZATION_NVP( sc );
            ar & BOOST_SERIALIZATION_NVP( pc );
        }
    };
    typedef MassSpectraData             celldata_t;
    typedef array2d<celldata_t>         data_matrix_t;
    typedef size_t                      size_type;
    typedef std::vector<size_t>         hit_counts_type;

private:
    typedef OPAProbe*                   probe_ptr_t;
    typedef OPAAssay*                   assay_ptr_t;
    typedef OPAObject*                  object_ptr_t;
    typedef boost::dynamic_bitset<>     hit_mask_type;

    typedef boost::ptr_vector<OPAProbe>                             probe_vector_t;
    typedef boost::ptr_vector<OPAObject>                            object_vector_t;
    typedef std::vector<OPAAssay>                                   assay_vector_t;
    typedef boost::unordered_map<probe_label_t, probe_ptr_t>        probe_by_label_map_t;
    typedef boost::unordered_map<object_label_t, object_ptr_t>      object_by_label_map_t;
    typedef boost::unordered_map<assay_index_t, probe_index_t>      assay_to_probe_map_t;
    typedef boost::unordered_map<assay_label_t, assay_index_t>      assay_label_to_index_map_t;
    typedef boost::unordered_multimap<object_index_t, probe_index_t>   bait_to_probe_mmap_t;

    probe_vector_t              _probes;
    object_vector_t             _objects;

    assay_vector_t              _assays;

    data_matrix_t               _matrix;

    mutable bool                _hitsDirty;
    /**
     *  For each object a binary vector of assays size -- 1 if this object is seen in given assay.
     */
    mutable std::vector<hit_mask_type>  _objectHits;
    /**
     *  For each probe a binary vector of object size -- 1 if this object is seen in all probe's assays.
     */
    mutable std::vector<hit_mask_type>  _probeHits;

    probe_by_label_map_t        _probeLabelMap;
    object_by_label_map_t       _objectLabelMap;
    assay_label_to_index_map_t  _assaysLabelMap;
    bait_to_probe_mmap_t        _baitToProbes;

protected:
    object_ptr_t object( const object_label_t& objLabel ) {
        object_by_label_map_t::const_iterator it = _objectLabelMap.find( objLabel );
        return ( it != _objectLabelMap.end() ? it->second : NULL );
    }
    probe_ptr_t probe( const probe_label_t& probeLabel ) {
        probe_by_label_map_t::const_iterator it = _probeLabelMap.find( probeLabel );
        return ( it != _probeLabelMap.end() ? it->second : NULL );
    }

    void initMatrixDataContainers();
    void resetHits() const;
    void resetIndexes();

    friend class boost::serialization::access;

    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> boost::serialization::make_nvp( "probes", _probes );
        ar >> boost::serialization::make_nvp( "objects", _objects );
        ar >> boost::serialization::make_nvp( "assays", _assays );
        ar >> boost::serialization::make_nvp( "matrix", _matrix );
        resetIndexes();
        resetHits();
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << boost::serialization::make_nvp( "probes", _probes );
        ar << boost::serialization::make_nvp( "objects", _objects );
        ar << boost::serialization::make_nvp( "assays", _assays );
        ar << boost::serialization::make_nvp( "matrix", _matrix );
    }

public:
    typedef OPAProbe const*                 const_probe_ptr_t;
    typedef OPAObject const*                const_object_ptr_t;
    typedef bait_to_probe_mmap_t::const_iterator const_bait_to_probe_iterator;

    OPAData();

    size_type objectsCount() const {
        return ( _objects.size() );
    }

    const OPAObject& object( object_index_t objIx ) const {
        return ( _objects[ objIx ] );
    }
    const_object_ptr_t object( const object_label_t& objLabel ) const {
        object_by_label_map_t::const_iterator it = _objectLabelMap.find( objLabel );
        return ( it != _objectLabelMap.end() ? const_object_ptr_t( it->second ) : const_object_ptr_t() );
    }
    const_object_ptr_t addObject( const object_label_t& label, size_t seqLength );

    size_type probesCount() const {
        return ( _probes.size() );
    }

    const OPAProbe& probe( probe_index_t probeIx ) const {
        return ( _probes[ probeIx ] );
    }
    const_probe_ptr_t probe( const probe_label_t& probeLabel ) const {
        probe_by_label_map_t::const_iterator it = _probeLabelMap.find( probeLabel );
        return ( it != _probeLabelMap.end() ? const_probe_ptr_t( it->second ) : const_probe_ptr_t() );
    }

    std::pair<const_bait_to_probe_iterator, const_bait_to_probe_iterator>
    probesOfBait( object_index_t baitIx ) const {
        return ( _baitToProbes.equal_range( baitIx ) );
    }
    const_probe_ptr_t addProbe( const probe_label_t& label, const object_label_t& baitLabel );

    assay_bitset_t probesToAssays( const probe_bitset_t& probes ) const;

    size_type assaysCount() const {
        return ( _assays.size() );
    }
    size_type assaysCount( probe_index_t probeIx ) const;
    size_type assaysCount( const probe_bitset_t& probes ) const;

    assay_index_t assayIndex( const assay_label_t& label ) const {
        assay_label_to_index_map_t::const_iterator ait = _assaysLabelMap.find( label );
        return ( ait != _assaysLabelMap.end() ? ait->second : ASSAY_NA );
    }
    const OPAAssay& assay( assay_index_t assayIx ) const {
        return ( _assays[ assayIx ] );
    }
    assay_index_t addAssay( const assay_label_t& label, const probe_label_t& probeLabel, prob_t multiplier = 1.0 );

    const celldata_t& measurement( object_index_t objIx, assay_index_t assayIx ) const {
        return ( _matrix( objIx, assayIx ) );
    }

    /**
     *  Pointer to the first measurement of given object in given probe.
     */
    const celldata_t* measurements( object_index_t objIx, probe_index_t probeIx ) const {
        return ( _matrix.data( objIx, _probes[ probeIx ]._assayIndices.front() ) );
    }

    void addMeasurement( const object_label_t& objectLabel, const assay_label_t& assayLabel, celldata_t measurement );

    const celldata_t& maxMeasurement() const;

    object_set_t probesToBaits( const probe_bitset_t& probes ) const;

    hit_counts_type getObjectsHitCounts( const object_set_t& objects ) const;
    hit_counts_type getProbesHitCounts( const probe_bitset_t& probes ) const;

    objects_label_map_type objectsLabelMap() const;
    probes_label_map_type probesLabelMap() const;

    static OPAData load( const char* filename );
    void save( const char* filename ) const;
};

BOOST_CLASS_IMPLEMENTATION( OPAData::MassSpectraData, object_serializable )
BOOST_CLASS_IMPLEMENTATION( OPAData, object_serializable )

