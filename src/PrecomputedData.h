#pragma once

#include "BasicTypedefs.h"

#include <boost/unordered_map.hpp>
#include <boost/dynamic_bitset.hpp>

#include "statically_tracked.h"

#include "OPAData.h"
#include "Signals.h"
#include "CoOccurrenceGraph.h"

#include "CellSignalLLHMaximizer.h"
#include "OrderedDistanceMatrix.h"
#include "CoOccurenceDistanceMatrix.h"

#include "math/DistributionCache.h"

struct ObjectsContext
{
    const object_set_t& objects;

    ObjectsContext( const object_set_t& objects )
    : objects( objects )
    {}

    bool operator()( object_index_t objIx ) const {
        return ( objects.find( objIx ) != objects.end() );
    }
};

struct ProbesContext
{
    const probe_bitset_t& probes;

    ProbesContext( const probe_bitset_t& probes )
    : probes( probes )
    {}

    bool operator()( probe_index_t probeIx ) const {
        return ( probes.test( probeIx ) );
    }
};

/**
 *  Parameters for precomputing data.
 *
 *  @see PrecomputedData
 */
struct PrecomputedDataParams {
    prob_t    objectFreqThreshold;      /** max frequency of object (for probes neighbourhood) */
    prob_t    probeFreqThreshold;       /** max frequency of probe (for objects neighbourhood) */

    prob_t    objectNeighbourThreshold; /** max p-value for objects to be neighbours */
    prob_t    probeNeighbourThreshold;  /** max p-value for probes to be neighbours */

    PrecomputedDataParams();

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( objectFreqThreshold );
        ar & BOOST_SERIALIZATION_NVP( probeFreqThreshold );
        ar & BOOST_SERIALIZATION_NVP( objectNeighbourThreshold );
        ar & BOOST_SERIALIZATION_NVP( probeNeighbourThreshold );
    }
};

/**
 *  Preliminary cell signals to aid the signals sampling.
 *  These signals are generated individually for each object x probe pair,
 *  regardless the absence of spectral counts.
 */
class PrecomputedData {
public:
    typedef array2d<signal_t> signals_matrix_type;

    typedef OrderedDistanceMatrix<object_index_t, log_prob_t> obj_cosignal_matrix_t;
    typedef OrderedDistanceMatrix<probe_index_t, log_prob_t>  probe_cosignal_matrix_t;
    typedef CoOccurenceDistanceMatrix<object_index_t> obj_cohits_matrix_t;
    typedef CoOccurenceDistanceMatrix<probe_index_t>  probe_cohits_matrix_t;
    typedef obj_cosignal_matrix_t::dist_to_elm_type dist_to_obj_t;
    typedef probe_cosignal_matrix_t::dist_to_elm_type dist_to_probe_t;
    typedef obj_cohits_matrix_t::obs_mask_type observations_mask_t;
    typedef ObjectsClusterSignal::distrib_cache_type sc_cache_type;

    PrecomputedData( const OPAData& data, const PrecomputedDataParams& params,
                     const CellSignalParams& signalParams );

    const OPAData& data() const {
        return ( _data );
    }

    const CellSignalParams& signalParams() const {
        return ( _llhMaximizer.signalParams() );
    }

    const signals_matrix_type& preliminarySignals() const {
        return ( _llhMaximizer.preliminarySignals() );
    }

    const obj_cosignal_matrix_t& objectCoSignalDistances() const {
        return ( _objectCoSignalDistances );
    }
    const probe_cosignal_matrix_t& probeCoSignalDistances() const {
        return ( _probeCoSignalDistances );
    }
    const obj_cohits_matrix_t& objectCoHitsDistances() const {
        return ( _objectCoHitsDistances );
    }
    const probe_cohits_matrix_t& probeCoHitsDistances() const {
        return ( _probeCoHitsDistances );
    }

    const co_occurrence_graph_t& objectsNeighbourhood() const {
        return ( _objectsNeighbourhood );
    }

    const co_occurrence_graph_t& probesNeighbourhood() const {
        return ( _probesNeighbourhood );
    }

    typedef CellSignalLLHMaximizer::multiple_map_t multiple_map_t;

    signal_t clusterSignal( const object_set_t& objects,
                            const probe_bitset_t& probes,
                            const multiple_map_t& multiples ) const
    {
        return ( _llhMaximizer.clusterSignalS( objects, probes, multiples ) );
    }

    log_prob_t maximizeSignalLLH( const object_set_t& objects,
                                  const probe_bitset_t& probes,
                                  const multiple_map_t& multiples,
                                  signal_t& signal0 ) const
    {
        return ( _llhMaximizer.maximizeSignalLLHS( objects, probes, multiples, signal0 ) );
    }

    const observations_mask_t& objectsMask( probe_index_t probeIx ) const {
        return ( _probeCoHitsDistances.observations( probeIx ) );
    }
    const observations_mask_t& probesMask( object_index_t objIx ) const {
        return ( _objectCoHitsDistances.observations( objIx ) );
    }

    observations_mask_t objectsUnion( const probe_bitset_t& probes ) const;
    observations_mask_t probesUnion( const object_set_t& objects ) const;

    const observations_mask_t& specificObjectsMask() const {
        return ( _probeCoHitsDistances.specificObservations() );
    }
    const observations_mask_t& specificProbesMask() const {
        return ( _objectCoHitsDistances.specificObservations() );
    }
    dist_to_obj_t nearestExternalObjectHitsMatch( const object_set_t& objects ) const;
    dist_to_probe_t nearestExternalProbeHitsMatch( const probe_bitset_t& probes ) const;

    dist_to_obj_t rankedObject( object_index_t objIx, const object_set_t& context,
                                bool withinContext, size_t rank, bool ascending,
                                bool signalDistance ) const
    {
        return ( ( signalDistance ? _objectCoSignalDistances : _objectCoHitsDistances )
                 .rankedDistance( objIx, ObjectsContext( context ),
                                  withinContext, rank, ascending ) );
    }
    dist_to_probe_t rankedProbe( probe_index_t probeIx, const probe_bitset_t& context,
                                bool withinContext, size_t rank, bool ascending,
                                bool signalDistance ) const
    {
        return ( ( signalDistance ? _probeCoSignalDistances : _probeCoHitsDistances )
                 .rankedDistance( probeIx, ProbesContext( context ),
                                  withinContext, rank, ascending ) );
    }

    sc_cache_type& scDistribCache() const {
        return ( _scDistribCache );
    }

    static std::vector<observations_mask_t> ProbesObservations( const OPAData& data );
    static std::vector<observations_mask_t> ObjectsObservations( const OPAData& data );

private:
    const OPAData&          _data;
    mutable sc_cache_type   _scDistribCache;
    CellSignalLLHMaximizer  _llhMaximizer;

    obj_cohits_matrix_t     _objectCoHitsDistances;
    probe_cohits_matrix_t   _probeCoHitsDistances;
    obj_cosignal_matrix_t   _objectCoSignalDistances;
    probe_cosignal_matrix_t _probeCoSignalDistances;

    co_occurrence_graph_t   _objectsNeighbourhood;  /** objects sampling neighbourhood */
    co_occurrence_graph_t   _probesNeighbourhood;   /** probes sampling neighbourhood */

    friend class boost::serialization::access;

    PrecomputedData( const OPAData& data, const CellSignalParams& signalParams );

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( _llhMaximizer );
        ar & BOOST_SERIALIZATION_NVP( _objectCoHitsDistances );
        ar & BOOST_SERIALIZATION_NVP( _probeCoHitsDistances );
        ar & BOOST_SERIALIZATION_NVP( _objectCoSignalDistances );
        ar & BOOST_SERIALIZATION_NVP( _probeCoSignalDistances );
        ar & BOOST_SERIALIZATION_NVP( _objectsNeighbourhood );
        ar & BOOST_SERIALIZATION_NVP( _probesNeighbourhood );
    }

    template<class Archive>
    friend inline void load_construct_data(
        Archive & ar, PrecomputedData* t, const unsigned int file_version
    ){
        const OPAData* pData;
        statically_tracked<OPAData> tData( "data", pData );
        ar >> boost::serialization::make_nvp( "data", tData );

        const CellSignalParams* pSignalParams;
        statically_tracked<CellSignalParams> tSignalParams( "signalParams", pSignalParams );
        ar >> boost::serialization::make_nvp( "signalParams", tSignalParams );

        // invoke inplace constructor to initialize instance of my_class
        ::new(t)PrecomputedData( *pData, *pSignalParams );
    }

    template<class Archive>
    friend inline void save_construct_data(
        Archive & ar, const PrecomputedData* t, const unsigned int file_version
    ){
        // save data required to construct instance
        const OPAData* pData = &t->data();
        statically_tracked<OPAData> tData( "data", pData );
        ar << boost::serialization::make_nvp( "data", tData );

        const CellSignalParams* pSignalParams = &t->signalParams();
        statically_tracked<CellSignalParams> tSignalParams( "signalParams", pSignalParams );
        ar << boost::serialization::make_nvp( "signalParams", tSignalParams );
    }
};

BOOST_CLASS_IMPLEMENTATION( PrecomputedData, object_serializable )
