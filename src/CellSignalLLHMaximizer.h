#pragma once

#include "BasicTypedefs.h"

#include <gsl/gsl_min.h>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>

#include "symmetric_array2d.h"

#include "OPAData.h"
#include "Signals.h"

class CellSignalLLHMaximizer {
public:
    typedef array2d<signal_t> signals_matrix_type;
    typedef std::vector<size_t>    multiple_map_t;

    CellSignalLLHMaximizer( const OPAData& data,
                            const CellSignalParams& signalParams,
                            ObjectsClusterSignal::distrib_cache_type& scDistribCache,
                            bool  evalPreliminarySignals = true );
    ~CellSignalLLHMaximizer();

    const OPAData& data() const {
        return ( _data );
    }

    const CellSignalParams& signalParams() const {
        return ( _signalParams );
    }
    const signals_matrix_type& preliminarySignals() const {
        return ( _preliminarySignals );
    }

    log_prob_t evalObjectsPrelimDistance( object_index_t obj1ix,
                                      object_index_t obj2ix ) const;
    log_prob_t evalProbesPrelimDistance( probe_index_t probe1ix,
                                      probe_index_t probe2ix ) const;
    symmetric_array2d<log_prob_t> evalObjectCoSignalDistances() const;
    symmetric_array2d<log_prob_t> evalProbeCoSignalDistances() const;

    signal_t clusterSignalA( const object_set_t& objects,
                             const assay_bitset_t& assays,
                             const multiple_map_t& multiples ) const;
    signal_t clusterSignalS( const object_set_t& objects,
                             const probe_bitset_t& probes,
                             const multiple_map_t& multiples ) const
    {
        return ( clusterSignalA( objects, _data.probesToAssays( probes ),  multiples ) );
    }
    log_prob_t maximizeSignalLLHA( const object_set_t& objects,
                                   const assay_bitset_t& assays,
                                   const multiple_map_t& multiples,
                                   signal_t& signal ) const;
    log_prob_t maximizeSignalLLHS( const object_set_t& objects,
                                   const probe_bitset_t& probes,
                                   const multiple_map_t& multiples,
                                   signal_t& signal ) const
    {
        return ( maximizeSignalLLHA( objects, _data.probesToAssays( probes ),  multiples, signal ) );
    }

    ObjectsClusterSignal::distrib_cache_type& scDistribCache() const {
        return ( *_scDistribCache );
    }

    static signals_matrix_type PreliminarySignals( const OPAData& data, const CellSignalParams& signalParams );
    static signals_matrix_type PreliminaryAssaySignals( const OPAData& data, const CellSignalParams& signalParams );

private:
    const OPAData&                  _data;
    const CellSignalParams&         _signalParams;  /** params for preliminary signals */
    ObjectsClusterSignal::distrib_cache_type* _scDistribCache;

    mutable multiple_map_t  _multCache;                 /** object multiples cache */
    signals_matrix_type     _preliminarySignals;        /** preliminary signal estimation for each object-probe pair */
    signals_matrix_type     _preliminaryAssaySignals;   /** preliminary signal estimation for each object-assay pair */
    gsl_min_fminimizer*     _llhMaximizer;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( _preliminarySignals );
        ar & BOOST_SERIALIZATION_NVP( _preliminaryAssaySignals );
    }
};

BOOST_CLASS_IMPLEMENTATION( CellSignalLLHMaximizer, object_serializable )
