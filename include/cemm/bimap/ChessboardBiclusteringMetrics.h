#pragma once

#include "BasicTypedefs.h"

namespace cemm { namespace bimap {

struct LLHPartitionWeights {
    log_prob_t  topo;
    log_prob_t  conf;

    LLHPartitionWeights()
    : topo( 1.0 ), conf( 1.0 )
    {}

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( topo );
        ar & BOOST_SERIALIZATION_NVP( conf );
    }
};

struct LLHWeights {
    LLHPartitionWeights objects;
    LLHPartitionWeights probes;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( objects );
        ar & BOOST_SERIALIZATION_NVP( probes );
    }
};

/**
 *  Vector of log-likelihood metrics.
 */
struct LLHMetrics {
    log_prob_t  quant;   /** quantitative log likelihood */
    log_prob_t  topo;    /** clusters topological log likelihood */
    log_prob_t  conf;    /** conformance log likelihood */

    LLHMetrics( log_prob_t iniValue = cemm::math::unset<log_prob_t>() )
    : quant( iniValue ), topo( iniValue ), conf( iniValue )
    {}

    void unset() {
        quant = topo = conf = cemm::math::unset<log_prob_t>();
    }
    bool is_unset() const {
        return (    cemm::math::is_unset( quant )
                 || cemm::math::is_unset( topo )
                 || cemm::math::is_unset( conf ) );
    }
    double total( const LLHPartitionWeights& w ) const {
        BOOST_ASSERT( !is_unset() );
        return ( quant + w.topo * topo + w.conf * conf );
    }
    LLHMetrics& operator+=( const LLHMetrics& b ) {
        quant += b.quant;
        topo += b.topo;
        conf += b.conf;
        return ( *this );
    }
    LLHMetrics& operator-=( const LLHMetrics& b ) {
        quant -= b.quant;
        topo -= b.topo;
        conf -= b.conf;
        return ( *this );
    }
    log_prob_t operator()( const LLHPartitionWeights& w ) const {
        return ( total( w ) );
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( quant );
        ar & BOOST_SERIALIZATION_NVP( topo );
        ar & BOOST_SERIALIZATION_NVP( conf );
    }
};

/**
 *  Vector of statistical metrics.
 */
struct StatsMetrics {
    log_prob_t  lppBlocks;      /** blocks log prior probability */
    log_prob_t  lppObjsPtn;     /** object's partition log prior prob. */
    log_prob_t  lppProbesPtn;   /** probe's partition log prior prob. */
    LLHMetrics  llhObjs;        /** object's likelihood metrics + quantitative */
    LLHMetrics  llhProbes;      /** probe's likelihood metrics */

    StatsMetrics()
    : lppBlocks( cemm::math::unset<log_prob_t>() )
    , lppObjsPtn( cemm::math::unset<log_prob_t>() )
    , lppProbesPtn( cemm::math::unset<log_prob_t>() )
    {}

    void unset( bool objects = true, bool probes = true ) {
        // any modification always affect quantitative component
        llhObjs.quant = llhProbes.quant = cemm::math::unset<log_prob_t>();
        // any modification affects blocks LPP
        lppBlocks = cemm::math::unset<log_prob_t>();
        if ( objects ) {
            lppObjsPtn = cemm::math::unset<log_prob_t>();
            llhObjs.topo = llhObjs.conf = cemm::math::unset<log_prob_t>();
            llhProbes.conf = cemm::math::unset<log_prob_t>();
            // probes topo component is not affected
        }
        if ( probes ) {
            lppProbesPtn = cemm::math::unset<log_prob_t>();
            llhProbes.topo = llhProbes.conf = cemm::math::unset<log_prob_t>();
            llhObjs.conf = cemm::math::unset<log_prob_t>();
            // objects topo component is not affected
        }
    }

    bool is_unset() const {
        return ( cemm::math::is_unset( lppBlocks )
                || cemm::math::is_unset( lppObjsPtn )
                || cemm::math::is_unset( lppProbesPtn )
                || llhObjs.is_unset() || llhProbes.is_unset() );
    }

    log_prob_t lpp() const {
        return ( lppBlocks + lppObjsPtn + lppProbesPtn );
    }

    log_prob_t llh( const LLHWeights& w ) const {
        return (   llhObjs.quant
                 + llhObjs.topo * w.objects.topo
                 + llhObjs.conf * w.objects.conf 
                 + llhProbes.topo * w.probes.topo
                 + llhProbes.conf * w.probes.conf );
    }

    log_prob_t totalLnP( const LLHWeights& w ) const {
        BOOST_ASSERT( !is_unset() );
        return ( lpp() + llh( w ) );
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( lppBlocks );
        ar & BOOST_SERIALIZATION_NVP( lppObjsPtn );
        ar & BOOST_SERIALIZATION_NVP( lppProbesPtn );
        ar & BOOST_SERIALIZATION_NVP( llhObjs );
        ar & BOOST_SERIALIZATION_NVP( llhProbes );
    }
};

} }