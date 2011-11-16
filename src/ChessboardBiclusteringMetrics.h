#pragma once

#include "BasicTypedefs.h"

/**
 *  Vector of log-likelihood metrics.
 */
struct LLHMetrics {
    log_prob_t  quant;   /** quantitative log likelihood */
    log_prob_t  topo;    /** clusters topological log likelihood */
    log_prob_t  conf;    /** conformance log likelihood */

    LLHMetrics( log_prob_t iniValue = ::unset() )
    : quant( iniValue ), topo( iniValue ), conf( iniValue )
    {}

    void unset() {
        quant = topo = conf = ::unset();
    }
    bool is_unset() const {
        return (    ::is_unset( quant )
                 || ::is_unset( topo )
                 || ::is_unset( conf ) );
    }
    double total() const {
        BOOST_ASSERT( !is_unset() );
        return ( quant + topo + conf );
    }
    LLHMetrics& operator+=( const LLHMetrics& b ) {
        quant += b.quant;
        topo += b.topo;
        conf += b.conf;
    }
    LLHMetrics& operator-=( const LLHMetrics& b ) {
        quant -= b.quant;
        topo -= b.topo;
        conf -= b.conf;
    }
    operator log_prob_t() const {
        return ( total() );
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
    log_prob_t  lppObjClu;      /** object's clustering log prior prob. */
    log_prob_t  lppProbesClu;   /** probe's clustering log prior prob. */
    LLHMetrics  llhObjs;        /** object's likelihood metrics + quantitative */
    LLHMetrics  llhProbes;      /** probe's likelihood metrics */

    StatsMetrics()
    : lppBlocks( ::unset() ), lppObjClu( ::unset() ), lppProbesClu( ::unset() )
    {}

    void unset( bool objects = true, bool probes = true ) {
        // any modification always affect quantitative component
        llhObjs.quant = llhProbes.quant = ::unset();
        // any modification affects blocks LPP
        lppBlocks = ::unset();
        if ( objects ) {
            lppObjClu = ::unset();
            llhObjs.topo = llhObjs.conf = ::unset();
            llhProbes.conf = ::unset();
            // probes topo component is not affected
        }
        if ( probes ) {
            lppProbesClu = ::unset();
            llhProbes.topo = llhProbes.conf = ::unset();
            llhObjs.conf = ::unset();
            // objects topo component is not affected
        }
    }

    bool is_unset() const {
        return ( ::is_unset( lppBlocks ) || ::is_unset( lppObjClu ) || ::is_unset( lppProbesClu )
                    || llhObjs.is_unset() || llhProbes.is_unset() );
    }

    log_prob_t lpp() const {
        return ( lppBlocks + lppObjClu + lppProbesClu );
    }

    log_prob_t llh() const {
        return ( llhObjs.total() + llhProbes.conf + llhProbes.topo );
    }

    log_prob_t totalLnP() const {
        BOOST_ASSERT( !is_unset() );
        return ( lpp() + llh() );
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( lppBlocks );
        ar & BOOST_SERIALIZATION_NVP( lppObjClu );
        ar & BOOST_SERIALIZATION_NVP( lppProbesClu );
        ar & BOOST_SERIALIZATION_NVP( llhObjs );
        ar & BOOST_SERIALIZATION_NVP( llhProbes );
    }
};

