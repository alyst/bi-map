#pragma once

#include "BasicTypedefs.h"

#include <boost/serialization/serialization.hpp>

#include "cemm/math/LagrangianPoisson.h"
#include "cemm/math/Distributions.h"
#include "cemm/math/DistributionCache.h"

#include "OPAAssay.h"
#include "OPAObject.h"

namespace cemm { namespace bimap {

/**
    Parameters of data affecting signal inference.
 */
class CellSignalParams {
public:
    typedef GammaDistribution preliminary_signal_prior_type;
    typedef ZeroInflatedLagrangianPoissonDistribution sc_distrib_type;
    typedef sc_distrib_type distribution_type;

    double                          sequenceLengthFactor;   /// the power to raise sequence length to for normalizing
    double                          scShape;                /// shape parameter for object's spectral counts distribution
    signal_t                        scRateStep;             /// discretization step for spectral counts distribution parameter
    signal_t                        zeroInflationFactor;    /// zero-inflation coefficient
    preliminary_signal_prior_type   preliminarySignalPrior; /// prior distribution of preliminary signal rate

    CellSignalParams( log_prob_t sequenceLengthFactor = 0.5,
                      log_prob_t scShape = 0,
                      signal_t   scRateStep = 1E-2,
                      signal_t   zeroInflationFactor = 0.0,
                      log_prob_t preliminarySignalShape = 1.5,
                      log_prob_t preliminarySignalScale = 2 )
    : sequenceLengthFactor( sequenceLengthFactor )
    , scShape( scShape )
    , scRateStep( scRateStep )
    , zeroInflationFactor( zeroInflationFactor )
    , preliminarySignalPrior( preliminarySignalShape, preliminarySignalScale )
    {
    }

    distribution_type operator()( signal_t lnRate ) const {
        return ( distribution_type( lnRate, scShape, zeroInflationFactor ) );
    }

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( sequenceLengthFactor );
        ar & BOOST_SERIALIZATION_NVP( scShape );
        ar & BOOST_SERIALIZATION_NVP( scRateStep );
        ar & BOOST_SERIALIZATION_NVP( zeroInflationFactor );
        ar & BOOST_SERIALIZATION_NVP( preliminarySignalPrior );
    }
};

/**
    Parameters of the distribution of object-probe measurements,
    according to the model.
    @todo rename to CellSignal
 */
struct ObjectsClusterSignal {
    typedef CellSignalParams  signal_params_type;
    typedef CellSignalParams::sc_distrib_type sc_distrib_type;

    signal_t    _lnScRate;
    signal_t    _scShape;

    ObjectsClusterSignal( signal_t lnScRate = unset<signal_t>(), signal_t scShape = 0 )
    : _lnScRate( lnScRate ), _scShape( scShape )
    {}

    ObjectsClusterSignal( const signal_params_type& signalParams,
                          const OPAObject& object,
                          size_t objMultiple = 1 ) 
    : _lnScRate( signalParams.sequenceLengthFactor * log_int( object.sequenceLength() )
                 + log_int( objMultiple ) )
    , _scShape( signalParams.scShape )
    {
    }

    ObjectsClusterSignal( const ObjectsClusterSignal& baseSignal,
                          const OPAAssay& assay,
                          signal_t lnSignal )
    : _lnScRate( baseSignal._lnScRate + lnSignal + assay.lnMultiplier() )
    , _scShape( baseSignal._scShape )
    {
    }

    ObjectsClusterSignal operator+( signal_t lnSignal ) const {
        return ( ObjectsClusterSignal( _lnScRate + lnSignal ) );
    }

    ObjectsClusterSignal& operator+=(signal_t lnSignal) {
        _lnScRate += lnSignal;
        return ( *this );
    }

    ObjectsClusterSignal& operator=(signal_t lnScRate ) {
        _lnScRate = lnScRate;
        return ( *this );
    }

    double lnScRate() const {
        return ( _lnScRate );
    }
    double scShape() const {
        return ( _scShape );
    }

    bool operator==(const ObjectsClusterSignal& that) const {
        return ( _lnScRate == that._lnScRate && _scShape == that._scShape );
    }

    bool operator!=(const ObjectsClusterSignal& that) const {
        return ( !operator==(that) );
    }

    template<class Stream>
    friend Stream& operator<<( Stream& out, const ObjectsClusterSignal& a )
    {
        return ( out << "lnScRate=" << a._lnScRate << ", scShape=" << a._scShape );
    }

    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & boost::serialization::make_nvp( "lnScRate", _lnScRate );
        ar & boost::serialization::make_nvp( "scShape", _scShape );
    }

    typedef DistributionCache<CellSignalParams, signal_t> distrib_cache_type;
    typedef distrib_cache_type::const_distrib_ref_type distrib_ref_type;

    distrib_ref_type distribTable( distrib_cache_type& cache ) const
    {
        /// @TODO: check scShape
        return ( cache.distrib( _lnScRate ) );
    }
};

} }

BOOST_CLASS_IMPLEMENTATION( cemm::bimap::CellSignalParams, object_serializable )
BOOST_IS_BITWISE_SERIALIZABLE( cemm::bimap::CellSignalParams )

BOOST_CLASS_IMPLEMENTATION( cemm::bimap::ObjectsClusterSignal, object_serializable )
BOOST_IS_BITWISE_SERIALIZABLE( cemm::bimap::ObjectsClusterSignal )
