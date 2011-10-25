#include "TurbineCascadeUnit.h"

tcascade_structure_type createTurbineCascadeStructure(
    const std::vector<size_t>& turbinesPerLevel, size_t units
){
    size_t turbinesCount = 0;
    for ( size_t level = 0; level < turbinesPerLevel.size(); level++ ) {
        turbinesCount += turbinesPerLevel[ level ];
    }

    tcascadeunit_ix_t unitIx = 0;
    turbine_ix_t curIx = 0;
    tcascade_structure_type res( turbinesCount );
    for ( size_t level = 0; level < turbinesPerLevel.size(); level++ ) {
        size_t turbinesOnLevel = turbinesPerLevel[ level ];
        for ( size_t i = 0; i < turbinesOnLevel; i++ ) {
            TurbineNode  turbineProps;
            turbineProps.level = level;
            turbineProps.turbineIx = curIx++;
            turbineProps.unitIx = unitIx;

            if ( (turbineProps.turbineIx+1) * units >= turbinesCount * (unitIx+1) ) {
                unitIx++;
            }
            res[ turbineProps.turbineIx ] = turbineProps;
            LOG_INFO( "Cascade: putting turbine #" << turbineProps.turbineIx
                      << " of level " << turbineProps.level
                      << " to unit #" << turbineProps.unitIx );
        }
    };
    if ( unitIx < units ) {
        LOG_WARN( "Only " << unitIx << " units of " << units << " are used for the cascade" );
    }
    if ( curIx != turbinesCount ) {
        THROW_RUNTIME_ERROR( "Count of turbines in the cascade (" << curIx 
            << ") does not match the requested count (" << turbinesCount << ")" );
    }
    BOOST_ASSERT( unitIx <= units );
    return ( res );
}

tcascade_structure_type createTurbineCascadeStructure(
    size_t levels, size_t turbines, size_t units
){
    if ( levels > turbines ) THROW_EXCEPTION( std::invalid_argument, "More levels (" << levels << ") than turbines (" << turbines << ")" );

    std::vector<size_t> turbinesPerLevel( levels );
    for ( size_t level = 0; level < turbinesPerLevel.size(); level++ ) {
        turbinesPerLevel[ level ] = level == 0 ? ( turbines + 1 - levels ) : 1;
    }
    return ( createTurbineCascadeStructure( turbinesPerLevel, units ) );
}

TurbineCascadeParams::TurbineCascadeParams()
    : turbinesCount( 2 )
    , levelsCount( 2 )
    , temperatureMultiplier( 1.5 )
    , burnInIterations( 1000 )
    , ladderAdjustPeriod( 100 )
    , broadcastStatusPeriod( 200 )
{
}

std::vector<energy_type> eval_even_segments(
    const std::vector<energy_type>&  sortedValues,
    size_t                      segmentsCount,
    size_t                      maxIx = std::numeric_limits<size_t>::max(),
    bool                        includeEnd = false,
    energy_type                 eps = std::numeric_limits<energy_type>::epsilon()
){
    typedef std::vector<energy_type> val_vector_t;

    if ( maxIx == std::numeric_limits<size_t>::max() ) maxIx = sortedValues.size();

//    std::sort( values.begin(), values.end() );
    val_vector_t    res( segmentsCount );
    if ( sortedValues.empty() ) {
        throw std::runtime_error( "No values provided" );
    }
    #if 0
    if ( sortedValues.size() < segmentsCount + ( includeEnd ? 1 : 0 ) ) {
        throw std::runtime_error( "Not enough values provided" );
    }
    #endif

    bool ascending = sortedValues.front() < sortedValues[ maxIx ];

#if 1
    res[ 0 ] = sortedValues[ 0 ];
    for ( size_t i = 1; i < segmentsCount; i++ ) {
        const energy_type prevValue = res[ i - 1 ];
        const double pos = i * maxIx / segmentsCount;
        size_t ix1 = (int)pos;
        size_t ix2 = std::min( maxIx - 1, (size_t)(pos + 0.999) );
        res[ i ] = ix1 == ix2 ? sortedValues[ ix1 ] 
                 : (pos - ix1) * sortedValues[ ix1 ] + (ix2 - pos) * sortedValues[ ix2 ];
        if ( ascending && res[ i ] <= prevValue ) res[ i ] = prevValue + eps;
        else if ( !ascending && res[ i ] >= prevValue ) res[ i ] = prevValue - eps;
    }
    if ( includeEnd ) {
        energy_type endVal = sortedValues[ maxIx ];
        const energy_type prevValue = res.back();
        if ( ascending && endVal <= prevValue ) endVal = prevValue + eps;
        else if ( !ascending && endVal >= prevValue ) endVal = prevValue - eps;
        res.push_back( endVal );
    }
#else
    bool start = sortedValues.front();
    bool delta = ( sortedValues[ maxIx ] - sortedValues.front() ) / segmentsCount;

    for ( size_t i = 0; i < segmentsCount; i++ ) {
        res.push_back( start + i * delta );
    }

    if ( includeEnd ) {
        energy_type endVal = sortedValues.back();
        const energy_type prevValue = res.back();
        if ( ascending && endVal <= prevValue ) endVal = prevValue + eps;
        else if ( !ascending && endVal >= prevValue ) endVal = prevValue - eps;
        res.push_back( endVal );
    }
#endif

    return ( res );
}

TurbineEnergyLadder::TurbineEnergyLadder(
    const TurbineCascadeParams& params,
    const TurbineEnergyLadder::energy_landscape_type& energyLandscape
) : steps( params.levelsCount )
{
    // adjust energy bounds of existing turbines
    for ( size_t stepIx = 0; stepIx < steps.size(); stepIx++ ) {
        EnergyTransform& step = steps[ stepIx ];
        if ( stepIx > 0 ) {
            energy_landscape_type::const_iterator lit = energyLandscape.find( stepIx - 1 );
            if ( lit == energyLandscape.end() || lit->second.empty() ) {
                step = EnergyTransform( std::numeric_limits<energy_type>::quiet_NaN(), 1.0 );
                EELOG_WARN( "Step #" << stepIx << ": energy landscape empty, skipping" );
                continue;
            } else {
                size_t quantileParticleIx = lit->second.size() * params.turbineParams.eeLadderEnergyQuantile;
                step = EnergyTransform( lit->second[ quantileParticleIx ], 
                                        gsl_pow_int( params.temperatureMultiplier, stepIx ) );
            }
        } else {
            step = EnergyTransform( -std::numeric_limits<energy_type>::infinity(), 1.0 );
        }
        EELOG_DEBUG0( "Step #" << stepIx << " params: E_min=" << step.minEnergyThreshold
                       << " T=" << step.temperature );
    }
}
