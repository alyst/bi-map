#include "Turbine.h"

TurbineParams::TurbineParams()
    : eeJumpRate( 0.1 )
    , eeJumpEnergyTolerance( 0.03 )
    , eeMaxEnergyQuantile( 0.9 )
    , eeLadderEnergyQuantile( 0.5 )
    , maxParticles( 500 )
    , particleSnapshotPeriod( 13 )
    , particleSamplingPeriod( 67 )
    , generateRate( 0.01 )
    , prisonerOnWalkSwitchRate( 0.2 )
    , detentionIterations( 19 )
{
}

bool is_equienergy_jump_accepted(
    const gsl_rng*  rng,
    energy_type     curEnergy,
    const EnergyTransform& curTransform,
    energy_type     jumpEnergy,
    const EnergyTransform& jumpTransform
){
    const energy_type jumpjump = -std::max( jumpEnergy, jumpTransform.minEnergyThreshold ) / jumpTransform.temperature;
    const energy_type jumpcur = -std::max( jumpEnergy, curTransform.minEnergyThreshold ) / curTransform.temperature;
    const energy_type curcur = -std::max( curEnergy, curTransform.minEnergyThreshold) / curTransform.temperature;
    const energy_type curjump = -std::max( curEnergy, jumpTransform.minEnergyThreshold ) / jumpTransform.temperature;
    const energy_type logRatio = curjump + jumpcur - jumpjump - curcur;

    return ( logRatio >= 0 || gsl_rng_uniform( rng ) <= exp( logRatio ) );
}

