#include <stdarg.h>

#include "ConsolePTCExecutionMonitor.h"

ConsolePTCExecutionMonitor::ConsolePTCExecutionMonitor(
    TurbineCascadeExecutionMonitor::log_level_t verbosity,
    const char* particleLog,
    const char* eeJumpLog
) : TurbineCascadeExecutionMonitor( verbosity )
{
    if ( particleLog ) {
        _particlesLogStream.open( particleLog );
        _particlesLogStream << "disk\tenergy\tstart\tend\treleased\n";
    }
    if ( eeJumpLog ) {
        _eeJumpLogStream.open( eeJumpLog );
        _eeJumpLogStream << "iteration\tdisk\tenergy\tnew_energy\treleased\taccepted\n";
    }
}

ConsolePTCExecutionMonitor::~ConsolePTCExecutionMonitor() {
    if ( _particlesLogStream.is_open() ) {
        _particlesLogStream.close();
    }
    if ( _eeJumpLogStream.is_open() ) {
        _eeJumpLogStream.close();
    }
}

void ConsolePTCExecutionMonitor::notifyCascadeCreated()
{
    log_printf( "EE turbine cascade created\n" );
    _iteration = -1;
}

void ConsolePTCExecutionMonitor::notifyCascadeInitialized()
{
    log_printf( "EE turbine cascade initialized\n" );
    _iteration = -1;
    _timer.restart();
}

void ConsolePTCExecutionMonitor::notifyIteration(size_t iteration) {
    _iteration = iteration;
    if ( _iteration % 1000 == 0 ) {
        log_printf( "EESamper it=%d (%.0f ms)...\n", _iteration, _timer.elapsed() * 1000 );
    }
}

void ConsolePTCExecutionMonitor::resizeCounters(size_t turbineIx)
{
    if ( _turbineStats.size() < turbineIx + 1 ) {
        _turbineStats.resize( turbineIx + 1 );
    }
}

void ConsolePTCExecutionMonitor::notifyEquiEnergyJump(
    size_t turbineIx, bool released, double prevEnergy, double nextEnergy, bool accepted
) {
    resizeCounters( turbineIx );
    TurbineStats& stats = _turbineStats[ turbineIx ];
    stats.eeJumps++;
    if ( accepted ) stats.eeSuccessfulJumps++;
    if ( released ) stats.eeReleasedJumps++;
    if ( accepted & released ) stats.eeReleasedSuccJumps++;

    if ( verbosity() >= 2 ) {
        LogMessage msgOutput = logOutput( 2 );
        msgOutput << "EE jump at t#" << turbineIx
                  << ": E(a)=" << prevEnergy << "<->E(d)=" << nextEnergy 
                  << ( released ? " (released)" : "" );
    }
    if ( _eeJumpLogStream.is_open() ) {
        _eeJumpLogStream << _iteration << '\t' << (int)turbineIx << '\t' 
                            << prevEnergy << '\t' << nextEnergy << '\t' << '\t' << released << '\t' << accepted << '\n';
    }
}

void ConsolePTCExecutionMonitor::notifyTurbineCreated(bool calibrating, size_t turbineIx, const EnergyTransform& params) {
    LogMessage msgOutput = logOutput( 0 );
    if ( calibrating ) {
        msgOutput << "started calibration turbine";
    } else {
        msgOutput << "started turbine #" << turbineIx;
    }
    msgOutput << ": E_threshold=" << params.minEnergyThreshold << " T=" << params.temperature;
}

void ConsolePTCExecutionMonitor::notifyRingsLadderAdjusted(size_t turbineIx, const std::vector< double >& energyThresholds, const std::vector< size_t >& particleCounts) {
    if ( verbosity() < 1 ) return;

    LogMessage msgOutput = logOutput( 1 );

    msgOutput << "energy rings of turbine #" << turbineIx << " adjusted:";
    for ( size_t i = 0; i < energyThresholds.size(); i++ ) {
        msgOutput << '\n' << "Ring #" << i << ": E_min=" << energyThresholds[ i ] << " size=" << particleCounts[ i ];
    }
}

void ConsolePTCExecutionMonitor::notifyFinished() {
    log_printf( "EE sampler finished, %d iterations, %.0f ms\n", _iteration + 1, _timer.elapsed() * 1000 );
    for ( size_t turbineIx = 0; turbineIx < _turbineStats.size(); turbineIx++ ) {
        const TurbineStats& stats = _turbineStats[ turbineIx ];
        log_printf( "Turbine %d: %d(%d)/%d(%d) EE jumps, %d/%d particles released/generated\n", 
                    turbineIx, 
                    stats.eeSuccessfulJumps, stats.eeReleasedSuccJumps, 
                    stats.eeJumps, stats.eeReleasedJumps,
                    stats.particlesReleased, stats.particlesGenerated );
    }
}

void ConsolePTCExecutionMonitor::notifyParticleReleased(size_t turbineIx, double energy)
{
    resizeCounters( turbineIx );
    _turbineStats[ turbineIx ].particlesReleased++;
}

void ConsolePTCExecutionMonitor::notifyParticlesGenerated(size_t turbineIx, size_t particlesCount)
{
    resizeCounters( turbineIx );
    _turbineStats[ turbineIx ].particlesGenerated += particlesCount;
}

void ConsolePTCExecutionMonitor::notifyParticlePushed(size_t turbineIx, double energy)
{

}

void ConsolePTCExecutionMonitor::notifyParticlePopped(size_t turbineIx, double energy, size_t pushedIteration, bool released)
{
    if ( _particlesLogStream.is_open() ) {
        _particlesLogStream << (int)turbineIx << '\t' 
                            << energy << '\t' << pushedIteration << '\t' << _iteration << '\t' << released << '\n';
    }
}

void StdOutPTCExecutionMonitor::log_printf( const char* format, ... )
{
    va_list ap;
    va_start( ap, format );
    vprintf( format, ap );
    va_end( ap );
}
