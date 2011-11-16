#pragma once

#include "BasicTypedefs.h"

#include <boost/timer.hpp>
#include <fstream>
#include <sstream>

#include "eesampler/ExecutionMonitor.h"

/**
    TurbineCascadeExecutionMonitor writing to textual console.
 */
class ConsolePTCExecutionMonitor: public TurbineCascadeExecutionMonitor
{
private:
    struct TurbineStats {
        size_t eeJumps; /// EE jumps
        size_t eeSuccessfulJumps; /// succ. EE jumps
        size_t eeReleasedJumps; // jumps to released particles
        size_t eeReleasedSuccJumps; // succ. jumps to released particles
        size_t particlesGenerated; /// generated (from this turbine) particles
        size_t particlesReleased; /// particles released to the cache

        TurbineStats()
        : eeJumps( 0 ), eeSuccessfulJumps( 0 )
        , eeReleasedJumps( 0 ), eeReleasedSuccJumps( 0 )
        , particlesGenerated( 0 ), particlesReleased( 0 )
        {}
    };
    boost::timer    _timer;
    size_t _iteration; /// current iteration of EE sampler
    std::vector<TurbineStats> _turbineStats; /// per turbine statistics
    std::ofstream   _particlesLogStream;
    std::ofstream   _eeJumpLogStream;

protected:
    virtual void log_printf( const char* format, ... ) = 0;
    void resizeCounters( size_t turbineIx );

public:

    ConsolePTCExecutionMonitor( log_level_t verbosity = 0, const char* particleLog = NULL, const char* eeJumpLog = NULL );

    ~ConsolePTCExecutionMonitor();

    virtual void notifyCascadeCreated();
    virtual void notifyCascadeInitialized();
    virtual void notifyIteration( size_t iteration );

    virtual void notifyEquiEnergyJump( size_t turbineIx, bool released, double prevEnergy, double nextEnergy, bool accepted );
    
    virtual void log( log_level_t level, const std::string& message )
    {
        if ( level <= verbosity() ) {
            log_printf( "EEItn=#%d: %s\n", _iteration, message.c_str() );
        }
    }

    virtual void notifyTurbineCreated( bool calibrating, size_t turbineIx, const EnergyTransform& params );

    virtual void notifyRingsLadderAdjusted( size_t turbineIx, const std::vector<double>& energyThresholds, const std::vector<size_t>& particleCounts );
    
    virtual void notifyFinished();
    
    virtual void notifyParticlesGenerated( size_t turbineIx, size_t particlesCount );
    virtual void notifyParticleReleased( size_t turbineIx, double energy );

    virtual void notifyParticlePushed( size_t turbineIx, double energy );
    virtual void notifyParticlePopped( size_t turbineIx, double energy, size_t pushedIteration, bool released );
};

/**
    TurbineCascadeExecutionMonitor printing to STDOUT of the process.
 */
class StdOutPTCExecutionMonitor: public ConsolePTCExecutionMonitor {
protected:
    virtual void log_printf( const char* format, ... );

public:
    StdOutPTCExecutionMonitor( log_level_t verbosity = 0 )
    : ConsolePTCExecutionMonitor( verbosity )
    {}
    virtual ~StdOutPTCExecutionMonitor()
    {}
};
