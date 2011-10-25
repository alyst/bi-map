#pragma once

#include "../BasicTypedefs.h"

#include <sstream>
#include <vector>

#include "Particle.h"
#include "ITurbineConnection.h"

#define TURBINE_LOG_REPORT_PERIOD 1000

/**
    ParticleTurbineCascade execution monitor.
 */
class TurbineCascadeExecutionMonitor {
public:
    typedef int log_level_t;

    /**
        Proxy class to buffer message string
        and ouput it to execution monitor upon
        destruction.
     */
    class LogMessage {
    protected:
        friend class TurbineCascadeExecutionMonitor;

        std::ostringstream                  _msg;
        TurbineCascadeExecutionMonitor&     _parent;
        log_level_t                         _level;
        bool                                _enabled;

        LogMessage( TurbineCascadeExecutionMonitor& parent, 
                    log_level_t level, bool enabled )
        : _parent( parent ), _level( level ), _enabled( enabled )
        {}

    public:
        LogMessage( const LogMessage& logMsg )
        : _parent( logMsg._parent ), _level( logMsg._level ), _enabled( logMsg._enabled )
        {
            _msg << logMsg._msg.str();
        }

        ~LogMessage()
        {
            if ( _enabled && _msg.tellp() > 0 ) {
                _parent.log( _level, _msg.str() );
            }
        }

        bool enabled() const {
            return ( _enabled );
        }

        operator std::ostringstream&() {
            return ( _msg );
        }
        
        template<class Printable>
        friend LogMessage& operator<<( LogMessage& out, const Printable& a )
        {
            if ( out._enabled ) out._msg << a;
            return ( out );
        }
    };

private:
    log_level_t _verbosity;

public:
    TurbineCascadeExecutionMonitor( log_level_t verbosity = 0 )
    : _verbosity( verbosity )
    {}

    virtual void log( log_level_t level, const std::string& message ) = 0;

    LogMessage logOutput( log_level_t level ) 
    {
        return ( LogMessage( *this, verbosity(), level <= verbosity() ) );
    }

    virtual void notifyCascadeCreated() = 0;
    virtual void notifyCascadeInitialized() = 0;
    virtual void notifyIteration( size_t iteration ) = 0;
    virtual void notifyEquiEnergyJump( size_t turbineIx, bool released, double prevEnergy, double nextEnergy, bool accepted ) = 0;
    virtual void notifyTurbineCreated( bool calibration, size_t turbineIx, const EnergyTransform& params ) = 0;
    virtual void notifyRingsLadderAdjusted( size_t turbineIx, const std::vector<double>& energyThresholds, const std::vector<size_t>& particleCounts ) = 0;
    virtual void notifyParticlesGenerated( size_t turbineIx, size_t particlesCount ) = 0;
    virtual void notifyParticleReleased( size_t turbineIx, double energy ) = 0;
    virtual void notifyParticlePushed( size_t turbineIx, double energy ) = 0;
    virtual void notifyParticlePopped( size_t turbineIx, double energy, size_t pushedIteration, bool released ) = 0;
    virtual void notifyFinished() = 0;

    log_level_t verbosity() const
    {
        return ( _verbosity );
    }
};

#if 0
#define EELOG_DEBUG( pMonitor, level, message ) \
if ( (pMonitor) != NULL && (pMonitor)->verbosity() >= (level) ) { \
    TurbineCascadeExecutionMonitor::LogMessage monMsgProxy = (pMonitor)->logOutput( level ); \
    if ( monMsgProxy.enabled() ) monMsgProxy << message; \
}
#else
#define EELOG_DEBUG( pMonitor, level, message )
#endif

#define EELOG_DEBUG0( message ) EELOG_DEBUG( _pMonitor, 0, message )
#define EELOG_WARN( message )   EELOG_DEBUG( _pMonitor, 0, message )
#define EELOG_DEBUG1( message ) EELOG_DEBUG( _pMonitor, 1, message )
#define EELOG_DEBUG2( message ) EELOG_DEBUG( _pMonitor, 2, message )
