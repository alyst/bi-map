#pragma once

#include <stdexcept>
#include <sstream>

#if defined(NDEBUG)
#define BOOST_DISABLE_ASSERTS
#else
#define _DEBUG
#endif

#define BOOST_ENABLE_ASSERT_HANDLER

#include <boost/assert.hpp>

namespace boost
{

inline void assertion_failed(char const * expr, char const * function, char const * file, long line)
{
    std::ostringstream out;
    out << file << "(" << line << "): assertion failed: " << expr;
    throw std::runtime_error( out.str() );
}

} 

#define LOG_STDERR( msg ) ( std::cerr << msg << '\n' )

#define LOG_INFO( msg ) LOG_STDERR( msg )
#define LOG_INFO_IF( condition, msg ) if ( condition ) { LOG_STDERR( msg ); }

#define LOG_WARN( msg ) LOG_STDERR( msg )
#define LOG_WARN_IF( condition, msg ) if ( condition ) { LOG_STDERR( msg ); }

#if DEBUG_LEVEL >= 0
#define LOG_DEBUG0( msg ) LOG_STDERR( msg )
#define LOG_DEBUG0_IF( condition, msg ) if ( condition ) { LOG_STDERR( msg ); }
#else
#define LOG_DEBUG0( msg )
#define LOG_DEBUG0_IF( condition, msg )
#endif

#if DEBUG_LEVEL >= 1
#define LOG_DEBUG1( msg ) LOG_STDERR( msg )
#define LOG_DEBUG1_IF( condition, msg ) if ( condition ) { LOG_STDERR( msg ); }
#else
#define LOG_DEBUG1( msg )
#define LOG_DEBUG1_IF( condition, msg )
#endif

#if DEBUG_LEVEL >= 2
#define LOG_DEBUG2( msg ) LOG_STDERR( msg )
#define LOG_DEBUG2_IF( condition, msg ) if ( condition ) { LOG_STDERR( msg ); }
#else
#define LOG_DEBUG2( msg )
#define LOG_DEBUG2_IF( condition, msg )
#endif

#if DEBUG_LEVEL >= 3
#define LOG_DEBUG3( msg ) LOG_STDERR( msg )
#define LOG_DEBUG3_IF( condition, msg ) if ( condition ) { LOG_STDERR( msg ); }
#else
#define LOG_DEBUG3( msg )
#define LOG_DEBUG3_IF( condition, msg )
#endif

#define THROW_EXCEPTION( excp_class, msg ) { std::ostringstream __excp__msg__; __excp__msg__ << msg; throw excp_class( __excp__msg__.str() ); }
#define THROW_RUNTIME_ERROR( msg ) THROW_EXCEPTION( std::runtime_error, msg )
