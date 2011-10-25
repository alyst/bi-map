#pragma once

#include <iostream>
#include <set>

#include <boost/unordered_map.hpp>

#include <boost/math/special_functions/fpclassify.hpp>

#include <limits>

#if defined(NDEBUG)
#define DEBUG_LEVEL -1
#else
#define DEBUG_LEVEL 1
#endif

#include "../include/debug.h"

/// index of system's probe
typedef std::size_t probe_index_t;

#define PROBE_NA ( (probe_index_t)-1 )

/// index of system's probe assay
typedef std::size_t assay_index_t;

#define ASSAY_NA ( (assay_index_t)-1 )

/// index of system's object
typedef std::size_t object_index_t;

#define OBJECT_NA ( (object_index_t)-1 )

typedef std::string assay_label_t;
typedef std::string probe_label_t;
typedef std::string object_label_t;

typedef boost::unordered_map<object_index_t, object_label_t> objects_label_map_type;
typedef boost::unordered_map<probe_index_t, probe_label_t> probes_label_map_type;

/// floating-point type for probability
typedef double prob_t;
/// floating-point type for probability log
typedef double log_prob_t;
/// floating-point type for signal (cluster abundance)
typedef double signal_t;

inline log_prob_t unset()
{
    return ( std::numeric_limits<log_prob_t>::quiet_NaN() );
}

template<typename T>
inline bool is_unset( T val )
{
    return ( boost::math::isnan( val ) );
}

template<typename T>
inline bool is_finite( T val )
{
    return ( boost::math::isfinite( val ) );
}

template<typename T>
inline bool is_infinite( T val )
{
    return ( boost::math::isinf( val ) );
}
