#pragma once

#if defined(NDEBUG)
#define DEBUG_LEVEL -1
#else
#define DEBUG_LEVEL 1
#endif

#include <cemm/debug.h>

#include <boost/unordered_map.hpp>

#include <cemm/math/common.h>

namespace cemm { namespace bimap {

using namespace cemm::math;

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

typedef value_t signal_t;

} }