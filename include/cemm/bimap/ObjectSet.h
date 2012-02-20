#pragma once

#include <cemm/containers/set_misc.h>

namespace cemm { namespace bimap {

typedef std::set<object_index_t>     object_set_t;

typedef cemm::containers::SetDistance<object_index_t> ObjectSetDistance;

} }