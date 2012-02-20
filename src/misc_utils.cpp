#include "cemm/bimap/misc_utils.h"

/**
 *  Test if string ends with given substring.
 */
bool endsWith(
    const char*     str,
    const char*     ending
){
    std::string oStr( str );
    std::string oEnd( ending );

    return ( oStr.length() > oEnd.length()
             && 0 == oStr.compare( oStr.length() - oEnd.length(), oEnd.length(), oEnd ) );
}
