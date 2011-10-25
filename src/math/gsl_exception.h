#pragma once

#include <gsl/gsl_errno.h>
#include <exception>

/**
 *  Wrapper for GSL errors reported by GSL error handler.
 */
class gsl_exception: public std::exception {
public:
    gsl_exception( const char * reason,
                   const char * file,
                   int line,
                   int gsl_errno )
    : std::exception()
    {
    }
};

/**
 *  Handles GSL error by throwing gsl_exception.
 */
void raise_gsl_exception( const char * reason,
                          const char * file,
                          int line,
                          int gsl_errno
){
    throw gsl_exception( reason, file, line, gsl_errno );
}
