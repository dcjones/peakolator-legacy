
#include "miscmath.hpp"
#include <cmath>
#include <algorithm>

using namespace std;

#ifdef HAVE_LIBGSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#endif



double logaddexp( double x, double y )
{
    double u = x - y;
    if( u > 0.0 ) {
        return x + log1p( exp( -u ) );
    }
    else if( u <= 0.0 ) {
        return y + log1p( exp( u ) );
    }
    else {
        return x + y;
    }
}



double logsubexp( double x, double y )
{
    double u = x - y;
    if( u > 0.0 ) {
        return x + log1p( -exp( -u ) );
    }
    else if( u <= 0.0 ) {
        return y + log1p( -exp( u ) );
    }
    else {
        return x - y;
    }
}


#ifdef HAVE_LIBGSL
double log_beta_inc( double a, double b, double x )
{
    /* Using Soper's Method, which possibly converges more slowly than the
     * continued fraction, but is easier to do in log space. */

    if( x > (a + 1) / (a + b +2) || 1 - x < (b + 1) / (a + b + 2) ) {
        return logsubexp( 0, log_beta_inc( b, a, 1 - x ) );
    }

    size_t maxiter = 200;
    double result0;
    double result = NAN;
    double r;

    while( maxiter-- ) {
        result0 = result;

        r  = a * log(x) + (b - 1) * log( 1 - x );
        r -= log(a) + gsl_sf_lnbeta( a, b );

        if( result == NAN ) result = r;
        else                result = logaddexp( result, r );

        a += 1;
        b = max( 0.0, b - 1 );

        if( gsl_fcmp( result, result0, 1e-16 ) == 0 ) break;
    }

    return result;
}


double lpnbinom( unsigned int q, double r, double p )
{
    return logsubexp( 0, log_beta_inc( (double)(q + 1), r, p ) );
}
#endif



double pgev( double q, double loc, double scale, double shape )
{
    double p;
    q = (q - loc) / scale;
    if( shape == 0.0 ) p = exp( -exp( -q ) );
    else               p = exp( -pow( max( 1 + shape * q, 0.0 ), -1/shape ) );

    return 1.0 - p;
}


double ldgev( double x, double loc, double scale, double shape )
{
    double hx;
    if( scale < 0.0 ) return -HUGE_VAL;
    if( shape == 0.0 ) {
        hx = (loc - x) / scale;
        return -log( scale ) + hx - exp( hx );
    }
    else {
        hx = 1.0 + shape * (x - loc) / scale;
        if( hx <= 0 ) return -HUGE_VAL;
        return -log( scale ) - (1 + 1/shape) * log(hx) - pow( hx, -1/shape );
    }
}



void colcpy( double* dest, const double* src, size_t j, size_t n, size_t m )
{
    size_t i;
    for( i = 0; i < n; i++ ) dest[i] = src[ i * m + j ];
}

void vecadd( double* u, double* v, size_t n )
{
    size_t i;
    for( i = 0; i < n; i++ ) {
        u[i] += v[i];
    }
}

void vecsub( double* u, double* v, size_t n )
{
    size_t i;
    for( i = 0; i < n; i++ ) {
        u[i] -= v[i];
    }
}

void vecaddcol( double* u, double* V, size_t n, size_t m, size_t j )
{
    size_t i;
    for( i = 0; i < n; i++ ) {
        u[i] += V[ i * m + j ];
    }
}

void vecsubcol( double* u, double* V, size_t n, size_t m, size_t j )
{
    size_t i;
    for( i = 0; i < n; i++ ) {
        u[i] -= V[ i * m + j ];
    }
}

void matsetcol( double* U, double* v, size_t n, size_t m, size_t j )
{
    size_t i;
    for( i = 0 ; i < n; i++ )  {
        U[ i * m + j ] = v[i];
    }
}


