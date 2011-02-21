
#include "miscmath.hpp"
#include <cmath>
#include <cstdio>
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
    double u = y - x;
    if( gsl_finite(u) ) return x + log1p( -exp( u ) );
    else                return x - y;
}


#ifdef HAVE_LIBGSL

double binco( double n, double k )
{
   return gsl_sf_gamma( n + 1.0 ) /
          ( gsl_sf_gamma( n - k + 1.0 ) * gsl_sf_gamma( k + 1.0 ) );
}


double lbinco( double n, double k )
{
   return gsl_sf_lngamma( n + 1.0 )
          - gsl_sf_lngamma( n - k + 1.0 )
          - gsl_sf_lngamma( k + 1.0 );
}


double log_beta_inc( double a, double b, double x )
{
    /* Using Soper's Method, which possibly converges more slowly than the
     * continued fraction, but is easier to do in log space. */

    if( a < (a+b)*x ) {
        return logsubexp( 0, log_beta_inc( b, a, 1.0 - x ) );
    }

    size_t maxiter = 200;
    double result = GSL_NAN;
    double r;

    int i;
    int s = (int)(b + (1.0 - x) * (a + b));

    for( i = 0; i < maxiter; i++ ) {

        if( i < s && b > 1.0 ) {
            r  = a * log(x) + (b-1.0) * log(1.0-x);
            r -= log(a) + gsl_sf_lnbeta(a,b);
            a += 1;
            b -= 1;
        }
        else {
            r  = a * log(x) + (b) * log(1.0-x);
            r -= log(a) + gsl_sf_lnbeta(a,b);
            a += 1;
        }

        if( gsl_isnan( result ) ) result = r;
        else                      result = logaddexp( result, r );

        if( !gsl_finite(result) || !gsl_finite(r) || abs(r) < 1e-12 ) break;
    }

    return result;
}


double ldbinom( unsigned int x, double p, unsigned int n )
{
    return lbinco( n, x ) + (double)x * log(p) + (double)(n-x) * log(1.0 - p);
}


double ldnbinom( unsigned int x_, double r, double p )
{
    double x = (double)x_;
    double ans;
    ans = r * log(p) + x * log( 1.0 - p );
    ans += gsl_sf_lngamma( r + x );
    ans -= (gsl_sf_lnfact( x_ ) + gsl_sf_lngamma( r ));

    return ans;
}


double lpnbinom( unsigned int q, double r, double p, bool lower_tail )
{
    if( lower_tail ) return log_beta_inc( r, (double)q + 1.0, p );
    else             return log_beta_inc( (double)q + 1.0, r, 1.0 - p );
}



double lddnbinom( unsigned int x, double r, double p )
{
    if( x == 0 ) return GSL_NEGINF;
    return lbinco( (double)x + r - 1, x )
         + x * log(p)
         - logsubexp( -r * log(1.0 - p), 0.0 );
}



double ddnbsum( unsigned int x_, double r, double p, unsigned int d_ )
{
    if( x_ < d_ ) return 0.0;

    int sgn = -1;
    double a, ans = 0.0;

    double x = (double)x_;
    double d = (double)d_;

    double u1 = pow( p, x );
    double u2 = pow( pow( 1.0 - p, -r ) - 1.0, d );

    unsigned int k;
    for( k = d; k > 0; k-- ) {
        sgn *= -1;

        a =  binco( d, k );
        a *= binco( (double)x + k * r - 1, x );
        a *= u1;
        a /= u2;

        if( sgn < 0 ) ans -= a;
        else          ans += a;
    }

    return ans;
}



double pdnbsum( unsigned int x_, double r, double p,
                unsigned int d_, bool lower_tail )
{
    if( x_ < d_ ) return 0.0;

    double x = (double)x_;
    double d = (double)d_;

    int sgn = -1;
    double a, ans = 0.0;


    double k;
    for( k = d; k > 0; k-- ) {
        sgn *= -1;
        a = binco( d, k );
        a *= pow( 1.0 - p, - (double)k * r );
        a *= gsl_sf_beta_inc( x + 1.0, k * r, p );

        if( sgn < 0 ) ans -= a;
        else          ans += a;
    }

    ans *= pow( pow( 1.0 - p, -r ), -d );

    if( lower_tail ) ans = 1.0 - ans;

    return ans;
}


double lddnbsum( unsigned int x_, double r, double p, unsigned int d_ )
{
    if( x_ < d_ ) return GSL_NEGINF;

    double x = (double)x_;
    double d = (double)d_;

    double lp = log( p );
    double lq = log( 1.0 - p );

    double a, ans = GSL_NAN;
    int sgn = -1;

    double k;
    for( k = d; k > 0; k-- ) {
        sgn *= -1;

        a =  lbinco( d, k );
        a += lbinco( x + k * r - 1.0, x );
        a += x * lp;
        a -= d * logsubexp( -r * lq, 0.0 );

        if( k == d )         ans = a;
        else if( sgn < 0.0 ) ans = logsubexp( ans, a );
        else                 ans = logaddexp( ans, a );
    }

    return ans;
}


double lpdnbsum( unsigned int x_, double r, double p,
                 unsigned int d_, bool lower_tail )
{
    if( x_ < d_ ) return GSL_NEGINF;

    double x = (double)x_;
    double d = (double)d_;

    double lp = log( p );
    double lq = log( 1.0 - p );

    double a, ans = GSL_NAN;
    int sgn = -1;

    double k;
    for( k = d; k > 0; k-- ) {
        sgn *= -1;

        a = lbinco( d, k );
        a += -k * r * lq;
        a += log_beta_inc( x + 1.0, k * r, p );

        if( k == d )       ans = a;
        else if( sgn < 0 ) ans = logsubexp( ans, a );
        else               ans = logaddexp( ans, a );
    }

    ans -= d * logsubexp( -r * lq, 0.0 );

    if( lower_tail ) ans = logsubexp( 0.0, ans );

    return ans;
}


double ldzinb( unsigned int x, double r, double p, double a )
{
    if( x == 0.0 ) return log(a);
    else           return log(1.0 - a) + lddnbinom( x, r, p );
}



#endif



double pgev( double q, double loc, double scale, double shape, bool lower_tail )
{
    double p;
    q = (q - loc) / scale;
    if( shape == 0.0 ) p = exp( -exp( -q ) );
    else               p = exp( -pow( max( 1 + shape * q, 0.0 ), -1/shape ) );

    if( lower_tail ) return 1.0 - p;
    else             return p;
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


