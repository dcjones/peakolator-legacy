
#ifndef PEAKOLATOR_MISCMATH
#define PEAKOLATOR_MISCMATH

#include <cstdlib>

/* Miscelaneous math functions. */

/* Compute log( exp(a) + exp(b) ) avoiding underflow. */
double logaddexp( double a, double b );


/* Compute log( exp(a) - exp(b) ) avoiding underflow. */
double logsubexp( double a, double b );


#ifdef HAVE_LIBGSL
/* Compute: log( I_x( a, b ) ), where I is the regularized incomplete beta
 * function.  Equivalent to log( gsl_sf_beta_inc( a, b, x ) ), but computed in
 * such a way as to avoid underflow.
 */
double log_beta_inc( double a, double b, double x );


/* The binomial log density function */
double ldbinom( unsigned int x, double p, unsigned int n );

/* The binomial log distribution function */
double lpbinom( unsigned int x, double p, unsigned int n, bool lower_tail = true );


/* The negative binomial log density function */
double ldnbinom( unsigned int x, double r, double p );

/* The negative binomial log comulative distribution function. */
double lpnbinom( unsigned int q, double r, double p, bool lower_tail = false );



/* The truncated negative binomial log density function */
double ldtnbinom( unsigned int x, double r, double p,
                 unsigned int k );

/* The truncated negative binomial log distribution function. */
double lptnbinom( unsigned int q, double r, double p,
                 unsigned int k, bool lower_tail = false );



/* binomial coefficient */
double binco( double n, double k );

/* log binomial coefficient */
double lbinco( double n, double k );


/* The log density function of the decapitated negative binomial distribution.
 */
double lddnbinom( unsigned x, double r, double p );


/* The density function of the the distribution over the sum of d i.i.d.
 * decapitated negative binomial variables. */
double ddnbsum( unsigned int x, double r, double p, unsigned int d );

/* The distribution function of the the distribution over the sum of d i.i.d.
 * decapitated negative binomial variables. */
double pdnbsum( unsigned int x, double r, double p,
                unsigned int d, bool lower_tail = false );


/* The log density function of the the distribution over the sum of d i.i.d.
 * decapitated negative binomial variables. */
double lddnbsum( unsigned int x, double r, double p, unsigned int d );

/* The log distribution function of the the distribution over the sum of d
 * i.i.d.  decapitated negative binomial variables. */
double lpdnbsum( unsigned int x, double r, double p,
                 unsigned int d, bool lower_tail = false );


/* log density of the zero inflated negative binomial distribution Where f( x;
 * r, p ) is the decapitated negative binomial density function, then density
 * function for the zinb distribution is,
 *
 *  g( x; r, p, a ) = a                           if x = 0
 *                    (1.0 - a) * f( x; r, p )    otherwise
 */
double ldzinb( size_t x, double r, double p, double a );




/* The log distribution function over the joint distribution of
 * the 'joint summed zero-inflated negative binomial distribution'.
 * If x_1, ..., x_n ~ zinb( r, p, a ), this function computes
 *
 * Pr( y, z | r, p, a )
 *
 * Where y = z_1 + ... + z_n and
 *       z = #{ z_i | z_i = 0 }
 *
*/
//double lpjzinbsum( size_t x, size_t z, double r, double p, double a, size_t d );


/* TODO */
/*
 * -> log density function for zero inflated negative binomial i.e.  log( P( x | d ) )
 * -> log distribution function for P( x, z |  d ) = P( x | z, d ) P( z | d )
 *
 * -> What about lookup tables? The scan is going to be horrendously slow if we
 *    must compute the incomplete beta function several million times to test
 *    one interval.
 *
 *    One possible scheme would be to simple build a huge (static?) lookup
 *    table for the incomplete beta function, but would this really help when we
 *    are computing the sum decabitated nb over a very large region?
 *
 *    Another possibility would be to see how the terms of the sum decapitated
 *    nb vary and cut it short when they are tiny (assuming they progress
 *    monotonically).
 */



#endif


/* Generalived extreme value cumulative distribution function. */
double pgev( double q, double loc, double scale, double shape, bool lower_tail = false );


/* Generalized extreme value log-density function. */
double ldgev( double x, double loc, double scale, double shape );



/* (very) simple vector/matrix arithmatic */
void colcpy( double* dest, const double* src, size_t j, size_t n, size_t m );
void vecadd( double* u, double* v, size_t n );
void vecsub( double* u, double* v, size_t n );
void vecaddcol( double* u, double* V, size_t n, size_t m, size_t j );
void vecsubcol( double* u, double* V, size_t n, size_t m, size_t j );
void matsetcol( double* U, double* v, size_t n, size_t m, size_t j );

#endif


