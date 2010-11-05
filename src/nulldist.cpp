
#include "nulldist.hpp"
#include "logger.h"

#include <ctime>
#include <gsl/gsl_randist.h>
#include <boost/math/distributions/negative_binomial.hpp>


nulldist::nulldist()
    : p(0.0), m(2000), n(100), A(NULL)
{
    /* initialize rng */
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,time(NULL));
}

bool nulldist::ready() const
{
    return A != NULL;
}


void nulldist::build( double r, double p, int m, int n )
{
    if( A ) delete[] A;

    if( m != 0 ) this->m = m; else m = this->m;
    if( n != 0 ) this->n = n; else n = this->n;

    this->p = p;
    this->r = r;
    A = new mpfr_class[m*n];

    boost::math::negative_binomial_distribution<mpfr_class> dist( 1.0, 1.0 );

    int i,j;
    A[0] = 1.0;
    for( i = 1; i < m; i++ ) {
        dist = boost::math::negative_binomial_distribution<mpfr_class>( (mpfr_class)i, (mpfr_class)p );

        A[i*n+0] = 1.0;
        for( j = 1; j < n; j++ ) {
            A[i*n+j] = boost::math::cdf( boost::math::complement( dist, j-1 ) );
        }
    }
}


nulldist::nulldist( const nulldist& y )
{
    r = y.r;
    p = y.p;

    m = y.m;
    n = y.n;
    A = new mpfr_class[m*n];

    int i,j;
    for( i = 0; i < m; i++ )  {
        for( j = 0; j < n; j++ ) {
            A[i*n+j] = y.A[i*n+j];
        }
    }
}



void nulldist::operator=( const nulldist& y )
{
    if( A ) delete[] A;

    r = y.r;
    p = y.p;

    m = y.m;
    n = y.n;
    A = new mpfr_class[m*n];

    int i,j;
    for( i = 0; i < m; i++ )  {
        for( j = 0; j < n; j++ ) {
            A[i*n+j] = y.A[i*n+j];
        }
    }
}


nulldist::~nulldist()
{
    if(A) delete[] A;
    gsl_rng_free(rng);
}



mpfr_class nulldist::QX( double r_i, rcount x_i )
{
    if( x_i == 0 || r*r_i <= 0.0 ) return 1.0;

    long int rd = (long int)ceil(r*r_i);

    if( A == NULL || rd >= m || x_i >= (rcount)n ) {
        boost::math::negative_binomial_distribution<mpfr_class> dist( r*r_i, p );
        return boost::math::cdf( boost::math::complement( dist, x_i-1 ) );
    }
    else {
        /* linear interpolation between ceil(r) and floor(r) */
        double z =  (r*r_i) - (rd-1);
        return A[rd*n+x_i]*z + A[(rd-1)*n+x_i]*(1.0-z);
    }
}



rcount nulldist::rand( double r_i )
{
    return (rcount)gsl_ran_negative_binomial( rng, p, r*r_i );
}


