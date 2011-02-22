
#include "nulldist.hpp"
#include "miscmath.hpp"
#include "logger.h"

#include <cmath>
#include <ctime>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>


nulldist::nulldist()
    : p(0.0), m(2000), n(100), A(NULL)
{
    /* initialize rng */
    rng = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng,time(NULL));
}

bool nulldist::ready() const
{
    return A != NULL;
}


void nulldist::build( double r, double p, double a, int m, int n )
{
    if( A ) delete[] A;

    if( m != 0 ) this->m = m; else m = this->m;
    if( n != 0 ) this->n = n; else n = this->n;

    this->p = p;
    this->r = r;
    this->a = a;
    A = new double[m*n];

    int i,j;
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ ) {
            A[i*n+j-1] = lpdnbsum( j+1, r, p, i+1 );
        }
    }
}


nulldist::nulldist( const nulldist& y )
{
    r = y.r;
    p = y.p;
    a = y.a;

    m = y.m;
    n = y.n;
    A = new double[m*n];

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
    a = y.a;

    m = y.m;
    n = y.n;
    A = new double[m*n];

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



double nulldist::QX( rcount x, unsigned int z, unsigned int d )
{
    if( x == 0 || z >= d ) return 0.0;

    double ans;

    if( A == NULL || x >= (rcount) n || d - z >= m ) {
        ans = lpdnbsum( x, r, p, d - z );
    }
    else {
        ans = A[ (d-z-1)*n + (x-1) ];
    }

    ans += lpbinom( z, a, d, true );

    return ans;
}


rcount nulldist::rand()
{
    if( gsl_ran_bernoulli( rng, a ) ) return 0;

    /* This uses a rejection sampling scheme worked out by Charles Geyer,
     * detailed in his notes "Lower-Truncated Poisson and Negative Binomial
     * Distributions".
     */
    rcount x;
    double accp;
    while( true ) {
        x = gsl_ran_negative_binomial( rng, p, r + 1.0 ) + 1;

        accp = gsl_sf_lnfact( x - 1 ) - gsl_sf_lnfact( x );

        if( gsl_ran_bernoulli( rng, exp( accp ) ) ) break;
    }

    return x;
}




