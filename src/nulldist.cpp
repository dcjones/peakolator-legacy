
#include "nulldist.hpp"
#include "miscmath.hpp"
#include "logger.h"

#include <cmath>
#include <ctime>
#include <gsl/gsl_randist.h>


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


void nulldist::build( double r, double p, int m, int n )
{
    if( A ) delete[] A;

    if( m != 0 ) this->m = m; else m = this->m;
    if( n != 0 ) this->n = n; else n = this->n;

    this->p = p;
    this->r = r;
    A = new double[m*n];

    int i,j;
    A[0] = 1.0;
    for( i = 1; i < m; i++ ) {

        A[i*n+0] = 1.0;
        for( j = 1; j < n; j++ ) {
            A[i*n+j] = lpnbinom( j-1, r * (double)i, p );
        }
    }
}


nulldist::nulldist( const nulldist& y )
{
    r = y.r;
    p = y.p;

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



double nulldist::QX( double r_i, rcount x_i )
{
    if( x_i == 0 || r*r_i <= 0.0 ) return 0.0;

    long int rd = (long int)ceil(r*r_i);

    if( A == NULL || rd >= m || x_i >= (rcount)n ) {
        return lpnbinom( x_i - 1, r*r_i, p );
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


