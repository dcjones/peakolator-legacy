
/*
 * Pvalue adjustments using Monte Carlo methods.
 */



#ifndef PEAKOLATOR_PVAL
#define PEAKOLATOR_PVAL

#include "common.hpp"

#include <nlopt.h>

class parameters;

class emppval
{
    public:
        emppval( parameters* params );
        emppval( const emppval& );
        ~emppval();
        mpfr_class adjust( const mpfr_class&, pos len ) const;

    private:
        bool fit_gev_to_emperical( double* mu, double* sigma, double* xi );

        parameters* params;

        int n;
        pos spacing;

        /* results from monte carlo trials */
        double*  qx_mc;

        /* generalized extreme value distribution parameters */
        /* linear models for the dependence of the extreme value distribution
         * parameters on the length of the sequence being scanned. */
        mpfr_class c_mu[2];
        mpfr_class c_sigma[2];
        mpfr_class c_xi[2];

        friend double gev_objf( unsigned int, const double*, double*, void* );
};

#endif

