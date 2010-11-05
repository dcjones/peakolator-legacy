
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
        mpfr_class* gumbel_locs;
        mpfr_class* gumbel_scales;

        /* results from monte carlo trials */
        double*  qx_mc;
        double   qx_mc_mean;

        /* generalized extreme value distribution parameters */
        mpfr_class* mu;
        mpfr_class* sigma;
        mpfr_class* xi;
        
        friend double gev_objf( unsigned int, const double*, double*, void* );
};

#endif

