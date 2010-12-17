
/*
 * Pvalue adjustments using Monte Carlo methods.
 *
 * If a interval of length l is scanned, discovering the most significant
 * scoring region with score s. A p-value is produced according to a model,
 *
 * log(s) ~ GEV( mu, sigma, xi )
 *
 * Where, GEV is a three-parameter generalized extreme value distribution,
 * and the parameters are linear in the log length. That is,
 *
 *      Position: mu    = mu0    + mu1    * log(l)
 *      Scale:    sigma = sigma0 + sigma1 * log(l)
 *      Shape:    xi    = xi0    + xi1    * log(l)
 *
 * This class handles computing these p-values given the model, as well as
 * performing maximimum likelihood fitting of the model using monte-carlo
 * trials.
 *
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

        mpfr_class c_mu[2];
        mpfr_class c_sigma[2];
        mpfr_class c_xi[2];

        friend double gev_objf( unsigned int, const double*, double*, void* );
};

#endif

