

#ifndef PEAKOLATOR_DIST
#define PEAKOLATOR_DIST

#include "common.hpp"
#include <gsl/gsl_rng.h>


class peakolator_parameters;

/* efficiently evaluate a probability distribution / density function using a
 * lookup table and interpolation. */
class peakolator_dist
{
    public:
        peakolator_dist();
        peakolator_dist( const peakolator_dist& );
        ~peakolator_dist();

        void build( double r, double p, int m = 0, int n = 0 );

        bool ready() const;

        void operator=( const peakolator_dist& );

        mpfr_class QX( double r_i, rcount x_i );
        rcount rand( double r_i );


    private:
        double r;
        double p;
        peakolator_parameters* param;
        int m; /* maximum rate */
        int n; /* maximum count */
        mpfr_class* A; /* rate by count */

        /* random number generator */
        gsl_rng* rng;
};


#endif


