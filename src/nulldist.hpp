

#ifndef PEAKOLATOR_DIST
#define PEAKOLATOR_DIST

#include "common.hpp"
#include "mpfr_bindings.hpp"
#include <gsl/gsl_rng.h>


class parameters;

/* efficiently evaluate a probability distribution / density function using a
 * lookup table and interpolation. */
class nulldist
{
    public:
        nulldist();
        nulldist( const nulldist& );
        ~nulldist();

        void build( double r, double p, int m = 0, int n = 0 );

        bool ready() const;

        void operator=( const nulldist& );

        mpfr_class QX( double r_i, rcount x_i );
        rcount rand( double r_i );


    private:
        double r;
        double p;
        parameters* param;
        int m; /* maximum rate */
        int n; /* maximum count */
        mpfr_class* A; /* rate by count */

        /* random number generator */
        gsl_rng* rng;
};


#endif


