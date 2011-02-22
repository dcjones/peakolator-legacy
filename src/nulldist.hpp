

#ifndef PEAKOLATOR_DIST
#define PEAKOLATOR_DIST

#include "common.hpp"
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

        void build( double r, double p, double a, int m = 0, int n = 0 );

        bool ready() const;

        void operator=( const nulldist& );

        double QX( rcount x, unsigned int z, unsigned int d );
        rcount rand();


    private:
        double r;
        double p;
        double a;
        parameters* param;
        int m;     /* maximum rate */
        int n;     /* maximum count */
        double* A; /* rate by count */

        /* random number generator */
        gsl_rng* rng;
};


#endif


