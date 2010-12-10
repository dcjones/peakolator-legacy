
#ifndef PEAKOLATOR_PEAKOLATOR_SCANNER
#define PEAKOLATOR_PEAKOLATOR_SCANNER

#include "common.hpp"
#include "context.hpp"
#include "intervals.hpp"
#include "parameters.hpp"

#include <cstdlib>

#include "samtools/sam.h"



class emppval;


class scanner
{
    public:
        scanner( 
                parameters* params,
                context*    context );

        interval_stack* run();
        subinterval least_likely_interval( pos l, pos r, double alpha );

        ~scanner();


    private:
        /* raw p-value for subinterval of duration d with read count x */
        mpfr_class QX( double r, rcount x );

        /* push a new bound onto the priority queue if it is valid
         * (i.e. within [d_min,d_max] and lower bound does not rule out.) */
        void conditional_push_copy(
                 subinterval_bound_pqueue& q,
                 subinterval_bound& x,
                 const mpfr_class& p_max );

        context*    ctx;
        parameters* params;

        void     trim_candidate( subinterval& I );
        emppval* padj;
};




#endif

