#ifndef PEAKOLATOR_PARAMETERS
#define PEAKOLATOR_PARAMETERS

#include "common.hpp"
#include "nulldist.hpp"


class emppval;

class parameters
{
    public:
        parameters();
        parameters( const parameters& );
        ~parameters();

        parameters* copy() const;

        void rebuild_lookup( int m = 0, int n = 0 );
        void build_padj();

        /* adjusted p-value at which we consider a sub-interval significant */
        double alpha;

        /* zero-inflated negative binomial parameters */
        double r,p,a;

        /* Hard limits on the duration of segments. */
        /* IMPORTANT: for technical reasons, d_min > 2, must be true.
         * There is no upper limit for d_max, but d_min <= d_max must be true. */
        pos d_min; 
        pos d_max;

        /* distribution lookup */
        nulldist dist;

        /* pvalue model */
        emppval* padj;

    private:
        void init_fsolve();
};



#endif



