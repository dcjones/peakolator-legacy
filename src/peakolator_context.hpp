#ifndef PEAKOLATOR_CONTEXT
#define PEAKOLATOR_CONTEXT

#include "common.hpp"
#include "peakolator_dataset.hpp"
#include "peakolator_dist.hpp"

#include <cstdio>
#include "samtools/bam.h"


/* constructed by peakolator_dataset */
class peakolator_context
{
    public:
        peakolator_context();
        void set( peakolator_dataset* dataset, const char* chrom,
                  pos start, pos end, int strand );
        void set( peakolator_dataset* dataset, const interval& );
        void set_noise( peakolator_dist&, pos len );

        void clear();


        void print_adjusted_unadjusted_bias( FILE* out = stdout );

        /* get total readcount for subinterval, coordinates are relative to the start */
        rcount count() const;
        rcount count( pos i ) const;
        rcount count( pos i, pos j ) const;
        rcount count( pos i, pos j, int strand ) const;

        double rate() const;
        double rate( pos i ) const;
        double rate( pos i, pos j ) const;
        double rate( pos i, pos j, int strand ) const;

        double min_rate( const subinterval_bound& B, pos d_min ) const;
        double min_rate( const subinterval_bound& B, pos d_min, int strand ) const;

        pos length() const;

        pos   start, end;
        char* seqname;
        int   strand;

        /* read counts for both strands */
        rcount* xs[2];

        /* sequencing biases for both strands */
        double* rs[2]; 

        friend int bam_fetch_callback( const bam1_t*, void*);
};


#endif


