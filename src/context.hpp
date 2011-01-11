#ifndef PEAKOLATOR_CONTEXT
#define PEAKOLATOR_CONTEXT

#include "common.hpp"
#include "dataset.hpp"
#include "nulldist.hpp"

#include <cstdio>
#include "samtools/bam.h"


/* constructed by dataset */
class context
{
    public:
        context();
        void set( dataset* ds, const char* chrom,
                  pos start, pos end, int strand );
        void set( dataset* ds, const interval& );
        void set_noise( nulldist&, pos len );

        void clear();

        void adjust_interval_by_coverage( interval& ) const;
        //void adjust_subinterval_by_coverage( subinterval& ) const;

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

        const char* get_seqname() const { return seqname; }
        pos         get_start() const   { return start; }
        pos         get_end() const     { return end; }
        int         get_strand() const  { return strand; }

    private:

        dataset* ds;

        pos   start, end;
        char* seqname;
        int   strand;

        /* read counts for both strands */
        rcount* xs[2];

        /* temporary coverage vector */
        void get_coverage( pos u, pos v, int s ) const;
        rcount* cs;

        /* sequencing biases for both strands */
        double* rs[2]; 

        friend int bam_fetch_callback( const bam1_t*, void*);
};


#endif


