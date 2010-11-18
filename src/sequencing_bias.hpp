
#ifndef PEAKOLATOR_BIAS_CORRECTION
#define PEAKOLATOR_BIAS_CORRECTION


#include "common.hpp"
#include "kmers.hpp"
#include "hash.h"

#include "samtools/sam.h"
#include "samtools/faidx.h"

class sequencing_bias
{
    public:
        sequencing_bias( const char* ref_fn,
                         const char* reads_fn,
                         size_t n, pos L, pos R );

        ~sequencing_bias();

        double* get_bias( const char* seqname, pos start, pos end, int strand );

        sequencing_bias* copy() const;
        void clear();

    private:
        sequencing_bias();
        void build( const char* ref_fn,
                    const char* reads_fn,
                    size_t n, pos L, pos R );

        void back_step( char** seqs );
        double ll( char** seqs );

        void hash_reads( table* T, samfile_t* reads_fn,
                         size_t limit = 0 ) const;
        void sample_foreground( char* seq, size_t seqlen,
                                struct hashed_value* v,
                                bool count_dups = true );
        void sample_background( char* seq, size_t seqlen,
                                struct hashed_value* v,
                                bool count_dup = true );
        void compute_ws();
        void markov_normalize( double* g, double* g_markov );

        static const double pseudocount;

        unsigned int n;    /* number observations to train on */
        pos L, R; /* left and right sequence context */

        char*   local_seq; /* used while sampling */

        faidx_t* ref_f;

        char* ref_fn;
        char* reads_fn;
};


#endif
