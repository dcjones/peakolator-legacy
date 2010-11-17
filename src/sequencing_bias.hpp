
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
                         pos L, pos R, unsigned int k,
                         bool count_dups = true, double q = 0.1,
                         const char* training_seqname = NULL );

        ~sequencing_bias();

        double* get_bias( const char* seqname, pos start, pos end, int strand );

        sequencing_bias* copy() const;
        void clear();

    private:
        sequencing_bias();
        void build( const char* ref_fn,
                    const char* reads_fn,
                    pos L, pos R, unsigned int k,
                    bool count_dups = true, double q = 0.1,
                    const char* training_seqname = NULL );

        void back_step( char** seqs );
        double ll( char** seqs );

        void hash_reads( table* T, samfile_t* reads_fn,
                         const char* training_seqname = NULL ) const;
        void sample_foreground( char* seq, size_t seqlen,
                                struct hashed_value* v,
                                bool count_dups = true );
        void sample_background( char* seq, size_t seqlen,
                                struct hashed_value* v,
                                bool count_dup = true );
        void compute_ws();
        void markov_normalize( double* g, double* g_markov );

        static const double pseudocount;

        unsigned int n;    /* number of positions being considered */
        unsigned int m;    /* size of the weight vector: 4^k * (L+R+1) */
        unsigned int k;    /* kmer size */
        pos L, R; /* left and right sequence context */

        /* posterior frequencing probabilities */
        kmer_matrix* ws;

        /* kmer frequencies surrounding the read start */
        kmer_matrix* fg;

        /* background kmer frequencies */
        kmer_matrix* bg;

        /* kmer bitmask */
        kmer kmer_mask;

        /* 4^k, precomputed for conveniance */
        size_t four_to_k;


        /* background is sampled by considering two regions of size bg_len.  One
         * 'bg_left' to the left of the read start and the other 'bg_right' to
         * the right. */
        pos bg_left;
        pos bg_right;
        pos bg_len;

        char*   local_seq; /* used while sampling */

        faidx_t* ref_f;

        char* ref_fn;
        char* reads_fn;
};


#endif
