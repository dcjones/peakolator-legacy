
#ifndef PEAKOLATOR_BIAS_CORRECTION
#define PEAKOLATOR_BIAS_CORRECTION


#include "common.hpp"
#include "hash.h"

#include "samtools/sam.h"
#include "samtools/faidx.h"
#include <boost/cstdint.hpp>


/* kmers are encoded in 16 bits, allowing for k <= 8 */
typedef boost::uint16_t kmer;

class sequencing_bias
{
    public:
        sequencing_bias();
        sequencing_bias( const char* ref_fn,
                         const char* reads_fn,
                         pos L, pos R, unsigned int k );

        sequencing_bias* copy() const;

        void clear();

        void build( const char* ref_fn,
                    const char* reads_fn,
                    pos L, pos R, unsigned int k );

        void print_kmer_frequencies( pos L, pos R, unsigned int k,
                                     bool adjusted = true, FILE* fout = stdout ) const;

        double* get_bias( const char* seqname, pos start, pos end, int strand );

        ~sequencing_bias();

    private:
        void hash_reads( table* T, samfile_t* reads_fn ) const;
        void sample_foreground( char* seq, size_t seqlen,
                                struct hashed_value* v );
        void sample_background( char* seq, size_t seqlen,
                                struct hashed_value* v );
        void compute_ws();
        void markov_normalize( double* g, double* g_markov );

        static const double pseudocount;

        unsigned int n;    /* number of positions being considered */
        unsigned int m;    /* size of the weight vector: 4^k * (L+R+1) */
        unsigned int k;    /* kmer size */
        pos L, R; /* left and right sequence context */

        /* posterior frequencing probabilities */
        double* ws;

        /* kmer frequencies surrounding the read start */
        double* fg;

        /* background kmer frequencies */
        double* bg;

        /* kmer bitmask */
        kmer kmer_mask;

        /* 4^k, precomputed for conveniance */
        size_t four_to_k;

        /* background distribution to which we are calibrating */
        double* bgs;    /* kmer frequencies across a window */
        double* bg_div; /* KL divergance sample */

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
