
#include "sequencing_bias.hpp"
#include "logger.h"
#include "common.hpp"

#include <cmath>
#include <cctype>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "samtools/faidx.h"

#include <algorithm>
using namespace std;



double sq( double x ) { return x*x; }


/* pseudocount used when sampling foreground and background nucleotide * frequencies */
const double sequencing_bias::pseudocount = 1;


/* compute the Kullback-Leibler divergance for two multinomial distributions */
double kl_div( const double* a, const double* b, size_t n )
{
    double kl = 0.0;
    uint32_t i;
    for( i = 0; i < n; i++ ) {
        kl += a[i] * log( a[i] / b[i] );
    }
    return kl;
}


sequencing_bias::sequencing_bias()
    : ref_f(NULL)
    , ref_fn(NULL)
    , reads_fn(NULL)
{}



sequencing_bias::sequencing_bias( const char* ref_fn,
                                  const char* reads_fn,
                                  size_t n, pos L, pos R )
    : ref_f(NULL)
    , ref_fn(NULL)
    , reads_fn(NULL)
{
    build( ref_fn, reads_fn, n, L, R );
}

sequencing_bias* sequencing_bias::copy() const
{
    sequencing_bias* sb = new sequencing_bias();
    sb->n = n;
    sb->L = L;
    sb->R = R;

    sb->ref_fn   = ref_fn   ? strdup(ref_fn)   : NULL;
    sb->reads_fn = reads_fn ? strdup(reads_fn) : NULL;

    if( ref_fn ) {
        sb->ref_f = fai_load(ref_fn);
        if( sb->ref_f == NULL ) {
            log_printf( LOG_ERROR, "Can't open fasta file '%s'\n", ref_fn );
            exit(1);
        }
    }
    else sb->ref_f = NULL;

    return sb;
}



void sequencing_bias::clear()
{
    if( ref_f ) {
        fai_destroy(ref_f);
        ref_f = NULL;
    }
    free(ref_fn);  ref_fn   = NULL;
    free(reads_fn);reads_fn = NULL;
}


void sequencing_bias::build( const char* ref_fn,
                             const char* reads_fn,
                             size_t n, pos L, pos R )
{
    log_puts( LOG_MSG, "Determining sequencing bias...\n" );
    log_indent();

    clear();
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);

    this->ref_fn   = strdup(ref_fn);
    this->reads_fn = strdup(reads_fn);
    
    this->L = L;
    this->R = R;

    unsigned int i;

    samfile_t* reads_f = samopen( reads_fn, "rb", NULL );
    if( reads_f == NULL ) {
        log_printf( LOG_ERROR, "Can't open bam file '%s'.\n", reads_fn );
        exit(1);
    }

    table T;
    hash_reads( &T, reads_f, n );

    /* resort the remaining (1-q)*n by position */
    log_puts( LOG_MSG, "sorting by position ... " );
    struct hashed_value** S;
    table_sort_by_position( &T, &S );
    log_puts( LOG_MSG, "done.\n" );


    /* sample foreground and background kmer frequencies */
    log_puts( LOG_MSG, "sampling sequence bias ...\n" );
    log_indent();


    ref_f = fai_load(ref_fn);
    if( ref_f == NULL ) {
        log_printf( LOG_ERROR, "Can't open fasta file '%s'\n", ref_fn );
        exit(1);
    }

    std::deque<sequence*> fg;
    std::deque<sequence*> bg;


    pos bg_offset;
    char*          seqname   = NULL;
    int            seqlen    = 0;
    int            curr_tid  = -1;
    char*          seq       = NULL;

    local_seq = (char*)safe_malloc(sizeof(char)*(L+R+2));
    local_seq[L+R+1] = '\0';


    for( i = 0; i < n; i++ ) {

        /* Load/switch sequences (chromosomes) as they are encountered in the
         * read stream. The idea here is to avoid thrashing by loading a large
         * sequence, but also avoid overloading memory by only loading one
         * chromosome at a time. */
        if( S[i]->pos.tid != curr_tid ) {
            seqname = reads_f->header->target_name[S[i]->pos.tid];
            seqlen  = reads_f->header->target_len[S[i]->pos.tid];
            if( seq ) free(seq); 

            log_printf( LOG_MSG, "reading sequence '%s' ... ", seqname );

            seq = faidx_fetch_seq( ref_f, seqname, 0, seqlen-1, &seqlen );

            if( seq == NULL ) {
                log_puts( LOG_WARN, "warning: reference sequence not found, skipping.\n" );
            }
            else {
                for( char* c = seq; *c; c++ ) *c = tolower(*c);
            }

            curr_tid = S[i]->pos.tid;

            log_puts( LOG_MSG, "done.\n" );
        }

        if( seq == NULL ) continue;
        
        /* add a foreground sequence */
        if( S[i]->pos.strand ) {
            if( S[i]->pos.pos < R ) continue;
            memcpy( local_seq, seq + S[i]->pos.pos - R, (L+1+R)*sizeof(char) );
            seqrc( local_seq, L+1+R );
        }
        else {
            if( S[i]->pos.pos < L ) continue;
            memcpy( local_seq, seq + (S[i]->pos.pos-L), (L+1+R)*sizeof(char) );
        }

        fg.push_back( new sequence( local_seq ) );


        /* add a background sequence */
        /* adjust the current read position randomly, and sample */
        bg_offset = (pos)gsl_ran_flat( rng, 50.0, 100.0 );
        if( gsl_ran_flat( rng, -1.0, 1.0 ) < 0.0 ) bg_offset = -bg_offset;

        S[i]->pos.pos += bg_offset;

        /* add a foreground sequence */
        if( S[i]->pos.strand ) {
            if( S[i]->pos.pos < R ) continue;
            memcpy( local_seq, seq + S[i]->pos.pos - R, (L+1+R)*sizeof(char) );
            seqrc( local_seq, L+1+R );
        }
        else {
            if( S[i]->pos.pos < L ) continue;
            memcpy( local_seq, seq + (S[i]->pos.pos-L), (L+1+R)*sizeof(char) );
        }

        bg.push_back( new sequence( local_seq ) );
    }


    size_t max_k = 5;
    motif M0( L+1+R, max_k, &bg );
    motif M1( L+1+R, max_k, &fg );

    train_motifs( M0, M1 );


    /* TODO: keep these, clean them up with the destructor */
    std::deque<sequence*>::iterator i_seq;
    for( i_seq = fg.begin(); i_seq != fg.end(); i_seq++ ) {
        delete *i_seq;
    }

    for( i_seq = bg.begin(); i_seq != bg.end(); i_seq++ ) {
        delete *i_seq;
    }


    log_unindent();

    free(seq);
    free(local_seq);
    samclose(reads_f);
    table_destroy(&T);
    free(S);
    gsl_rng_free(rng);

    log_unindent();
}


sequencing_bias::~sequencing_bias()
{
    clear();
}



void sequencing_bias::hash_reads( table* T, samfile_t* reads_f, size_t limit ) const
{
    log_puts( LOG_MSG, "hashing read positions..." );

    table_create(T);

    bam1_t* read = bam_init1();

    while( (limit == 0 || T->m < limit) && samread( reads_f, read ) > 0 ) {
        table_inc( T, read );
    }

    bam_destroy1(read);
}


double* sequencing_bias::get_bias( const char* seqname, pos start, pos end, int strand )
{
    /* TODO */
    return NULL;
}

