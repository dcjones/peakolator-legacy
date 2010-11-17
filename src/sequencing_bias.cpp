
#include "sequencing_bias.hpp"
#include "logger.h"
#include "common.hpp"

#include <cmath>
#include <cctype>
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
    : ws(NULL)
    , fg(NULL)
    , bg(NULL)
    , bgs(NULL)
    , bg_div(NULL)
    , ref_f(NULL)
    , ref_fn(NULL)
    , reads_fn(NULL)
{}



sequencing_bias::sequencing_bias( const char* ref_fn,
                                  const char* reads_fn,
                                  pos L, pos R, unsigned int k,
                                  bool count_dups, double q,
                                  const char* training_seqname )
    : ws(NULL)
    , fg(NULL)
    , bg(NULL)
    , bgs(NULL)
    , bg_div(NULL)
    , ref_f(NULL)
    , ref_fn(NULL)
    , reads_fn(NULL)
{
    build( ref_fn, reads_fn, L, R , k, count_dups, q, training_seqname );
}

sequencing_bias* sequencing_bias::copy() const
{
    sequencing_bias* sb = new sequencing_bias();
    sb->n = n;
    sb->m = m;
    sb->k = k;
    sb->L = L;
    sb->R = R;
    sb->kmer_mask = kmer_mask;
    sb->four_to_k = four_to_k;
    sb->bg_left   = bg_left;
    sb->bg_right  = bg_right;
    sb->bg_len    = bg_len;

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

    sb->ws = (double*)safe_malloc( m*sizeof(double) );
    memcpy( sb->ws, ws, m*sizeof(double) );

    sb->fg = (double*)safe_malloc( m*sizeof(double) );
    memcpy( sb->fg, fg, m*sizeof(double) );

    sb->fg = (double*)safe_malloc( four_to_k*sizeof(double) );
    memcpy( sb->fg, fg, four_to_k*sizeof(double) );

    sb->bgs = (double*)safe_malloc( four_to_k*2*bg_len*sizeof(double) );
    memcpy( sb->bgs, bgs, four_to_k*2*bg_len*sizeof(double) );

    return sb;
}



void sequencing_bias::clear()
{
    free(ws);     ws     = NULL;
    free(fg);     fg     = NULL;
    free(bg);     bg     = NULL;
    free(bgs);    bgs    = NULL;
    free(bg_div); bg_div = NULL;
    if( ref_f ) {
        fai_destroy(ref_f);
        ref_f = NULL;
    }
    free(ref_fn);  ref_fn   = NULL;
    free(reads_fn);reads_fn = NULL;
}


void sequencing_bias::build( const char* ref_fn,
                             const char* reads_fn,
                             pos L, pos R, unsigned int k,
                             bool count_dups, double q,
                             const char* training_seqname )
{
    log_puts( LOG_MSG, "Determining sequencing bias...\n" );
    log_indent();

    clear();

    this->ref_fn   = strdup(ref_fn);
    this->reads_fn = strdup(reads_fn);
    
    this->bg_len = 50;
    this->bg_left  = 50;
    this->bg_right = 50;

    this->k = k;
    this->L = L;
    this->R = R;

    unsigned int i, j;

    kmer_mask = 0;
    for( i = 0; i < k; i++ ) kmer_mask = (kmer_mask<<2) | 0x3;

    four_to_k = 1<<(2*k);

    m = four_to_k*(L+1+R);

    /* currently, I am packing kmers into 16bits, using two bit encoding, thus
     * k <= 8 must be true. */
    if( k > 4*sizeof(kmer) ) {
        log_printf( LOG_ERROR, "Assertion failed: k <= %d\n", 4*sizeof(kmer) );
        exit(1);
    }

    samfile_t* reads_f = samopen( reads_fn, "rb", NULL );
    if( reads_f == NULL ) {
        log_printf( LOG_ERROR, "Can't open bam file '%s'.\n", reads_fn );
        exit(1);
    }

    table T;
    hash_reads( &T, reads_f, training_seqname );





    /* determine top percentile to avoid the sample being thrown off by a few
     * extremely abundant reads */
    struct hashed_value** S;

    log_puts( LOG_MSG, "sorting by count ... " );
    table_sort_by_count( &T, &S );
    log_puts( LOG_MSG, "done.\n" );


    /* ignore the top q*n */
    unsigned int n = T.m - (unsigned int)((double)T.m * q);

    /* resort the remaining (1-q)*n by position */
    log_puts( LOG_MSG, "sorting by position ... " );
    qsort( S, n, sizeof(struct hashed_value*), compare_hashed_pos );
    log_puts( LOG_MSG, "done.\n" );


    /* sample foreground and background kmer frequencies */
    log_puts( LOG_MSG, "sampling sequence bias ...\n" );
    log_indent();


    ref_f = fai_load(ref_fn);
    if( ref_f == NULL ) {
        log_printf( LOG_ERROR, "Can't open fasta file '%s'\n", ref_fn );
        exit(1);
    }


    fg  = (double*)safe_malloc( m*sizeof(double) );
    memset( fg, 0, m*sizeof(double) );

    bgs = (double*)safe_malloc(sizeof(double)*four_to_k*2*bg_len);
    memset( bgs, 0, sizeof(double)*four_to_k*2*bg_len );


    /* initialize using pseudocounts */
    for( i = 0; i < m; i++ ) fg[i] = pseudocount;
    for( i = 0; i < four_to_k*2*bg_len; i++ ) bgs[i] = pseudocount;




    char*          seqname   = NULL;
    int            seqlen    = 0;
    int            curr_tid  = -1;
    char*          seq       = NULL;

    local_seq = (char*)safe_malloc(sizeof(char)*(L+R+bg_len+k+1));


    for( i = 0; i < n; i++ ) {

        /* Load/switch sequences (chromosomes) as they are encountered in the
         * read stream. The idea here is to avoid thrashing by loading a large
         * sequence, but also avoid overloading memory by only loading one
         * chromosome at a time. */
        if( S[i]->pos.tid != curr_tid ) {
            seqname = reads_f->header->target_name[S[i]->pos.tid];
            seqlen  = reads_f->header->target_len[S[i]->pos.tid];
            if( seq ) free(seq); 

            log_printf( LOG_MSG, "reading sequence '%s'...\n", seqname );

            seq = faidx_fetch_seq( ref_f, seqname, 0, seqlen-1, &seqlen );

            if( seq == NULL ) {
                log_puts( LOG_WARN, "warning: reference sequence not found, skipping.\n" );
            }
            else {
                for( char* c = seq; *c; c++ ) *c = tolower(*c);
            }

            curr_tid = S[i]->pos.tid;
        }

        if( seq == NULL ) continue;


        sample_foreground( seq, seqlen, S[i], count_dups );
        sample_background( seq, seqlen, S[i], count_dups );
    }

    free(seq);
    free(local_seq);
    samclose(reads_f);
    table_destroy(&T);
    free(S);


    /* normalize background kmer frequencies */
    bg = (double*)safe_malloc(sizeof(double)*four_to_k);
    memset( bg, 0, sizeof(double)*four_to_k );

    for( i = 0; i < (unsigned int)2*bg_len; i++ ) {
        for( j = 0; j < four_to_k; j++ ) {
            bg[j] += bgs[i*four_to_k+j];
        }
    }

    double z = 0.0;


    /* normalize background */
    for( j = 0; j < four_to_k; j++ ) z += bg[j];
    for( j = 0; j < four_to_k; j++ ) bg[j] /= z;

    /* normalize foreground */
    for( i = 0; i < (unsigned int)(L+R+1); i++ ) {
        z = 0.0;
        for( j = 0; j < four_to_k; j++ ) z += fg[i*four_to_k+j];
        for( j = 0; j < four_to_k; j++ ) fg[i*four_to_k+j] /= z;
    }

    /* normalize positional background  */
    for( i = 0; i < (unsigned int)2*bg_len; i++ ) {
        z = 0.0;
        for( j = 0; j < four_to_k; j++ ) z += bgs[i*four_to_k+j];
        for( j = 0; j < four_to_k; j++ ) bgs[i*four_to_k+j] /= z;
    }



    /* compute Kullback-Leibler distance for each background position */
    bg_div = (double*)safe_malloc(sizeof(double)*2*bg_len);
    for( i = 0; i < (unsigned int)2*bg_len; i++ ) {
        bg_div[i] = kl_div( bgs+(i*four_to_k), bg, four_to_k );
    }

    log_unindent();

    compute_ws();

    log_unindent();


    /*
     * TODO:
     * Here's the plan. We do a step back procedure, doing a greedy reduction of
     * parameters, evaluating liklihood on several thousand randomly chosen
     * reads, and minimizing the AIC.
     *
     * 1. sample 10000 sequences
     * 2. write a function to evaluate liklihood
     * 3. think about how to enumerate parameters and zero them out.
     *
     */
}


sequencing_bias::~sequencing_bias()
{
    clear();
}



void sequencing_bias::sample_foreground( char* seq, size_t seqlen,
                                         struct hashed_value* v,
                                         bool count_dups )
{
    pos j;

    if( v->pos.strand ) {
        if( v->pos.pos < R ) return;
        memcpy( local_seq, seq + v->pos.pos - R, (L+R+k)*sizeof(char) );
        seqrc( local_seq, L+R+k );
    }
    else {
        if( v->pos.pos < L+(pos)(k-1) ) return;
        memcpy( local_seq, seq + (v->pos.pos-L-(k-1)), (L+R+k)*sizeof(char) );
    }

    /* sample kmers */
    kmer K = 0;
    for( j = 0; j < L+R+(pos)k; j++ ) {
        K = ((K<<2) | nt2num(local_seq[j])) & kmer_mask;
        if( j >= (pos)(k-1) )
            fg[ (j-(pos)(k-1))*four_to_k + K ] += count_dups ? v->count : 1;
    }
}


void sequencing_bias::sample_background( char* seq, size_t seqlen,
                                         struct hashed_value* v,
                                         bool count_dups )
{
    kmer K;
    unsigned int j;

    /* sample left of read start */

    local_seq[0] = '\0';
    if( v->pos.strand ) {
        if( (size_t)(v->pos.pos + bg_left) < seqlen )  {
            memcpy( local_seq,
                    seq + v->pos.pos + bg_right,
                    (bg_len+k-1)*sizeof(char) );
            seqrc( local_seq, bg_len+k-1 );
        }
    }
    else {
        if( v->pos.pos > bg_len+bg_left-(pos)(k-1) ) {
            memcpy( local_seq,
                    seq + v->pos.pos - bg_len - bg_left - (pos)(k-1),
                    (bg_len+k-1)*sizeof(char) );
        }
    }

    if( local_seq ) {
        K = 0;
        for( j = 0; j < (unsigned int)bg_len+(k-1); j++ ) {
            K = ((K<<2) | nt2num(local_seq[j])) & kmer_mask;
            if( j >= k-1 )
                bgs[four_to_k*(j-(k-1))+K] += count_dups ? v->count : 1;
        }
    }



    /* sample right of read start */

    local_seq[0] = '\0';
    if( v->pos.strand ) {
        if( v->pos.pos > bg_right + bg_len ) {
            memcpy( local_seq,
                    seq + v->pos.pos - (bg_right + bg_len),
                    (bg_len+k-1)*sizeof(char) );
            seqrc( local_seq, bg_len+k-1 );
        }
    }
    else {
        if( (size_t)(v->pos.pos + bg_right - (k-1)) < seqlen ) {
            memcpy( local_seq,
                    seq + v->pos.pos + bg_right,
                    (bg_len+k-1)*sizeof(char) );
        }
    }

    if( local_seq[0] ) {
        K = 0;
        for( j = 0; j < bg_len+(k-1); j++ ) {
            K = ((K<<2) | nt2num(local_seq[j])) & kmer_mask;
            if( j >= k-1 )
                bgs[four_to_k*(j+bg_len-(k-1))+K] += count_dups ? v->count : 1;
        }
    }
}


void sequencing_bias::hash_reads( table* T, samfile_t* reads_f,
                                  const char* training_seqname ) const
{
    log_puts( LOG_MSG, "hashing read positions..." );


    /* Restrict the table to reads aligned to the given seqname. */
    /* Find the corresponding tid. */
    int i, training_tid = -1;
    if( training_seqname != NULL ) {
        for( i = 0; i < reads_f->header->n_targets; i++ ) {
            if( strcmp( reads_f->header->target_name[i], training_seqname ) == 0 ) {
                training_tid = i;
                break;
            }
        }
    }

    table_create(T);

    bam1_t* read = bam_init1();

    while( samread( reads_f, read ) > 0 ) {
        if( training_tid == -1 || read->core.tid == training_tid ) {
            table_inc( T, read );
        }
    }
    bam_destroy1(read);
}


/* normalize kmer frequencies to form a k-1 order markov * chain */
void sequencing_bias::markov_normalize( double* g, double* g_markov )
{
    unsigned int i,j;
    double z;
    kmer K;

    for( i = 0; i < four_to_k>>2; i++ ) {
        z = 0.0;
        for( j = 0; j < 4; j++ ) {
            K = (i<<2)|j;
            z += g[K];
        }
        for( j = 0; j < 4; j++ ) {
            K = (i<<2)|j;
            g_markov[K] = g[K] /= z;
        }
    }
}


void sequencing_bias::compute_ws()
{
    log_puts( LOG_MSG, "computing posterior probabilities..." );

    double* bg_markov = (double*)safe_malloc( four_to_k*sizeof(double) );
    markov_normalize( bg, bg_markov );

    ws = (double*)safe_malloc( m*sizeof(double) );
    memset( ws, 0, m*sizeof(double) );

    double* fg_markov = (double*)safe_malloc( four_to_k*sizeof(double) );

    pos i;
    unsigned int j;
    for( i = 0; i < L+1+R; i++ ) {
        markov_normalize( fg+i*four_to_k, fg_markov );

        for( j = 0; j < four_to_k; j++ ) {
            ws[i*four_to_k+j] = fg_markov[j] / bg_markov[j];
        }
    }

    free(fg_markov);
    free(bg_markov);
}





double* sequencing_bias::get_bias( const char* seqname, pos start, pos end, int strand )
{
    if( strand < 0 || ref_f == NULL || ws == NULL ) return NULL;

    int seqlen = end-start+1;
    kmer K;
    pos i, j;

    pos left_offset = L + (pos)(k-1);
    double* bias = (double*)safe_malloc( seqlen*sizeof(double) );

    for( i = 0; i < seqlen; i++ ) bias[i] = 1.0;

    if( strand == 1 ) {
        local_seq = faidx_fetch_seq_forced_lower( ref_f, seqname,
                                                  start-R, end+L+(k-1) );
        if( local_seq ) seqrc( local_seq, seqlen );
    }
    else {
        local_seq = faidx_fetch_seq_forced_lower( ref_f, seqname,
                                                  start-left_offset, end+R );
    }

    if( local_seq ) {
        /* i,j are relative to start */
        K = 0;
        for( i = -left_offset; i < (end-start+1)+R; i++ ) {

            K = ( (K<<2) | nt2num( local_seq[i+left_offset] ) ) & kmer_mask;

            for( j = max((pos)0,i-R); j <= min((pos)seqlen-1,i+L); j++ ) {
                bias[j] *= ws[ (L+i-j)*four_to_k + K ];
            }
        }

        free(local_seq);

        if( strand == 1 ) {
            rev(bias,seqlen);
        }
    }

    return bias;
}

