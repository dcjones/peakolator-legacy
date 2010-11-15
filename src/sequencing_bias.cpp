
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
    build( ref_fn, reads_fn, L, R , k, training_seqname );
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

    /* sort by position so we can read the reference sequence one chromosome at
     * a time */
    log_puts( LOG_MSG, "sorting read positions ... " );
    table_sort_by_position( &T, &S );
    log_puts( LOG_MSG, "done.\n" );

    n = T.m;

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


        sample_foreground( seq, seqlen, S[i] );
        sample_background( seq, seqlen, S[i] );
    }

    free(seq);
    free(local_seq);
    samclose(reads_f);
    table_destroy(&T);


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
}


sequencing_bias::~sequencing_bias()
{
    clear();
}



/* sample and print unadjusted and adjusted kmer frequencies. There's a fair
 * amount of code duplication with 'build', unfortuenately... */
void sequencing_bias::print_kmer_frequencies( pos L_, pos R_, unsigned int k_,
                                              bool adjusted, FILE* fout ) const
{
    unsigned int i;
    pos j;

    /* get a few constants out of the way */
    unsigned int four_to_k_ = 1<<(2*k_);
    unsigned int m_ = four_to_k_*(L_+R_+1);
    kmer kmer_mask_ = 0;
    for( i = 0; i < k_; i++ ) kmer_mask_ = (kmer_mask_<<2) | 0x3;


    /* allocation */
    double* freq = (double*)safe_malloc( m_*sizeof(double) );
    memset( freq, 0, m_*sizeof(double) );

    double* freq_adjusted = NULL;
    if( adjusted ){
        freq_adjusted = (double*)safe_malloc( m_*sizeof(double) );
        memset( freq_adjusted, 0, m_*sizeof(double) );
    }

    samfile_t* reads_f = samopen( reads_fn, "rb", NULL );
    if( reads_f == NULL ) {
        log_printf( LOG_ERROR, "Can't open bam file '%s'.", reads_fn );
        exit(1);
    }

    /* hash read locations for efficiency, remove the top percentile to avoid
     * overfitting */
    table T;
    hash_reads( &T, reads_f );
    struct hashed_value** S;
    table_sort_by_position( &T, &S );

    unsigned int n = T.m;


    char*          seqname   = NULL;
    int            seqlen    = 0;
    int            curr_tid  = -1;
    char*          seq       = NULL;
    kmer K;

    double* kl_bg     = (double*)safe_malloc(sizeof(double)*four_to_k_);
    memset( kl_bg, 0, four_to_k_*sizeof(double) );
    char*   kl_bg_seq = (char*)safe_malloc(sizeof(char)*bg_len);

    char* local_seq = (char*)safe_malloc(sizeof(char)*(L_+R_+k_+1));
    char* bias_seq  = (char*)safe_malloc(sizeof(char)*(L+R+k+1));

    for( i = 0; i < n; i++ ) {
        if( S[i]->pos.tid != curr_tid ) {
            seqname = reads_f->header->target_name[S[i]->pos.tid];
            seqlen  = reads_f->header->target_len[S[i]->pos.tid];
            if( seq ) free(seq); 

            log_printf( LOG_MSG, "reading sequence '%s' ...\n", seqname );

            seq = faidx_fetch_seq( ref_f, seqname, 0, seqlen-1, &seqlen );
            for( char* c = seq; *c; c++ ) *c = tolower(*c);

            if( seq == NULL ) {
                log_puts( LOG_WARN, "warning: reference sequence not found, skipping." );
            }

            curr_tid = S[i]->pos.tid;
        }

        if( seq == NULL ) continue;


        hashed_value* v = S[i];

        /* sample background (left) */
        kl_bg_seq[0] = '\0';
        if( v->pos.strand ) {
            if( v->pos.pos + bg_left < seqlen )  {
                memcpy( kl_bg_seq,
                        seq + v->pos.pos + bg_right,
                        (bg_len+k_-1)*sizeof(char) );
                seqrc( kl_bg_seq, bg_len+k_-1 );
            }
        }
        else {
            if( v->pos.pos > bg_len+bg_left-(pos)(k_-1) ) {
                memcpy( kl_bg_seq,
                        seq + v->pos.pos - bg_len - bg_left - (pos)(k_-1),
                        (bg_len+k_-1)*sizeof(char) );
            }
        }

        if( kl_bg_seq[0] ) {
            K = 0;
            for( j = 0; j < bg_len+((pos)k_-1); j++ ) {
                K = ((K<<2) | nt2num(kl_bg_seq[j])) & kmer_mask_;
                if( j >= (pos)k_-1 ) kl_bg[K] += 1;
            }
        }



        /* sample background (right) */
        kl_bg_seq[0] = '\0';
        if( v->pos.strand ) {
            if( v->pos.pos > bg_right + bg_len ) {
                memcpy( kl_bg_seq,
                        seq + v->pos.pos - (bg_right + bg_len),
                        (bg_len+k_-1)*sizeof(char) );
                seqrc( kl_bg_seq, bg_len+k_-1 );
            }
        }
        else {
            if( v->pos.pos + bg_right - ((pos)k_-1) < seqlen ) {
                memcpy( kl_bg_seq,
                        seq + v->pos.pos + bg_right,
                        (bg_len+k_-1)*sizeof(char) );
            }
        }

        if( kl_bg_seq[0] ) {
            K = 0;
            for( j = 0; j < bg_len+((pos)k_-1); j++ ) {
                K = ((K<<2) | nt2num(kl_bg_seq[j])) & kmer_mask_;
                if( j >= (pos)k_-1 ) kl_bg[K] += 1;
            }
        }



        /* sample foreground */
        if( v->pos.strand ) {
            if( v->pos.pos < R ||  v->pos.pos < R_ ) continue;
            memcpy( local_seq, seq + v->pos.pos - R_, (L_+R_+k_)*sizeof(char) );
            seqrc( local_seq, L_+R_+k_ );

            memcpy( bias_seq, seq + v->pos.pos - R, (L+R+k)*sizeof(char) );
            seqrc( bias_seq,  L+R+k );
        }
        else {
            if( v->pos.pos < L_+(pos)(k_-1) ) continue;
            memcpy( local_seq, seq + (v->pos.pos-L_-(k_-1)), (L_+R_+k_)*sizeof(char) );
            memcpy( bias_seq,  seq + (v->pos.pos-L-(k-1)), (L+R+k)*sizeof(char) );
        }

        

        /* sample kmers */
        K = 0;
        for( j = 0; j < L_+R_+(pos)k_; j++ ) {
            K = ((K<<2) | nt2num(local_seq[j])) & kmer_mask_;
            if( j >= (pos)(k_-1) ) {
                freq[ (j-(pos)(k_-1))*four_to_k_ + K ] += 1;
            }
        }

        /* sample adjusted kmer frequencies */
        if( adjusted ) {
            double w = 1.0;
            K = 0;
            for( j = 0; j < L+R+(pos)k; j++ ) {
                K = ((K<<2) | nt2num(bias_seq[j])) & kmer_mask;
                if( j >= (pos)(k-1) ) {
                    w *= ws[ (j-(pos)(k-1))*four_to_k + K ];
                }
            }

            K = 0;
            for( j = 0; j < L_+R_+(pos)k_; j++ ) {
                K = ((K<<2) | nt2num(local_seq[j])) & kmer_mask_;
                if( j >= (pos)(k_-1) ) {
                    freq_adjusted[ (j-(pos)(k_-1))*four_to_k_ + K ] += 1.0/w;
                }
            }
        }
    }

    /* normalize */
    double z;
    for( j = 0; j < L_+R_+1; j++ ) {
        z = 0.0;
        for( K = 0; K < four_to_k_; K++ ) {
            z += freq[ j*four_to_k_+K ];
        }

        for( K = 0; K < four_to_k_; K++ ) {
            freq[ j*four_to_k_+K ] /= z;
        }
    }

    if( adjusted ) {
        for( j = 0; j < L_+R_+1; j++ ) {
            z = 0.0;
            for( K = 0; K < four_to_k_; K++ ) {
                z += freq_adjusted[ j*four_to_k_+K ];
            }

            for( K = 0; K < four_to_k_; K++ ) {
                freq_adjusted[ j*four_to_k_+K ] /= z;
            }
        }
    }

    z = 0.0;
    for( K = 0; K < four_to_k_; K++ ) z += kl_bg[K];
    for( K = 0; K < four_to_k_; K++ ) kl_bg[K] /= z;
    

    /* compute kl divergance for each position */
    double* kl = (double*)safe_malloc( (L_+R_+1)*sizeof(double) );
    double* kl_adjusted = NULL;
    if( adjusted ) kl_adjusted = (double*)safe_malloc( (L_+R_+1)*sizeof(double) );

    for( j = 0; j < L_+R_+1; j++ ) {
        kl[j] = kl_div( kl_bg, freq+(j*four_to_k_), four_to_k_ );
        if( adjusted ) {
            kl_adjusted[j] = kl_div( kl_bg, freq_adjusted+(j*four_to_k_), four_to_k_ );
        }
    }




    /* print */
    char* ntstr = (char*)safe_malloc( (k_+1)*sizeof(char) );
    for( j = 0; j < L_+R_+1; j++ ) {
        fprintf( fout, "%ld\tNA\t%e\tunadjusted\n", j - L_, kl[j] );
        if( adjusted ) {
            fprintf( fout, "%ld\tNA\t%e\tadjusted\n", j - L_, kl_adjusted[j] );
        }

        for( K = 0; K < four_to_k_; K++ ) {
            num2nt( K, ntstr, k_ );
            fprintf( fout, "%ld\t%s\t%e\tunadjusted\n", j - L_, ntstr, freq[j*four_to_k_+K] );
            if( adjusted ) {
                fprintf( fout, "%ld\t%s\t%e\tadjusted\n", j - L_,
                        ntstr, freq_adjusted[j*four_to_k_+K] );
            }
        }
    }

    free(ntstr);


    free(S);
    free(kl_bg);
    free(kl_bg_seq);
    free(local_seq);
    free(bias_seq);
    samclose(reads_f);

    free(kl);
    free(kl_adjusted);
    free(freq);
    free(freq_adjusted);
}





void sequencing_bias::sample_foreground( char* seq, size_t seqlen,
                                         struct hashed_value* v )
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
        if( j >= (pos)(k-1) ) fg[ (j-(pos)(k-1))*four_to_k + K ] += 1;
    }
}


void sequencing_bias::sample_background( char* seq, size_t seqlen,
                                         struct hashed_value* v )
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
            if( j >= k-1 ) bgs[four_to_k*(j-(k-1))+K] += 1;
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
            if( j >= k-1 ) bgs[four_to_k*(j+bg_len-(k-1))+K] += 1;
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

