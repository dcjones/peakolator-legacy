

#include "context.hpp"
#include "logger.h"
#include <gsl/gsl_math.h>
#include <algorithm>

using namespace std;


context::context()
    : start(-1), end(-1), seqname(NULL), strand(-1), ds(NULL)
{
    xs[0] = NULL; xs[1] = NULL;
    rs[0] = NULL; rs[1] = NULL;
    cs = NULL;
}

context::~context()
{
    clear();
}


void context::clear()
{
    start = -1;
    end   = -1;
    free(seqname);
    seqname = NULL;
    strand  = -1;
    ds      = NULL;

    delete[] xs[0]; delete[] xs[1];
    xs[0] = xs[1] = NULL;

    delete[] xs[0]; delete[] xs[1];
    rs[0] = rs[1] = NULL;

    delete[] cs;
    cs = NULL;
}

void context::set( dataset* dataset, const interval& i )
{
    this->set( dataset, i.seqname, i.start, i.end, i.strand );
}

void context::set( dataset* ds, const char* seqname,
                   pos start, pos end, int strand )
{
    clear();
    this->seqname = strdup(seqname);
    this->start   = start;
    this->end     = end;
    this->strand  = strand;
    this->ds      = ds;


    /* count reads */
    bam_iter_t it;
    bam1_t* read;
    pos i;

    int tid = bam_get_tid( ds->reads_f->header, seqname );
    if( tid < 0 ) {
        failf( "Sequence '%s' is not present in BAM header.\n", seqname );
    }


    cs = new rcount[ length() ];
    memset( cs, 0, length()*sizeof(rcount) );

    if( strand == 0 || strand == -1 ) {
        xs[0] = new rcount[ length() ];
        memset( xs[0], 0, length()*sizeof(rcount) );
    }

    if( strand == 1 || strand == -1 ) {
        xs[1] = new rcount[ length() ];
        memset( xs[1], 0, length()*sizeof(rcount) );
    }

    it   = bam_iter_query( ds->reads_index, tid, start, end );
    read = bam_init1();

    while( bam_iter_read( ds->reads_f->x.bam, it, read ) >= 0 ) {
        if( strand != -1 && bam1_strand(read) != strand ) continue;

        if( bam1_strand(read) == 0 ) {
            i = read->core.pos;
            if( start <= i && i <= end ) {
                xs[0][i - start]++;
            }
        }
        else {
            i = bam_calend( &read->core, bam1_cigar(read) ) - 1;
            if( start <= i && i <= end ) {
                xs[1][i - start]++;
            }
        }
    }

    bam_destroy1(read);
    bam_iter_destroy(it);


    /* get sequencing bias */
    if( ds->bias ) {
        if( strand == -1 || strand == 0 ) {
            rs[0] = ds->bias->get_bias( seqname, start, end, 0 );
        }

        if( strand == -1 || strand == 1 ) {
            rs[1] = ds->bias->get_bias( seqname, start, end, 1 );
        }
    }
    else rs[0] = rs[1] = NULL;

    this->strand = strand;
}


void context::get_coverage( pos u, pos v, int s ) const
{
    if( ds == NULL || cs == NULL ) return;

    bam1_t* read = bam_init1();
    bam_iter_t it;
    uint32_t* cigar;
    uint8_t  c_op;
    pos      c_end;
    pos i;
    size_t j;

    int tid = bam_get_tid( ds->reads_f->header, seqname );
    if( tid < 0 ) {
        failf( "Sequence '%s' is not present in BAM header.\n", seqname );
    }

    memset( cs, 0, length() * sizeof(rcount) );

    it = bam_iter_query( ds->reads_index, tid, u, v );

    while( bam_iter_read( ds->reads_f->x.bam, it, read ) >= 0 ) {
        if( bam1_strand(read) != s ) continue;

        cigar = bam1_cigar(read);
        i = read->core.pos;

        for( j = 0; j < read->core.n_cigar && i <= end; j++ ) {

            c_op  = cigar[j] &  BAM_CIGAR_MASK;
            c_end = i + (pos)(cigar[j] >> BAM_CIGAR_SHIFT);

            if( c_op != BAM_CMATCH || c_end <= start ) {
                i = c_end;
            }
            else {
                for( i = max( i, start ); i < min( c_end, end + 1 ); i++ ) {
                    cs[i - start]++;
                }
            }
        }
    }

    bam_destroy1(read);
}


void context::adjust_interval_by_coverage( interval& I ) const
{
    rcount limit;

    /* extend end */
    if( strand == -1 || strand == 0 ) {
        limit = xs[0][ I.end - start ];
        if( limit > 0 ) {
            get_coverage( I.end, I.end+1, 0 );
            while( I.end <= end && cs[I.end - start] >= limit ) I.end++;
        }
    }

    /* extend start */
    if( strand == -1 || strand == 1 ) {
        limit = xs[1][ I.start - start ];
        if( limit > 0 ) {
            get_coverage( I.start, I.start+1, 1 );
            while( I.start >= start && cs[I.start - start] >= limit ) I.start--;
        }
    }
}


void context::set_noise( nulldist& dist, pos len )
{
    if( strand != 0 || length() != len ){
        clear();
        this->seqname = strdup("NOISE");
        this->start    = 0;
        this->end      = len-1;
        this->strand   = 0;
        this->rs[0]    = new double[ len ];
        this->xs[0]    = new rcount[ len ];
    }

    pos i;
    for( i = 0; i < len; i++ ) {
        rs[0][i] = 1.0;
        xs[0][i] = dist.rand();
    }
}


rcount context::count() const
{
    return count( 0, length()-1 );
}

rcount context::count( pos i ) const
{
    if( strand == -1 ) {
        return (xs[0] ? xs[0][i] : 0) + (xs[1] ? xs[1][i] : 0);
    }
    else {
        return xs[strand] ? xs[strand][i] : 0;
    }
}


rcount context::count( pos i, pos j ) const
{
    return count( i, j, strand );
}

rcount context::count( pos i, pos j, int strand ) const
{
    rcount total = 0;
    if( strand == -1 ) {
        while( i <= j ) {
            total += xs[0] ? xs[0][i] : 0;
            total += xs[1] ? xs[1][i] : 0;
            i++;
        }
    }
    else {
        while( i <= j ) {
            total += xs[strand] ? xs[strand][i] : 0;
            i++;
        }
    }

    return total;
}

double context::rate() const
{
    return rate( 0, length()-1 );
}

double context::rate( pos i ) const
{
    if( strand == -1 ) {
        return (rs[0] ? rs[0][i] : 1.0) + (rs[1] ? rs[1][i] : 1.0);
    }
    else {
        return rs[strand] ? rs[strand][i] : 1.0;
    }
}


double context::rate( pos i, pos j ) const
{
    return rate( i, j, strand );
}

double context::rate( pos i, pos j, int strand ) const
{
    double total = 0.0;
    if( strand == -1 ) {
        while( i <= j ) {
            total += rs[0] ? rs[0][i] : 1.0;
            total += rs[1] ? rs[1][i] : 1.0;
            i++;
        }
    }
    else while( i <= j ) {
        total += rs[strand] ? rs[strand][i] : 1.0;
        i++;
    }

    return total;
}

void context::count_zeros( rcount* c, unsigned int* z )
{
    count_zeros( c, z, 0, length()-1 );
}

void context::count_zeros( rcount* c, unsigned int* z, pos i )
{
    rcount c_i;
    *c = 0;
    *z = 0;

    if( strand == -1 ) {
        c_i = xs[0] ? xs[0][i] : 0;
        if( c_i == 0 ) *z += 1;
        else           *c += c_i;

        c_i = xs[1] ? xs[1][i] : 0;
        if( c_i == 0 ) *z += 1;
        else           *c += c_i;
    }
    else {
        c_i = xs[strand] ? xs[strand][i] : 0;
        if( c_i == 0 ) *z += 1;
        else           *c += c_i;
    }
}

void context::count_zeros( rcount* c, unsigned int* z,
                           pos i, pos j )
{
    count_zeros( c, z, i, j, strand );
}

void context::count_zeros( rcount* c, unsigned int* z,
                           pos i, pos j, int strand )
{
    rcount c_i;
    *c = 0;
    *z = 0;

    if( strand == -1 ) {
        while( i <= j ) {
            c_i = xs[0] ? xs[0][i] : 0;
            if( c_i == 0 ) *z += 1;
            else           *c += c_i;

            c_i = xs[1] ? xs[1][i] : 0;
            if( c_i == 0 ) *z += 1;
            else           *c += c_i;
            i++;
        }
    }
    else {
        while( i <= j ) {
            c_i = xs[strand] ? xs[strand][i] : 0;
            if( c_i == 0 ) *z += 1;
            else           *c += c_i;
            i++;
        }
    }
}


unsigned int context::min_zeros( const subinterval_bound& B, pos d_min ) const
{
    unsigned int mind = (unsigned int)min_duration( B, d_min );
    unsigned int maxd = (unsigned int)B.max_length();

    if( B.zeros <= maxd - mind ) return 0;
    else                         return B.zeros - (maxd - mind);
}

pos context::min_duration( const subinterval_bound& B, pos d_min ) const
{
    if( B.disjoint_bounds() ) {
        return max( B.J_min - B.I_max + 1, d_min );
    }
    else {
        return min( B.J_max - B.I_min + 1, d_min );
    }
}


double context::min_rate( const subinterval_bound& B, pos d_min ) const
{
    return min_rate( B, d_min, strand );
}

double context::min_rate( const subinterval_bound& B, pos d_min, int strand ) const
{
    if( B.disjoint_bounds() ) {
        return rate( B.I_max, B.J_min, strand );
        /* NOTE: this is incorrect when B.J_min - B.I_max + 1 < d_min, but
         * it will only underestimate min rate, thus the bound will still be
         * correct, just not as tight.
         */
    }

    else if( B.equal_bounds() ) {
        double rate_min = GSL_POSINF;
        double r = 0.0;
        pos i;

        for( i = B.I_min; i <= B.I_max - d_min + 1; i++ ) {
            r += rate(i);

            if( i - B.I_min + 1 >= d_min ) {
                r -= rate(i-1);
                if( r < rate_min ) rate_min = r;
            }
        }

        /* this, of course, should never happen */
        if( gsl_isinf(rate_min) || gsl_isnan(rate_min) ) rate_min = (double)d_min;

        return rate_min;
    }

    else return 0.0;
}


pos context::length() const
{
    return end - start + 1;
}


void context::print_adjusted_unadjusted_bias( FILE* out_f )
{
    int s, s0, s1;
    if( strand == -1 ) { s0 = 0; s1 = 1; }
    else s0 = s1 = strand;
    int i;
    for( s = s0; s <= s1; s++ ) {
        for( i = 0; i < length(); i++ ) {
            fprintf( out_f, "unadjusted\t%ld\t%u\n",  start+i, xs[s][i] );
            fprintf( out_f, "adjusted\t%ld\t%.5e\n" , start+i, rs[s][i]*(double)xs[s][i] );
            fprintf( out_f, "bias\t%ld\t%.5e\n",      start+i, rs[s][i] );
        }
    }

}


