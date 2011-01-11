

#include "context.hpp"
#include <gsl/gsl_math.h>
#include <algorithm>

using namespace std;



int bam_fetch_callback( const bam1_t* b, void* a_ )
{
    context* a = (context*)a_;
    if( (a->strand == -1 || bam1_strand(b) == a->strand) && b->core.pos >= a->start ) {
        a->xs[bam1_strand(b)][b->core.pos - a->start]++;
    }
    return 0;
}


context::context()
    : start(-1), end(-1), seqname(NULL), strand(-1)
{
    xs[0] = NULL; xs[1] = NULL;
    cs[0] = NULL; cs[1] = NULL;
    rs[0] = NULL; rs[1] = NULL;
}


void context::clear()
{
    start = -1;
    end   = -1;
    free(seqname);
    seqname = NULL;
    strand  = -1;

    delete[] xs[0]; delete[] xs[1];
    xs[0] = xs[1] = NULL;

    delete[] cs[0]; delete[] cs[1];
    cs[0] = cs[1] = NULL;

    delete[] xs[0]; delete[] xs[1];
    rs[0] = rs[1] = NULL;
}

void context::set( dataset* dataset, const interval& i )
{
    this->set( dataset, i.seqname, i.start, i.end, i.strand );
}

void context::set( dataset* ds, const char* seqname,
                   pos start, pos end, int strand )
{
    /* determine wether to examine both strand or just one */
    int s, s0, s1;
    if( strand >= 0 ) {
        s0 = strand;
        s1 = strand;
    }
    else {
        s0 = 0;
        s1 = 1;
    }

    clear();
    this->seqname = strdup(seqname);
    this->start   = start;
    this->end     = end;
    this->strand  = strand;


    /* count reads */
    bam_iter_t it;
    bam1_t* read;
    pos i;
    size_t j;
    pos read_end;
    uint32_t* cigar;
    uint8_t   cigar_op;
    pos       cigar_opend;


    for( s = s0; s <= s1; s++ ) {

        /* so that bam_fetch_callback counts strand seperately */
        this->strand = s;

        char* region;
        int bam_ref_id, bam_start, bam_end, region_error;
        int c = asprintf( &region, "%s:%ld-%ld", seqname, start, end );
        if( c <= 0 ) continue;
        region_error = bam_parse_region( ds->reads_f->header, region,
                                         &bam_ref_id, &bam_start, &bam_end );
        free(region);
        if( region_error != 0 || bam_ref_id < 0 ) {
            xs[s] = NULL;
            cs[s] = NULL;
            continue;
        }

        xs[s] = new rcount[ length() ];
        cs[s] = new rcount[ length() ];
        memset( xs[s], 0, length()*sizeof(rcount) );
        memset( cs[s], 0, length()*sizeof(rcount) );

        it   = bam_iter_query( ds->reads_index, bam_ref_id, bam_start, bam_end );
        read = bam_init1();

        while( bam_iter_read( ds->reads_f->x.bam, it, read ) >= 0 ) {
            if( bam1_strand(read) != s ) continue;

            cigar = bam1_cigar(read);

            i = read->core.pos;
            for( i = read->core.pos, j = 0; i <= end, j < read->core.n_cigar; i++, j++ ) {
                cigar_op    = cigar[j] &  BAM_CIGAR_MASK;
                cigar_opend = i + (pos)(cigar[j] >> BAM_CIGAR_SHIFT);

                if( cigar_op == BAM_CMATCH && cigar_opend >= start ) {
                    for( ; i < cigar_opend, i <= end; i++ ) {
                        if( i >= start ) cs[s][i - start]++;
                    }
                }
                else i = cigar_opend;
            }

            /* set xs */
            if( s == 0 && start <= read->core.pos && read->core.pos <= end ) {
                xs[s][read->core.pos - start]++;
            }
            else {
                read_end = bam_calend( &read->core, cigar ) - 1;
                if( s == 1 && start <= read_end && read_end <= end ) {
                    xs[s][read_end - start]++;
                }
            }
        }

        bam_destroy1(read);
    }




        /* OLD OLD OLD */

/*
 *        while( bam_iter_read( ds->reads_f->x.bam, it, read ) >= 0 ) {
 *            if( strand != -1 && bam1_strand(read) != strand ) continue;
 *
 *            [> positive strand <]
 *            if( bam1_strand(read) == 0 ) {
 *                i = read->core.pos;
 *                read_end = bam_calend( &read->core, bam1_cigar(read) );
 *                if( start <= i && i <= end ) xs[0][i-start]++;
 *                for( i = max(start,i); i < read_end && i <= end; i++ ) {
 *                    cs[0][i-start]++;
 *                }
 *            }
 *            [> negative strand <]
 *            else {
 *                i = bam_calend( &read->core, bam1_cigar(read) ) - 1;
 *                read_end = read->core.pos;
 *                if( start <= i && i <- end ) xs[1][i-start]++;
 *                for( i = min(end,i); i >= read_end && i >= start; i-- ) {
 *                    cs[1][i-start]++;
 *                }
 *            }
 *        }
 *
 *        bam_destroy1(read);
 *    }
 */


    /* get sequencing bias */
    if( ds->bias ) {
        for( s = s0; s <= s1; s++ ) {
            rs[s] = ds->bias->get_bias( seqname, start, end, s );
        }
    }
    else rs[0] = rs[1] = NULL;

    this->strand = strand;
}


void context::adjust_interval_by_coverage( interval& I ) const
{
    rcount limit;
    if( I.strand == 0 || I.strand == -1 ) {
        if( cs[0] != NULL && xs[0] != NULL && start <= I.end && I.end <= end ) {
            limit = max( xs[0][I.end - start], cs[0][I.end - start] );
            while( I.end <= end && cs[0][I.end - start] >= limit ) I.end++;
        }
    }

    if( I.strand == 1 || I.strand == -1 ) {
        if( cs[1] != NULL && xs[1] != NULL && start <= I.start && I.start <= end ) {
            limit = max( xs[1][I.start - start], cs[0][I.end - start] );
            while( I.start >= start && cs[1][I.start - start] >= limit ) I.start--;
        }
    }
}


/*
 *void context::adjust_subinterval_by_coverage( subinterval& I ) const
 *{
 *    rcount limit;
 *    if( I.strand == 0 || I.strand == -1 ) {
 *        if( cs[0] != NULL && xs[0] != NULL && 0 <= I.end && I.end < length() ) {
 *            limit = xs[0][I.end];
 *            while( I.end < length() && cs[0][I.end] >= limit ) I.end++;
 *        }
 *    }
 *
 *    if( I.strand == 1 || I.strand == -1 ) {
 *        if( cs[1] != NULL && xs[1] != NULL && 0 <= I.start && I.start < length() ) {
 *            limit = xs[1][I.start];
 *            while( I.start >= 0 && cs[1][I.start] >= limit ) I.start--;
 *        }
 *    }
 *}
 */

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
        xs[0][i] = dist.rand( 1.0 );
        /* don't bother generating coverage */
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


double context::min_rate( const subinterval_bound& B, pos d_min ) const
{
    return min_rate( B, d_min, strand );
}

double context::min_rate( const subinterval_bound& B, pos d_min, int strand ) const
{
    if( B.disjoint_bounds() ) {
        return rate( B.I_max, B.J_min, strand );
        /* TODO: incorrect when B.J_min - B.I_max + 1 < d_min */
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
            fprintf( out_f, "unadjusted\t%ld\t%ld\n", start+i, xs[s][i] );
            fprintf( out_f, "adjusted\t%ld\t%.5e\n", start+i, rs[s][i]*(double)xs[s][i] );
            fprintf( out_f, "bias\t%ld\t%.5e\n",     start+i, rs[s][i] );
        }
    }

}


