

#include "context.hpp"
#include <gsl/gsl_math.h>



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
    rs[0] = NULL; rs[1] = NULL;
}


void context::clear()
{
    start = -1;
    end   = -1;
    free(seqname);
    seqname = NULL;
    strand  = -1;

    free(xs[0]); free(xs[1]);
    xs[0] = xs[1] = NULL;

    free(rs[0]); free(rs[1]);
    rs[0] = rs[1] = NULL;
}

void context::set( dataset* dataset, const interval& i )
{
    this->set( dataset, i.seqname, i.start, i.end, i.strand );
}

void context::set( dataset* dataset, const char* seqname,
                              pos start, pos end, int strand )
{
    /* determine wether to examine both strand or just one */
    int s, s0, s1;
    if( strand >= 0 ) s0 = s1 = strand;
    else              s0 = 0; s1 = 1;

    clear();
    this->seqname = strdup(seqname);
    this->start   = start;
    this->end     = end;
    this->strand  = strand;


    /* count reads */
    for( s = s0; s <= s1; s++ ) {

        /* so that bam_fetch_callback counts strand seperately */
        this->strand = s;

        char* region;
        int bam_ref_id, bam_start, bam_end, region_error;
        int c = asprintf( &region, "%s:%ld-%ld", seqname, start, end );
        if( c <= 0 ) continue;
        region_error = bam_parse_region( dataset->reads_f->header, region,
                                         &bam_ref_id, &bam_start, &bam_end );
        free(region);
        if( region_error != 0 || bam_ref_id < 0 ) {
            xs[s] = NULL;
            continue;
        }

        xs[s] = (rcount*)safe_malloc( length()*sizeof(rcount) );
        memset( xs[s], 0, length()*sizeof(rcount) );
        bam_fetch( dataset->reads_f->x.bam, dataset->reads_index,
                   bam_ref_id, bam_start, bam_end,
                   (void*)this, bam_fetch_callback );
    }


    /* get sequencing bias */
    if( dataset->bias ) {
        for( s = s0; s <= s1; s++ ) {
            rs[s] = dataset->bias->get_bias( seqname, start, end, s );
        }
    }
    else rs[0] = rs[1] = NULL;

    this->strand = strand;
}

void context::set_noise( nulldist& dist, pos len )
{
    if( strand != 0 || length() != len ){
        clear();
        this->seqname = strdup("NOISE");
        this->start    = 0;
        this->end      = len-1;
        this->strand   = 0;
        this->rs[0]    = (double*)safe_malloc( len*sizeof(double) );
        this->xs[0]    = (rcount*)safe_malloc( len*sizeof(double) );
    }

    pos i;
    for( i = 0; i < len; i++ ) {
        rs[0][i] = 1.0;
        xs[0][i] = dist.rand( 1.0 );
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
        total += rs[strand] ? rs[strand][i++] : 1.0;
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
    }

    else if( B.equal_bounds() ) {
        double rate_min = GSL_POSINF;
        double r = 0.0;
        pos i;

        for( i = B.I_min; i <= B.J_min; i++ ) {
            r += rate(i);;

            if( i - start + 1 >= d_min ) {
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


