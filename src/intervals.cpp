
#include "intervals.hpp"
#include "logger.h"

#include <cstring>
#include <algorithm>

using namespace std;


bool subinterval_ordering::operator()( const subinterval& A, const subinterval& B )
{
    if( A.start < B.start )       return true;
    else if( A.start == B.start ) return A.end < B.end;
    else                          return false;
}


subinterval::subinterval() : start(-1), end(-1), score(0.0), count(0), rate(0) {}

subinterval::subinterval( pos start, pos end, double score,
                          rcount count, unsigned int zeros, double rate )
    : start(start)
    , end(end)
    , score(score)
    , count(count)
    , zeros(zeros)
    , rate(rate)
{}

subinterval::subinterval( const subinterval& y )
    : start(y.start)
    , end(y.end)
    , score(y.score)
    , count(y.count)
    , zeros(y.zeros)
    , rate(y.rate) {}

void subinterval::operator=( const subinterval& y ) { 
    start  = y.start;
    end    = y.end;
    score  = y.score;
    count  = y.count;
    zeros  = y.zeros;
    rate   = y.rate;
}

bool subinterval::operator< ( const subinterval& y ) const
{
    return score < y.score;
}

bool subinterval::operator> ( const subinterval& y ) const
{
    return score > y.score;
}


void subinterval::set( pos start, pos end, double score,
                       rcount count, unsigned int zeros, double rate ) {
    this->start  = start;
    this->end    = end;
    this->score  = score;
    this->count  = count;
    this->zeros  = zeros;
    this->rate   = rate;
}

pos subinterval::length() const
{
    return end - start + 1;
}



void subinterval_bound::check() const
{
    if( I_min > I_max || J_min > J_max || rate < 0.0 )  {
        failf( "Fatal error: nonsense interval bounds encountered\n"
               "bounds = (%ld,%ld,%ld,%ld) rate = %0.4e\n"
               "Report/investigate.\n",
               I_min, I_max, J_min, J_max, rate );
    }
}



subinterval_bound::subinterval_bound()
    : I_min(0)
    , I_max(0)
    , J_min(0)
    , J_max(0)
    , rate(1.0)
    , count(0)
    , zeros(0)
    , score(0.0)
{
}


subinterval_bound::subinterval_bound( pos min, pos max,
                                      double       rate,
                                      rcount       count,
                                      unsigned int zeros,
                                      double       score )
    : I_min(min), I_max(max)
    , J_min(min), J_max(max)
    , rate(rate), count(count), zeros(zeros), score(score)
{
    check();
}

subinterval_bound::subinterval_bound( pos I_min, pos I_max,
                                      pos J_min, pos J_max,
                                      double       rate,
                                      rcount       count,
                                      unsigned int zeros,
                                      double       score )
    : I_min(I_min), I_max(I_max)
    , J_min(J_min), J_max(J_max)
    , rate(rate), count(count), zeros(zeros), score(score)
{
    check();
}

subinterval_bound::subinterval_bound( const subinterval_bound& y )
    : I_min(y.I_min), I_max(y.I_max)
    , J_min(y.J_min), J_max(y.J_max)
    , rate(y.rate), count(y.count), zeros(y.zeros), score(y.score)
{
    check();
}


void subinterval_bound::operator=( const subinterval_bound& y )
{
    this->I_min = y.I_min;
    this->I_max = y.I_max;
    this->J_min = y.J_min;
    this->J_max = y.J_max;
    this->rate  = y.rate;
    this->count = y.count;
    this->zeros = y.zeros;
    this->score  = y.score;
    check();
}


void subinterval_bound::set( pos min, pos max,
          double       rate,
          rcount       count,
          unsigned int zeros,
          double       score )
{
    this->I_min  = min;
    this->I_max  = max;
    this->J_min  = min;
    this->J_max  = max;
    this->rate   = rate;
    this->count  = count;
    this->zeros  = zeros;
    this->score  = score;
    check();
}


void subinterval_bound::set( pos I_min, pos I_max, pos J_min, pos J_max,
          double       rate,
          rcount       count,
          unsigned int zeros,
          double       score )
{
    this->I_min  = I_min;
    this->I_max  = I_max;
    this->J_min  = J_min;
    this->J_max  = J_max;
    this->rate   = rate;
    this->count  = count;
    this->zeros  = zeros;
    this->score  = score;
    check();
}


bool subinterval_bound::disjoint_bounds() const
{
    return I_max < J_min || J_max < I_min;
}


bool subinterval_bound::equal_bounds() const
{
    return I_min == J_min && I_max == J_max;
}

/* Count the number of subintervals contained within [i,j], given a min and max
 * duration of d_min_, d_max_, respectively. */
double subinterval_bound::subinterval_count_equal( pos i, pos j, pos d_min_, pos d_max_ ) const
{
    /* note, this is done with doubles to avoid wrap-around on overflow */
    double d_min = (double)d_min_;
    double d_max = (double)d_max_;

    if( i > j ) return 0;
    d_min_ = min<double>( d_min, 1 );
    double n = (double)(j - i + 1);

    if( n < d_min )     return 0;
    if( d_max < d_min ) return 0;

    double m = min( n, d_max );

    return (1+n)*(1+m-d_min) - (m*(m+1) - (d_min*(d_min-1))) / 2;
}


/* Count the number of subintervals contined bounded by this object, given a
 * min and max duration of d_min and d_max, respectively. */
double subinterval_bound::subinterval_count( pos d_min, pos d_max ) const
{
    if( disjoint_bounds() ) {
        double T, L, R, C;
        T = subinterval_count_equal( I_min, J_max, d_min, d_max );
        L = subinterval_count_equal( I_min, J_min-1, d_min, d_max );
        R = subinterval_count_equal( I_max+1, J_max, d_min, d_max );
        C = subinterval_count_equal( I_max+1, J_min-1, d_min, d_max );

        return (size_t)(T + C - L - R);
    }
    else if( equal_bounds() ) {
        return subinterval_count_equal( I_min, J_max, d_min, d_max );
    }
    else return 0;
}


pos subinterval_bound::subinterval_bound::max_length() const
{
    return J_max - I_min + 1;
}

pos subinterval_bound::subinterval_bound::min_length() const
{

    return max( J_min - I_max + 1, (pos)0 );
}


bool subinterval_bound_priority::operator()( subinterval_bound*& a,
                                             subinterval_bound*& b )
{
    return a->score > b->score;
}


subinterval_bound_pqueue::subinterval_bound_pqueue()
{
}

subinterval_bound_pqueue::~subinterval_bound_pqueue()
{
    iterator i;
    for( i = begin(); i != end(); i++ ) {
        delete *i;
    }
}


void subinterval_bound_pqueue::push( subinterval_bound* x )
{
    push_back( x );
    push_heap( begin(), end(), subinterval_bound_priority() );
}


subinterval_bound* subinterval_bound_pqueue::pop()
{
    pop_heap( begin(), end(), subinterval_bound_priority() );
    subinterval_bound* B = back();
    pop_back();
    return B;
}



interval::interval()
    : seqname(NULL)
    , start(-1)
    , end(-1)
    , strand(-1)
{
}

interval::interval( const char* seqname, pos start, pos end, int strand )
    : seqname(strdup(seqname))
    , start(start)
    , end(end)
    , strand(strand)
{
}

interval::interval( const subinterval& i, const char* seqname, pos start, int strand )
    : seqname(strdup(seqname))
    , start(start+i.start)
    , end(start+i.end)
    , strand(strand)
    , score(i.score)
{
}

interval::interval( const interval& i )
    : seqname(strdup(i.seqname))
    , start(i.start)
    , end(i.end)
    , strand(i.strand)
    , score(i.score)
{
}

interval::~interval()
{
    free(seqname);
}

void interval::set( const char* seqname, pos start, pos end, int strand )
{
    free(this->seqname);
    this->seqname = strdup(seqname);
    this->start  = start;
    this->end    = end;
    this->strand = strand;
}

void interval::set( const subinterval& i, const char* seqname, pos start, int strand )
{
    free(this->seqname);
    this->seqname  = strdup(seqname);
    this->start    = start+i.start;
    this->end      = start+i.end;
    this->strand   = strand;
    this->score    = i.score;
}

pos interval::length() const
{
    return end - start + 1;
}



void interval_stack_push( interval_stack* is, char* seqname, pos start, pos end, int strand )
{
    is->push_back( interval( seqname, start, end, strand ) );
}


interval_stack::interval_stack()
{
}

interval_stack::interval_stack( const subinterval_stack& sis,
                                const char* chrom,pos start, int strand )
{
    subinterval_stack::const_iterator i;
    for( i = sis.begin(); i != sis.end(); i++ ) {
        push_back( interval( *i, chrom, start, strand ) );
    }
}

interval_stack::~interval_stack()
{
}


subinterval_stack::subinterval_stack()
{

}

/*
void subinterval_stack::print( const char* seqname, char strand, pos offset,
                            const char* gene_name )
{
    if( gene_name == NULL ) gene_name = "";

    subinterval_stack::iterator i;
    for( i = this->begin(); i != this->end(); i++ ) {
        printf("%s\t%d\t%d\t%s\t%0.8e\t%c\n",
               seqname, i->start+offset-1, i->end+offset-1,
               gene_name, i->score, strand );
        fflush(stdout);
    }
}
*/





