
#include "peakolator_interval.hpp"
#include "peakolator_model.hpp"
#include "logger.h"

#include <algorithm>

using namespace std;


bool subinterval_ordering::operator()( const subinterval& A, const subinterval& B )
{
    if( A.start < B.start ) return true;
    else if( A.start == B.start ) return A.end < B.end;
    else                          return false;
}


subinterval::subinterval() : start(-1), end(-1), pval(1.0), count(0), rate(0) {}

subinterval::subinterval( pos start, pos end, const mpfr_class& pval, rcount count, double rate )
    : start(start), end(end), pval(pval), count(count), rate(rate) {}

subinterval::subinterval( const subinterval& y )
    : start(y.start), end(y.end), pval(y.pval), count(y.count), rate(y.rate) {}

void subinterval::operator=( const subinterval& y ) { 
    start = y.start;
    end   = y.end;
    pval  = y.pval;
    count = y.count;
    rate  = y.rate;
}

bool subinterval::operator< ( const subinterval& y ) const { return pval < y.pval; }
bool subinterval::operator> ( const subinterval& y ) const { return pval > y.pval; }


//bool subinterval::operator< ( const subinterval& y ) const
//{ 
    //return ((double)count / rate) < ((double)y.count / y.rate);
//}

//bool subinterval::operator> ( const subinterval& y ) const
//{ 
    //return ((double)count / rate) > ((double)y.count / y.rate);
//}

void subinterval::set( pos start, pos end, const mpfr_class& pval, rcount count, double rate ) {
    this->start = start;
    this->end   = end;
    this->pval  = pval;
    this->count = count;
    this->rate  = rate;
}

pos subinterval::length() const
{
    return end - start + 1;
}



void subinterval_bound::check() const
{
    if( I_min > I_max || J_min > J_max || rate < 0.0 )  {
        log_puts( LOG_ERROR,
                "Fatal error: nonsense interval bounds encountered. Report/investigate." );
        exit(1);
    }
}



subinterval_bound::subinterval_bound()
    : I_min(0)
    , I_max(0)
    , J_min(0)
    , J_max(0)
    , rate(1.0)
    , count(0)
    , pval(1.0)
{
}


subinterval_bound::subinterval_bound( pos min, pos max,
                                      double rate,
                                      rcount count,
                                      const mpfr_class& pval )
    : I_min(min), I_max(max)
    , J_min(min), J_max(max)
    , rate(rate), count(count)
    , pval(pval)
{
    check();
}

subinterval_bound::subinterval_bound( pos I_min, pos I_max,
                                      pos J_min, pos J_max,
                                      double rate,
                                      rcount count,
                                      const mpfr_class& pval )
    : I_min(I_min), I_max(I_max)
    , J_min(J_min), J_max(J_max)
    , rate(rate), count(count)
    , pval(pval)
{
    check();
}

subinterval_bound::subinterval_bound( const subinterval_bound& y )
    : I_min(y.I_min), I_max(y.I_max)
    , J_min(y.J_min), J_max(y.J_max)
    , rate(y.rate), count(y.count)
    , pval(y.pval)
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
    this->pval  = y.pval;
    check();
}


void subinterval_bound::set( pos min, pos max,
          double rate,
          rcount count,
          const mpfr_class& pval )
{
    this->I_min = min;
    this->I_max = max;
    this->J_min = min;
    this->J_max = max;
    this->rate  = rate;
    this->count = count;
    this->pval  = pval;
    check();
}


void subinterval_bound::set( pos I_min, pos I_max, pos J_min, pos J_max,
          double rate,
          rcount count,
          const mpfr_class& pval )
{
    this->I_min = I_min;
    this->I_max = I_max;
    this->J_min = J_min;
    this->J_max = J_max;
    this->rate  = rate;
    this->count = count;
    this->pval  = pval;
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

uint64_t subinterval_bound::subinterval_count_equal( pos i, pos j, pos d_min_, pos d_max_ ) const
{
    uint64_t d_min = (uint64_t)d_min_;
    uint64_t d_max = (uint64_t)d_max_;

    if( i > j ) return 0;
    d_min_ = min<uint64_t>( d_min, 1 );
    uint64_t n = (uint64_t)(j - i + 1);

    if( n < d_min )     return 0;
    if( d_max < d_min ) return 0;

    uint64_t m = min( n, d_max );

    return (1+n)*(1+m-d_min) - (m*(m+1) - (d_min*(d_min-1))) / 2;
}


size_t subinterval_bound::subinterval_count( pos d_min, pos d_max ) const
{
    if( disjoint_bounds() ) {
        uint64_t T, L, R, C;
        T = subinterval_count_equal( I_min, J_max, d_min, d_max );
        L = subinterval_count_equal( I_min, J_min-1, d_min, d_max );
        R = subinterval_count_equal( I_max+1, J_max, d_min, d_max );
        C = subinterval_count_equal( I_max+1, J_min-1, d_min, d_max );

        return (size_t)(T + C - L - R);
    }
    else if( equal_bounds() ) {
        return (size_t)subinterval_count_equal( I_min, J_max, d_min, d_max );
    }
    else return 0;
}


pos subinterval_bound::subinterval_bound::max_length() const
{
    return J_max - I_min + 1;
}

pos subinterval_bound::subinterval_bound::min_length() const
{
    return J_min - I_max + 1;
}



#ifdef PEAKOLATOR_PVAL_HEURISTIC
bool subinterval_bound_priority::operator()( subinterval_bound*& a,
                                             subinterval_bound*& b )
{
    return a->pval > b->pval;
}
#else
bool subinterval_bound_priority::operator()( subinterval_bound*& a,
                                             subinterval_bound*& b )
{
    return a->count/a->rate < b->count/b->rate;

}
#endif


subinterval_bound_pqueue::subinterval_bound_pqueue( peakolator_model* model )
    : model(model)
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

void subinterval_bound_pqueue::conditional_push_copy(
        subinterval_bound& x,
        const mpfr_class& p_max )
{
    if( x.count > 0 &&
        x.min_length() <= model->params->d_max &&
        x.max_length() >= model->params->d_min )
    {
        x.pval = model->QX(
                model->context->min_rate( x, model->params->d_min ),
                x.count );
        if( x.pval < p_max ) push( new subinterval_bound(x) );
    }
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
    , pval(i.pval)
{
}

interval::interval( const interval& i )
    : seqname(strdup(i.seqname))
    , start(i.start)
    , end(i.end)
    , strand(i.strand)
    , pval(i.pval)
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
    this->seqname = strdup(seqname);
    this->start   = start+i.start;
    this->end     = start+i.end;
    this->strand  = strand;
    this->pval    = i.pval;
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

void subinterval_stack::print( const char* seqname, char strand, pos offset,
                            const char* gene_name )
{
    if( gene_name == NULL ) gene_name = "";

    subinterval_stack::iterator i;
    for( i = this->begin(); i != this->end(); i++ ) {
        mpfr_printf("%s\t%d\t%d\t%s\t%0.8Re\t%c\n",
               seqname, i->start+offset-1, i->end+offset-1,
               gene_name, i->pval.get_mpfr_t(), strand );
        fflush(stdout);
    }
}




