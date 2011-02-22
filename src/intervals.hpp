

#ifndef PEAKOLATOR_INTERVAL
#define PEAKOLATOR_INTERVAL

#include "common.hpp"

#include <deque>
#include <queue>



/* subinterval: an interval with relative coordinates and no strand or
 * chromosome */
class subinterval
{
    public:

        subinterval();
        subinterval( pos start, pos end,
                  double       score = 0.0,
                  rcount       count = 0,
                  unsigned int zeros = 0,
                  double       rate  = 1.0 );

        subinterval( const subinterval& y );

        void operator=( const subinterval& y );
        bool operator<( const subinterval& y ) const;
        bool operator>( const subinterval& y ) const;

        void set( pos start, pos end,
                  double       score = 1.0,
                  rcount       count = 0,
                  unsigned int zeros = 0,
                  double       rate  = 1.0 );

        pos  length() const;

        pos start;
        pos end;
        double       score;
        rcount       count;
        unsigned int zeros;
        double       rate;
};


class subinterval_ordering
{
    public:
        bool operator()( const subinterval&, const subinterval& );
};


class subinterval_stack : public std::deque<subinterval>
{
    public:
        subinterval_stack();
        void print( const char* seqname, char strand, pos offset,
                    const char* gene_name = NULL );
};


/* a bound on the start and end position of a subinterval */
class subinterval_bound
{
    public:

        subinterval_bound();

        subinterval_bound( pos min, pos max,
                           double       rate  = 1.0,
                           rcount       count = 0,
                           unsigned int zeros = 0,
                           double       score = 0.0 );

        subinterval_bound( pos I_min, pos I_max, pos J_min, pos J_max,
                           double       rate  = 1.0,
                           rcount       count = 0,
                           unsigned int zeros = 0,
                           double       score = 0.0 );

        subinterval_bound( const subinterval_bound& );

        void operator=( const subinterval_bound& );

        void set( pos min, pos max,
                  double       rate  = 1.0,
                  rcount       count = 0,
                  unsigned int zeros = 0,
                  double       score = 0.0 );

        void set( pos I_min, pos I_max, pos J_min, pos J_max,
                  double       rate  = 1.0,
                  rcount       count = 0,
                  unsigned int zeros = 0,
                  double       score = 0.0 );

        bool   disjoint_bounds()   const;
        bool   equal_bounds()      const;
        double subinterval_count( pos d_min, pos d_max ) const;

        pos max_length() const;
        pos min_length() const;

        pos I_min;
        pos I_max;
        pos J_min;
        pos J_max;
        double       rate;
        rcount       count;
        unsigned int zeros;
        double       score;

    private:
        void check() const;

        double subinterval_count_equal( pos i, pos j, pos d_min, pos d_max ) const;
};


class subinterval_bound_priority
{
    public:
        bool operator()( subinterval_bound*&, subinterval_bound*& );
};


class subinterval_bound_pqueue : private std::vector<subinterval_bound*>
{
    public:
        subinterval_bound_pqueue();
        ~subinterval_bound_pqueue();

        void push( subinterval_bound* );
        subinterval_bound* pop();
        using std::vector<subinterval_bound*>::empty;

        pos total_length() const;
};


/* interval: an absolute interval without extra count, rate data */
class interval
{
    public:
        interval();
        interval( const char* seqname, pos start, pos end, int strand );
        interval( const subinterval&, const char* seqname, pos start, int strand );
        interval( const interval& );
        ~interval();

        void set( const char* seqname, pos start, pos end, int strand );
        void set( const subinterval&, const char* seqname, pos start, int strand );

        pos length() const;

        char* seqname;
        pos start, end;
        int strand;
        double score;
};

class interval_stack : public std::deque<interval>
{
    public:
        interval_stack();
        interval_stack( const subinterval_stack& sis, const char* chrom, pos start, int strand );
        ~interval_stack();
};

/* this makes things easier for cython */
void interval_stack_push( interval_stack*, char* seqname, pos start, pos end, int strand );

#endif
