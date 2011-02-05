
#include "table.h"



/* From: http://planetmath.org/encyclopedia/GoodHashTablePrimes.html */
static const size_t num_primes = 26;
static const uint32_t primes[] = {
/*                                 53,         97, */
/*          193,       389,       769,       1543, */
         3079,      6151,     12289,      24593,
        49157,     98317,    196613,     393241,
       786433,   1572869,   3145739,    6291469,
     12582917,  25165843,  50331653,  100663319,
    201326611, 402653189, 805306457, 1610612741 };


static const double max_load = 0.75;




/* From Thomas Wang (http://www.cris.com/~Ttwang/tech/inthash.htm) */
uint32_t hash( uint32_t a)
{
    a = (a ^ 61) ^ (a >> 16);
    a = a + (a << 3);
    a = a ^ (a >> 4);
    a = a * 0x27d4eb2d;
    a = a ^ (a >> 15);
    return a;
}



/* simple quadratic probing */
uint32_t probe( uint32_t h, uint32_t i )
{
    const double c1 = 0.5;
    const double c2 = 0.5;

    return h
           + (uint32_t)(c1 * (double)i)
           + (uint32_t)(c2 * (double)(i*i));
}


void subtable_create( struct subtable* T )
{
    T->n = 0;
    T->m = 0;
    T->A = malloc(sizeof(struct hashed_value)*primes[T->n]);
    size_t i;
    for( i = 0; i < primes[T->n]; i++ ) {
        T->A[i].pos = -1;
        T->A[i].count = 0;
    }
    T->max_m = (size_t)(((double)primes[T->n]) * max_load);
}


void subtable_destroy( struct subtable* T )
{
    free( T->A );
    T->A = NULL;
}


void subtable_rehash( struct subtable* T, size_t new_n );


bool subtable_inc( struct subtable* T, int32_t pos )
{
    if( T->m == T->max_m ) subtable_rehash( T, T->n + 1 );

    uint32_t h = hash( (uint32_t)pos );
    uint32_t i = 1;
    uint32_t j = h % primes[T->n];

    while( T->A[j].pos != -1 && T->A[j].pos != pos ) {
        j = probe( h, ++i ) % primes[T->n];
    }

    if( T->A[j].pos == -1 ) {
        T->A[j].pos = pos;
        T->A[j].count = 1;
        T->m++;
        return true;
    }
    else {
        T->A[j].count++;
        return false;
    }
}


void subtable_set( struct subtable* T, int32_t pos, uint32_t count )
{
    if( T->m == T->max_m ) subtable_rehash( T, T->n + 1 );

    uint32_t h = hash( (uint32_t)pos );
    uint32_t i = 1;
    uint32_t j = h % primes[T->n];

    while( T->A[j].pos != -1 && T->A[j].pos != pos ) {
        j = probe( h, ++i ) % primes[T->n];
    }

    if( T->A[j].pos == -1 ) {
        T->A[j].pos = pos;
        T->A[j].count = count;
    }
    else {
        T->A[j].count = count;
    }
}



void subtable_rehash( struct subtable* T, size_t new_n )
{
    if( new_n >= num_primes ) {
        fail( "a table has grown too large" );
    }

    struct subtable U;
    U.n = new_n;
    U.A = malloc( sizeof(struct hashed_value) * primes[U.n] );
    size_t i;
    for( i = 0; i < primes[U.n]; i++ ) {
        U.A[i].pos = -1;
        U.A[i].count = 0;
    }

    U.max_m = (size_t)(((double)primes[U.n]) * max_load);


    for( i = 0; i < primes[T->n]; i++ ) {
        if( T->A[i].pos == -1 ) continue;
        subtable_set( &U, T->A[i].pos, T->A[i].count );
    }

    free(T->A);
    T->A = U.A;
    T->n = U.n;
    T->max_m = U.max_m;
}



uint32_t subtable_count( struct subtable* T, int32_t pos )
{
    uint32_t h = hash( (uint32_t)pos );
    uint32_t i = 1;
    uint32_t j = h % primes[T->n];

    while( T->A[j].pos != -1 && T->A[j].pos != pos ) {
        j = probe( h, ++i ) % primes[T->n];
    }

    if( T->A[j].pos == pos ) return T->A[j].count;
    else                     return 0;
}


void table_create( struct table* T, size_t n )
{
    T->seq_names = NULL;
    T->n = n;
    T->m = 0;

    T->ts[0] = malloc( n * sizeof(struct subtable) );
    T->ts[1] = malloc( n * sizeof(struct subtable) );

    size_t i, j;
    for( i = 0; i <= 1; i++ ) {
        for( j = 0; j < n; j++ ) {
            subtable_create( &T->ts[i][j] );
        }
    }
}


void table_destroy( struct table* T )
{
    size_t i, j;
    for( i = 0; i <= 1; i++ ) {
        for( j = 0; j < T->n; j++ ) {
            subtable_destroy( &T->ts[i][j] );
        }
    }

    free( T->ts[0] );
    free( T->ts[1] );

    T->n = 0;
}



void table_inc( struct table* T, bam1_t* read )
{
    int32_t pos;
    if( bam1_strand(read) ) pos = bam_calend( &read->core, bam1_cigar(read) ) - 1;
    else                    pos = read->core.pos;

    table_inc_pos( T, read->core.tid, pos, bam1_strand(read) );
}



void table_inc_pos( struct table* T, int32_t tid, int32_t pos, uint32_t strand )
{
    if( tid < 0 || tid >= T->n ) return;
    if( subtable_inc( &T->ts[strand][tid], pos ) ) T->m++;
}


uint32_t table_count( struct table* T, bam1_t* read )
{
    int32_t pos;
    if( bam1_strand(read) ) pos = bam_calend( &read->core, bam1_cigar(read) ) - 1;
    else                    pos = read->core.pos;

    return table_count_pos( T, read->core.tid, pos, bam1_strand(read) );
}


uint32_t table_count_pos( struct table* T, int32_t tid, int32_t pos, uint32_t strand )
{
    if( tid < 0 || tid >= T->n ) return 0;
    return subtable_count( &T->ts[strand][tid], pos );
}


void table_dump( struct table* T, struct read_pos** A_, size_t* N_ )
{
    struct read_pos* A;
    size_t N = 0;
    size_t i, j;


    for( i = 0; i <= 1; i++ ) {
        for( j = 0; j < T->n; j++ ) {
            N += T->ts[i][j].m;
        }
    }

    A = malloc( N * sizeof(struct read_pos) );


    size_t u = 0;
    size_t v;
    for( i = 0; i <= 1; i++ ) {
        for( j = 0; j < T->n; j++ ) {
            for( v = 0; v < primes[T->ts[i][j].n]; v++ ) {
                if( T->ts[i][j].A[v].pos != -1 ) {
                    A[u].tid = j;
                    A[u].strand = i;
                    A[u].pos    = T->ts[i][j].A[v].pos;
                    A[u].count  = T->ts[i][j].A[v].count;
                    u++;
                }
            }
        }
    }

    *A_ = A;
    *N_ = N;
}






