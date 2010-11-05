/*
 *           hash
 *           A quick and dirty hash table implementation.
 *
 *           Daniel Jones <dcjones@cs.washington.edu>
 *           July 2010
 *
 */

#include "hash.h"


#define INITIAL_TABLE_SIZE 128
#define MAX_LOAD 0.75
#define MIN_LOAD 0.05 /* make sure this is less than MAX_LOAD/2 */



/* This is Jenkin's hash. The implementation is stolen from the Linux kernel,
 * modified only slightly.
 */

#define __jhash_mix(a, b, c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

/* The golden ration: an arbitrary value */
#define JHASH_GOLDEN_RATIO	0x9e3779b9

/* The most generic version, hashes an arbitrary sequence
 * of bytes.  No alignment or length assumptions are made about
 * the input key.
 */
uint32_t hash( const void *key, uint32_t length )
{
	uint32_t a, b, c, len;
	const uint8_t *k = key;

	len = length;
	a = b = JHASH_GOLDEN_RATIO;

    c = 0;

	while (len >= 12) {
		a += (k[0] +((uint32_t)k[1]<<8) +((uint32_t)k[2]<<16) +((uint32_t)k[3]<<24));
		b += (k[4] +((uint32_t)k[5]<<8) +((uint32_t)k[6]<<16) +((uint32_t)k[7]<<24));
		c += (k[8] +((uint32_t)k[9]<<8) +((uint32_t)k[10]<<16)+((uint32_t)k[11]<<24));

		__jhash_mix(a,b,c);

		k += 12;
		len -= 12;
	}

	c += length;
	switch (len) {
	case 11: c += ((uint32_t)k[10]<<24);
	case 10: c += ((uint32_t)k[9]<<16);
	case 9 : c += ((uint32_t)k[8]<<8);
	case 8 : b += ((uint32_t)k[7]<<24);
	case 7 : b += ((uint32_t)k[6]<<16);
	case 6 : b += ((uint32_t)k[5]<<8);
	case 5 : b += k[4];
	case 4 : a += ((uint32_t)k[3]<<24);
	case 3 : a += ((uint32_t)k[2]<<16);
	case 2 : a += ((uint32_t)k[1]<<8);
	case 1 : a += k[0];
	};

	__jhash_mix(a,b,c);

	return c;
}




void rehash( struct table* T, size_t new_n );



/* Create an empty hash table. */
void table_create( struct table* T )
{
    T->A = malloc(sizeof(struct hashed_value*)*INITIAL_TABLE_SIZE);
    memset( T->A, 0, sizeof(struct hashed_value*)*INITIAL_TABLE_SIZE );
    T->n = INITIAL_TABLE_SIZE;
    T->m = 0;
    T->max_m = T->n * MAX_LOAD;
    T->min_m = T->n * MIN_LOAD;
}



/* Remove all elements in the table. */
void table_clear( struct table* T )
{
    struct hashed_value* u;
    size_t i;
    for( i = 0; i < T->n; i++ ) {
        while( T->A[i] ){
            u = T->A[i]->next;
            free(T->A[i]);
            T->A[i] = u;
        }
    }
    T->m = 0;
}



/* Free all memory associated with a table. */
void table_destroy( struct table* T )
{
    table_clear(T);
    free(T->A);
}



void table_inc( struct table* T, bam1_t* read )
{
    if( T->m == T->max_m ) rehash( T, T->n*2 );

    struct read_pos pos;
    pos.tid = read->core.tid;
    /*                                      XXX: I don't know why I must
     *                                      subtract one here. It bothers me
     *                                      that I don't, but the results do not
     *                                      come out correctly otherwise.
     *                                      Possibly, the function gives the
     *                                      nucleotide immediately after the
     *                                      read. */
    if( bam1_strand(read) ) pos.pos = bam_calend( &read->core, bam1_cigar(read) ) - 1;
    else                    pos.pos = read->core.pos;
    pos.strand = bam1_strand(read);

    uint32_t h = hash((void*)&pos, sizeof(struct read_pos)) % T->n;

    struct hashed_value* u = T->A[h];

    while(u) {
        if( memcmp( &u->pos, &pos, sizeof(struct read_pos) ) == 0 ) {
            u->count++;
            return;
        }

        u = u->next;
    }

    u = malloc(sizeof(struct hashed_value));
    memcpy( &u->pos, &pos, sizeof(struct read_pos) );

    u->count = 1;

    u->next = T->A[h];
    T->A[h] = u;

    T->m++;
}



/* Insert existing entries without copying sequences. Used for rehashing. */
bool table_insert_without_copy( struct table* T, struct hashed_value* V )
{
    if( T->m == T->max_m ) rehash( T, T->n*2 );

    uint32_t h = hash((void*)&V->pos,sizeof(struct read_pos)) % T->n;

    V->next = T->A[h];
    T->A[h] = V;

    T->m++;

    return true;
}


/* Rezise the table T to new_n. */
void rehash( struct table* T, size_t new_n )
{
    struct table U;
    U.n = new_n;
    U.m = 0;
    U.A = malloc( sizeof(struct hashed_value*) * U.n );
    memset( U.A, 0, sizeof(struct hashed_value*) * U.n );


    struct hashed_value *j,*k;
    size_t i;
    for( i = 0; i < T->n; i++ ) {
        j = T->A[i];
        while( j ) {
            k = j->next;
            table_insert_without_copy( &U, j );
            j = k;
        }
        T->A[i] = NULL;
    }

    free(T->A);
    T->A = U.A;
    T->n = U.n;
    T->max_m = T->n*MAX_LOAD;
    T->min_m = T->n*MIN_LOAD;
}

int compare_count( const void* x, const void* y )
{
    struct hashed_value* const * a = x;
    struct hashed_value* const * b = y;

    return (*a)->count - (*b)->count;
}

int compare_pos( const void* x, const void* y )
{
    struct hashed_value* const * a = x;
    struct hashed_value* const * b = y;


    int c = (*a)->pos.tid - (*b)->pos.tid;
    if( c == 0 ) {
        if( (*a)->pos.pos == (*b)->pos.pos ) return 0;
        else return (*a)->pos.pos > (*b)->pos.pos ? 1 : -1;
    }
    else return c;
}


void sort_table( struct table* T,
                 struct hashed_value*** S_,
                 int(*compar)(const void *, const void *) )
{
    struct hashed_value** S = malloc( sizeof(struct hashed_value*) * T->m );
    memset( S, 0, sizeof(struct hashed_value*) * T->m );

    struct hashed_value* j;
    size_t i,k;
    for( i=0, k=0; i < T->n; i++ ) {
        j = T->A[i];
        while( j ) {
            S[k] = j;
            k++;
            j = j->next;
        }
    }

    qsort( S, T->m, sizeof(struct hashed_value*), compar );

    *S_ = S;
}

void table_sort_by_count( struct table* T,
                    struct hashed_value*** S_ )
{
    sort_table( T, S_, compare_count );
}

void table_sort_by_position( struct table* T,
                    struct hashed_value*** S_ )
{
    sort_table( T, S_, compare_pos );
}


bool table_member( struct table* T, bam1_t* read )
{
    struct read_pos pos;
    pos.tid = read->core.tid;
    pos.pos = read->core.pos;
    pos.strand = bam1_strand(read);

    uint32_t h = hash((void*)&pos, sizeof(struct read_pos)) % T->n;

    struct hashed_value* u = T->A[h];

    while(u) {
        if( memcmp( &u->pos, &pos, sizeof(struct read_pos) ) == 0 ) {
            return true;
        }

        u = u->next;
    }

    return false;
}



void rehash_tail( struct table* T, int32_t q1, int32_t q2 )
{
    struct table U;
    U.n = T->n;
    U.m = 0;
    U.A = malloc( sizeof(struct hashed_value*) * U.n );
    memset( U.A, 0, sizeof(struct hashed_value*) * U.n );

    struct hashed_value** S;
    table_sort_by_count( T, &S );

    int32_t i;
    for( i = q1; i < q2; i++ ) {
        table_insert_without_copy( &U, S[i] );
    }

    /* free the rest */
    for( ; i < T->m; i++ ) free(S[i]);

    free(S);
    free(T->A);
    T->A = U.A;
    T->m = U.m;
}

