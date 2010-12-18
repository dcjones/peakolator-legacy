
#include "annotations.hpp"
#include "logger.h"

#include <cstdio>
#include <cstring>
#include <set>
#include <algorithm>

using namespace std;

gtf_row::gtf_row()
    : start(-1)
    , end(-1)
    , strand('.')
    , frame(-1)
{}    

gtf_row::gtf_row( const gtf_row& r )
    : seqname(r.seqname)
    , source(r.source)
    , feature(r.feature)
    , start(r.start)
    , end(r.end)
    , score(r.score)
    , strand(r.strand)
    , frame(r.frame)
    , attributes(r.attributes)
{

}

gtf_row::~gtf_row()
{
}


rows::rows()
{
}

rows::~rows()
{
    rows::iterator i;
    for( i = begin(); i != end(); i++ ) {
        delete *i;
        *i = NULL;
    }
}


void rows::build_gene_id_index()
{
    rows::iterator i;
    for( i = begin(); i != end(); i++ ) {
        // skip unlabelled annotations
        if( (*i)->attributes.find( "gene_id" ) ==
                (*i)->attributes.end() ) continue;

        gene_id_index.insert( pair<string,gtf_row*>(
                                (*i)->attributes["gene_id"], *i ) );
    }
}

void rows::build_transcript_id_index()
{
    rows::iterator i;
    for( i = begin(); i != end(); i++ ) {
        // skip unlabelled annotations
        if( (*i)->attributes.find( "transcript_id" ) ==
                (*i)->attributes.end() ) continue;

        transcript_id_index.insert( pair<string,gtf_row*>(
                                (*i)->attributes["transcript_id"], *i ) );
    }
}


gtf_row* new_gtf_row_from_gtf_line( char* line )
{
    if( line == NULL || line[0] == '\0' || line[0] == '\n' ) return NULL;

    gtf_row* r = new gtf_row();
    int k = 0;
    char* saveptr;
    char *s, *c;

    char *key_start, *key_end, *val_start, *val_end;
    string key, val;

    s = strtok_r( line, "\t\n", &saveptr );
    while( s ) {

        switch( k ) {
            case 0: // seqname
                r->seqname = s;
                break;
            case 1:  // source
                r->source = s;
                break;
            case 2: // feature
                r->feature = s;
                break;
            case 3: // start
                r->start = (pos)atoi(s) - 1; // make 0-based
                break;
            case 4: //end
                r->end   = (pos)atoi(s) - 1; // make 0-based
                break;
            case 5: // score
                r->score = s;
                break;
            case 6: // strand
                r->strand = (s[0] == '+' || s[0] == '-') ? s[0] : '.';
                break;
            case 7: // frame
                r->frame = s[0] == '.' ? -1 : atoi(s);
                break;
            default: // (k >= 8) attributes
                key_start = key_end = val_start = val_end = NULL;
                for( c = s; *c; c++ ) {
                    if( key_start == NULL ) {
                        if( !isspace(*c) ) key_start = c;
                    }
                    else if( key_end == NULL ) {
                        if( isspace(*c) ) key_end = c;
                    }
                    else if( val_start == NULL ) {
                        if( !isspace(*c) && *c != '\"' ) val_start = c;
                    }
                    else if( val_end == NULL ) {
                        if( isspace(*c) || *c == '\"' ) val_end = c;
                    }
                    else break;
                }

                if( key_start == NULL ||
                    key_end   == NULL ||
                    val_start == NULL ||
                    val_end   == NULL )
                {
                    log_printf( LOG_WARN, "Malformed GTF attribute: %s\n", s );
                    break;
                }

                key.assign( key_start, key_end - key_start );
                val.assign( val_start, val_end - val_start );

                r->attributes[key] = val;
        }

        k++;
        if( k < 8 ) s = strtok_r( NULL, "\t\n", &saveptr );
        else        s = strtok_r( NULL, ";\n", &saveptr );
    }

    if( k < 9 ) {
        log_printf( LOG_WARN, "Malformed GTF line: %s\n", line );
        delete r;
        return NULL;
    }

    return r;
}



gtf_row* new_gtf_row_from_bed_line( char* line )
{
    if( line == NULL || line[0] == '\0' || line[0] == '\n' ) return NULL;

    gtf_row* r = new gtf_row();
    int k = 0;
    char* saveptr;
    char *s;


    s = strtok_r( line, "\t\n", &saveptr );
    while( s && k <= 5 ) {

        switch( k ) { 
            case 0: // seqname
                r->seqname = s;
                break;
            case 1: // start
                r->start = atoi(s);
                break;
            case 2: // end
                r->end = atoi(s) - 1; // make end inclusive
                break;
            case 3: // name/feature
                r->feature = s;
                break;
            case 4: // score
                r->score = s;
                break;
            case 5:
                r->strand = (s[0] == '+' || s[0] == '-') ? s[0] : '.';
                break;
        }

        k++;
        s = strtok_r( NULL, "\t\n", &saveptr );
    }

    if( k < 3 ) {
        log_printf( LOG_WARN, "Malformed BED line: %s\n", line );
        delete r;
        return NULL;
    }

    return r;
}




rows* parse_file( const char* fn, gtf_row* (*parse_line)( char* ) )
{
    log_printf( LOG_MSG, "parsing %s ... ", fn );
    rows* rs = new rows;
    gtf_row* r;
    const size_t max_line = 4096;
    char* line = new char[max_line];

    FILE* f = fopen( fn, "r" );
    if( f == NULL ) {
        log_printf( LOG_ERROR, "Can\'t open file %s\n", fn );
        exit(1);
    }

    while( fgets( line, max_line, f ) ) {
        r = parse_line( line );
        if( r ) {
            if( r->feature == "exon" ) rs->push_back( r );
            else delete r;
        }
    }

    fclose( f );

    delete[] line;
    log_printf( LOG_MSG, "done. (%zu annotations)\n", rs->size() );
    return rs;
}


rows* parse_gtf( const char* fn ) {
    return parse_file( fn, new_gtf_row_from_gtf_line );
}


rows* parse_bed( const char* fn ) {
    return parse_file( fn, new_gtf_row_from_bed_line );
}


bool cmp_gtf_row_stranded( const gtf_row* a, const gtf_row* b )
{
    if( a->seqname < b->seqname ) return true;
    else if( a->strand < b->strand ) return true;
    else return a->start < b->start;
}

bool cmp_gtf_row_unstranded( const gtf_row* a, const gtf_row* b )
{
    if( a->seqname < b->seqname ) return true;
    else return a->start < b->start;
}

rows* gtf_intersect_rows( rows* r1, rows* r2, bool stranded )
{
    rows* rs = new rows;
    gtf_row *r, *ri, *rj;

    bool (*cmp_gtf_row)( const gtf_row*, const gtf_row* );
    cmp_gtf_row = stranded ? cmp_gtf_row_stranded : cmp_gtf_row_unstranded;

    sort( r1->begin(), r1->end(), cmp_gtf_row );
    sort( r2->begin(), r2->end(), cmp_gtf_row );

    size_t i, j;
    i = j = 0;

    while( i < r1->size() && j < r2->size() ) {
        ri = (*r1)[i];
        rj = (*r2)[j];

        if(      stranded && ri->strand < rj->strand ) i++;
        else if( stranded && ri->strand > rj->strand ) j++;
        else if( ri->seqname < rj->seqname ) i++;
        else if( ri->seqname > rj->seqname ) j++;
        else if( ri->end   < rj->start ) i++;
        else if( ri->start > rj->end )   j++;

        /* intersection case */
        else {
            r = new gtf_row;
            r->seqname    = ri->seqname;
            r->source     = ri->source;
            r->feature    = ri->feature;
            r->start      = max( ri->start, rj->start );
            r->end        = min( ri->end, rj->end );
            r->score      = ri->score;
            r->strand     = stranded ? ri->strand : '.';
            r->frame      = ri->frame;
            r->attributes = ri->attributes;

            if     ( ri->end < rj->end ) i++;
            else if( ri->end > rj->end ) j++;
            else {
                i++;
                j++;
            }

            rs->push_back(r);
        }
    }

    return rs;
}


rows* get_constitutive_exons( rows* rs, bool stranded )
{
    log_puts( LOG_MSG, "getting genes and transcripts ... " );

    /* build an index mapping gene_ids to sets of transcript_ids */
    map< string, set<string> > gene_transcripts;
    rows::iterator i;
    for( i = rs->begin(); i != rs->end(); i++ ) {
        if( (*i)->attributes.find( "gene_id" ) ==
                (*i)->attributes.end() ) continue;
        if( (*i)->attributes.find( "transcript_id" ) ==
                (*i)->attributes.end() ) continue;

        gene_transcripts[ (*i)->attributes["gene_id"] ].insert(
                (*i)->attributes["transcript_id"] );
    }

    /* index rows by transcript_id */
    rs->build_transcript_id_index();

    /* count the number of unique transcripts */
    size_t trans_count = 0;
    map< string, set<string> >::iterator j;
    for( j = gene_transcripts.begin(); j != gene_transcripts.end(); j++ ) {
        trans_count += j->second.size();
    }

    log_printf( LOG_MSG, "done. (%zu genes, %zu transcripts)\n",
                gene_transcripts.size(),
                trans_count );


    rows gene_exons;
    rows* exons = new rows;
    log_puts( LOG_MSG, "getting constitutive exons ... " );

    for( j = gene_transcripts.begin(); j != gene_transcripts.end(); j++ ) {

        /* Intersect everything in j->second.
         * (The next fifty lines of horendous code
         * replaces a python one liner. Gah.) */

        gene_exons.clear();
        set<string>::iterator trans = j->second.begin();

        rows *u, *v, *w;

        w = new rows;

        multimap<string,gtf_row*>::iterator k, k_last;

        k = rs->transcript_id_index.lower_bound( *trans );
        k_last = rs->transcript_id_index.upper_bound( *trans );


        for( ; k != k_last; k++ ) {
            w->push_back( new gtf_row( *(k->second) ) );
        }


        trans++;
        while( trans != j->second.end() ) {
            u = w;

            /* load next v from transcript */
            k = rs->transcript_id_index.lower_bound( *trans );
            k_last = rs->transcript_id_index.upper_bound( *trans );

            v = new rows;

            for( ; k != k_last; k++ ) {
                v->push_back( new gtf_row( *(k->second) ) );
            }

            /* intersect into w */
            w = gtf_intersect_rows( u, v );

            delete u;
            delete v;

            trans++;
        }


        /* add everything remaining in w to exons */
        for( i = w->begin(); i != w->end(); i++ ) {
            exons->push_back(*i);
        }
        w->clear();
        delete w;
    }

    exons->build_gene_id_index();
    log_printf( LOG_MSG, "done. (%zu constitutive exons)", exons->size() );

    return exons;
}



