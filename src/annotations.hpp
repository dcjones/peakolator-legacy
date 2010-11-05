
#ifndef PEAKOLATOR_ANNOTATIONS
#define PEAKOLATOR_ANNOTATIONS

#include "common.hpp"

#include <string>
#include <map>
#include <deque>

class gtf_row
{
    public:
        gtf_row();
        gtf_row( const gtf_row& );
        ~gtf_row();

        std::string seqname;
        std::string source;
        std::string feature;
        pos start;
        pos end;
        std::string score;
        char strand;
        int frame;
        std::map<std::string,std::string> attributes;
};


class rows : public std::deque<gtf_row*>
{
    public:
        rows();
        ~rows();

        void build_gene_id_index();
        void build_transcript_id_index();

        std::multimap<std::string,gtf_row*> gene_id_index;
        std::multimap<std::string,gtf_row*> transcript_id_index;
};


rows* parse_gtf( const char* fn );
rows* parse_bed( const char* fn );

rows* gtf_intersect_rows( rows* r1, rows* r2, bool stranded = true );
rows* get_constitutive_exons( rows* rs, bool stranded );
//rows* get_transcript_extents( rows* r1 );
//rows* get_interval_gaps( rows* r1 );

// TODO: and so on



#endif
