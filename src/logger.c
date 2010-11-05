

#include "logger.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

const char* INDENT = "  ";

struct logger
{
    FILE* f;
    int indent;
    int verbosity;
};


struct logger* g_log = NULL;


void log_init(void)
{
    if( !g_log ) {
        g_log = malloc(sizeof(struct logger));
        g_log->indent = 1;
        g_log->verbosity = LOG_MSG;
        g_log->f = stderr;
    }
}


void log_indent()
{
    log_init();
    g_log->indent++;
}


void log_unindent()
{
    log_init();
    g_log->indent--;
}


void log_verbosity( int v )
{
    log_init();
    g_log->verbosity = v;
}


void log_puts( int vl, const char* msg )
{
    log_printf( vl, "%s",  msg );
}

void log_printf( int vl, const char* fmt, ... )
{
    log_init();

    if( vl > g_log->verbosity ) return;

    va_list ap;
    va_start(ap,fmt);

    char* fmt_ = malloc(sizeof(char) * (strlen(fmt)+g_log->indent*strlen(INDENT)+1));
    size_t m = strlen(INDENT);
    int i = g_log->indent;
    int j = 0;
    while( i-- ) {
        memcpy( fmt_+j, INDENT, m*sizeof(char) );
        j += m;
    }

    m = strlen(fmt);
    memcpy( fmt_+j, fmt, m*sizeof(char) );
    j += m;
    fmt_[j] = '\0';

    vfprintf( g_log->f, fmt_, ap );

    va_end(ap);
}



