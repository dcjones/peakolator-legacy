
SUBDIRS=gmpfrxx samtools yaml-cpp

dist_bin_SCRIPTS = peakolate peakolate-seqbias
bin_PROGRAMS = peakolate-seqbias
pyexec_LTLIBRARIES = peakolator.la

peakolator_la_SOURCES = hash.h              hash.c \
						logger.h            logger.c \
						superfasthash.h     superfasthash.c \
						miscmath.hpp        miscmath.cpp \
						annotations.hpp     annotations.cpp \
						common.hpp          common.cpp \
						context.hpp         context.cpp \
						dataset.hpp         dataset.cpp \
						emppval.hpp         emppval.cpp \
						intervals.hpp       intervals.cpp \
						kmers.hpp           kmers.cpp \
						nulldist.hpp        nulldist.cpp \
						parameters.hpp      parameters.cpp \
						scanner.hpp         scanner.cpp \
						sequencing_bias.hpp sequencing_bias.cpp \
						mpfr_bindings.hpp \
						peakolator.cpp      peakolator.pyx \
						peakolator.pxd      sequencing_bias.pxd \
						sequencing_bias.pxi 


peakolator_la_CPPFLAGS = $(PYTHON_CPPFLAGS)
peakolator_la_LDFLAGS = -module -avoid-version
peakolator_la_LIBADD = gmpfrxx/libgmpfrxx.la \
					   samtools/libbam.la \
					   yaml-cpp/libyaml-cpp.la \
					   -lz \
					   -lgsl \
					   -lgslcblas \
					   -lgmp \
					   -lmpfr \
					   -lnlopt 
peakolator.cpp : peakolator.pyx
	cython --cplus -o $@ $^


peakolate_seqbias_SOURCES = peakolate_seqbias.pyx peakolate-seqbias.cpp \
							sequencing_bias.pxi sequencing_bias.pxd \
							sequencing_bias.hpp sequencing_bias.cpp \
							common.hpp          common.cpp \
							miscmath.hpp        miscmath.cpp \
							logger.h            logger.c \
							superfasthash.h     superfasthash.c \
							hash.h              hash.c \
							kmers.hpp           kmers.cpp


peakolate_seqbias_LDADD = samtools/libbam.la \
						  yaml-cpp/libyaml-cpp.la \
						  -lgsl \
						  -lgslcblas \
						  $(PYTHON_LDFLAGS) \
						  $(PYTHON_EXTRA_LDFLAGS)

peakolate_seqbias_CPPFLAGS = $(PYTHON_CPPFLAGS)


peakolate-seqbias.cpp : peakolate_seqbias.pyx sequencing_bias.pxi sequencing_bias.pxd
	cython --embed --cplus -o $@ $<

CLEANFILES = peakolator.cpp peakolate-seqbias.cpp

