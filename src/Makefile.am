
SUBDIRS=samtools yaml-cpp

dist_bin_SCRIPTS = peakolate
bin_PROGRAMS = peakolate-seqbias
pyexec_LTLIBRARIES = peakolator.la

peakolator_la_SOURCES = table.h             table.c \
						logger.h            logger.c \
						ushuffle.h          ushuffle.c \
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
						peakolator.pyx      peakolator.pxd \
						sequencing_bias.pxd sequencing_bias.pxi 

nodist_peakolator_la_SOURCES = peakolator.cpp

peakolator_la_CPPFLAGS = $(PYTHON_CPPFLAGS)
peakolator_la_LDFLAGS = -module -avoid-version
peakolator_la_LIBADD = samtools/libbam.la \
					   yaml-cpp/libyaml-cpp.la \
					   -lz \
					   -lgsl \
					   -lgslcblas \
					   -lnlopt 

peakolator.cpp : peakolator.pyx sequencing_bias.pxi sequencing_bias.pxd
	cython --cplus -o $@ $<


peakolate_seqbias_SOURCES = peakolate_seqbias.pyx \
							sequencing_bias.pxi sequencing_bias.pxd \
							sequencing_bias.hpp sequencing_bias.cpp \
							common.hpp          common.cpp \
							miscmath.hpp        miscmath.cpp \
							logger.h            logger.c \
							table.h             table.c \
							ushuffle.h          ushuffle.c \
							kmers.hpp           kmers.cpp \
							intervals.hpp       intervals.cpp

nodist_peakolate_seqbias_SOURCES = peakolate-seqbias.cpp


peakolate_seqbias_LDADD = samtools/libbam.la \
						  yaml-cpp/libyaml-cpp.la \
						  -lz \
						  -lgsl \
						  -lgslcblas \
						  $(PYTHON_LDFLAGS) \
						  $(PYTHON_EXTRA_LDFLAGS)

peakolate_seqbias_CPPFLAGS = $(PYTHON_CPPFLAGS)


peakolate-seqbias.cpp : peakolate_seqbias.pyx sequencing_bias.pxi sequencing_bias.pxd
	cython --embed --cplus -o $@ $<

CLEANFILES = peakolator.cpp peakolate-seqbias.cpp

