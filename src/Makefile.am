
SUBDIRS=gmpfrxx samtools yaml-cpp

dist_bin_SCRIPTS = peakolate peakolate-seqbias
pyexec_LTLIBRARIES = peakolator.la

peakolator_la_SOURCES = hash.h              hash.c \
						logger.h            logger.c \
						superfasthash.h     superfasthash.c \
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
						peakolator.cpp      peakolator.pyx \
						mpfr_bindings.hpp

peakolator_la_CPPFLAGS = $(PYTHON_INCLUDES)
peakolator_la_LDFLAGS = -module -avoid-version
peakolator_la_LIBADD = gmpfrxx/libgmpfrxx.la \
					   samtools/libbam.la \
					   yaml-cpp/libyaml-cpp.la \
					   -lz \
					   -lgmp \
					   -lmpfr \
					   -lgsl \
					   -lgslcblas \
					   -lnlopt


peakolator.cpp : peakolator.pyx
	cython -cplus -o $@ $^

CLEANFILES = peakolator.cpp

