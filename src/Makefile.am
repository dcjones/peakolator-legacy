
SUBDIRS=gmpfrxx samtools yaml-cpp

pyexec_LTLIBRARIES = peakolator.la

peakolator_la_SOURCES = hash.c \
						logger.c \
						superfasthash.c \
						annotations.cpp \
						common.cpp \
						context.cpp \
						dataset.cpp \
						emppval.cpp \
						intervals.cpp \
						kmers.cpp \
						nulldist.cpp \
						parameters.cpp \
						peakolator.cpp \
						scanner.cpp \
						sequencing_bias.cpp

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

