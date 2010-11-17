
include ../config.inc

OBJ=common.o \
	logger.o \
	annotations.o \
	nulldist.o \
	intervals.o \
	parameters.o \
	model.o \
	dataset.o \
	context.o \
	emppval.o \
	sequencing_bias.o \
	kmers.o \
	superfasthash.o \
	hash.o \
	gmpfrxx/libgmpfrxx.a \
	samtools/libbam.a

all : peakolator.so

clean :
	rm -f ${OBJ}
	rm -f peakolator.so
	rm -f peakolator_module.cpp
	rm -f plot_aggregate_bias
	rm -f plot_bias_tracks
	cd samtools ; ${MAKE} clean
	cd gmpfrxx ; ${MAKE} clean

.c.o :
	gcc ${CFLAGS} ${INC} -c $<

.cpp.o :
	g++ ${CPPFLAGS} ${INC} -c $<


# python module
peakolator.cpp : peakolator.pyx
	cython --cplus -o $@ $< 

peakolator.so : peakolator.cpp ${OBJ}
	g++ ${CPPFLAGS} ${INC} -shared -o $@ $^ ${LIBS}

samtools/libbam.a :
	cd samtools ; ${MAKE} libbam.a

gmpfrxx/libgmpfrxx.a :
	cd gmpfrxx ; ${MAKE} libgmpfrxx.a

