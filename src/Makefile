
include ../config.inc

OBJ=common.o \
	logger.o \
	annotations.o \
	peakolator_dist.o \
	peakolator_interval.o \
	peakolator_model.o \
	peakolator_dataset.o \
	peakolator_context.o \
	peakolator_pval.o \
	bias_correction.o \
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

