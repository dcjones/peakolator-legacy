
include ../config.inc

OBJ= \
	peakolator.o \
	common.o \
	logger.o \
	annotations.o \
	nulldist.o \
	intervals.o \
	parameters.o \
	scanner.o \
	dataset.o \
	context.o \
	emppval.o \
	sequencing_bias.o \
	kmers.o \
	superfasthash.o \
	hash.o \
	gmpfrxx/libgmpfrxx.a \
	samtools/libbam.a \
	yaml-cpp/libyaml-cpp.a

SEQBIAS_OBJ= \
	peakolate-seqbias.o \
	sequencing_bias.o \
	common.o \
	logger.o \
	superfasthash.o \
	hash.o \
	kmers.o \
	samtools/libbam.a \
	yaml-cpp/libyaml-cpp.a


all : peakolator.so peakolate-seqbias

clean :
	rm -f ${OBJ}
	rm -f peakolator.so
	rm -f peakolator.cpp
	cd samtools ; ${MAKE} clean
	cd gmpfrxx ; ${MAKE} clean
	cd yaml-cpp ; ${MAKE} clean

.c.o :
	gcc ${CFLAGS} ${INC} -c $<

.cpp.o :
	g++ ${CPPFLAGS} ${INC} -c $<

peakolator.o : peakolator.cpp
	g++ ${CPPFLAGS} ${INC} -c $<

# python module
peakolator.cpp : peakolator.pyx
	cython --cplus -o $@ $< 

peakolate-seqbias.cpp : peakolate_seqbias.pyx sequencing_bias.pxd sequencing_bias.pxi
	cython --embed --cplus -o $@ $<

peakolator.so : ${OBJ}
	g++ ${CPPFLAGS} ${INC} -shared -o $@ $^ ${GSL_LIBS} ${PY_LIBS} ${MATH_LIBS} -lz

peakolate-seqbias : ${SEQBIAS_OBJ}
	g++ ${CPPFLAGS} ${INC} -o $@ $^ ${PY_LIBS} ${GSL_LIBS}

samtools/libbam.a :
	cd samtools ; ${MAKE} libbam.a

gmpfrxx/libgmpfrxx.a :
	cd gmpfrxx ; ${MAKE} libgmpfrxx.a

yaml-cpp/libyaml-cpp.a :
	cd yaml-cpp ; ${MAKE} libyaml-cpp.a

