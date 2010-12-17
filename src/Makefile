
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

all : peakolator.so

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


# python module
peakolator.cpp : peakolator.pyx
	cython --cplus -o $@ $< 

peakolator.so : ${OBJ}
	g++ ${CPPFLAGS} ${INC} -shared -o $@ $^ ${LIBS}

samtools/libbam.a :
	cd samtools ; ${MAKE} libbam.a

gmpfrxx/libgmpfrxx.a :
	cd gmpfrxx ; ${MAKE} libgmpfrxx.a

yaml-cpp/libyaml-cpp.a :
	cd yaml-cpp ; ${MAKE} libyaml-cpp.a
