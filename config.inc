
//BASE_FLAGS=-fPIC -Wall -g -O0
BASE_FLAGS=-fPIC -Wall -g -O3 -DGSL_RANGE_CHECK_OFF -DHAVE_INLINE

CFLAGS=$(shell gsl-config --cflags) \
       $(shell python-config --cflags) \
       --std=c99 $(BASE_FLAGS)

CPPFLAGS=$(shell gsl-config --cflags) \
         $(shell python-config --cflags | sed s/-Wstrict-prototypes// ) \
         $(BASE_FLAGS)

LIBS=$(shell gsl-config --libs) $(shell python-config --libs) -lz -lgmp -lmpfr -lnlopt
INC=$(shell python-config --includes)


